/* Copyright (c) 2006, National ICT Australia
 * All rights reserved.
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * Authors      : Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 * Created      : 17/11/2006
 * Last Updated :
 */

#ifndef _SPARSEVEC_CPP_
#define _SPARSEVEC_CPP_

#include "sparsevec.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <fstream>

using namespace std;


CSparseVec::CSparseVec():len(0),fidx(0),fval(0)
{
}


CSparseVec::CSparseVec(int _len):len(_len),
                                 fidx(new FIDX[len]),
                                 fval(new FVAL[len]) 
{
}


CSparseVec::CSparseVec(int _len, FIDX *_fidx, FVAL *_fval):len(_len),
                                                           fidx(0),
                                                           fval(0) 
{
  fidx = new FIDX[len];
  fval = new FVAL[len];
  memcpy(fidx, _fidx, sizeof(FIDX)*len);
  memcpy(fval, _fval, sizeof(FVAL)*len);
}


CSparseVec::~CSparseVec() 
{
  len = 0;
  if(fidx) {delete [] fidx; fidx = 0;}
  if(fval) {delete [] fval; fval = 0;}
}




/** Assignment. CSparseVec to CSparseVec
 */
void CSparseVec::Assign(const CSparseVec *rhs)
{
  if(rhs == this) return;

  if(len > 0) 
  {
    assert(fidx && fval);
    delete [] fidx;
    delete [] fval;   
  }
  
  len = rhs->len;
  fidx = new FIDX[len];
  fval = new FVAL[len];   
  memcpy(fidx, rhs->fidx, sizeof(FIDX)*len);
  memcpy(fval, rhs->fval, sizeof(FVAL)*len);

  return;
}

/** Assignment. Array to CSparseVec
 */
void CSparseVec::Assign(const double *array, int dim)
{
  assert(dim > 0);
  register int i, j;

  for(i=0; i<dim; i++) 
    if(fabs(array[i]) > 1e-30) len++;

  fidx = new FIDX[len];
  fval = new FVAL[len];

  for(i=0, j=0; i<dim; i++)
    if(fabs(array[i]) > 1e-30)
    {
      fidx[j] = i+1;
      fval[j] = array[i];
      j++;
    }
}


/** Count number of unique fidx.
 */
inline int CSparseVec::UniqueFidx(FIDX *a, int alen, FIDX *b, int blen)
{
  register int a_idx, b_idx;
  int uniqueLen = 0;

  a_idx = b_idx = 0;
  while(a_idx < alen && b_idx < blen)
  {
    if(a[a_idx] == b[b_idx]) {a_idx++; b_idx++;}
    else if(a[a_idx] < b[b_idx]) a_idx++;
    else if(a[a_idx] > b[b_idx]) b_idx++;
    uniqueLen++;
  }

  return uniqueLen + (alen-a_idx) + (blen-b_idx);
}


/** Addition.
 */
void CSparseVec::Add(const CSparseVec *rhs)
{
  if(len == 0)
  {
    this->Assign(rhs);
    return;
  }

  if(rhs == this)
  {
    for(int i=0; i<len; i++) 
      fval[i] += fval[i];
    return;
  }

  assert(len > 0);
  assert(fidx);
  assert(fval);

  int n = UniqueFidx(fidx, len, rhs->fidx, rhs->len);

  FIDX *newfidx = new FIDX[n];
  FVAL *newfval = new FVAL[n];

  assert(newfidx && newfval);

  int i=0, j=0, k=0;
  while(i < rhs->len && j < len) 
  {
    if(rhs->fidx[i] == fidx[j])
    {
      newfidx[k] = fidx[j];
      newfval[k] = fval[j] + rhs->fval[i];
      i++; j++; k++;
      
    }
    else if(rhs->fidx[i] < fidx[j]) 
    {
      newfidx[k] = rhs->fidx[i];
      newfval[k] = rhs->fval[i];
      i++; k++;
    }
    else
    {
      newfidx[k] = fidx[j];
      newfval[k] = fval[j];
      j++; k++;
    }
  }

  while(i < rhs->len) 
  {
    newfidx[k] = rhs->fidx[i];
    newfval[k] = rhs->fval[i];
    i++; k++;
  }

  while(j < len) 
  {
    newfidx[k] = fidx[j];
    newfval[k] = fval[j];
    j++; k++;
  }
  
  delete [] fidx;
  delete [] fval;
  
  fidx = newfidx;
  fval = newfval;
  len  = n;

  return;
}


/** Addition with weighted CSparseVec.
 */
void CSparseVec::Add(const CSparseVec *rhs, double weight)
{
  if(len == 0)
  {
    this->Assign(rhs);
    this->Mult(weight);
    return;
  }

  if(rhs == this)
  {
    for(int i=0; i<len; i++) 
      fval[i] += weight*fval[i];
    return;
  }
  
  if(ABS(weight) < ZERO_EPS) return;

  int n = UniqueFidx(fidx, len, rhs->fidx, rhs->len);
  FIDX *newfidx = new FIDX[n];
  FVAL *newfval = new FVAL[n];

  int i=0, j=0, k=0;
  while(i < rhs->len && j < len) {
    if(rhs->fidx[i] == fidx[j])
    {
      newfidx[k] = fidx[j];
      newfval[k] = fval[j] + weight*rhs->fval[i];
      i++; j++; k++;
      
    }
    else if(rhs->fidx[i] < fidx[j]) 
    {
      newfidx[k] = rhs->fidx[i];
      newfval[k] = weight*rhs->fval[i];
      i++; k++;
    }
    else
    {
      newfidx[k] = fidx[j];
      newfval[k] = fval[j];
      j++; k++;
    }
  }

  while(i < rhs->len) {
    newfidx[k] = rhs->fidx[i];
    newfval[k] = weight*rhs->fval[i];
    i++; k++;
  }

  while(j < len) {
    newfidx[k] = fidx[j];
    newfval[k] = fval[j];
    j++; k++;
  }

  delete [] fidx;
  delete [] fval;
  
  fidx = newfidx;
  fval = newfval;
  len  = n;

  return;
}


/** Multiply with scalar.
 */
void CSparseVec::Mult(const double rhs)
{
  for(int i=0; i<len; i++)
    fval[i] *= rhs;

  return;
}


/** Compute dot-product with another CSparseVec.
 */
double CSparseVec::Dot(const CSparseVec *rhs)
{
  double res = 0.0;
  register int i=0, j=0;
  FIDX *rhs_fidx = rhs->fidx;
  FVAL *rhs_fval = rhs->fval;
  int rhs_len = rhs->len;

  while(i < rhs_len && j < len) 
  {
    if (rhs_fidx[i] == fidx[j]) 
    {
      res += rhs_fval[i]*fval[j];
      i++; j++;
    }
    else if(rhs_fidx[i] > fidx[j]) j++;
    else if(rhs_fidx[i] < fidx[j]) i++;
  }

  return res;
}


/** Compute norm square.
 */
double CSparseVec::NormSq()
{
  double norm = 0.0;
  
  for(int i=0; i<len; i++)
    norm += fval[i]*fval[i];

  return norm;
}


/** Compute L2 norm.
 */
double CSparseVec::L2Norm()
{
  double norm = 0.0;
  
  for(int i=0; i<len; i++)
    norm += fval[i]*fval[i];
  
  return sqrt(norm);
}


/** Compute L_infty norm of vector (*this - rhs).
 */
double CSparseVec::LInftyNorm(const CSparseVec *rhs) 
{
  register int i=0, j=0;
  double norm = -1;
  int rhs_len = rhs->len;

  FIDX *rhs_fidx = rhs->fidx;
  FVAL *rhs_fval = rhs->fval;


  while(i < rhs_len && j < len) 
  {
    if (rhs_fidx[i] == fidx[j]) 
    {
      norm = MAX(norm, fabs(rhs_fval[i] - fval[j]));
      i++; j++;
    }
    else if(rhs_fidx[i] > fidx[j]) j++;
    else if(rhs_fidx[i] < fidx[j]) i++;
  }

  return norm;
}


/** Dump CSparseVec into a file.
 */
void CSparseVec::ToFile(string filename)
{
  ofstream ofp(filename.c_str());
  
  if(!ofp.good())
  {
    cout << "ERROR: can not open file <" << filename << "> !" << endl;
    return;
  }

  ofp.setf(ios::fixed);

  for(int i=0; i<len; i++)
    ofp << fidx[i] << ":" << fval[i] << endl;
  
  ofp.close();
}


/**  Set *this to a zero vector.
 */
void CSparseVec::Zero()
{
  // chteo: when len==0, fidx and fval should be null too.
  if(len==0 || this==0) return;

  len = 0;
  if(fidx) {delete [] fidx; fidx = 0;}
  if(fval) {delete [] fval; fval = 0;}
}


/** Print vector to stdout
 */
void CSparseVec::Print()
{
  for(int i=0; i<len; i++)
  {
    printf("%3ld:%.4f\n",fidx[i], fval[i]);
  }
  printf("\n");
}

/** Convert CSparseVec to array.
 *
 *  \param array [write] Dense vector
 *  \param dim [read] Dimension of dense vector
 */
void CSparseVec::ToArray(double *array, int dim)
{
  if(dim < len)
    cout << "dim : " << dim << "   len :" << len << endl;

  assert(dim >= len);
  register int i;

  memset(array, 0, sizeof(double)*dim);
  
  for(i=0; i<len; i++)
    array[fidx[i]-1] = fval[i];    
}

#endif
