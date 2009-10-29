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
 * Created      : 05/01/2007
 * Last Updated :
 */

#ifndef _DENSEVEC_CPP_
#define _DENSEVEC_CPP_

#include "densevec.hpp"
#include "sparsevec.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iomanip>

using namespace std;

// CDenseVec::CDenseVec():len(0),array(0)
// {
// }


CDenseVec::CDenseVec(int _len):len(_len),
                               array(0)
{
  // chteo: malloc + memset is slower than calloc
  //array = (double*)malloc(len*sizeof(double));
  //memset(array, 0, len*sizeof(double));
  array = (double*)calloc(len, sizeof(double));
}

CDenseVec::CDenseVec(int _len, double *rhs):len(_len),
                                          array(rhs)
{
  array = (double*)malloc(len*sizeof(double));
  memcpy(array, rhs, len*sizeof(double));
}

CDenseVec::CDenseVec(int _len, CSparseVec *svec):len(_len),
                                                 array(0)
{
  assert(svec->fidx[svec->len-1] <= len);

  // chteo: malloc + memset is slower than calloc
  //array = (double*)malloc(len*sizeof(double));
  //memset(array, 0, len*sizeof(double));
  array = (double*)calloc(len, sizeof(double));

  for(register int i=0; i<svec->len; i++)
    array[svec->fidx[i]-1] = svec->fval[i];
}



CDenseVec::~CDenseVec() 
{
  if(array) {free(array); array = 0;}
}



/** Assignment. CDenseVec to CDenseVec
 *
 *  \param rhs [read] Dense vector to be assigned to *this
 */
void CDenseVec::Assign(const CDenseVec *rhs)
{
  assert(len == rhs->len);

  if(rhs == this) return;
  
  memcpy(array, rhs->array, len*sizeof(double));
}


/** Assignment. CSparseVec to CDenseVec
 *  
 *  \param svec [read] Sparse vector to be made dense
 */
void CDenseVec::Assign(const CSparseVec *svec)
{
  // zero sparse vector
  if(svec->len == 0)
  {
    memset(array, 0, len*sizeof(double));
    return;
  }
   
  assert(svec->fidx[svec->len-1] <= len);

  for(register int i=0; i<svec->len; i++)
    array[svec->fidx[i]-1] = svec->fval[i];
}


/** Assignment. Array to CDenseVec
 *
 *  \param rhs [read] Array to be assigned to this->array
 *  \param dim [read] Dimension of rhs
 */
void CDenseVec::Assign(const double *rhs, int dim)
{
  assert(dim ==  len);
 
  if(rhs == array) return;

  memcpy(array, rhs, len*sizeof(double));
}


/** Addition. Add CSparseVec to CDenseVec.
 *
 *  \param svec [read] Sparse vector to be added to *this
 */
void CDenseVec::Add(const CSparseVec *svec)
{
  if(svec->len == 0) return;

  assert(len > 0);
  assert(svec->len > 0);
  assert(svec->fidx[svec->len-1] <= len);

  for(register int i=0; i<svec->len; i++)
    array[svec->fidx[i]-1] += svec->fval[i];
}


/** Addition with weighted CSparseVec.
 *
 *  \param svec [read] Sparse vector to be added to *this
 *  \param weight [read] Coefficient of rhs
 */
void CDenseVec::Add(const CSparseVec *svec, double weight)
{
  if(svec->len == 0) return;
  if(ABS(weight) < ZERO_EPS) return;
  /*if(weight == 1.0) {
    Add(svec);
    return;
    }*/

  assert(len > 0);
  assert(svec->len > 0);
  assert(svec->fidx[svec->len-1] <= len);

  for(register int i=0; i<svec->len; i++)
    array[svec->fidx[i]-1] += weight * svec->fval[i];
}


/** Addition with CDenseVec.
 *
 *  \param rhs [read] Dense vector to be added to *this
 */
void CDenseVec::Add(const CDenseVec *rhs)
{
  assert(len == rhs->len);

  for(register int i=0; i<len; i++)
    array[i] += rhs->array[i];
}



/** Addition with weighted CDenseVec.
 *
 *  \param rhs [read] Dense vector to be added to *this
 *  \param weight [read] Coefficient of rhs
 */
void CDenseVec::Add(const CDenseVec *rhs, double weight)
{
  assert(len == rhs->len);
  if(ABS(weight) < ZERO_EPS) return;

  for(register int i=0; i<len; i++)
    array[i] += weight*rhs->array[i];
}


/** Addition with weighted array.
 *
 *  \param n      [read] Dimension of rhs
 *  \param rhs    [read] Array to be added to *this
 *  \param weight [read] Coefficient of rhs
 */
void CDenseVec::Add(int n, double *rhs, double weight)
{
  assert(len == n);
  if(ABS(weight) < ZERO_EPS) return;

  for(register int i=0; i<len; i++)
    array[i] += weight*rhs[i];
}


/** Multiply with scalar.
 *
 *  \param weight [read] Multiplication operand
 */
void CDenseVec::Mult(const double weight)
{
  for(register int i=0; i<len; i++)
    array[i] *= weight;
}


/**  Set *this to a zero vector.
 */
void CDenseVec::Zero()
{
  assert(len > 0);
  memset(array, 0, len*sizeof(double));
}


/** Compute dot-product with an array
 *
 *  \param n   [read] dimension of rhs
 *  \param rhs [read] array to be dotted with *this
 */
double CDenseVec::Dot(int n, double *rhs)
{
  assert(n == len);

  double res = 0.0;
  
  for(register int i=0; i<len; i++)
    res += array[i]*rhs[i]; 
  
  return res;
}


/** Compute dot-product with a CDenseVec.
 */
double CDenseVec::Dot(const CDenseVec *rhs)
{
  assert(rhs->len == len);

  double res = 0.0;
  
  for(register int i=0; i<len; i++)
    res += array[i]*rhs->array[i]; 
  
  return res;
}


/** Compute dot-product with a CSparseVec.
 */
double CDenseVec::Dot(const CSparseVec *svec)
{
  double res = 0.0;
  
  if(svec->len == 0) return res;
  
  assert(svec->len > 0);

  for(register int i=0; i<svec->len; i++)
    res += array[svec->fidx[i]-1] * svec->fval[i];

  return res;
}


/** Compute norm square.
 */
double CDenseVec::NormSq()
{
  double norm = 0.0;
  
  for(int i=0; i<len; i++)
    norm += array[i]*array[i];

  return norm;
}


/** Compute the sum of all entries
 */
double CDenseVec::Sum()
{
  double sum=0;
  
  for(int i=0; i<len; i++)
    sum += array[i];

  return sum;
}


/** Compute L1 norm.
 */
double CDenseVec::L1Norm()
{
  double norm = 0;

  for(int i=0; i<len; i++)
    norm += ABS(array[i]);

  return norm;
}


/** Compute L2 norm.
 */
double CDenseVec::L2Norm()
{
  double norm = 0.0;
  
  for(int i=0; i<len; i++)
    norm += array[i]*array[i];
  
  return sqrt(norm);
}



/** Compute L_infty norm of vector (*this - rhs).
 * 
 *  \param svec [read] Sparse vector
 */
double CDenseVec::LInftyNorm(const CSparseVec *svec) 
{
  double norm = -1;

  if(svec->len == 0)
  {
    for(register int i=0; i<len; i++)
      norm = MAX(norm, ABS(array[i]));
  }
  else
  {
    assert(svec->len > 0);
    for(register int i=0; i<svec->len; i++)
      norm = MAX(norm, ABS(svec->fval[i] - array[svec->fidx[i]-1]));
  }

  return norm;
}


/** Compute L_infty norm of vector (*this - rhs).
 *
 *  \param rhs [read] Dense vector
 */
double CDenseVec::LInftyNorm(const CDenseVec *rhs) 
{
  double norm = -1;

  for(register int i=0; i<len; i++)
    norm = MAX(norm, ABS(rhs->array[i] - array[i]));

  return norm;
}


/** Dump CDenseVec into a file.
 */
void CDenseVec::ToFile(string filename)
{
  ofstream ofp(filename.c_str());
  
  if(!ofp.good())
  {
    printf("ERROR: can not open file <%s> !\n", filename.c_str());
    return;
  }

  ofp.setf(ios::fixed);

  for(int i=0; i<len; i++)
    ofp << setprecision(12) << array[i] << endl;
  
  ofp.close();
}


/** Read dense vector from file
 */  
void CDenseVec::FromFile(string filename)
{
  ifstream ifp(filename.c_str());
  
  if(!ifp.good())
  {
    printf("ERROR: cannot read dense vector from file <%s> !\n",filename.c_str());
    exit(0);
  }

  for(int i=0; i<len && !ifp.eof(); i++)
    ifp >> array[i];
}

/** Print vector to stdout
 */
void CDenseVec::Print()
{
  for(int i=0; i<len; i++)
    printf("%3d:%.4f\n",i, array[i]);
  printf("\n");
}


/** Convert CDenseVec to CSparseVec.
 */
CSparseVec* CDenseVec::ToCSparseVec()
{
  int  n = 0;
  FIDX *fidx = 0;
  FVAL *fval = 0;

  for(int i=0; i<len; i++)
    if(ABS(array[i]) > 1e-15)
      n++;

  fidx = new FIDX[n];
  fval = new FVAL[n];
  
  int j=0;
  for(int i=0; i<len; i++)
    if(ABS(array[i]) > 1e-15)
    {
      fidx[j] = i+1;
      fval[j] = array[i];
      j++;
    }
  
  CSparseVec *svec = new CSparseVec(n, fidx, fval);

  return svec;
}


/** Randomly initialize the entries
 */
void CDenseVec::RandInit()
{
  for(int i=0; i<len; i++)
    array[i] = (double)rand()/RAND_MAX;
}

#endif


/*   NOTE: 
 *   1. [050107:2256] CSparseVec index starts from 1 !
 *
 *
 */
