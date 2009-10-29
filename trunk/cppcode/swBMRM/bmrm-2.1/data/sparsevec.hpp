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
 * Last Updated : 29/01/2008
 */

#ifndef _SPARSEVEC_HPP_
#define _SPARSEVEC_HPP_


/** Sparse vector class.
 */
class CSparseVec 
{
protected:
        int UniqueFidx(FIDX *a, int alen, FIDX *b, int blen);
        
public:

  int  len;
  FIDX *fidx;
  FVAL *fval;

  /// Constructors
  CSparseVec();
  CSparseVec(int _len);
  CSparseVec(int _len, FIDX *_fidx, FVAL *_fval);

  /// Destructor
  ~CSparseVec();

  /// Methods
  
  void Assign(const CSparseVec *rhs);
  void Assign(const double *array, int dim);
  void Add(const CSparseVec *rhs);
  void Add(const CSparseVec *rhs, double weight);
  void Mult(const double rhs);
  void Zero();

  double Dot(const CSparseVec *rhs);
  double Norm();
  double LInftyNorm(const CSparseVec *rhs); 

  void ToFile(string filename);
  void ToArray(double *array, int dim);
  void Print();

};

#endif
