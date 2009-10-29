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

#ifndef _DENSEVEC_HPP_
#define _DENSEVEC_HPP_

#include "errorcode.hpp"
#include "common.hpp"
#include "dataset.hpp"

class CSparseVec;

/** Dense vector class.
 */
class CDenseVec {
public:

  int  len;
  double *array;

  /// Constructors
  CDenseVec(int _len);
  CDenseVec(int _len, double *rhs);
  CDenseVec(int _len, CSparseVec *svec);

  /// Destructor
  ~CDenseVec();

  /// Methods
  void Assign(const CDenseVec *rhs);
  void Assign(const CSparseVec *svec);
  void Assign(const double *rhs, int dim);
  void Add(const CDenseVec *rhs);
  void Add(const CSparseVec *svec);
  void Add(const CDenseVec *rhs, double weight);
  void Add(const CSparseVec *svec, double weight);
  void Add(int n, double *rhs, double weight);
  void Mult(const double weight);
  void Zero();
  void RandInit();

  double Dot(const CDenseVec *rhs);
  double Dot(const CSparseVec *rhs);
  double Dot(int n, double *rhs);
  double NormSq();
  double Sum();
  double L1Norm();
  double L2Norm();
  double LInftyNorm(const CDenseVec *rhs); 
  double LInftyNorm(const CSparseVec *rhs); 

  void ToFile(string filename);
  void FromFile(string filename);
  CSparseVec* ToCSparseVec();
  void Print();

};

#endif
