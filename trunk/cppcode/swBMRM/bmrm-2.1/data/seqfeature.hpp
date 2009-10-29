/* Copyright (c) 2006, National ICT Australia 
 * All rights reserved. 
 * 
 * The contents of this file are subject to the Mozilla Public License 
 * Version 1.1 (the "License"); you may not use this file except in 
 * compliance with the License. You may obtain a copy of the License at 
 * http://www.mozilla.org/MPL/ 
 * 
 * Software distributed under the License is distributed on an "AS IS" 
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the 
 * License for the specific language governing rights and limitations 
 * under the License. 
 * 
 * Authors: Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *
 * Created: (29/01/2008) 
 *
 * Last Updated:
 */

#ifndef _SEQFEATURE_HPP_
#define _SEQFEATURE_HPP_

#include <vector>
#include <iostream>

#include "common.hpp"
#include "sml.hpp"


/** Container for sequence feature vectors
 *  Mainly for Automatic Paragraph Segmentation 
 *  
 *   Example format: 
 *  
 *   maxDuration:<P>
 *   minDuration:<N>
 *   globalFeatureDim:<P>
 *   sequence:<N> 
 *   phi:1 
 *   pos:<N> <sparse_vector>
 *   ...
 *   phi:2 
 *   pos:<N>,<N> <sparse_vector>
 *   ...
 *   sequence:<N>
 *   phi:1 
 *   pos:<N> <sparse_vector>
 *   ...
 *   phi:2 
 *   pos:<N>,<N> <sparse_vector>
 *   ...
 *
 *   where
 *   <P>              .=. positive integer
 *   <N>              .=. natural number
 *   <sparse_vector>  .=. <index_1>:<value_1> [<index_2>:<value_2> ... <index_k>:<value_k>]
 *   <index_j>        .=. index (positive integer) of j-th non-zero element of a particular sparse feature vector
 *                        Note that index_i > index_j for j>i
 *   <value_j>        .=. value (scalar) of j-th non-zero element of a particular feature vector
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 *   2.  For centralized dataset and serial computation mode only!
 */
class CSeqFeature
{
protected:               
        /** verbosity level
         */
        int seqfeature_verbosity;

        /** Feature vectors
         */
        vector<TheMatrix> phi_1;
        vector<vector<TheMatrix> > phi_2;
       
        /** Maximum duration of a segment in a sequence
         */
        unsigned int maxDuration;
        
        /** Minimum duration of a segment in a sequence
         */
        unsigned int minDuration;
        
        /** number of examples in this sub-dataset
         */
        unsigned int numOfSeq;

        /** number of examples in the whole dataset
         */
        unsigned int numOfAllSeq;

        /** dimensionality of the feature vector
         */
        unsigned int featureDimension;        

        /** Name of the file containing feature vectors
         */
        std::string featureFile;
      
        /** A string template for sparse vector element
         */
        std::string svec_feature_index_and_value_format;

        /** A string template for dense vector element
         */
        std::string scalar_value_format;
      
        /** number of nonzero features
         */
        unsigned int nnz;

        /** numOfNonzero / total number of entries in feature matrix
         */
        double density;
       
        
        /** Allocate data matrix and load fetures from data file
         */
        virtual void LoadFeatures();     
        
public:
        CSeqFeature();
        virtual ~CSeqFeature(){}
      
};

#endif
