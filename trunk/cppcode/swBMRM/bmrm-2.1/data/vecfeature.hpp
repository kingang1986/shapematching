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
 * Created: (26/01/2008) 
 *
 * Last Updated:
 */

#ifndef _VECFEATURE_HPP_
#define _VECFEATURE_HPP_

#include <vector>
#include <iostream>

#include "common.hpp"
#include "sml.hpp"


/** Container for feature vectors
 *    
 *   Feature vector file format: rows of <feature_vector>
 *   where
 *   <feature_vector> .=. [qid:<id>] <dense_vector> | <sparse_vector>
 *   <dense_vector>   .=. <value_1> [<value_2> ... <value_k>]
 *   <sparse_vector>  .=. <index_1>:<value_1> [<index_2>:<value_2> ... <index_k>:<value_k>]
 *   <id>             .=. a number (examples with same qid must be adjacent to each other in the data file)
 *   <index_j>        .=. index (positive integer) of j-th non-zero element of a particular sparse feature vector
 *                        Note that index_i > index_j for j>i
 *   <value_j>        .=. value (scalar) of j-th non-zero element of a particular feature vector
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 *   2.  For centralized dataset and serial computation mode only!
 */
class CVecFeature
{
protected:               
        /** verbosity level
         */
        int vecfeature_verbosity;

        /** Feature matrix (each row is an example)
         */
        TheMatrix* X;
      
        /** A list of examples as (explicit) row vectors
         */
        TheMatrix** x;
      
        /** Types of feature vector
         */
        enum FEATURE_TYPE {SPARSE_FEATURE, DENSE_FEATURE};
        
        /** whether feature vector comes with artificial bias feature
         */
        bool biasFlag;
        
        /** number of examples in this sub-dataset
         */
        unsigned int numOfExample;

        /** number of examples in the whole dataset
         */
        unsigned int numOfAllExample;

        /** dimensionality of the feature vector
         */
        unsigned int featureDimension;        
 
        /** Number of subsets in this (sub-)dataset
         */
        unsigned int numOfSubset; 
      
        /** Number of subsets in the dataset
         */
        unsigned int numOfAllSubset;

        /** Number of examples in each subset
         */
        std::vector<int>* subsetSizes; 
              
        /** To create matrix row view for feature matrix
         */
        bool featureMatrixRowView;
            
        /** The first example of the dataset (skip #startExample# examples before reading the rest)
         */
        unsigned int startExample;

        /** The first subset of the dataset (skip #startSubsets# subsets before reading the rest)
         */
        unsigned int startSubset;
      
        /** Flag for type of feature vector e.g., dense or sparse
         */
        unsigned int featureType;

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

        /** Structure to keep some info. of subsets of examples in the dataset
         */
        struct Subset {
                int ID;
                int startIndex;
                int size;
        };
        
        /** Scan the feature file and determine some feature set properties such as dimension of example        
         */
        virtual void ScanFeatureFile();
        
        /** Allocate data matrix and load fetures from data file
         */
        virtual void LoadFeatures();     
        
public:
        CVecFeature();
        virtual ~CVecFeature();

        /** Given a weight vector w return the prediction f = Xw. It is the
         *  duty of the caller to ensure that w and f have conforming
         *  dimensions.  
         * 
         *  @param w [read]  the weight vector
         *  @param f [write] the prediction 
         */
        virtual void ComputeF(const TheMatrix& w, TheMatrix& f);
      
       
        /** Given a weight vector w return the prediction f_i = x_i*w. 
         *  It is the duty of the caller to ensure that w and f have conforming
         *  dimensions.  
         * 
         *  @param w [read]  the weight vector (i.e. a matrix)
         *  @param f [write] the prediction 
         */
        virtual void ComputeFi(const TheMatrix& w, TheMatrix& f, const unsigned int i);
      
        /** Given a weight vector w return the prediction f = x_i*w. 
         *  It is the duty of the caller to ensure that w has conforming
         *  dimensions.  
         * 
         *  @param w [read]  the weight vector
         */
        virtual void ComputeFi(const TheMatrix& w, Scalar & f, const unsigned int i);
      
      
        /** f = X*w, where X is a data matrix (with rows of feature vectors) 
         */
        virtual void XMultW(const TheMatrix& w, TheMatrix& f)
        {
                X->Dot(w, f);
        }


        /** f = X^T*w, where X is a data matrix (with rows of feature vectors) 
         */
        virtual void XTMultW(const TheMatrix& w, TheMatrix& f)
        {
                X->TransposeDot(w, f);
        }

        
        /** f = w*X, where X is a data matrix (with rows of feature vectors) 
         */
        virtual void WMultX(const TheMatrix& w, TheMatrix& f)
        {
                w.Dot(*X, f);
        }
        
        
        /** f = w*X, where X is a data matrix (with rows of feature vectors) 
         */
        virtual void WTMultX(const TheMatrix& w, TheMatrix& f)
        {
                w.TransposeDot(*X, f);
        }
        
        /** Adds scale times x_i to w. 
         *  It is the duty of the caller to ensure that w has conforming
         *  dimensions.
         *
         *  @param w [write] weight vector
         *  @param i [read] position of the x_i
         *  @param scale [read] scaling factor    
         */
        virtual void AddElement(TheMatrix& w, const unsigned int &i,Scalar scale);


        /** Given the loss vector l return the gradient grad = lX. It is the
         *  duty of the caller to ensure that l and grad have conforming
         *  dimensions.
         * 
         *  @param l    [read]  the loss vector
         *  @param grad [write] the gradient
         */
        virtual void Grad(const TheMatrix& l, TheMatrix& grad);
 

        /** Create explicit row views of feature matrix
         */
        virtual void CreateFeatureMatrixRowViews();  // create matrix row view for every example vector

        /** return the number of subsets of this sub-dataset
         */
        int NumOfSubset() {return numOfSubset;}
        
        /** return the number of subsets of the WHOLE dataset when used in distributed environment
         *  for centralized dataset, this is equivalent to NumOfSubset()
         */
        int NumOfAllSubset() {return numOfAllSubset;}
        
               
        /** subset information
         */
        Subset *subset;              
};

#endif
