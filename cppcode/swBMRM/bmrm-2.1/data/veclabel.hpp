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

#ifndef _VECLABEL_HPP_
#define _VECLABEL_HPP_

#include <vector>
#include <iostream>

#include "common.hpp"
#include "sml.hpp"


/** Container for vector labels (each label must be same length!)
 *    
 *   Label file format: rows of <label>
 *   where
 *   <label>   .=. <value_1> [value_2> ... <value_k>]
 *   <value_i> .=. i-th label (scalar value) of a particular example
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 *   2.  For centralized dataset and serial computation mode only!
 */
class CVecLabel
{        
protected:
        /** Verbosity level
         */
        int veclabel_verbosity;
        
        /** Label matrix (each row is a label)
         */
        TheMatrix* Y;
            
        /** A list of labels as (explicit) row vectors
         */
        TheMatrix** y;
                    
        /** Types of label
         */
        enum LABEL_TYPE {SINGLE_LABEL, VECTOR_LABEL};
        
        /** dimensionality of the label
         */
        unsigned int labelDimension;

        
        /** Smallest label (needed in determining the number of classes in multi-class classification)
         */
        Scalar minLabel;

        /** Largest label (needed in determining the number of classes in multi-class classification)
         */
        Scalar maxLabel;
      
        /** Number of labels in the sub-dataset
         */
        unsigned int numOfLabel;
        
        /** Number of labels in the WHOLE dataset
         *  This is the same as numOfLabel in centralised dataset, serial computation mode
         */
        unsigned int numOfAllLabel;

        /** The first example of the dataset (skip #startExample# examples before reading the rest)
         */
        unsigned int startExample;
      
        /** Flag for type of label e.g., single or multiple
         */
        int labelType;
      
        /** Name of the file containing label vectors
         */
        std::string labelFile;

        /** A string template for floating point value
         */
        std::string scalar_value_format;
        
          /** To create matrix row view for label matrix
         */
        bool labelMatrixRowView;
        
        /** Scan label file to determine some label set properties such as type of label
         */
        virtual void ScanLabelFile();
        
        /** Allocate label vector (or matrix) and loat label from label file
         */
        virtual void LoadLabels(); 
        
public:
        CVecLabel();
        virtual ~CVecLabel();

        /** Return the labels. 
         */
        virtual const TheMatrix& labels(void){ return *Y; }
        
        /** 
         * Return the maximum of labels
         */
        Scalar const MaxLabel() {return maxLabel;}

        /** 
         * Return the minimum of labels
         */
        Scalar const MinLabel() {return minLabel;}
  
        /** Create matrix row view for every label vector 
         */
        virtual void CreateLabelMatrixRowViews();  

};

#endif
