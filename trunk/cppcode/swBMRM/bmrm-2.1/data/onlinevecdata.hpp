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
 * Created: (02/11/2007) 
 *
 * Last Updated: (06/11/2007)   
 */

#ifndef _ONLINEVECDATA_HPP_
#define _ONLINEVECDATA_HPP_

#include <vector>
#include <iostream>

#include "common.hpp"
#include "sml.hpp"
#include "data.hpp"
#include "vecdata.hpp"


/**
 * Container for dataset of type vector label and vector feature
 * and used in a online  setting
 * 
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 *   2.  For centralized dataset and serial computation mode only!
 */
class COnlineVecData: public CVecData
{
protected:
        unsigned int batchsize;  // size of minibatch
        unsigned int usenum; // number of times minibatch will be accessed
        unsigned int position; // start position in the whole data of the current minibatch
        unsigned int called; // number of times current minibatch was accessed
        unsigned int batchnum; // number of minibatches
        unsigned int iternum; // number of pass through data set
        TheMatrix * currentlabels; // labels of current minibatch
        
        virtual void newcalc();

public:
        COnlineVecData(unsigned int batchsize, unsigned int usenum);
        virtual ~COnlineVecData();
    

        /** 
         * Given a weight vector w return the prediction f = Xw. It is the
         * duty of the caller to ensure that w and f have conforming
         * dimensions.  
         * 
         * @param w [read]  the weight vector
         * @param f [write] the prediction 
         */
        virtual void ComputeF(const TheMatrix& w, TheMatrix& f);


        /** 
         * Given the loss vector l return the gradient grad = lX. It is the
         * duty of the caller to ensure that l and grad have conforming
         * dimensions.
         * 
         * @param l    [read]  the loss vector
         * @param grad [write] the gradient
         */
        virtual void Grad(const TheMatrix& l, TheMatrix& grad);

    
        /** f = X*w, where X is a data matrix (with rows of feature vectors) 
         */
        virtual void XMultW(const TheMatrix& w, TheMatrix& f)
        {
                ComputeF(w,f);
        }
        

        /** f = X^T*w, where X is a data matrix (with rows of feature vectors) 
         */
        virtual void XTMultW(const TheMatrix& w, TheMatrix& f)
        {        
                ComputeF(w,f);
        }
    
    
        /** f = w*X, where X is a data matrix (with rows of feature vectors) 
         */
        virtual void WMultX(const TheMatrix& w, TheMatrix& f)
        {
                Grad(w,f);
        }
        
    
        /** f = w*X, where X is a data matrix (with rows of feature vectors) 
         */
        virtual void WTMultX(const TheMatrix& w, TheMatrix& f)
        {
                Grad(w,f);
        }
    
        
        /** returns labels of current minibatch */
        virtual const TheMatrix & labels() 
        {
                return *currentlabels;
        }
      
        
        /** 
         * Return the number of examples in MINIBATCH
         */    
        virtual unsigned int slice_size(void) const { return batchsize; }
            
};

#endif
