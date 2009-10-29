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
 * Authors: Simon Guenter (simon.guenter@nicta.com.au)
 *
 * Created: (03/01/2008) 
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "onlinevecdata.hpp"
//#include "vecdata.hpp"
//#include "info.hpp"


/**  Constructor
 */
COnlineVecData::COnlineVecData(unsigned int batchsize, unsigned int usenum)
    : CVecData(),
      position(0),
      called(0)
{
        this->batchsize = batchsize;
        this->usenum = usenum;
    
        if (!matrixRowView) 
        {
                // always use matrixRowView
                CreateFeatureMatrixRowViews(); 
                CreateLabelMatrixRowViews();  
        }
        
        currentlabels = new TheMatrix(batchsize,1);
        for(unsigned int i=0;i<batchsize;i++) 
        {
                Scalar label;
                Y->Get(i,label);
                currentlabels->Set(i,label);
        }
}


// destructor
COnlineVecData::~COnlineVecData()
{
        delete currentlabels;
}


/** changes mini batches if necessary */
inline void COnlineVecData::newcalc() 
{
        called++;
   
        if(usenum<called) 
        {
                batchnum++;
                position = position + batchsize;
                if(position > numOfAllExample) 
                {
                        position -= numOfAllExample;
                        iternum++;
                }
                called = 1;
                for(unsigned int i=0;i<batchsize;i++) 
                {
                        // start at the beginning if mini wrap would
                        // exceed size of data set
                        unsigned int pos = (position + i) % numOfAllExample;        
                        Scalar label;
                        Y->Get(pos,label);
                        currentlabels->Set(i,label);
                }
        }    
}


void COnlineVecData::ComputeF(const TheMatrix& w, TheMatrix& f)
{
        newcalc();
        f.Zero();
        for(unsigned int i=0;i<batchsize;i++) 
        {
                unsigned int pos = (position + i) % numOfAllExample;
                Scalar value;        
                x[pos]->Dot(w, value);  
                f.Set(i,value);
        } 
}



void COnlineVecData::Grad(const TheMatrix& l, TheMatrix& grad)
{
        grad.Zero();
        for(unsigned int i=0;i<batchsize;i++) 
        {
                unsigned int pos = (position + i) % numOfAllExample;
                Scalar scale;
                l.Get(i,scale);
                grad.ScaleAdd(scale,*x[pos]);
        }
}
