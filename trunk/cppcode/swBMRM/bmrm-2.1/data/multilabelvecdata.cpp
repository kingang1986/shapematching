/* Copyright (c) 2006, National ICT Austtralia 
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

#ifndef _SEQDATA_CPP_
#define _SEQDATA_CPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "multilabelvecdata.hpp"


   
/**  Constructor
 */
CMultilabelVecData::CMultilabelVecData()
   : CData(), 
     CVecFeature(), 
     CVarLenVecLabel()
{
        // sanity check
        if(numOfAllExample != numOfAllLabel)
        {
                std::ostringstream msg;
                msg << "Number of examples (" << numOfAllExample << ") does not match number of labels (" << numOfAllLabel << ")";
                throw CBMRMException(msg.str(),"CMultilabelVecData::CMultilabelVecData()");
        }
        
        // dataset statistics
        if(verbosity)
        {
                std::cout << "Dataset properties:"  << std::endl;
                std::cout << "1.  Feature file              : " << featureFile << std::endl;
                std::cout << "2.  Label file                : " << labelFile << std::endl;
                std::cout << "3.  Number of subsets         : " << numOfAllSubset << std::endl;
                std::cout << "4.  Number of examples        : " << numOfAllExample << std::endl;    
                std::cout << "5.  Number of nonzero features: " << nnz; 
                if(biasFlag)
                        std::cout << "(+ no. of examples for shifted hyperplane)";
                std::cout << std::endl;
                std::cout << "7.  Average label dimension   : " << avgLabelDimension << std::endl;
                std::cout << "8.  Dim. of feature vector    : " << featureDimension;
                if(biasFlag)
                        std::cout << "(+ 1 for shifted hyperplane)";
                std::cout << std::endl; 
                std::cout << "9.  Dataset density           : " <<  density*100.0 << std::endl;
                std::cout << "10. Avg. nnz(feature vec.)    : " << ((double)nnz)/numOfAllExample << std::endl;
                std::cout << std::endl;
        }
}

#endif
