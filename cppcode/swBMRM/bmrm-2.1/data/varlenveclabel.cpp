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

#ifndef _VARLENVECLABEL_CPP_
#define _VARLENVECLABEL_CPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "varlenveclabel.hpp"

   
/**  Constructor
 */
CVarLenVecLabel::CVarLenVecLabel()
   : varlenveclabel_verbosity(0),
     Y(0),
     avgLabelDimension(0.0),
     minLabel(SML::INFTY),
     maxLabel(-SML::INFTY),
     numOfAllLabel(0),
     startExample(0),
     labelFile("")     
{ 
        CTimer scanlabeltime;
        CTimer loadlabeltime;
   
        // decide the format string to use
        if(sizeof(Scalar) == sizeof(double))       
                scalar_value_format = "%lf";
        else 
                scalar_value_format = "%f";
        
        // get configurations
        Configuration &config = Configuration::GetInstance();
        
        labelFile = config.GetString("Data.labelFile");
        if(config.IsSet("Data.verbosity"))
                varlenveclabel_verbosity = config.GetInt("Data.verbosity");
        
        // collect some properties of the dataset
        if(varlenveclabel_verbosity >= 1) 
                std::cout << "Scanning label file... "<< std::endl;

        scanlabeltime.Start();
        ScanLabelFile();
        scanlabeltime.Stop();            
            
        // for serial computation, we don't split dataset
        numOfLabel = numOfAllLabel;
        
        // read labels into memory
        if(varlenveclabel_verbosity >= 1)
                std::cout << "Loading label file... "<< std::endl;
      
        loadlabeltime.Start();
        LoadLabels();
        loadlabeltime.Stop();
  
        if(varlenveclabel_verbosity >= 2)
        {           
                std::cout << "scanlabeltime   : " << scanlabeltime.CPUTotal() << std::endl;
                std::cout << "loadlabeltime   : " << loadlabeltime.CPUTotal() << std::endl;
        }
}


// destructor
CVarLenVecLabel::~CVarLenVecLabel()
{}



/**  Determine the following properties of the label file:
 *   1.  type of label
 *   2.  number of labels
 *   3.  dimensionality of label
 *   4.  largest label
 */
void CVarLenVecLabel::ScanLabelFile()
{  
        Scalar tmpLabel = 0;
        std::string line = "";
        std::string token = "";
        std::ifstream labelFp;
        
        numOfAllLabel = 0;

        labelFp.open(labelFile.c_str());
        if(!labelFp.good()) 
        {
                string msg = "Cannot open label file <" + labelFile + ">!";
                throw CBMRMException(msg, "CVarLenVecLabel::ScanLabelFile()");
        }
   
        while(!labelFp.eof()) 
        {
                getline(labelFp, line);
                trim(line);
                if(IsBlankLine(line)) continue;  // blank line
                if(line[0] == '#') continue;  // comment line
                istringstream iss(line);
     
                while(!iss.eof()) 
                {       
                        iss >> token;
                        if(token[0] == '#') break;
                        if(sscanf(token.c_str(), scalar_value_format.c_str(), &tmpLabel) == 1) 
                        {
                                minLabel = std::min(minLabel, tmpLabel);
                                maxLabel = std::max(maxLabel, tmpLabel);
                        } 
                        else 
                        {
                                ostringstream ostr;
                                ostr << "Label file contains invalid label at line " << numOfAllLabel+1 << "!\n" 
                                     << "Token: " << token << "\n";
                                throw CBMRMException(ostr.str(), "CVarLenVecLabel::ScanLabelFile()");
                        }
                }     
                numOfAllLabel++;
        }
   
        if(numOfAllLabel <= 0) 
                throw CBMRMException("Label set is empty!","CVarLenVecLabel::ScanLabelFile()");
   
        labelFp.close();
}



/**  Load labels (the Y-part of the dataset) from file into memory.
 */
void CVarLenVecLabel::LoadLabels()
{
        unsigned int lineCnt = 0;
        unsigned int labelCnt = 0;        // number of labels per example
        std::string line;
        std::ifstream labelFp;
        Scalar tmpLabel;
        vector<Scalar> tmpLabelset;
   
        // open label file
        labelFp.open(labelFile.c_str());

        if(!labelFp.good()) 
        {
                string msg = "Cannot open label file <" + labelFile + ">!";
                throw CBMRMException(msg, "CVarLenVecLabel::LoadLabels()");
        }
   
        // read labels  
        for(lineCnt=0; lineCnt < numOfLabel and !labelFp.eof() ;) 
        {
                // skip lines with comments
                do {
                        getline(labelFp, line);
                        trim(line);
                } while((IsBlankLine(line) || line[0] == '#') && !labelFp.eof());
      
                if(labelFp.eof() && IsBlankLine(line)) break;
               
                istringstream isslabel(line);
      
                for(labelCnt=0; !isslabel.eof(); labelCnt++)  
                {
                        isslabel >> tmpLabel;
                        tmpLabelset.push_back(tmpLabel);
                        avgLabelDimension += 1;
                        
                }
                assert(tmpLabelset.size() != 0);
                Y.push_back(tmpLabelset);
                tmpLabelset.clear();
                lineCnt++;
        }
        assert((unsigned int)Y.size() == numOfAllLabel);
        
        avgLabelDimension /= numOfAllLabel;
        labelFp.close();
}

#endif
