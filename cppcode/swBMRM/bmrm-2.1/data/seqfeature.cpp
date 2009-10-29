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
 * Created: (29/01/2008) 
 *
 * Last Updated:
 */

#ifndef _SEQFEATURE_CPP_
#define _SEQFEATURE_CPP_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "seqfeature.hpp"

   
/**  Constructor
 */
CSeqFeature::CSeqFeature()   
        :seqfeature_verbosity(0),
         maxDuration(0),
         minDuration(0),
         numOfSeq(0),
         numOfAllSeq(0),                
         featureDimension(0),                
         featureFile(""),
         nnz(0),
         density(0)
{
        CTimer loadfeaturetime;

        // decide the format string to use
        if(sizeof(Scalar) == sizeof(double)) 
        {
                svec_feature_index_and_value_format = "%d:%lf";
                scalar_value_format = "%lf";
        } 
        else 
        {
                svec_feature_index_and_value_format = "%d:%f";
                scalar_value_format = "%f";
        }
   
        // get configurations
        Configuration &config = Configuration::GetInstance();
   
        featureFile = config.GetString("Data.featureFile");
        if(config.IsSet("Data.verbosity"))
                seqfeature_verbosity = config.GetInt("Data.verbosity");   
   
        // read dataset into memory
        if(seqfeature_verbosity >= 1)
                std::cout << "Loading feature file... "<< std::endl;
   
        loadfeaturetime.Start();
        LoadFeatures();
        loadfeaturetime.Stop();
          
        // in centralized dataset, serial computation mode
        numOfAllSeq = numOfSeq;
        
        if(seqfeature_verbosity >= 2)
        {           
                std::cout << "loadfeaturetime : " << loadfeaturetime.CPUTotal() << std::endl;
        }
}



/** Read examples into memory
 */
void CSeqFeature::LoadFeatures()
{ 
        unsigned int tmpFidx = 0;
        Scalar tmpFval = 0;
        unsigned int featureCnt = 0;
        unsigned int seqNum = 0;
        unsigned int phiNum = 0;
        unsigned int posNum1 = 0, posNum2 = 0;
        std::string line = "";
        std::string token = "";
        std::ifstream featureFp;
   
        featureFp.open(featureFile.c_str());   
        if(!featureFp.good()) 
        {
                string msg = "Cannot open feature file <" + featureFile + ">!";
                throw CBMRMException(msg, "CSeqFeature::ScanFeatureFile()");
        }
   
        // read header information
        int headerInfoCnt = 3; // min duration, max duration, feature dimension
        do {
                getline(featureFp, line);
                trim(line);
                if(IsBlankLine(line)) continue;  // blank line
                if(line[0] == '#') continue;  // comment line
                if(sscanf(line.c_str(),"maxDuration:%d",&maxDuration)==1) headerInfoCnt--;
                if(sscanf(line.c_str(),"minDuration:%d",&minDuration)==1) headerInfoCnt--;
                if(sscanf(line.c_str(),"globalFeatureDim:%d",&featureDimension)==1) headerInfoCnt--;
        } while(!featureFp.eof() && (headerInfoCnt != 0));
        
        assert(maxDuration >= minDuration);
        assert(featureDimension < (1<<30));  // featureDimension is normally less then 1 billion
                
        if(featureFp.eof())
                throw CBMRMException("Feature file does not contain valid examples","CSeqFeature::LoadFeatures()");
        
        // read sequences
        nnz = 0;
        while(!featureFp.eof()) 
        {
                // read sequence number
                do {
                        getline(featureFp, line);
                        trim(line);
                        if(IsBlankLine(line)) continue;  // blank line
                        if(line[0] == '#') continue;  // comment line
                        if(sscanf(line.c_str(),"sequence:%d",&seqNum)==1) break;
                } while(!featureFp.eof());
                
                if(featureFp.eof())
                        throw CBMRMException("Feature file does not contain valid phi:*","CSeqFeature::LoadFeatures()");
                
                
                // read phi:1 tag
                phiNum = 0;
                do {
                        getline(featureFp, line);
                        trim(line);
                        if(IsBlankLine(line)) continue;  // blank line
                        if(line[0] == '#') continue;  // comment line
                        if(sscanf(line.c_str(),"phi:%d",&phiNum)==1) break;
                } while(!featureFp.eof());
                
                if(featureFp.eof() || (phiNum != 1))
                        throw CBMRMException("Feature file does not contain valid phi:1","CSeqFeature::LoadFeatures()");
                
                // read phi:1 sparse vectors
                do {
                        getline(featureFp, line);
                        trim(line);
                        if(IsBlankLine(line)) continue;  // blank line
                        if(line[0] == '#') continue;  // comment line
                        
                        if(sscanf(line.c_str(),"phi:%d",&phiNum) == 1)
                                break;
                        
                        istringstream iss(line);
                        iss >> token;
                        if((sscanf(token.c_str(),"pos:%d",&posNum1) != 1))
                                throw CBMRMException("Feature file does not contain valid pos tag in phi:1","CSeqFeature::LoadFeatures()");
                        
                        TheMatrix svec(1,featureDimension,SML::SPARSE);
                        featureCnt = 0;
                        while(!iss.eof())
                        {
                                iss >> token;
                                if(sscanf(token.c_str(),svec_feature_index_and_value_format.c_str(),&tmpFidx, &tmpFval) != 2)
                                {
                                        ostringstream msg;
                                        msg << "Invalid #" << featureCnt + 1 << " sparse vector element in phi:"<< phiNum << " seq:" << seqNum << " pos:" << posNum1;
                                        throw CBMRMException(msg.str(),"CSeqFeature::LoadFeatures()");
                                }
                                svec.Set(0,tmpFidx,tmpFval);       
                                nnz++;
                        }
                        
                        if(featureCnt == 0)
                                throw CBMRMException("Feature file does not contain valid phi:2 sparse vector","CSeqFeature::LoadFeatures()");
                        
                        phi_1.push_back(svec);
                } while(!featureFp.eof());
                
                if(phi_1.size() < 1)
                        throw CBMRMException("Feature file does not contain valid phi:1","CSeqFeature::LoadFeatures()");
                
                numOfSeq = phi_1.size();
                
                if(featureFp.eof() || (phiNum != 2))
                        throw CBMRMException("Feature file does not contain valid phi:2","CSeqFeature::LoadFeatures()");
                
                // read phi:2 sparse vectors
                unsigned int prevPosNum1 = 0, prevPosNum2 = 0; 
                vector<TheMatrix> tmp_phi_2_svecs;
                featureCnt = 0;
                do {
                        getline(featureFp, line);
                        trim(line);
                        if(IsBlankLine(line)) continue;  // blank line
                        if(line[0] == '#') continue;  // comment line
                        
                        if((sscanf(line.c_str(),"phi:%d",&phiNum) == 1))
                                break;
                        
                        istringstream iss(line);
                        iss >> token;
                        if((sscanf(token.c_str(),"pos:%d,%d",&posNum1,&posNum2) != 2))
                                throw CBMRMException("Feature file does not containt valid pos tag in phi:2","CSeqFeature::LoadFeatures()");
                        
                        if(prevPosNum2 >= posNum2)
                        {
                                ostringstream msg;
                                msg << "previous posNum2 must be > current posNum2 in phi:2 (phi:2 pos:" << posNum1 << "," << posNum2;
                                throw CBMRMException(msg.str(),"CSeqFeature::LoadFeatures()");
                        }
                        
                        if(prevPosNum1 >= posNum1)
                        {
                                ostringstream msg;
                                msg << "previous posNum1 must be > current posNum1 in phi:2 (phi:2 pos:" << posNum1 << "," << posNum2;
                                throw CBMRMException(msg.str(),"CSeqFeature::LoadFeatures()");
                        }
                        
                        if(posNum1 != prevPosNum1)
                        {
                                phi_2.push_back(tmp_phi_2_svecs);
                                tmp_phi_2_svecs.clear();                                
                        }
                        
                        TheMatrix svec(1,featureDimension,SML::SPARSE);
                        featureCnt = 0;
                        while(!iss.eof())
                        {
                                iss >> token;
                                if(sscanf(token.c_str(),svec_feature_index_and_value_format.c_str(),&tmpFidx, &tmpFval) != 2)
                                {
                                        ostringstream msg;
                                        msg << "Invalid #" << featureCnt + 1 << " sparse vector element in phi:"<< phiNum << " seq:" << seqNum << " pos:" << posNum1;
                                        throw CBMRMException(msg.str(),"CSeqFeature::LoadFeatures()");
                                }
                                svec.Set(0,tmpFidx,tmpFval);    
                                nnz++;
                        }
                        
                        if(featureCnt == 0)
                                throw CBMRMException("Feature file does not containt valid phi:2 sparse vector","CSeqFeature::LoadFeatures()");
                        
                        tmp_phi_2_svecs.push_back(svec);

                } while(!featureFp.eof());
                
                if(phi_2.size() < 1)
                        throw CBMRMException("Feature file does not contain phi:2","CSeqFeature::LoadFeatures()");
        }
        
        // data matrix density
        density = ((double)nnz/featureDimension)/numOfSeq;
   
        featureFp.close();
}


#endif
