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
 * Authors: S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
 *          Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *
 * Created: (02/11/2007) 
 *
 * Last Updated: (07/11/2007)   
 */

#ifndef _DATAFACTORY_HPP_
#define _DATAFACTORY_HPP_

#include <string>

#include "data.hpp"
#include "configuration.hpp"
#include "vecdata.hpp"
#include "multilabelvecdata.hpp"
#include "bmrmexception.hpp"
#include "genericdata.hpp"

/**  
 * Factory class for creating new Data instances. 
 *
 * When you subclass the CData class you must also add an entry into
 * this factory class (grep YOUR_DATA_FORMAT for an example). 
 */
class CDataFactory
{
      
   public:
      
      /**  Return a data object based on user's argument in configuration file
       *
       *   @return data object
       */
      static CData* GetData(void)
      {         
         CData *ds = 0;
         Configuration &config = Configuration::GetInstance();
         
         // default to this format
         std::string dataFormat = "VECTOR_LABEL_VECTOR_FEATURE";
         
         // unless the user specifies otherwise in the config 
         if(config.IsSet("Data.format"))
            dataFormat = config.GetString("Data.format");
         
         // svnvish: BUGBUG
         // Want to do this with a switch 
         if(dataFormat == "VECTOR_LABEL_VECTOR_FEATURE")
         {
            ds = new CVecData();            
         }
         else if(dataFormat == "VARIABLE_LENGTH_VECTOR_LABEL_VECTOR_FEATURE")
         {
            ds = new CMultilabelVecData();            
         }
         else if(dataFormat == "GENERIC")
         {
            ds = new CGenericData();            
         }
         // else if(dataFormat == "YOUR_DATA_FORMAT")
         //{
         //   ds = new CYourDataFormat();
         //}
         else
         {
            throw CBMRMException("ERROR: unrecognised data format ("+dataFormat+")\n", "CDataFactory::GetData()");
         }
         return ds;
      }
      
};

#endif
