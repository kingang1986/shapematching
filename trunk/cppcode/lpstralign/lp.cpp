#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Data.h"
#include "SW.h"
#include "lp.h"
#include "math.h"
#include <glpk.h>

/////////////////////////////////////////////////////////////
//class CSample
/////////////////////////////////////////////////////////////


CSample::CSample()
{
	strcpy(m_szFolder, ".");
	m_nShapes = 0;
	m_nPattern = 0;
}

CSample::~CSample()
{

	Release();

}

int CSample::Release()
{
	for (int i = 0; i < (int) m_vFileNames.size(); i ++)
	{
		delete [] m_vFileNames[i];
	}

	for (int i = 0; i < (int) m_vSS.size(); i ++)
	{
		delete m_vSS[i];
	}
	m_vFileNames.clear();
	m_vSS.clear();
	return 1;
}

CSetOfSeq* CSample::LoadSoS(const char* strSoSFile)
{
	m_nShapes ++;
	for (int i = 0; i < (int)m_vFileNames.size(); i ++)
	{
		if (strcmp(strSoSFile, m_vFileNames[i]) == 0)
		{
			return m_vSS[i];
		}
	}

	CSetOfSeq* pSS = new CSetOfSeq();
	char tmp[500];
	sprintf(tmp, "%s/%s", m_szFolder, strSoSFile);
	//fprintf(stderr, "load: %s from %s in folder %s...", strSoSFile, tmp, m_szFolder);
	if (m_bBinaryData)
	{
		pSS->LoadSSBinary(tmp);
	}
	else
	{
		pSS->LoadSS(tmp);
	}
	char* str = new char[strlen(strSoSFile) + 2];
	strcpy(str, strSoSFile);
	m_vFileNames.push_back(str);
	m_vSS.push_back(pSS);

	fprintf(stderr, "\r %d Shapes (%d uniq) loaded, in %d patterns ", m_nShapes, (int) m_vSS.size(), m_nPattern + 1);

	return pSS;
}

int CSample::LoadSample(const char* strFile)
{

	FILE* fp = fopen(strFile, "r");
	char line[327680];
	char tmp[327680];
	if (fp == NULL)
	{
		fprintf(stderr, "can't open file %s for loading samples\n", strFile);
		exit(1);
	}
	fprintf(stderr, "Loading samples from file ... \n");
   	while(!feof(fp))
    {
    	line[0] = 0;
    	fgets(line, 327680, fp);
    	if (strlen(line) < 2)
    	{
    	    continue;
    	}
    	CPattern* pPattern = new CPattern();
    	chompstr(line);
   		memcpy(tmp, line, 327680);
		char tmpbuffer[512];
		char* pbuffer = tmp;
		int iColumn = count_column(line, ',');
		pbuffer = getnextstring(pbuffer, tmpbuffer, ',');
		pPattern->m_original = LoadSoS(tmpbuffer);
		pbuffer = getnextstring(pbuffer, tmpbuffer, ',');
		pPattern->m_iHomolog = atoi(tmpbuffer);
             //  fprintf(stderr, "%d homologs\n", pPattern->m_iHomolog);
		pPattern->m_iDecoy = iColumn - 2 - pPattern->m_iHomolog;
		for (int k = 0; k < pPattern->m_iHomolog; k ++)
	    {
    		pbuffer = getnextstring(pbuffer, tmpbuffer, ',');
        	pPattern->m_vHomolog.push_back(LoadSoS(tmpbuffer));
        	//fprintf(stderr, "%f ", atof(tmpbuffer));
    	}
		for (int k = 0; k < iColumn - 2 - pPattern->m_iHomolog; k ++)
	    {
    		pbuffer = getnextstring(pbuffer, tmpbuffer, ',');
        	pPattern->m_vDecoy.push_back(LoadSoS(tmpbuffer));
        	//fprintf(stderr, "%f ", atof(tmpbuffer));
    	}
    	m_vPatterns.push_back(pPattern);
    	m_nPattern ++;
    }
    fclose(fp);
    fprintf(stderr, " done\n");
    m_vActive.clear();
    for (int i = 0; i <m_nPattern; i ++)
    {
    	m_vActive.push_back(1);
    }
    return (int) m_vPatterns.size();

}

int CSample::SetAllActive()
{
	for (int i = 0; i < m_nPattern; i ++)
	{
		m_vActive[i] = 1;
	}
	return m_nPattern;
}
void CSample::UpdateHomologScore(CModel* pModel)
{
	int iParamDim = pModel->m_iParamDim;
	double* pw1 = new double[iParamDim];
	for (int i = 0; i < m_nPattern; i ++)
	{
		CPattern* pPattern = m_vPatterns[i];
		pPattern->m_fHomologScore = GetPhi(&(pPattern->m_homo_align), pw1, iParamDim, pModel);
	}
	delete[] pw1;
	return;
}
double CSample::AddConstraint(CConstraints* pCC, CModel* pModel, bool bActiveOnly)
{
	int iParamDim = pModel->m_iParamDim;
	double fmaxvio = 0;
	for (int i = 0; i < m_nPattern; i ++)
	{
		if (bActiveOnly && m_vActive[i] == 0)
		{
			continue;
		}
		double* pw1 = new double[iParamDim];
		double* pw2 = new double[iParamDim];
		CPattern* pPattern = m_vPatterns[i];
		double f1 = GetPhi(&(pPattern->m_homo_align), pw1, iParamDim, pModel);
		double f2 = GetPhi(&(pPattern->m_decoy_align), pw2, iParamDim, pModel);
		double vio =  f2 + pCC->m_fDistance - f1 - 0.001;
		fprintf(stderr, "%f - %f > %f - epsilon, violoate %f\n", f1, f2, pCC->m_fDistance, vio);
		for (int k = 0; k < iParamDim; k ++)
		{
			pw1[k] = pw1[k] - pw2[k];
		}
		//if (fmaxvio < vio) 			fmaxvio = vio;
		if (vio > 0 ) 			fmaxvio += vio;

		if (vio <= 0.0001)
		{
			m_vActive[i] = 0;
		}
		else
		{
			pCC->Add(pw1, i);
		}
	}
	return fmaxvio;
}

int CSample::AlignHomolog(CModel* pModel)
{
	for (int i = 0; i < m_nPattern; i ++)
	{

		CPattern* pPattern = m_vPatterns[i];
		double fmaxscore = -1.0;
		for (int j = 0; j < pPattern->m_iHomolog; j ++)
		{
			CSetOfSeq* pSS = pPattern->m_vHomolog[j];
			CAlignment align;
			double fscore = SmithWaterman_SetOfSeq(pPattern->m_original,  pSS, &align, pModel);
			if (fscore > fmaxscore)
			{
				pPattern->m_homo_align = align;
				pPattern->m_iBestHomolog = j;
				fmaxscore = fscore;
				pPattern->m_fHomologScore = fscore;

			}
		}
		fprintf(stderr, "\r %d of %d, %f", i + 1, m_nPattern, pPattern->m_fHomologScore);
	}
	return m_nPattern;

}
int CSample::AlignSamples(CModel* pModel, const char* outputfile)
{
        fprintf(stderr, "classifying samples ... ");
        FILE* fp = fopen(outputfile, "w");
        if (!fp)
        {
                fprintf(stderr, "cant' write result to %s. terminated. \n", outputfile);
                return -1;
        }
        int icorrect = 0;
        for (int i = 0; i < m_nPattern; i ++)
        {

                CPattern* pPattern = m_vPatterns[i];
                double fmaxscore = -1.0;
                for (int j = 0; j < pPattern->m_iHomolog; j ++)
                {
                       //fprintf(stderr, "%d homologs for testing \n", pPattern->m_iHomolog);
                        CSetOfSeq* pSS = pPattern->m_vHomolog[j];
                        CAlignment align;
                        double fscore = SmithWaterman_SetOfSeq(pPattern->m_original,  pSS, &align, pModel);
                        fprintf(fp, "%3.3f ", fscore);
                        if (fscore > fmaxscore)
                        {
                                pPattern->m_homo_align = align;
                                pPattern->m_iBestHomolog = j;
                                fmaxscore = fscore;
                                pPattern->m_fHomologScore = fscore;
                        }
                }
                //fprintf(fp, "\n"); 
                fmaxscore = -1.0;
                for (int j = 0; j < pPattern->m_iDecoy; j ++)
                {
                        CSetOfSeq* pSS = pPattern->m_vDecoy[j];
                        CAlignment align;
                        double fscore = SmithWaterman_SetOfSeq(pPattern->m_original,  pSS, &align, pModel);
                        if ( j % 20 == 0) fprintf(fp, "\n");
                        fprintf(fp, "%3.3f ", fscore);
                        
                        if (fscore > fmaxscore)
                        {
                                pPattern->m_decoy_align = align;
                                pPattern->m_iBestDecoy = j;
                                fmaxscore = fscore;
                                pPattern->m_fDecoyScore = fscore;
                        }
                }
                if (pPattern->m_fHomologScore >= pPattern->m_fDecoyScore) icorrect += 1;

                fprintf(stderr, "\r %d of %d, %f, %f, %d correct  ", i + 1, m_nPattern, pPattern->m_fHomologScore, pPattern->m_fDecoyScore, icorrect);
                //fprintf(fp, " %f %f\n", pPattern->m_fHomologScore, pPattern->m_fDecoyScore, icorrect);
               
        }
        //fprintf(fp, "%d of %d correct\n", icorrect, m_nPattern); 
        fprintf(fp, "\n"); 
        fprintf(stderr, "done\n");
        fclose(fp);
        return m_nPattern;

}


double CSample::AlignDecoy(CConstraints* pCC, CModel* pModel, bool bActiveOnly, int iMinNewConstraint, int& iNumVio)
{
	int iNewConstraint = 0;
	double fsumvio = 0.0;
	iNumVio = 0;
	for (int i = 0; i < m_nPattern; i ++)
	{
		if (bActiveOnly && m_vActive[i] == 0)
		   continue;
		CPattern* pPattern = m_vPatterns[i];

		double fmaxdecoy = -1;
		for (int j = 0; j < pPattern->m_iDecoy; j ++)
		{
			//fprintf(stderr, "%ok1\n");
			CSetOfSeq* pSS = pPattern->m_vDecoy[j];
			//fprintf(stderr, "%ok2\n");
			CAlignment align;
			//fprintf(stderr, "%ok3\n");
			double fscore = SmithWaterman_SetOfSeq(pPattern->m_original,  pSS, &align, pModel);
			//fprintf(stderr, "%d, %f\n", j, fscore);
			if (fscore > fmaxdecoy)
			{
				//fprintf(stderr, "%ok5\n");
				pPattern->m_decoy_align = align;
				pPattern->m_iBestDecoy = j;
				pPattern->m_fDecoyScore = fscore;
				fmaxdecoy = fscore;
				//fprintf(stderr, "%ok6\n");
			}
		}
		double vio = pPattern->m_fDecoyScore + pCC->m_fDistance - pPattern->m_fHomologScore - 0.001;
		fprintf(stderr, "\rPattern %d, homolog: %3.4f, decoy: %3.4f, %3.4f ", i + 1, pPattern->m_fHomologScore, pPattern->m_fDecoyScore, pPattern->m_fHomologScore - pPattern->m_fDecoyScore);
		if (vio <= 0.0001)
		{
			m_vActive[i] = 0;
			//fprintf(stderr, "\n");
		}
		else
		{
			int iParamDim = pModel->m_iParamDim;
			//fprintf(stderr, "%d ", iParamDim);
			double* pw1 = new double[iParamDim];
			double* pw2 = new double[iParamDim];
			CPattern* pPattern = m_vPatterns[i];
			GetPhi(&(pPattern->m_homo_align), pw1, iParamDim, pModel);
			GetPhi(&(pPattern->m_decoy_align), pw2, iParamDim, pModel);
                        double fsum = 0, fs2 = 0;
                        for (int k = 0; k < iParamDim; k ++)
                        {
  				fsum += pw1[k] * pModel->m_vWeight[k];
  				fs2 += pw2[k] * pModel->m_vWeight[k];
                         }

			//fprintf(stderr, "get phi ok, %f == %f ?, %f = %f?\n", fsum, pPattern->m_fHomologScore, fs2, pPattern->m_fDecoyScore);

			for (int k = 0; k < iParamDim; k ++)
			{
				pw1[k] = pw1[k] - pw2[k];
			}

			pCC->Add(pw1, i);
			iNewConstraint ++;
			fprintf(stderr, " ++\n");
			if (iNewConstraint >= iMinNewConstraint)
                        {
                             pCC->GLPK_lp(pModel);
                             UpdateHomologScore(pModel);
                             fprintf(stderr, "\n%d contraints, total %d violation, sum to %f \n", pCC->m_vWeights.size(), iNumVio, fsumvio);
                             iNewConstraint = 0;
                        }
		}
		if (vio > 0 )
		{  fsumvio += vio; iNumVio ++ ;}
	}


	return fsumvio;
}


int CSample::GetActiveNum()
{
	int count = 0;
	for (int i = 0; i < m_nPattern; i ++)
	{
		count += m_vActive[i];
	}
	return count;
}


//////////////////////////////////////////////////////
//class CConstraints
//////////////////////////////////////////////////////

CConstraints::CConstraints()
{
	m_fC = 100;
}
CConstraints::~CConstraints()
{
	Clear();
}
int CConstraints::Add(double* pWeight, int iPatternIndex)
{
	double * pw = new double[m_iWeightLength];
	memcpy(pw, pWeight, sizeof(double) * m_iWeightLength);
	m_vWeights.push_back(pw);
	m_vPatternIndex.push_back(iPatternIndex);
	return (int) m_vWeights.size();
}

void CConstraints::Clear()
{
	for (int i = 0; i < (int) m_vWeights.size(); i ++)
	{
		delete[] m_vWeights[i];
	}
	m_vWeights.clear();
	m_vPatternIndex.clear();
	return;
}
int CConstraints::Init(CModel* model, int iPatternNum, double fEpsilon)
{
	m_iWeightLength = model->m_iParamDim;
	m_fEpsilon = fEpsilon;
	m_iPatternNum = iPatternNum;
	return m_iWeightLength;
}

int CConstraints::PrintMathProg(const char* filename)
{
	FILE* fp = fopen(filename, "w");
	if (!fp)
	{
		fprintf(stderr, "Cant' write to file %s \n", filename);
		return -1;
	}
	//print variable definitions
	fprintf(fp, "/*Variable definitions*/\n");
	for (int i = 0 ;i < m_iWeightLength; i ++)
	{
		fprintf(fp, "var w%d >=0;\n", i + 1);
	}
	for (int i = 0 ;i < m_iPatternNum; i ++)
	{
		fprintf(fp, "var e%d >=0;\n", i + 1);
	}
	fprintf(fp, "/*Objective function*/\n");
	fprintf(fp, "minimize z: ");
	for (int i = 0 ;i < m_iWeightLength; i ++)
	{
		fprintf(fp, "+ w%d ", i + 1);
	}
	for (int i = 0 ;i < m_iPatternNum; i ++)
	{
		fprintf(fp, "+ %f*e%d ", m_fC, i + 1);
	}
	fprintf(fp, "; \n");

	fprintf(fp, "/*Constraints*/\n");

	for (int i = 0; i < (int) m_vWeights.size(); i ++)
	{
		fprintf(fp, "s.t. cc_%d: ", i + 1);
		double* pd = m_vWeights[i];
		for (int k = 0; k < m_iWeightLength; k ++)
		{
			fprintf(fp, "%+f*w%d ", pd[k], k + 1);
		}
		fprintf(fp, " + e%d ", m_vPatternIndex[i] + 1);
		fprintf(fp, " >= 1 - %f;\n", m_fEpsilon);
	}
	fprintf(fp, "end;\n");

	fclose(fp);
	return (int)m_vWeights.size();
	return 1;
}

int CConstraints::Save(const char* file)
{
	FILE* fp = fopen(file, "w");
	if (!fp)
	{
		fprintf(stderr, "can't open file %s for writing constraints.\n", file);
		return -1;
	}
	fprintf(fp, "%d %d %d\n", m_iWeightLength, m_iPatternNum, (int) m_vWeights.size());
	fprintf(fp, "%f %f\n", m_fC, m_fEpsilon);
	for (int i = 0; i < (int) m_vPatternIndex.size(); i ++)
	{
		fprintf(fp, "%d ", m_vPatternIndex[i]);
	}
	fprintf(fp, "\n");
	for (int i = 0; i < (int) m_vWeights.size(); i ++)
	{
		double * fw = m_vWeights[i];
		for (int k = 0; k < m_iWeightLength; k ++)
		{
			fprintf(fp, "%f ", fw[k]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return 0;
}
int CConstraints::Load(const char* file)
{
	Clear();
	FILE* fp = fopen(file, "r");
	if (!fp)
	{
		fprintf(stderr, "can't open file %s for reading constraints.\n", file);
		return -1;
	}
	int iConstraintsNum = 0;
	fscanf(fp, "%d %d %d", &m_iWeightLength, &m_iPatternNum, &iConstraintsNum);
	fscanf(fp, "%g %g\n", &m_fC, &m_fEpsilon);
	for (int i = 0; i < iConstraintsNum; i ++)
	{
		int patternindex = 0;
		fscanf(fp, "%d ", &patternindex);
		m_vPatternIndex.push_back(patternindex);
	}
	for (int i = 0; i < iConstraintsNum; i ++)
	{
		double * fw = new double[m_iWeightLength];
		float w = -1;
		for (int k = 0; k < m_iWeightLength; k ++)
		{
			fscanf(fp, "%f", &w);
			fw[k] = w;
		}
		m_vWeights.push_back(fw);
	}
	fclose(fp);
	return 0;


	return 0;
}

int CConstraints::GLPK_lp(CModel* pmodel)
{
	glp_prob* lp = glp_create_prob();


	int iRow = (int) m_vWeights.size();
	//int iCol = (int) m_iWeightLength + m_iPatternNum;
	int iSize = iRow * ( m_iWeightLength + 1);
    int ia[10 + iSize], ja[1 + iSize];
    double ar[1 + iSize];
    glp_set_prob_name(lp, "StrLP");
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_rows(lp, (int) m_vWeights.size());
    // setup rows
	for (int i = 0; i < (int) m_vWeights.size(); i ++)
	{
		char tmp[200];
		sprintf(tmp, "cc_%d", i + 1);
		glp_set_row_name(lp, i + 1, tmp);
		glp_set_row_bnds(lp, i + 1, GLP_LO, m_fDistance - m_fEpsilon, 0);
	}
	glp_add_cols(lp, m_iWeightLength + m_iPatternNum);
	for (int i = 0; i < m_iWeightLength; i ++)
	{
		char tmp[200];
		sprintf(tmp, "w%d", i + 1);
		glp_set_col_name(lp, i + 1, tmp);
		glp_set_col_bnds(lp, i + 1, GLP_LO, 0, 0.0);
		glp_set_obj_coef(lp, i + 1, 1.0);
	}

	for (int i = 0; i < m_iPatternNum; i ++)
	{
		char tmp[200];
		sprintf(tmp, "e%d", i + 1);
		glp_set_col_name(lp, m_iWeightLength + i + 1, tmp);
		glp_set_col_bnds(lp, m_iWeightLength + i + 1, GLP_LO, 0, 0.0);
		glp_set_obj_coef(lp, m_iWeightLength + i + 1, m_fC / m_iPatternNum);
	}
	int iIndex = 1;
	for (int i = 0; i < (int)m_vWeights.size(); i ++)
	{
		double* pd = m_vWeights[i];
		for (int j = 0; j < (int) m_iWeightLength; j ++)
		{
			ia[iIndex] = i + 1, ja[iIndex] = j + 1;
			if (pmodel->m_vSign[j] <= 0)
			{
				ar[iIndex] = -pd[j];
			}
			else
			{
				ar[iIndex] = pd[j];
			}
			iIndex ++;
		}
		ia[iIndex] = i + 1;
		ja[iIndex] = m_iWeightLength + m_vPatternIndex[i] + 1;
		ar[iIndex] = 1;
		iIndex ++;
	}
    glp_load_matrix(lp, iIndex - 1, ia, ja, ar);
    glp_simplex(lp, NULL);
    double z = glp_get_obj_val(lp);
    fprintf(stderr, "minimal value %f \n", z);
    for (int i = 0; i < m_iWeightLength; i ++)
    {
    	double x = glp_get_col_prim(lp, i + 1);
    	if (pmodel->m_vSign[i] <=0)
    	{
    		pmodel->m_vWeight[i] = -x;
    	}
    	else
    	{
    		pmodel->m_vWeight[i] = x;
    	}

    	if (x != 0) fprintf(stderr, "(w%d, %f)\t", i + 1, pmodel->m_vWeight[i]);
    }
    for (int i = 0; i < m_iPatternNum; i ++)
    {
    	double x = glp_get_col_prim(lp, m_iWeightLength + i + 1);
    	pmodel->m_vTheta[i] = x;
    	if (x != 0)  fprintf(stderr, "(e%d, %f)\t", i + 1, pmodel->m_vTheta[i]);
    }
    glp_delete_prob(lp);
	return 1;
}


//////////////////////////////////////////////////////
//class CStructureLearning
//////////////////////////////////////////////////////

CStructureLearning::CStructureLearning()
{
	m_fDistance = 1;
	m_iMaxIteration = 12;
	m_iMinNewConstraint = 10;
}
CStructureLearning::~CStructureLearning()
{
}
int CStructureLearning::FindMostViolate(CPattern* p)
{
	double fMaxScore = -1;
	CSetOfSeq* pOriginal = p->m_original;
	int iIndex = -1;
	for (int i = 0; i < p->m_iDecoy; i ++)
	{
		CSetOfSeq* pDecoy = p->m_vDecoy[i];
		CAlignment align;
		double fscore = SmithWaterman_SetOfSeq(pOriginal, pDecoy, &align, &m_Model);
		if (fscore > fMaxScore)
		{
			fMaxScore = fscore;
			p->m_decoy_align = align;
			iIndex = i;
		}
	}
	p->m_iBestDecoy = iIndex;
	return iIndex;
}

int CStructureLearning::Init(const char* datafile, const char* szpath, const char* modelfile, bool bBinaryData)
{
	strcpy(m_Sample.m_szFolder, szpath);
	m_Sample.LoadSample(datafile);
	m_Model.Read(modelfile);
	m_Model.Print();
	m_Model.InitTheta(m_Sample.m_nPattern);
	m_Sample.m_bBinaryData = bBinaryData;
	return 1;
}



int CStructureLearning::Learn(const char* outputfile)
{
	fprintf(stderr, "aligning homologs ...");
	m_Sample.AlignHomolog(&m_Model);
	fprintf(stderr, "done\n");
	CConstraints CC;
	CC.m_fC = m_fC;
	CC.m_fEpsilon =  m_fEpsilon;
	CC.m_fDistance = m_fDistance;
	CC.Init(&m_Model, m_Sample.m_nPattern, 0.001);
	double fepsilon = 1e9;
	int iIteration = 0;
        do
	{
		m_Sample.SetAllActive();
		fepsilon = 1e9;
		double pre = 1e10;
		int iStep = 0;
		int iNumVio = 0;
                do 
		{
			iStep ++;
			fprintf(stderr, "%d iteration %d step \n", iIteration, iStep);
			fprintf(stderr, "active patterns %d (out of %d)\n", m_Sample.GetActiveNum(), m_Sample.m_nPattern);
			pre = fepsilon;
			fepsilon =  m_Sample.AlignDecoy(&CC, &m_Model, true, m_iMinNewConstraint, iNumVio);
			//CC.PrintMathProg("tmp.mod");
			//CC.GLPK_lp(&m_Model);
			//m_Sample.UpdateHomologScore(&m_Model);
			fprintf(stderr, "\n%d contraints, total %d violation, sum to %f \n", CC.m_vWeights.size(), iNumVio, fepsilon);
		}while ((m_Sample.GetActiveNum() > 0) &&  (fepsilon > 0) && (fabs(pre - fepsilon) > 0.0001) && (iStep < m_iMaxStep) );
                iIteration ++;
	}while ((iIteration < m_iMaxIteration) && (fepsilon > 0));
	m_Model.Write(outputfile);

	return 1;

}

