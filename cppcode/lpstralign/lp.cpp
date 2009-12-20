#include <stdio.h>
#include <vector>
#include <algorithm>
#include <string.h>
#include <stdlib.h>
#include "sequence.h"
#include "model.h"
#include "alignment.h"
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
    m_iLastUpdated = -1;
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

CSetOfSeq* CSample::FindSoS(const char* strFile)
{
    for (int i = 0; i < (int)m_vFileNames.size(); i ++)
    {
        if (strcmp(strFile, m_vFileNames[i]) == 0)
        {
            return m_vSS[i];
        }
    }
    return NULL;
}

CSetOfSeq* CSample::LoadSoS(const char* strSoSFile)
{

    CSetOfSeq* pSS = new CSetOfSeq();
    char tmp[500];
    sprintf(tmp, "%s%s", m_szFolder, strSoSFile);
//    if (m_bBinaryData)
//    {
//        pSS->LoadSSBinary(tmp);
//    }
//    else
//    {
        pSS->Load(tmp);
//    }


    return pSS;
}


int CSample::LoadShape(const char* strFile)
{

    FILE* fp = fopen(strFile, "r");
    char line[32768];
    char tmp[32768];
    if (fp == NULL)
    {
        fprintf(stderr, "can't open file %s for loading samples\n", strFile);
        exit(1);
    }
    fprintf(stderr, "Loading shapes from file %s ... \n", strFile);
    while(!feof(fp))
    {
        line[0] = 0;
        fgets(line, 32768, fp);
        if (strlen(line) < 2)
        {
            continue;
        }
        chompstr(line);
        memcpy(tmp, line, 32768);
        char tmpbuffer[512];
        char* pbuffer = tmp;
        pbuffer = getnextstring(pbuffer, tmpbuffer, ',');
        CSetOfSeq* pSS = LoadSoS(tmpbuffer);
        char* str = new char[strlen(tmpbuffer) + 2];
        strcpy(str, tmpbuffer);
        m_vFileNames.push_back(str);
        pbuffer = getnextstring(pbuffer, tmpbuffer, ',');
        pSS->m_iShapeID = atoi(tmpbuffer);
        pbuffer = getnextstring(pbuffer, tmpbuffer, ',');
        pSS->m_iClassID = atoi(tmpbuffer);
        
//        fprintf(stderr, " %s, %d %d\n", str, pSS->m_iShapeID, pSS->m_iClassID);
        m_vSS.push_back(pSS);
    }
    fclose(fp);
    fprintf(stderr, " done! %d shapes loaded\n", int(m_vSS.size()));
    return (int) m_vSS.size();

}
int CSample::LoadPattern(const char* strFile)
{

    FILE* fp = fopen(strFile, "r");
    char line[32768];
    char tmp[32768];
    if (fp == NULL)
    {
        fprintf(stderr, "can't open file %s for loading samples\n", strFile);
        exit(1);
    }
    //fprintf(stderr, "Loading samples from file ... \n");
    while(!feof(fp))
    {
        line[0] = 0;
        fgets(line, 32768, fp);
        if (strlen(line) < 2)
        {
            continue;
        }
        CPattern* pPattern = new CPattern();
        chompstr(line);
        memcpy(tmp, line, 32768);
        char tmpbuffer[512];
        char* pbuffer = tmp;
        int iColumn = count_column(line, ',');
        pbuffer = getnextstring(pbuffer, tmpbuffer, ',');
        pPattern->m_original = FindSoS(tmpbuffer);
        for (int k = 0; k < iColumn - 1; k ++)
        {
            pbuffer = getnextstring(pbuffer, tmpbuffer, ',');
            CSetOfSeq* pSS = FindSoS(tmpbuffer);
            
 //           fprintf(stderr, " %s, %d, %d ", tmpbuffer,pSS->m_iShapeID,  pSS->m_iClassID);
            pPattern->m_vShapeClass[pSS->m_iClassID].push_back(pSS);
        }
        m_vPatterns.push_back(pPattern);
    }
    fclose(fp);
    fprintf(stderr, " done! %d patterns loaded\n", (int)m_vPatterns.size());
    m_vActive.clear();
    for (int i = 0; i < (int) m_vPatterns.size(); i ++)
    {
        m_vActive.push_back(1);
    }
    return (int) m_vPatterns.size();

}


int CSample::SetAllActive()
{
    for (int i = 0; i < (int)m_vPatterns.size(); i ++)
    {
        m_vActive[i] = 1;
    }
    return (int)m_vPatterns.size();
}


int CSample::AlignHomolog(CModel* pModel)
{
    fprintf(stderr, "\n");
    for (int i = 0; i < (int) m_vPatterns.size(); i ++)
    {
        CPattern* pPattern = m_vPatterns[i];
        fprintf(stderr, "\rPattern %d   ", i);
        pPattern->AlignHomolog(pModel);
    }
    return (int)m_vPatterns.size();

}


int CSample::Classify(CModel* pModel, const char* outputfile)
{
    fprintf(stderr, "classifying samples ... ");
    FILE* fp = fopen(outputfile, "w");
    if (!fp)
    {
        fprintf(stderr, "cant' write result to %s. terminated. \n", outputfile);
        return -1;
    }
    double fSumAcc = 0.0;
    for (int i = 0; i < (int) m_vPatterns.size(); i ++)
    {
        CPattern* pPattern = m_vPatterns[i];
        fprintf(stderr, "\n%d: ", i + 1);
        fSumAcc += pPattern->Classify(pModel);
    }
    fprintf(stderr, "sample # %d, avg accuracy %f \n", (int)m_vPatterns.size(), fSumAcc/(double) m_vPatterns.size());
    fclose(fp);
    return  (int) m_vPatterns.size();

}


//add constraints upto iMinNewConstraint, starting from last updated pattern

double CSample::UpdateConstraint(CConstraints* pCC, CModel* pModel, bool bActiveOnly, int iMinNewConstraint, int& iNumVio, double fEpsilon)
{
    int iNewConstraint = 0;
    double fsumvio = 0.0;
    int iTotalPattern = m_vPatterns.size();
    iNumVio = 0;
    double vio = 0.0;
    int lastupdate = m_iLastUpdated;
    for (int i = lastupdate + 1; i < iTotalPattern + lastupdate + 1 && iNewConstraint < iMinNewConstraint - 1; i ++)
    {
        int ii  = i;
        while (ii >= iTotalPattern) ii -= iTotalPattern;
        if (bActiveOnly && m_vActive[ii] == 0)
            continue;
        CPattern* pPattern = m_vPatterns[ii];
        pPattern->AlignDecoy(pModel);
        int iParamDim = pModel->m_iParamDim;
        double* pw1 = new double[iParamDim];
        double* pw2 = new double[iParamDim];
        double labelloss = pPattern->GetLabelLoss(pModel);
        pPattern->GetLabelPhi(pw1, iParamDim, pModel);
//        fprintf(stderr, " label loss %.3f ", labelloss);
        int iDecoyClass = pPattern->GetDecoyPhi(pw2, iParamDim, pModel);
        double fsum = 0, fs2 = 0;
        for (int k = 0; k < iParamDim; k ++)
        {
            fsum += pw1[k] * pModel->m_vWeight[k];
            fs2 += pw2[k] * pModel->m_vWeight[k];
        }
        vio = fsum - fs2 - labelloss;
        fprintf(stderr, "\n %d (class %d) %.3f, (%d %.3f %.3f) ", ii, pPattern->m_original->m_iClassID, fsum, iDecoyClass, fs2, -vio);
        if (vio > - fEpsilon)
        //if (labelloss < 0.00001)
        {
            m_vActive[ii] = 0;
        }
        else
        {
            for (int k = 0; k < iParamDim; k ++)
            {
                pw1[k] = pw1[k] - pw2[k];
            }

            pCC->Add(pw1, ii, labelloss);
            iNewConstraint ++;
            iNumVio ++;

            fsumvio += vio;
            fprintf(stderr, " ++  ");
        }
        m_iLastUpdated = ii;
    }
    return -fsumvio;
}


int CSample::GetActiveNum()
{
    int count = 0;
    for (int i = 0; i < (int) m_vPatterns.size(); i ++)
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


int CConstraints::Add(double* pWeight, int iPatternIndex, double loss)
{
    double * pw = new double[m_iWeightLength];
    memcpy(pw, pWeight, sizeof(double) * m_iWeightLength);
    m_vWeights.push_back(pw);
    m_vPatternIndex.push_back(iPatternIndex);
    m_vLoss.push_back(loss);
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
        fprintf(fp, " >= %f - %f;\n", m_vLoss[i], m_fEpsilon);
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
    fscanf(fp, "%lf %lf\n", &m_fC, &m_fEpsilon);
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
        //ar[iIndex] = 1;
        ar[iIndex] = m_vLoss[i];
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
//class CPattern
//////////////////////////////////////////////////////
int CPattern::m_iTopK = 3;
int CPattern::m_iTotalClass = 10;
CPattern::CPattern()
{
    m_pLabelPhi = NULL;
    for (int i = 0; i < m_iTotalClass; i ++)
    {
        vector<CSetOfSeq*> pSet;
        m_vShapeClass.push_back(pSet);
        vector<CAlignment> pAlign;
        m_vShapeAlign.push_back(pAlign);
        vector<CAlignment*> pSortedAlign;
        m_vSortedShapeAlign.push_back(pSortedAlign);
    }
}

CPattern::~CPattern()
{
    if (m_pLabelPhi != NULL)
    {
        delete[] m_pLabelPhi;
        m_pLabelPhi = NULL;
    }
}

bool cmpAlign(const CAlignment* a1, const CAlignment* a2)
{
    return a1->m_fScore > a2->m_fScore;
}


void CPattern::AlignDecoy(CModel* pModel)
{
    for (int i = 0; i < m_iTotalClass; i ++)
    {
        if (i == m_original->m_iClassID)
             continue;
        AlignClass(i, pModel);
    }
}

int CPattern::AlignClass(int iClass, CModel* pModel)
{
    vector<CAlignment>& aligns = m_vShapeAlign[iClass];   
    vector<CAlignment*>& sortedaligns = m_vSortedShapeAlign[iClass];   
    vector<CSetOfSeq*>& shapes = m_vShapeClass[iClass];   
    aligns.clear();
    for (int j = 0; j <(int) shapes.size(); j ++) 
    {
        CSWMatch * pMatch = new CSWMatch();
        CSetOfSeq* pSS = shapes[j];
        CAlignment align;
        pMatch->Match(m_original, pSS, &align, pModel);
//        fprintf(stderr, "[%d %d %f] ", iClass, j, align.m_fScore);
        aligns.push_back(align);
    }
    sortedaligns.clear();
    for (int i = 0; i < (int)aligns.size(); i ++)
        sortedaligns.push_back(&aligns[i]);
    sort(sortedaligns.begin(), sortedaligns.end(), cmpAlign);
/*
    fprintf(stderr, " (%d : ",  iClass);
    for (int i = 0; i < m_iTopK; i ++)
    {
        fprintf(stderr, " %.3f", sortedaligns[i]->m_fScore);
    }

    fprintf(stderr, " ) ");
*/
    return (int) aligns.size(); 
}

void CPattern::AlignHomolog(CModel* pModel)
{
    AlignClass(m_original->m_iClassID, pModel);
}

double CPattern::GetLabelLoss(CModel* pModel)
{
    return 1.0; 
}

double CPattern::GetLabelPhi(double* pw, int iParamDim, CModel* pModel)
{
    return GetClassPhi(pw, iParamDim, m_original->m_iClassID, pModel);
}

int CPattern::GetDecoyPhi(double* pw, int iParamDim, CModel* pModel)
{
    double* tmp = new double[iParamDim];
    double fmaxscore = -1;
    int maxClass = -1;
    for (int i = 0; i < m_iTotalClass; i ++)
    {
        if (i == m_original->m_iClassID) continue;
        double fscore = GetClassPhi(tmp, iParamDim, i, pModel);
        if (fscore > fmaxscore)
        {
            fmaxscore = fscore;
            memcpy(pw, tmp, sizeof(double)*iParamDim);
            maxClass = i;
        }
          
    }
    return maxClass;
}
  

double CPattern::GetClassPhi(double* pw, int iParamDim, int iClassID, CModel* pModel) 
{
    memset(pw, 0, sizeof(double) * iParamDim);
    double* tmp = new double[iParamDim];
    vector<CAlignment*> & sortedaligns = m_vSortedShapeAlign[iClassID];
    for (int i = 0; i < m_iTopK; i ++) 
    {
        sortedaligns[i]->GetPhi(tmp, iParamDim, pModel); 
        for (int j = 0; j < iParamDim; j ++) 
            pw[j] += tmp[j];
    } 
    double fsum = 0;
    for (int j = 0; j < iParamDim; j ++) 
    {
        pw[j] /= m_iTopK; 
        fsum += pw[j] * pModel->m_vWeight[j];
    }
    
    return fsum;
}

/*
void CPattern::UpdateHomologScore(int iParamDim, CModel* pModel)
{
    double * phi = new double[iParamDim];
    for (int i = 0; i <(int) m_vLabelHomoAlign.size(); i ++)
    {
        CAlignment* palign = &m_vLabelHomoAlign[i];
        palign->GetPhi(phi, iParamDim, pModel);
        palign->m_fScore = 0;
        for (int k = 0; k < iParamDim; k ++)
            palign->m_fScore += phi[k] * pModel->m_vWeight[k];
    }
    m_vSortedLabel.clear();
    for (int i = 0; i < (int)m_vLabelHomoAlign.size(); i ++)
        m_vSortedLabel.push_back(&m_vLabelHomoAlign[i]);
    sort(m_vSortedLabel.begin(), m_vSortedLabel.end(), cmpAlign);
    for (int i = 0; i < m_iTopK; i ++)
    {
        //fprintf(stderr, "%.3f ", m_vSortedLabel[i]->m_fScore);
    }
    return;
}
*/
double CPattern::Classify(CModel* pModel)
{ 
   return 1.0;
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


int CStructureLearning::Init(const char* datafile, const char* patternfile, const char* szpath, const char* modelfile, bool bBinaryData)
{
    m_Sample.m_bBinaryData = bBinaryData;
    strcpy(m_Sample.m_szFolder, szpath);
    
    m_Sample.LoadShape(datafile);
    m_Sample.LoadPattern(patternfile);
    m_Model.Read(modelfile);
    m_Model.InitTheta(m_Sample.m_vPatterns.size());

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
    CC.Init(&m_Model, m_Sample.m_vPatterns.size(), m_fEpsilon);
    double fepsilon = 1e9;
    int iIteration = 0;
    double falleps = 0.0;
    do
    {
        m_Sample.SetAllActive();
        m_Sample.m_iLastUpdated = -1;
        fepsilon = 1e9;
        double pre = 1e10;
        int iStep = 0;
        int iNumVio = 0;
        falleps = 0.0;
        do
        {
            iStep ++;
            fprintf(stderr, "%d iteration %d step \n", iIteration, iStep);
            fprintf(stderr, "active patterns %d (out of %d)\n", m_Sample.GetActiveNum(), m_Sample.m_vPatterns.size());
            pre = fepsilon;
            fepsilon =  m_Sample.UpdateConstraint(&CC, &m_Model, true, m_iMinNewConstraint, iNumVio, m_fEpsilon);
            if (fepsilon > 0)
            {
                CC.PrintMathProg("tmp.mod");
                CC.GLPK_lp(&m_Model);
            }
            falleps += fepsilon;
            fprintf(stderr, "\n%d contraints, total %d violation, sum to %f, all eps %f \n", CC.m_vWeights.size(), iNumVio, fepsilon, falleps);
        }while ((m_Sample.GetActiveNum() > 0) &&  (fepsilon > 0) && (fabs(pre - fepsilon) > 0.0001) && (iStep < m_iMaxStep) );
        iIteration ++;
    }while ((iIteration < m_iMaxIteration) && (falleps> 0));
    m_Model.Write(outputfile);

    return 1;

}
