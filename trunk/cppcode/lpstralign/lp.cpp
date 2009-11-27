#include <stdio.h>
#include <vector>
#include <algorithm>
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
    //  fprintf(stderr, "load: %s from %s in folder %s...", strSoSFile, tmp, m_szFolder);
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

    // fprintf(stderr, "\r %d Shapes (%d uniq) loaded, in %d patterns ", m_nShapes, (int) m_vSS.size(), m_nPattern + 1);

    return pSS;
}


int CSample::LoadSample(const char* strFile)
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
        pPattern->m_original = LoadSoS(tmpbuffer);
        pbuffer = getnextstring(pbuffer, tmpbuffer, ',');
        pPattern->m_iHomolog = atoi(tmpbuffer);
        //fprintf(stderr, "%d homologs\n", pPattern->m_iHomolog);
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


int CSample::AlignHomolog(CModel* pModel)
{
    for (int i = 0; i < m_nPattern; i ++)
    {
        CPattern* pPattern = m_vPatterns[i];
        fprintf(stderr, "\nPattern %d\n", i);
        pPattern->AlignHomolog(pModel);
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
    double fSumAcc = 0.0;
    for (int i = 0; i < m_nPattern; i ++)
    {
        CPattern* pPattern = m_vPatterns[i];
        fprintf(stderr, "\n%d: ", i + 1);
        fSumAcc += pPattern->Align(pModel);
    }
    fprintf(stderr, "sample # %d, avg accuracy %f \n", m_nPattern, fSumAcc/m_nPattern);
    fclose(fp);
    return m_nPattern;

}


//add constraints upto iMinNewConstraint, starting from last updated pattern

double CSample::UpdateConstraint(CConstraints* pCC, CModel* pModel, bool bActiveOnly, int iMinNewConstraint, int& iNumVio)
{
    int iNewConstraint = 0;
    double fsumvio = 0.0;
    iNumVio = 0;
    double vio = 0.0;
    int lastupdate = m_iLastUpdated;
    for (int i = lastupdate + 1; i < m_nPattern + lastupdate + 1 && iNewConstraint < iMinNewConstraint - 1; i ++)
    {
        int ii  = i;
        while (ii >= m_nPattern) ii -= m_nPattern;
        if (bActiveOnly && m_vActive[ii] == 0)
            continue;
        CPattern* pPattern = m_vPatterns[ii];
        fprintf(stderr, "\n %d", ii);
        pPattern->Align(pModel);
        int iParamDim = pModel->m_iParamDim;
        double* pw1 = new double[iParamDim];
        double* pw2 = new double[iParamDim];
        double labelloss = pPattern->GetLabeledPhi(pw1, iParamDim, pModel);
        fprintf(stderr, " label loss %.3f ", labelloss);
        pPattern->GetDecoyPhi(pw2, iParamDim, pModel);
        double fsum = 0, fs2 = 0;
        for (int k = 0; k < iParamDim; k ++)
        {
            fsum += pw1[k] * pModel->m_vWeight[k];
            fs2 += pw2[k] * pModel->m_vWeight[k];
        }
        vio = fsum - fs2 - labelloss;
        //if (vio > -0.00001)
        if (labelloss < 0.00001)
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
            fprintf(stderr, " ++ %.3f, %.3f %.3f ", fsum, fs2, -vio);
              
        }
        m_iLastUpdated = ii;
    }
    return -fsumvio;
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
CPattern::CPattern()
{
}


CPattern::~CPattern()
{
}


bool cmpAlign(const CAlignment* a1, const CAlignment* a2)
{
    return a1->m_fScore > a2->m_fScore;
}


//align all samples accoding to the pModel
double CPattern::Align(CModel* pModel)
{
    m_vAligns.clear();
    for (int j = 0; j < m_iHomolog; j ++)
    {
        CSWMatch * pMatch = new CSWMatch();
        CSetOfSeq* pSS = m_vHomolog[j];
        CAlignment align;
        pMatch->Match(m_original, pSS, &align, pModel);
        align.m_bSameClass = true;
        m_vAligns.push_back(align);

    }
    for (int j = 0; j < m_iDecoy; j ++)
    {
        CSetOfSeq* pSS = m_vDecoy[j];
        CAlignment align;
        CSWMatch* pMatch = new CSWMatch();
        pMatch->Match(m_original,  pSS, &align, pModel);
        align.m_fScore += 0.00001;
        align.m_bSameClass = false;
        m_vAligns.push_back(align);
    }
    m_vSortedAlign.clear();
    for (int i = 0; i < (int)m_vAligns.size(); i ++)
        m_vSortedAlign.push_back(&m_vAligns[i]);
    sort(m_vSortedAlign.begin(), m_vSortedAlign.end(), cmpAlign);
    fprintf(stderr, " top %d results : ", CPattern::m_iTopK);
    double iCorrect = 0;
    for (int i = 0; i < m_iTopK; i ++)
    {
        //fprintf(stderr, " %.3f", m_vSortedAlign[i]->m_fScore);
        if (m_vSortedAlign[i]->m_bSameClass)
        {
            fprintf(stderr, "*");
            iCorrect ++;
        }
    }

    return iCorrect / m_iTopK;
}


void CPattern::AlignHomolog(CModel* pModel)
{
    m_vLabelHomoAlign.clear();
    for (int j = 0; j < m_iHomolog; j ++)
    {
        CSetOfSeq* pSS = m_vHomolog[j];
        CSWMatch* pMatch = new CSWMatch();
        CAlignment align;
        pMatch->Match(m_original, pSS, &align, pModel);
        fprintf(stderr, "(%d, %f) ", j, align.m_fScore);
        m_vLabelHomoAlign.push_back(align);
    }
    
    m_vSortedLabel.clear();
    for (int i = 0; i < (int)m_vLabelHomoAlign.size(); i ++)
        m_vSortedLabel.push_back(&m_vLabelHomoAlign[i]);
    sort(m_vSortedLabel.begin(), m_vSortedLabel.end(), cmpAlign);
    return ;
}


double CPattern::GetLabelLoss()
{
    //sort the aligns;
    double iError = 0.0;
    for (int i = 0; i < m_iTopK; i ++)
    {
        if (!m_vSortedAlign[i]->m_bSameClass)
        {
            iError += m_iTopK - i;
        }
    }
    return iError /  (m_iTopK * m_iTopK);

}


void CPattern::GetDecoyPhi(double* pw, int iParamDim, CModel* pmodel)
{
    memset(pw, 0, sizeof(double) * iParamDim);
    int iIdx = 0;
    while (m_vSortedAlign[iIdx]->m_bSameClass) 
          iIdx ++;
    m_vSortedAlign[iIdx]->GetPhi(pw, iParamDim, pmodel);
    return;
}


double CPattern::GetLabeledPhi(double* pw, int iParamDim, CModel* pmodel)
{
    // find the first homolog that are not ranked in top k, return its phi
    for (int i = 0; i < m_iTopK; i ++) 
    {
       CSetOfSeq* pSS = m_vSortedLabel[i]->m_pSS2; 
       bool bFound = false;
       for (int j = 0; j < m_iTopK; j ++)
       {
            if (pSS == m_vSortedAlign[j]->m_pSS2)
            { 
                 bFound = true;
                 break;
            } 
       }
       if (! bFound) 
       {
          m_vSortedLabel[i]->GetPhi(pw, iParamDim, pmodel); 
          return (m_iTopK - i) / (double)m_iTopK;
       }
    } 
    return 0.0;
}


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


int CStructureLearning::Init(const char* datafile, const char* szpath, const char* modelfile, bool bBinaryData)
{
    m_Sample.m_bBinaryData = bBinaryData;
    strcpy(m_Sample.m_szFolder, szpath);
    m_Sample.LoadSample(datafile);
    m_Model.Read(modelfile);
    //    m_Model.Print();
    m_Model.InitTheta(m_Sample.m_nPattern);

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
    CC.Init(&m_Model, m_Sample.m_nPattern, m_fEpsilon);
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
            fprintf(stderr, "active patterns %d (out of %d)\n", m_Sample.GetActiveNum(), m_Sample.m_nPattern);
            pre = fepsilon;
            fepsilon =  m_Sample.UpdateConstraint(&CC, &m_Model, true, m_iMinNewConstraint, iNumVio);
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
