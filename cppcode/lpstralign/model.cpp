#include "model.h"
#include "stdio.h"
#include "stdlib.h"

/////////////////////////////////////////////////////////////
//class CModel
/////////////////////////////////////////////////////////////

CModel::CModel()
{
    m_vWeight = NULL;
    m_vFeatureIndex = NULL;
    m_vWeightIndex = NULL;
    m_vMapType = NULL;
    m_vMatchOrGap = NULL;
    m_vMatch = NULL;
    m_vGap = NULL;
    m_vSign = NULL;
    m_vTheta = NULL;

}

CModel::~CModel()
{
    Release();
}

void CModel::Print()
{
    /* Writes structural model sm to file file. */
    fprintf(stderr, "\nweight\tfeaidx\tw_idx\tmaptype\tm_or_g\n");
    for (int i = 0; i < m_iMapNum; i ++)
    {
        fprintf(stderr, "%d\t%d\t%d\t%d\n",   m_vFeatureIndex[i], m_vWeightIndex[i],m_vMapType[i], m_vMatchOrGap[i]);
    }
    for (int i = 0; i < m_iParamDim; i ++)
    {
        fprintf(stderr, "%f ", m_vWeight[i]);
    }
    fprintf(stderr, "\n");
    for (int i = 0; i < m_iParamDim; i ++)
    {
        fprintf(stderr, "%d ", m_vSign[i]);
    }

    fprintf(stderr, "\nmax weight index %d, mapping # %d\n", m_iParamDim, m_iMapNum);
}


int CModel::Write(const char* strFile)
{
    /* Writes structural model sm to file file. */
    FILE* fp = fopen(strFile, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "Can't open file %s for writing model.\n", strFile);
        exit(1);
    }
    fprintf(fp, "feaidx\tw_idx\tmaptype\tm_or_g\n");
    for (int i = 0; i < m_iMapNum; i ++)
    {
        fprintf(fp, "%d\t%d\t%d\t%d\n",  m_vFeatureIndex[i], m_vWeightIndex[i],m_vMapType[i], m_vMatchOrGap[i]);
    }
    for (int i = 0; i < m_iParamDim; i ++)
    {
        fprintf(fp, "%f ", m_vWeight[i]);
    }
    fprintf(fp, "\n");
    for (int i = 0; i < m_iParamDim; i ++)
    {
        fprintf(fp, "%d ", m_vSign[i]);
    }
    fprintf(fp, "\n");

    fclose(fp);
    return m_iParamDim;
}

int CModel::Read(const char* strFile)
{

    fprintf(stderr, "Reading struct model from %s ...", strFile);
    Release();

    int itotaline = count_line(strFile);
    char line[32768];
    FILE* fp = fopen(strFile, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "can't open model file %s\n", strFile);
        exit(0);
    }
    fgets(line, 32768, fp); //head line
    Init(itotaline - 1);
    int matchcount = 0;
    int gapcount = 0;
    int maxWeightIndex = -1;
    for (int i = 0; i < itotaline - 3; i ++)
    {
        int idx, maptype, matchgap,  wid;
        fscanf(fp, "%d %d %d %d", &idx, &wid, &maptype, &matchgap);
   //     fprintf(stderr, "%d %d %d %d",idx, wid, maptype, matchgap);
        if (wid > maxWeightIndex) maxWeightIndex = wid;
        m_vFeatureIndex[i] = idx;
        m_vWeightIndex[i] = wid;
        m_vMapType[i] = maptype;
        m_vMatchOrGap[i] = matchgap;
        if (matchgap == 101) matchcount ++;
        if (matchgap == 102) gapcount ++;
    }

    m_iParamDim = maxWeightIndex + 1;
    //fprintf(stderr, "param dim = %d\n",m_iParamDim);
    m_iMapNum = itotaline - 3;
    m_vWeight = new double[m_iParamDim];
    for (int i = 0; i < m_iParamDim; i ++)
    {
        float w;
        fscanf(fp, "%f", &w);
        //fprintf(stderr, "%f\t", w);
        m_vWeight[i] = w;
    }
    //fprintf(stderr, "\nsign of each weight\n");
    for (int i = 0; i < m_iParamDim; i ++)
    {
        int w;
        fscanf(fp, "%d", &w);
     //   fprintf(stderr, "%d\t", w);
        m_vSign[i] = w;
    }
    fclose(fp);
    fprintf(stderr, "done\n");

    //build
    m_vMatch = new int[matchcount];
    m_vGap = new int[gapcount];
    int idx1 = 0;
    int idx2 = 0;
    for (int i = 0; i < itotaline - 2; i ++)
    {
        if (m_vMatchOrGap[i] == MAPTYPE_MATCH)
        {
            m_vMatch[idx1] = i; idx1 ++ ;
        }
        else if (m_vMatchOrGap[i] == MAPTYPE_GAP)
        {
            m_vGap[idx2] = i; idx2 ++;
        }
    }
    m_iMatchCount = matchcount;
    m_iGapCount = gapcount;
    return m_iParamDim; 
} 
void CModel::Release()
{
    if (m_vWeight != NULL)
    {
        delete [] m_vWeight;
    }
    if (m_vFeatureIndex != NULL)
    {
        delete [] m_vFeatureIndex;
    }
    if (m_vWeightIndex != NULL)
    {
        delete [] m_vWeightIndex;
    }
    if (m_vMatchOrGap != NULL)
    {
        delete [] m_vMatchOrGap;
    }
    if (m_vTheta != NULL)
    {
        delete [] m_vTheta;
    }
    if (m_vMatch!= NULL)
    {
        delete [] m_vMatch;
    }
    if (m_vGap!= NULL)
    {
        delete [] m_vGap;
    }
    if (m_vSign != NULL)
    {
        delete [] m_vSign;
    }
    m_vWeight = NULL;
    m_vFeatureIndex = NULL;
    m_vWeightIndex = NULL;
    m_vMatchOrGap = NULL;
    m_vTheta = NULL;
    m_vMatch  = NULL;
    m_vGap = NULL;
    m_vSign  = NULL;
    m_vGap = NULL;
    m_vSign = NULL;
    
//    fprintf(stderr, "Released\n");
    return ;
}

void CModel::InitTheta(int iPatternNum)
{
    m_vTheta = new double[iPatternNum];
    m_iPatternNum = iPatternNum;
}

void CModel::Default(int iFeatureNum)
{
    m_iFeatureDim = iFeatureNum; 
    m_iParamDim = iFeatureNum + 2;// gap,  b - fx;
    Init(m_iParamDim);
    for (int i = 0;i < iFeatureNum; i ++)
    {
        m_vWeight[i] = -1.0;
        m_vMapType[i] = MAPTYPE_UNCHANGE;
        //m_vMapType[i] = MAPTYPE_ABS;
        //m_vMapType[i] = MAPTYPE_EU;
        m_vFeatureIndex[i] = i; 
        m_vWeightIndex[i] = i; 
        m_vSign[i] = -1; 
        m_vMatchOrGap[i] = MAPTYPE_MATCH; 
    }
    m_vWeight[iFeatureNum] = 0.2 * iFeatureNum;
    m_vMapType[iFeatureNum] = MAPTYPE_CONSTANT;
    m_vFeatureIndex[iFeatureNum] = -1;
    m_vWeightIndex[iFeatureNum] = iFeatureNum;
    m_vSign[iFeatureNum] = 1;
    m_vMatchOrGap[iFeatureNum] = MAPTYPE_MATCH;

    m_vWeight[iFeatureNum + 1] = -1;
    m_vMapType[iFeatureNum + 1] = MAPTYPE_CONSTANT;
    m_vFeatureIndex[iFeatureNum + 1] = -1;
    m_vWeightIndex[iFeatureNum + 1] = iFeatureNum + 1;
    m_vSign[iFeatureNum + 1] = -1;
    m_vMatchOrGap[iFeatureNum + 1] = MAPTYPE_GAP;
    //build
    m_vMatch = new int[iFeatureNum + 1];
    m_vGap = new int[1];
    int idx1 = 0;
    int idx2 = 0;
    for (int i = 0; i < m_iParamDim; i ++)
    {
        if (m_vMatchOrGap[i] == MAPTYPE_MATCH)
        {
            m_vMatch[idx1] = i; idx1 ++ ;
        }
        else if (m_vMatchOrGap[i] == MAPTYPE_GAP)
        {
            m_vGap[idx2] = i; idx2 ++;
        }
    }
    m_iMapNum = iFeatureNum + 2;
    m_iMatchCount = iFeatureNum + 1;
    m_iGapCount = 1;
    return;
}
void CModel::Init(int iParamDim)
{
    Release();
    m_vWeight = new double[iParamDim];
    m_vMapType = new int[iParamDim];
    m_vFeatureIndex = new int[iParamDim];
    m_vWeightIndex = new int[iParamDim];
    m_vSign = new int[iParamDim];
    m_vMatchOrGap = new int[iParamDim];
    m_iParamDim = iParamDim;
    return;
}

