#include "Data.h"
#include <cstdlib>
#include <string.h>
#include <stdio.h>
#include "sequence.h"


/////////////////////////////////////////////////////////
//class CAlignment
/////////////////////////////////////////////////////////
CAlignment::CAlignment()
{
    m_pSS1 = m_pSS2 = NULL;
    m_fScore = 0.0;
}

CAlignment::~CAlignment()
{

}

double CAlignment::AddAlignment(CAlignment& align)
{
    if ((m_pSS1 != NULL && m_pSS1 != align.m_pSS1) || (m_pSS2 != NULL && m_pSS2 != align.m_pSS2))
    {
        fprintf(stderr, "can't merge two alignment because the shapes are different\n");
        return -1;
    }
    m_pSS1 = align.m_pSS1;
    m_pSS2 = align.m_pSS2;
//    fprintf(stderr, "add alignment %d and %d \n", m_operation.size(), align.m_operation.size());
    for (int i = 0; i < (int) align.m_operation.size(); i ++)
    {
        m_operation.push_back(align.m_operation[i]);
        m_SeqIndex1.push_back(align.m_SeqIndex1[i]);
        m_SeqIndex2.push_back(align.m_SeqIndex2[i]);
        m_PointIndex1.push_back(align.m_PointIndex1[i]);
        m_PointIndex2.push_back(align.m_PointIndex2[i]);
    }
    m_fScore += align.m_fScore;
    return m_fScore;

}

double CAlignment::GetPhi(double* phi, int iParamDim, CModel* model)
{
    memset(phi, 0, sizeof(double) * iParamDim);

    int icount = 0;
    for (int i = 0; i < (int) m_operation.size(); i++)
    {
      DATATYPE* a;
      DATATYPE* b;

      try{
         a = m_pSS1->m_vSeqs[m_SeqIndex1[i]]->GetPoint(m_PointIndex1[i]);
         b = m_pSS2->m_vSeqs[m_SeqIndex2[i]]->GetPoint(m_PointIndex2[i]);
      }catch(...)
      {
          fprintf(stderr, "%d %d",i,  m_SeqIndex2[i]);
          exit(1); 
       }
        if (m_operation[i] == SUBST)
        {
            icount++;
            for (int j = 0; j < model->m_iMatchCount; j++)
            {

                int t = model->m_vMatch[j];
                int idx = model->m_vFeatureIndex[t];
                int wid = model->m_vWeightIndex[t];
                //fprintf(stderr, "match, %d, %d %d \n", t, idx, wid);
                if (model->m_vMapType[t] == MAPTYPE_CONSTANT)
                {
                    phi[wid] += 1;
                }
                else if (model->m_vMapType[t] == MAPTYPE_UNCHANGE)
                {
                    DATATYPE fa = a[idx];
                    DATATYPE fb = b[idx];
                    if (fabs(fa) + fabs(fb) > 0)
                    {
                        DATATYPE d = fa - fb;
                        DATATYPE f = d * d / (fabs(fa) + fabs(fb));
                        phi[wid] += f;
                        //fprintf(stderr, "/%f, %f/ ", model->m_vWeight[t],  d * d / (fabs(a[i]) + fabs(b[i])));
                    }
                }
                else if (model->m_vMapType[t] == MAPTYPE_ABS)
                {
                    DATATYPE fa = fabs(a[idx]);
                    DATATYPE fb = fabs(b[idx]);
                    if (fa + fb > 0)
                    {
                        DATATYPE d = fa - fb;
                        DATATYPE f = d * d / (fa + fb);
                        phi[wid] += f;
                        //fprintf(stderr, "/%f, %f/ ", model->m_vWeight[t],  d * d / (fabs(a[i]) + fabs(b[i])));
                    }
                }
                else if (model->m_vMapType[t] == MAPTYPE_EU)
                {
                    DATATYPE fa = a[idx];
                    DATATYPE fb = b[idx];
                    phi[wid] += fabs(fa - fb);
                }
                else if (model->m_vMapType[t] == MAPTYPE_ABSEU)
                {
                    DATATYPE fa = fabs(a[idx]);
                    DATATYPE fb = fabs(b[idx]);
                    phi[wid] += fabs(fa - fb);
                }

            }
        }
        else if (m_operation[i] == DELET || m_operation[i] == INSRT)
        {
            icount++;
            for (int j = 0; j < model->m_iGapCount; j++)
            {
                int t = model->m_vGap[j];
                int idx = model->m_vFeatureIndex[t];
                int wid = model->m_vWeightIndex[t];
                //fprintf(stderr, "gap, %d, %d %d \n", t, idx, wid);
                if (model->m_vMapType[t] == MAPTYPE_CONSTANT)
                {
                    phi[wid] += 1;
                }
                else if (model->m_vMapType[t] == MAPTYPE_UNCHANGE)
                {
                    DATATYPE fa = a[idx];
                    DATATYPE fb = b[idx];
                    if (m_operation[i] == DELET)
                        phi[wid] += fabs(fa);
                    else
                        phi[wid] += fabs(fb);
                }
            }
        }
        else
        {
            fprintf(stderr, "error operation code %d at %d, %d %d %d %d\n",
                    m_operation[i], icount, m_SeqIndex1[i],
                    m_SeqIndex2[i], m_PointIndex1[i],
                    m_PointIndex2[i]);
        }
    }
    double fsum = 0;
    for (int i = 0; i < iParamDim; i++)
    {
        fsum += model->m_vWeight[i] * phi[i];
    }
    //fprintf(stderr, "%d %d\n", icount, (int)m_operation.size());
    return fsum;
}

void CAlignment::GetBound(int& start1, int& end1, int& start2, int& end2)
{
/*    for (int i = 0; i < (int)m_PointIndex1.size(); i ++)
    {
        printf("%d ", m_PointIndex1[i]);
    }
    printf("\n");

    for (int i = 0; i < (int)m_PointIndex2.size(); i ++)
    {
        printf("%d ", m_PointIndex2[i]);
    }
    printf("\n");
*/
    end1    = m_PointIndex1[0];
    start1  = m_PointIndex1[m_PointIndex1.size() - 1];
    end2 = m_PointIndex2[0];
    start2  =  m_PointIndex2[m_PointIndex2.size() - 1]; 
}

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

    Release();

    int itotaline = count_line(strFile);
    char line[32768];
    //fprintf(stderr, "Reading struct model ...");
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
        //fprintf(stderr, "%f %d %d %d %d",idx, wid, maptype, matchgap);
        if (wid > maxWeightIndex) maxWeightIndex = wid;
        m_vFeatureIndex[i] = idx;
        m_vWeightIndex[i] = wid;
        m_vMapType[i] = maptype;
        m_vMatchOrGap[i] = matchgap;
        if (matchgap == 101) matchcount ++;
        if (matchgap == 102) gapcount ++;
    }

    m_iParamDim = maxWeightIndex + 1;
    fprintf(stderr, "param dim = %d\n",m_iParamDim);
    m_iMapNum = itotaline - 3;
    m_vWeight = new double[m_iParamDim];
    for (int i = 0; i < m_iParamDim; i ++)
    {
        float w;
        fscanf(fp, "%f", &w);
        fprintf(stderr, "%f\t", w);
        m_vWeight[i] = w;
    }
    fprintf(stderr, "\nsign of each weight\n");
    for (int i = 0; i < m_iParamDim; i ++)
    {
        int w;
        fscanf(fp, "%d", &w);
        fprintf(stderr, "%d\t", w);
        m_vSign[i] = w;
    }
    fclose(fp);
    fprintf(stderr, "\ndone\n");

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
    m_vMatch = m_vGap = NULL;
    m_vFeatureIndex = NULL;
    m_vWeightIndex = NULL;
    m_vMatchOrGap = NULL;
    m_vSign  = NULL;
    m_vTheta = NULL;
    return ;
}

void CModel::InitTheta(int iPatternNum)
{
    m_vTheta = new double[iPatternNum];
    m_iPatternNum = iPatternNum;
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

