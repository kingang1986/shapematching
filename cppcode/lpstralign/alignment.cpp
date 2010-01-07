#include "sequence.h"
#include "alignment.h"
#include <string.h>
#include "stdio.h"
#include "stdlib.h"
#include "model.h"
/////////////////////////////////////////////////////////
//class CAlignment
/////////////////////////////////////////////////////////
CAlignment::CAlignment()
{
    m_pSS1 = m_pSS2 = NULL;
    m_fScore = 0.0;
    m_iLayerCount = 0;
}

CAlignment::~CAlignment()
{

}
int CAlignment::GetOperNum() { return (int) m_operation.size();}
int CAlignment::GetOper(int iIndex, int& oper, int& seq1, int& seq2, int& pt1, int&pt2, int& layer)
{
    oper = seq1 = seq2 = pt1 = pt2 = -1;
    if (iIndex < 0 || iIndex >= (int)m_operation.size()) 
          return -1;
    oper = m_operation[iIndex];
    seq1 = m_SeqIndex1[iIndex];
    seq2 = m_SeqIndex2[iIndex];
    pt1 = m_PointIndex1[iIndex];
    pt2 = m_PointIndex2[iIndex];
    layer = m_layer[iIndex];
    return iIndex;

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
//   fprintf(stderr, "add alignment %d and %d \n", m_operation.size(), align.m_operation.size());
//    fprintf(stderr, "0: %d %d %d \n", align.m_operation[0], align.m_SeqIndex1[0], align.m_SeqIndex2[0]);
    for (int i = 0; i < (int) align.m_operation.size(); i ++)
    {
        m_operation.push_back(align.m_operation[i]);
        m_SeqIndex1.push_back(align.m_SeqIndex1[i]);
        m_SeqIndex2.push_back(align.m_SeqIndex2[i]);
        m_PointIndex1.push_back(align.m_PointIndex1[i]);
        m_PointIndex2.push_back(align.m_PointIndex2[i]);
        m_layer.push_back(m_iLayerCount);
    }
    m_fScore += align.m_fScore;
    m_iLayerCount ++;
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
         a = m_pSS1->GetSeq(m_SeqIndex1[i])->GetPointValue(m_PointIndex1[i]);
         b = m_pSS2->GetSeq(m_SeqIndex2[i])->GetPointValue(m_PointIndex2[i]);
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
/*
    for (int i = 0; i < (int)m_PointIndex1.size(); i ++)
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
    end1    = m_iEnd1; 
    start1  = m_iStart1; 
    end2 = m_iEnd2; 
    start2  =  m_iStart2; 
    if ( start1 > end1) { int t = end1; end1 = start1; start1 = t;}
    if ( start2 > end2) { int t = end2; end2 = start2; start2 = t;}
//    fprintf(stderr, "\n(bound %d %d %d %d)\n", start1, end1, start2, end2);
    return; 
}

void CAlignment::Clean()
{
    m_pSS1 = m_pSS2 = NULL;
    m_SeqIndex1.clear(); 
    m_SeqIndex2.clear(); 
    m_PointIndex1.clear();
    m_PointIndex2.clear();
    m_operation.clear();
    m_layer.clear();
    m_iLayerCount = 0;
    m_fScore = 0.0;
    m_iStart1 = m_iEnd1 = m_iStart2 = m_iEnd2 = 0;
    m_nLength1 = m_nLength2 = 0;
     
}

void CAlignment::Print()
{
    fprintf(stderr, " %d opers, %.4f scores\n", GetOperNum(), m_fScore);
    for (int i = 0; i < (int)m_SeqIndex1.size(); i ++)
    {
        //fprintf(stderr, "op %d: %d %d %d %d %d %d\n", i, m_SeqIndex1[i], m_SeqIndex2[i], m_PointIndex1[i], m_PointIndex2[i], m_operation[i], m_layer[i]);
        fprintf(stderr, "op %d: %d %d %d %d %d\n", i, m_SeqIndex1[i], m_SeqIndex2[i], m_PointIndex1[i], m_PointIndex2[i], m_operation[i]);
    }
}
