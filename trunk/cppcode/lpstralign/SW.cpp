#include "SW.h"
#include "math.h"
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>

using namespace std;

double CSWMatch::getSubstituteCost(DATATYPE* a, DATATYPE* b, CModel* model)
{
    double cost = 0;
    for (int i = 0; i < model->m_iMatchCount; i++)
    {
        int t = model->m_vMatch[i];
        int idx = model->m_vFeatureIndex[t];
        int wid = model->m_vWeightIndex[t];
        if (model->m_vMapType[t] == MAPTYPE_CONSTANT)
        {
            cost += model->m_vWeight[wid];
        }
        else if (model->m_vMapType[t] == MAPTYPE_UNCHANGE)
        {
            DATATYPE fa = a[idx];
            DATATYPE fb = b[idx];
            if (fabs(fa) + fabs(fb) > 0)
            {
                DATATYPE d = fa - fb;
                cost += model->m_vWeight[wid] * d * d / (fabs(fa) + fabs(fb));
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
                cost += model->m_vWeight[wid] * d * d / (fa + fb);
                //fprintf(stderr, "/%f, %f/ ", model->m_vWeight[t],  d * d / (fabs(a[i]) + fabs(b[i])));
            }
        }
        else if (model->m_vMapType[t] == MAPTYPE_EU)
        {
            DATATYPE fa = a[idx];
            DATATYPE fb = b[idx];
            cost += model->m_vWeight[wid] * fabs(fa - fb);
        }
        else if (model->m_vMapType[t] == MAPTYPE_ABSEU)
        {
            DATATYPE fa = fabs(a[idx]);
            DATATYPE fb = fabs(b[idx]);
            cost += model->m_vWeight[wid] * fabs(fa - fb);
        }
    }
    //double cost = model->m_vWeight[model->m_vMatchMap[model->m_iFeatureDim]] + f;
    //if (cost > 0)
    //fprintf(stderr, "(%f, %f, %d)", cost, model->m_vWeight[model->m_vMatchMap[model->m_iFeatureDim]], model->m_vMatchMap[model->m_iFeatureDim]);
    //fprintf(stderr, "(%f)", cost);
    return cost;
}

double CSWMatch::getGapCost(DATATYPE* a, CModel* model)
{
    //fprintf(stderr, "cal gap cost\n");
    double cost = 0;
    for (int i = 0; i < model->m_iGapCount; i++)
    {

        int t = model->m_vGap[i];
        //fprintf(stderr, "t %d\n", t);
        int idx = model->m_vFeatureIndex[t];
        int wid = model->m_vWeightIndex[t];

        if (model->m_vMapType[t] == MAPTYPE_CONSTANT)
        {
            cost += model->m_vWeight[wid];
        }
        else if (model->m_vMapType[t] == MAPTYPE_UNCHANGE)
        {
            //fprintf(stderr, "gap cost: %d %d %d\n", t, idx, wid);
            DATATYPE fa = a[idx];
            cost += model->m_vWeight[wid] * fabs(fa);
        }
    }
    //double cost = model->m_vWeight[model->m_vGapMap[model->m_iFeatureDim]] + f;
    //fprintf(stderr, "(%f, %f)", cost, model->m_vWeight[model->m_vGapMap[model->m_iFeatureDim]]);
    return cost;
    // return f; //debug
}

double CSWMatch::MatchSequence(CSequence* pSeqA, CSequence* pSeqB, CAlignment* pAlign, CModel* model)
{

    //allocate nodes
//    fprintf(stderr, "oper (%d,  %d), (%d, %d)\n",  pSeqA->m_iOriginalSeqId, pSeqA->m_iStartPos, pSeqB->m_iOriginalSeqId, pSeqB->m_iStartPos);

    CSWNode ** S;
    //fprintf(stderr, "ttok-1\n");
    int iLengthA = pSeqA->m_iPoint;
    int iLengthB = pSeqB->m_iPoint;
    S = new CSWNode*[iLengthA + 1];
    //fprintf(stderr, "ttok0 length %d %d\n", iLengthA, iLengthB);
    for (int i = 0; i < iLengthA + 1; i++)
    {
        S[i] = new CSWNode[iLengthB + 1];
    }
    for (int i = 0; i < iLengthA + 1; i++)
    {
        S[i][0].score = 0;
        S[i][0].prev = NULL;
        S[i][0].row_index = i;
        S[i][0].column_index = 0;
    }

    for (int i = 0; i < iLengthB + 1; i++)
    {
        S[0][i].score = 0;
        S[0][i].prev = NULL;
        S[0][i].row_index = 0;
        S[0][i].column_index = i;
    }
    //fprintf(stderr, "ttok1\n");
    int maxA = 0, maxB = 0;
    CSWNode* endSequence = &S[0][0];
    DATATYPE** A = pSeqA->m_vFeature;
    DATATYPE** B = pSeqB->m_vFeature;
    //fprintf(stderr, "ttok2\n");
    for (int i = 1; i < iLengthA + 1; i++)
    {
        for (int j = 1; j < iLengthB + 1; j++)
        {
            S[i][j].row_index = i;
            S[i][j].column_index = j;
            S[i][j].score = 0;
            S[i][j].prev = NULL;
            double subst_ij1 = getSubstituteCost(A[i - 1], B[j - 1], model) + S[i - 1][j - 1].score;
            //if (subst_ij1 > 0) fprintf(stderr, "sub cost %f\n", subst_ij1);
            if (subst_ij1 > S[i][j].score)
            {
                S[i][j].score = subst_ij1;
                S[i][j].prev = &S[i - 1][j - 1];
                S[i][j].operation = SUBST;
            }
            double subst_ij2 = getGapCost(A[i - 1], model) + S[i - 1][j].score;
            if (subst_ij2 > S[i][j].score)
            {
                S[i][j].score = subst_ij2;
                S[i][j].prev = &S[i - 1][j];
                S[i][j].operation = DELET;
            }
            double subst_ij3 = getGapCost(B[j - 1], model) + S[i][j - 1].score;
            //fprintf(stderr, "gap cost %f, %f\n", subst_ij2, subst_ij3);
            if (subst_ij3 > S[i][j].score)
            {
                S[i][j].score = subst_ij3;
                S[i][j].prev = &S[i][j - 1];
                S[i][j].operation = INSRT;
            }

            //record information if new maximum
            if (S[i][j].score > endSequence->score)
            {
                endSequence = &S[i][j];
                maxA = i;
                maxB = j;
            }
            //fprintf(stderr, "(%d, %d): %4.2f ", i, j, S[i][j].score);
        }
    }
    //fprintf(stderr, "ttok3\n");
    ///////////////////////////////////////////
    // Retrieve alignment
    ///////////////////////////////////////////
    /*get length of alignment*/
    CSWNode* temp_node = endSequence;
    pAlign->m_fScore = temp_node->score;
    int oper_count = 0;
    while (temp_node->prev)
    {
        pAlign->m_SeqIndex1.push_back(pSeqA->m_iOriginalSeqId);
        pAlign->m_SeqIndex2.push_back(pSeqB->m_iOriginalSeqId);
        pAlign->m_PointIndex1.push_back(temp_node->prev->row_index + pSeqA->m_iStartPos);
        pAlign->m_PointIndex2.push_back(temp_node->prev->column_index + pSeqB->m_iStartPos);
        pAlign->m_operation.push_back(temp_node->operation);
        oper_count++;
        temp_node = temp_node->prev;
    }
    //fprintf(stderr, "ttok4\n");

    /*free DP matrices */
    for (int i = 0; i < iLengthA + 1; i++)
    {
        free(S[i]);
    }
    free(S);
    return pAlign->m_fScore;

}
void CSWMatch::ReleaseMatchingTable()
{
    std::map<pair<int, int>, CAlignment*>::iterator itr; 
    for(itr = m_mapCachedMatching.begin(); itr != m_mapCachedMatching.end(); ++itr)
    {
         CAlignment* pAlign =(CAlignment*) (*itr).second;
         delete pAlign;
    }

}

void CSWMatch::RemoveTable(int shapei, int shapej)
{
        std::map<pair<int, int>, CAlignment*>::iterator itr; 
        std::map<pair<int, int>, CAlignment*>::iterator maxitr; 
        for(itr = m_mapCachedMatching.begin(); itr != m_mapCachedMatching.end(); ++itr)
        {
             
             int p1 = (*itr).first.first; 
             int p2 = (*itr).first.second;
             CAlignment* pAlign =(CAlignment*) (*itr).second;
             if (p1 == shapei || p2 == shapej)
             {
                 pAlign->m_fScore = -1; 
             }
        }
}

void CSWMatch::UpdateMatching()
{
    for (int i = 0; i < m_pSS1->m_iSeqNum; i ++)
    {
        CSequence* pS1 = m_pSS1->m_vSeqs[i];
        for (int j = 0; j < m_pSS2->m_iSeqNum; j ++)  
        {
            CSequence* pS2 = m_pSS2->m_vSeqs[j];
            CAlignment* pAlign = new CAlignment();
            pAlign->m_pSS1 = m_pOSS1;
            pAlign->m_pSS2 = m_pOSS2;
            pair<int, int> p ;
            p.first = pS1->m_iID;
            p.second = pS2->m_iID;
//            fprintf(stderr, "(%d, %d) ", p.first, p.second);
            if (m_mapCachedMatching.find(p) == m_mapCachedMatching.end())
            {
                MatchSequence(pS1, pS2, pAlign, m_model);
                m_mapCachedMatching[p] = pAlign;
            }
        }
    }  
    return; 
}

double CSWMatch::Match(CSetOfSeq* pSS1, CSetOfSeq* pSS2, CAlignment* pASet, CModel* model)
{
    m_pOSS1 = pSS1;
    m_pOSS2 = pSS2;
    m_pSS1 = new CSetOfSeq(*pSS1);
    m_pSS2 = new CSetOfSeq(*pSS2);
    m_model = model;

    m_pSS1->RemoveShortSeqs(3);
    m_pSS2->RemoveShortSeqs(3);

    UpdateMatching();
    double fTotalScore = 0;
    pASet->m_pSS1 = pSS1;
    pASet->m_pSS2 = pSS2;
    while (m_pSS1->m_iSeqNum > 0 && m_pSS2->m_iSeqNum > 0)
    {
//        fprintf(stderr, "\nmatching (%d seq  %d seq )", m_pSS1->m_iSeqNum, m_pSS2->m_iSeqNum);
 //      for (int kk = 0; kk < m_pSS1->m_iSeqNum; kk ++)
  //            fprintf(stderr, " (seq %d, id %d : %d ) ", kk, m_pSS1->m_vSeqs[kk]->m_iID, m_pSS1->m_vSeqs[kk]->m_iPoint);
   //    for (int kk = 0; kk < m_pSS2->m_iSeqNum; kk ++)
    //       fprintf(stderr, " (seq %d, id %d : %d ) ", kk, m_pSS2->m_vSeqs[kk]->m_iID, m_pSS2->m_vSeqs[kk]->m_iPoint);
 
        double fMaxScore = -1;
        int bestP1 = 0, bestP2 = 0;
        CAlignment* pBestAlign = NULL;
        std::map<pair<int, int>, CAlignment*>::iterator itr; 
        std::map<pair<int, int>, CAlignment*>::iterator maxitr; 
        for(itr = m_mapCachedMatching.begin(); itr != m_mapCachedMatching.end(); ++itr)
        {
             CAlignment* pAlign =(CAlignment*) (*itr).second;
             if (pAlign->m_fScore > fMaxScore) 
             {
                      
                 maxitr = itr;
                 fMaxScore = pAlign->m_fScore;
                 bestP1 = (*itr).first.first; 
                 bestP2 = (*itr).first.second;
     //           fprintf(stderr, "<< %d %d %.4f\n", bestP1, bestP2, fMaxScore);
                 pBestAlign = pAlign;
             }
        }
        pASet->AddAlignment(*pBestAlign);
        fTotalScore += fMaxScore;
        if (fMaxScore <= 0.000001) break;
       
        int start1, end1, start2, end2; 
        pBestAlign->GetBound(start1, end1, start2, end2);
        m_pSS1->SplitSeqByID(bestP1, start1, end1); 
        m_pSS2->SplitSeqByID(bestP2, start2, end2); 
        RemoveTable(bestP1, bestP2);
      //  fprintf(stderr, "(%d %d %d %d %f %f)\n ", m_pSS1->m_iSeqNum,m_pSS2->m_iSeqNum, bestP1, bestP2,fMaxScore, fTotalScore);
        m_pSS1->RemoveShortSeqs(3);
        m_pSS2->RemoveShortSeqs(3);
       // fprintf(stderr, "(%f %d %d  %d %d )\n", fMaxScore, start1, end1, start2, end2); 
        UpdateMatching();
    }
    ReleaseMatchingTable();
    delete m_pSS1;
    delete m_pSS2;
    return fTotalScore; 
}
