#include "SW.h"
#include "math.h"
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include "shapecontext.h"

using namespace std;

double CSWMatch::getSubstituteCost(DATATYPE* a, DATATYPE* b)
{
    double cost = 0;
    CModel* model = m_model;
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

double CSWMatch::getGapCost(DATATYPE* a)
{
    //fprintf(stderr, "cal gap cost\n");
    double cost = 0;
    CModel* model = m_model;
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

double CSWMatch::MatchSequence(CSequence* pSeqA, CSequence* pSeqB, CAlignment* pAlign)
{
    bool bTwoWay = m_bTwoWay;
    if (bTwoWay)   
    {
        CAlignment* pAl = new CAlignment(); 
        pAl-> m_pSS1 = pAlign->m_pSS1;
        pAl-> m_pSS2 = pAlign->m_pSS2;
        CSequence* pRevB = new CSequence(*pSeqB);
        pRevB->Reverse();
        double f = MatchSequenceOneWay(pSeqA, pRevB, pAl);
        double f1 =  MatchSequenceOneWay(pSeqA, pSeqB, pAlign);
//        fprintf(stderr, "(forward: %f <> reverse: %f)\n", f1, f);
        if (f > f1)
        {
             *pAlign = *pAl; // use default copy
             delete pRevB;
             return f;
        }
        delete pRevB;
        return f1;
        
    }
    else 
    {
        return MatchSequenceOneWay(pSeqA, pSeqB, pAlign);
    }
    return 0.0;
}

//oneway
double CSWMatch::MatchSequenceOneWay(CSequence* pSeqA, CSequence* pSeqB, CAlignment* pAlign)
{
    int iLengthA = pSeqA->GetPointNum();
    int iLengthB = pSeqB->GetPointNum();
    pAlign->m_nLength1 = iLengthA;
    pAlign->m_nLength2 = iLengthB;
    CSWNode** S = new CSWNode*[iLengthA + 1];
    if (!S) fprintf(stderr, "cant allocate memory matrix\n");
    //fprintf(stderr, "length %d %d\n", iLengthA, iLengthB);
    for (int i = 0; i < iLengthA + 1; i++)
    {
        S[i] = new CSWNode[iLengthB + 1];
        if (!S[i]) fprintf(stderr, "cant allocate memeory at %d \n", i );
    }
    for (int i = 0; i < iLengthA + 1; i++)
    {
        S[i][0].score = 0;
        S[i][0].prev = NULL;
        S[i][0].row_index = i;
        S[i][0].column_index = 0;
        S[i][0].operation = -1;
    }

    for (int i = 0; i < iLengthB + 1; i++)
    {
        S[0][i].score = 0;
        S[0][i].prev = NULL;
        S[0][i].row_index = 0;
        S[0][i].column_index = i;
        S[0][i].operation = -1;
    }
    int maxA = 0, maxB = 0;
    CSWNode* endSequence = &S[0][0];
    vector<DATATYPE*>& A = pSeqA->m_vFeature;
    vector<DATATYPE*>& B = pSeqB->m_vFeature;
    for (int i = 1; i < iLengthA + 1; i++)
    {
         
        for (int j = 1; j < iLengthB + 1; j++)
        {
            S[i][j].row_index = i;
            S[i][j].column_index = j;
            S[i][j].score = 0;
            S[i][j].prev = NULL;
            double subst_ij1 = getSubstituteCost(A[i - 1], B[j - 1]) + S[i - 1][j - 1].score;
            if (subst_ij1 > S[i][j].score)
            {
                S[i][j].score = subst_ij1;
                S[i][j].prev = &S[i - 1][j - 1];
                S[i][j].operation = SUBST;
            }
            double subst_ij2 = getGapCost(A[i - 1]) + S[i - 1][j].score;
            if (subst_ij2 > S[i][j].score)
            {
                S[i][j].score = subst_ij2;
                S[i][j].prev = &S[i - 1][j];
                S[i][j].operation = DELET;
            }
            double subst_ij3 = getGapCost(B[j - 1]) + S[i][j - 1].score;
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
        }
    }
    ///////////////////////////////////////////
    // Retrieve alignment
    ///////////////////////////////////////////
    //get length of alignment
    CSWNode* temp_node = endSequence;
    pAlign->m_fScore = temp_node->score;
    int oper_count = 0;
    pAlign->m_iStart1 = temp_node->row_index ; 
    pAlign->m_iStart2 = temp_node->column_index; 
    while (temp_node->prev)
    {
//        fprintf(stderr, "%d (%d, %d): %4.2f %d\n", oper_count, temp_node->row_index, temp_node->column_index, temp_node->score, temp_node->operation);
        if (temp_node->operation<0 || temp_node->operation > 10) fprintf(stderr, "[ERROR] wrong oper %d\n", temp_node->operation);
        pAlign->m_operation.push_back(temp_node->operation);
        if (temp_node->operation == SUBST || temp_node->operation == DELET)
        {
            pAlign->m_SeqIndex1.push_back(pSeqA->m_vPoints[temp_node->row_index-1]->m_iOriginalSeqIdx);
            pAlign->m_PointIndex1.push_back(pSeqA->m_vPoints[temp_node->row_index-1]->m_iOriginalPtIdx); 
        }
        else
        {
            pAlign->m_SeqIndex1.push_back(-1);
            pAlign->m_PointIndex1.push_back(-1);
        }
        
        if (temp_node->operation == SUBST || temp_node->operation == INSRT) 
        {
            pAlign->m_SeqIndex2.push_back(pSeqB->m_vPoints[temp_node->column_index-1]->m_iOriginalSeqIdx);
            pAlign->m_PointIndex2.push_back(pSeqB->m_vPoints[temp_node->column_index-1]->m_iOriginalPtIdx); 
        }
        else
        {
            pAlign->m_SeqIndex2.push_back(-1);
            pAlign->m_PointIndex2.push_back(-1);
        }

        oper_count++;
        pAlign->m_iEnd1 = temp_node->row_index - 1; 
        pAlign->m_iEnd2 = temp_node->column_index - 1; 
        temp_node = temp_node->prev;
  
    }
//    fprintf(stderr, "start_ends %d %d %d %d\n", pAlign->m_iStart1, pAlign->m_iEnd1, pAlign->m_iStart2, pAlign->m_iEnd2);
    /*free DP matrices */
    for (int i = 0; i < iLengthA + 1; i++)
    {
        delete[] S[i];
    }
    delete[] S;
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
    m_mapCachedMatching.clear();

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
    int iCount = 0;
    for (int i = 0; i < m_pSS1->GetSeqNum(); i ++)
    {
        //CSequence* pS1 = m_pSS1->m_vSeqs[i];
        CSequence* pS1 = m_pSS1->GetSeq(i);
        for (int j = 0; j < m_pSS2->GetSeqNum(); j ++)  
        {
            CSequence* pS2 = m_pSS2->GetSeq(j);
            CAlignment* pAlign = new CAlignment();
            pAlign->m_pSS1 = m_pOSS1;
            pAlign->m_pSS2 = m_pOSS2;
            pair<int, int> p ;
            p.first = pS1->m_iID;
            p.second = pS2->m_iID;
            if (m_mapCachedMatching.find(p) == m_mapCachedMatching.end())
            {
                MatchSequence(pS1, pS2, pAlign);
 //               if (pAlign->m_fScore > 100.0)
  //                  fprintf(stderr, "(%d (len %d) %d(len %d) -> %.4f)\n", p.first, pS1->GetPointNum(), p.second,pS2->GetPointNum(), pAlign->m_fScore);
                m_mapCachedMatching[p] = pAlign;
                iCount ++;
            }
        }
    }  
//    fprintf(stderr, "Table updated, %d matching are calcualted \n", iCount);
    return; 
}

void CSWMatch::showinfo(int bestP1, int bestP2, double fMaxScore)
{
   fprintf(stderr, "result: from (%d seq %d seq) %d %d score: %f total score%f)\n", m_pSS1->GetSeqNum(),m_pSS2->GetSeqNum(), bestP1, bestP2,fMaxScore, m_fTotalScore);
    for (int kk = 0; kk < m_pSS1->GetSeqNum(); kk ++)
        if (m_pSS1->GetSeq(kk)->m_iID == bestP1) 
            fprintf(stderr, "(seq %d, id %d : %d ) ", kk, m_pSS1->GetSeq(kk)->m_iID, m_pSS1->GetSeq(kk)->GetPointNum());
    
    fprintf(stderr, " and ");
    for (int kk = 0; kk < m_pSS2->GetSeqNum(); kk ++)
       if (m_pSS2->GetSeq(kk)->m_iID == bestP2) 
        fprintf(stderr, "(seq %d, id %d : %d )\n", kk, m_pSS2->GetSeq(kk)->m_iID, m_pSS2->GetSeq(kk)->GetPointNum());
}
double CSWMatch::Match(CSetOfSeq* pSS1, CSetOfSeq* pSS2, CAlignment* pASet, CModel* model, bool bTwoWay)
{
    //init
    m_model = model;
    m_bTwoWay = bTwoWay;
    m_mapCachedMatching.clear();
    m_pOSS1 = pSS1;
    m_pOSS2 = pSS2;
    m_pSS1 = new CSetOfSeq(*pSS1); // not so deep copy
    m_pSS2 = new CSetOfSeq(*pSS2); // not so deep copy
    m_pSS1->RemoveShortSeqs(3);
    m_pSS2->RemoveShortSeqs(3);
    m_pAlign = pASet;
    m_pAlign->m_pSS1 = pSS1;
    m_pAlign->m_pSS2 = pSS2;
    m_fTotalScore = 0;
    double fscore = 1;

    UpdateMatching();
    int iStepCount = 0;
    while (m_pSS1->GetSeqNum() > 0 && m_pSS2->GetSeqNum() > 0 && fscore > 0 && iStepCount < m_iMaxStep)
    {
        iStepCount += 1;
        double fMaxScore = -1;
        double fMaxProb = -1;
        int bestP1 = 0, bestP2 = 0;
        CAlignment* pBestAlign = NULL;
        std::map<pair<int, int>, CAlignment*>::iterator itr; 
        std::map<pair<int, int>, CAlignment*>::iterator maxitr; 
        for(itr = m_mapCachedMatching.begin(); itr != m_mapCachedMatching.end(); ++itr)
        {
             CAlignment* pAlign =(CAlignment*) (*itr).second;
             float fprob = pAlign->m_fScore / log(1 + pAlign->m_nLength1) / log(1 + pAlign->m_nLength2);
    //         float fprob = pAlign->m_fScore / pAlign->m_nLength1 /  pAlign->m_nLength2;
             if (fprob > fMaxProb)
             {
                 maxitr = itr;
                 fMaxScore = pAlign->m_fScore;
                 fMaxProb = fprob;
                 bestP1 = (*itr).first.first; 
                 bestP2 = (*itr).first.second;
     //           fprintf(stderr, "<< %d %d %.4f\n", bestP1, bestP2, fMaxScore);
                 pBestAlign = pAlign;
             }
        }
        //showinfo(bestP1, bestP2, fMaxScore);
        if (fMaxScore < 0.00001)
           return fMaxScore; 
        m_fTotalScore += fMaxScore;
        m_pAlign->AddAlignment(*pBestAlign);
       
        int start1, end1, start2, end2; 
        pBestAlign->GetBound(start1, end1, start2, end2);
        m_pSS1->SplitSeqByID(bestP1, start1, end1); 
        m_pSS2->SplitSeqByID(bestP2, start2, end2); 
        RemoveTable(bestP1, bestP2);
        m_pSS1->RemoveShortSeqs(3);
        m_pSS2->RemoveShortSeqs(3);
 //       fprintf(stderr, "(%f %d %d  %d %d )\n", fMaxScore, start1, end1, start2, end2); 
        UpdateMatching();
//        fprintf(stderr, "total score so far %f \n", m_fTotalScore);    
    }
    ReleaseMatchingTable();
    delete m_pSS1;
    delete m_pSS2;
    return m_fTotalScore;    
}

CAlignment* CSWMatch::FindBestMatch(int& bestP1, int& bestP2, double& fMaxScore)
{
    fMaxScore = -1;
    bestP1 = -1; bestP2 = -1;
    CAlignment* pBestAlign = NULL;
    
    for (int i = 0; i < m_pSS1->GetSeqNum(); i ++)
    {
        CSequence* pS1 = m_pSS1->GetSeq(i);
        for (int j = 0; j < m_pSS2->GetSeqNum(); j ++)  
        {
            CSequence* pS2 = m_pSS2->GetSeq(j);
            CAlignment* pAlign = new CAlignment();
            pAlign->m_pSS1 = m_pOSS1;
            pAlign->m_pSS2 = m_pOSS2;
            MatchSequence(pS1, pS2, pAlign);
//            fprintf(stderr, " * search : %d %d %.4f\n", i, j, pAlign->m_fScore); 
            if (fMaxScore < pAlign->m_fScore)
            {
                fMaxScore = pAlign->m_fScore;
                if (!pBestAlign) delete pBestAlign;
                pBestAlign = pAlign;
                bestP1 = pS1->m_iID; bestP2 = pS2->m_iID;
            }
            else
            {
                delete pAlign;
            }
        }
    }  
    return pBestAlign; 
}
/////////////////////////////////////////////////////////////////////////////////////////
// class CDynamicMatch
/////////////////////////////////////////////////////////////////////////////////////////

CDynamicMatch::CDynamicMatch()
{
}
CDynamicMatch::~CDynamicMatch()
{
}

double CDynamicMatch::InitCandidate(CAlignment* palign)
{
    if (m_pSS1 != NULL) delete m_pSS1;
    if (m_pSS2 != NULL) delete m_pSS2;
    m_pSS1 = new CSetOfSeq(*m_pOSS1);
    m_pSS2 = new CSetOfSeq(*m_pOSS2);
    m_fTotalScore = palign->m_fScore;
    int start1, end1, start2, end2; 
    palign->GetBound(start1, end1, start2, end2);
    int bestP1 = palign->m_SeqIndex1[0];
    int bestP2 = palign->m_SeqIndex2[0];
    m_pSS1->SplitSeqByID(bestP1, start1, end1); 
    m_pSS2->SplitSeqByID(bestP2, start2, end2); 
    fprintf(stderr, "(%f %d %d  %d %d )\n", palign->m_fScore, start1, end1, start2, end2); 
    return m_fTotalScore;
}

double CDynamicMatch::FindAllCandidate()
{
    for (int i = 0; i < m_pSS1->GetSeqNum(); i ++)
    {
        CSequence* pS1 = m_pSS1->GetSeq(i);
        for (int j = 0; j < m_pSS2->GetSeqNum(); j ++)  
        {
            CSequence* pS2 = m_pSS2->GetSeq(j);
            CAlignment* pAlign = new CAlignment();
            pAlign->m_pSS1 = m_pOSS1;
            pAlign->m_pSS2 = m_pOSS2;
            MatchSequence(pS1, pS2, pAlign);
            if (pAlign->m_fScore > 10.0)
            {
                m_vCandidate.push_back(pAlign);
                fprintf(stderr, "candidiate %d %d %f \n", i, j, pAlign->m_fScore); 
            }
        }
    }  
    fprintf(stderr, "%d candidates found\n", m_vCandidate.size());
    return 0;
}

double CDynamicMatch::DynaMatch(CSetOfSeq* pSS1, CSetOfSeq* pSS2, CAlignment* pASet, CModel* model, CModel* dynmodel,  bool bTwoWay)
{
    double fscore = 1;
    //init
    m_model = model;
    m_DynaModel = dynmodel;
    m_bTwoWay = bTwoWay;
    m_pOSS1 = pSS1;
    m_pOSS2 = pSS2;
    m_pSS1 = new CSetOfSeq(*pSS1); // not so deep copy
    m_pSS2 = new CSetOfSeq(*pSS2); // not so deep copy
    m_pSS1->RemoveShortSeqs(3);
    m_pSS2->RemoveShortSeqs(3);
    m_pAlign = pASet;
    m_pAlign->m_pSS1 = pSS1;
    m_pAlign->m_pSS2 = pSS2;
    m_fTotalScore = 0;

    FindAllCandidate();
    m_model = dynmodel;
    double fMaxScore = -1;
    for (int i =0; i < (int)m_vCandidate.size(); i ++)
    {
        fprintf(stderr, "\n=================\nStart with candidate %d\n", i);
        // use alignment i as the initial matching
        CAlignment* palign = m_vCandidate[i];
        InitCandidate(palign);
        while (m_pSS1->GetSeqNum() > 0 && m_pSS2->GetSeqNum() > 0 && fscore > 0)
        {
            vector<int> ref1;
            vector<int> ref2;
            GetRef(ref1, ref2);
            float dist1, dist2;
            GetMeanDist(dist1, dist2);
            CShapeContext::ExtractFeature(*m_pOSS1, ref1, (dist1 + dist2)/2.0);       
            CShapeContext::ExtractFeature(*m_pOSS2, ref2, (dist1 + dist2)/2.0);       
            int bestP1, bestP2;
            double fScore;
            CAlignment* pa = FindBestMatch(bestP1, bestP2, fScore);
            //pAlign->Print();
            palign->AddAlignment(*pa);
            int start1, end1, start2, end2; 
            pa->GetBound(start1, end1, start2, end2);
            m_pSS1->SplitSeqByID(bestP1, start1, end1); 
            m_pSS2->SplitSeqByID(bestP2, start2, end2); 
            m_pSS1->RemoveShortSeqs(3);
            m_pSS2->RemoveShortSeqs(3);
            //showinfo(bestP1, bestP2, fScore);
            delete pa;
        }
        fprintf(stderr, "total score %f \n", palign->m_fScore);    
        if (palign->m_fScore > fMaxScore) 
        {
            fMaxScore = palign->m_fScore;
            *m_pAlign = *palign;
        } 
    }
    delete m_pSS1;
    delete m_pSS2;
    return m_fTotalScore;    
}

void  CDynamicMatch::GetRef(vector<int>& ref1, vector<int>& ref2)
{
    ref1.clear();
    ref2.clear();
    for (int i = 0; i < m_pAlign->GetOperNum(); i ++)
    {
       int oper, seq1, seq2, pt1, pt2, layer; 
       m_pAlign->GetOper(i, oper, seq1, seq2, pt1, pt2, layer);
       if (seq1 >=0 && pt1 >= 0 && seq2 >= 0 && pt2 >=0)
       {
           ref1.push_back(m_pOSS1->GetPointIdx(seq1, pt1));
           ref2.push_back(m_pOSS2->GetPointIdx(seq2, pt2));
 //          fprintf(stderr, "oper %d: %d %d %d,  %d %d %d\n", i,  seq1, pt1, m_pOSS1->GetPointIdx(seq1, pt1),  seq2, pt2, m_pOSS2->GetPointIdx(seq2, pt2)); 
       }
    }
}
void  CDynamicMatch::GetMeanDist(float& dist1, float& dist2)
{
    float s1 = 0;
    float s2 = 0;
    vector<float> x1, y1; 
    vector<float> x2, y2; 
    float x, y;
    for (int i = 0; i < m_pAlign->GetOperNum(); i ++)
    {
       int oper, seq1, seq2, pt1, pt2, layer; 
       m_pAlign->GetOper(i, oper, seq1, seq2, pt1, pt2, layer);
       if (seq1 >=0 && pt1 >= 0)
       {
           m_pOSS1->GetXY(seq1, pt1, x, y);
           x1.push_back(x); y1.push_back(y);
       }
       if (seq2 >=0 && pt2 >= 0)
       {
           m_pOSS2->GetXY(seq2, pt2, x, y);
           x2.push_back(x); y2.push_back(y);
       }
    } 
    for (int i = 0; i <(int) x1.size(); i ++)
    {
        for (int j = 0; j < (int) x1.size(); j ++ )
        {
             s1 += sqrt((x1[i] - x1[j]) *(x1[i] - x1[j]) + (y1[i] - y1[j]) *(y1[i] - y1[j]));
//             fprintf(stderr, "%.3f ",  sqrt((x1[i] - x1[j]) *(x1[i] - x1[j]) + (y1[i] - y1[j]) *(y1[i] - y1[j])));
        }  
    }
    if (x1.size() > 0) s1 /= (x1.size() * x1.size());

 //   fprintf(stderr, "\n----\n ");
    for (int i = 0; i <(int) x2.size(); i ++)
    {
        for (int j = 0; j < (int) x2.size(); j ++ )
        {
             s2 += sqrt((x2[i] - x2[j]) *(x2[i] - x2[j]) + (y2[i] - y2[j]) *(y2[i] - y2[j]));
  //           fprintf(stderr, "%.3f ",  sqrt((x2[i] - x2[j]) *(x2[i] - x2[j]) + (y2[i] - y2[j]) *(y2[i] - y2[j])));
        }  
    }
    if (x2.size() > 0) s2 /= (x2.size() * x2.size());
//    fprintf(stderr, "ref sample %d %d\n", (int) x1.size(), x2.size());
    dist1 = s1; dist2 = s2;
}

