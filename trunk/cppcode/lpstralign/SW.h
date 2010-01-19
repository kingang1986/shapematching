#ifndef SW_H_
#define SW_H_
#include <vector>
#include "sequence.h"
#include "model.h"
#include "alignment.h"
#include <map>



using namespace std;


/* Useful structure for SW algorithm */
class CSWNode
{
public:
    double score;
    CSWNode* prev;
    int    operation;
    int    row_index;
    int    column_index;
};

//direct MSS matching

class CSWMatch
{
public:
    CSWMatch() 
    {
         m_pSS1 = m_pSS2 = NULL;
         m_model = NULL;
        m_iMaxStep = -1;
    }   
    ~CSWMatch(){};
    double Match(CSetOfSeq* pSS1, CSetOfSeq* pSS2, CAlignment* pASet, CModel* pModel, bool bTwoway = false);
    double MatchSequence(CSequence* pSeqA, CSequence* pSeqB, CAlignment* pAlign); 
    double MatchSequenceOneWay(CSequence* pSeqA, CSequence* pSeqB, CAlignment* pAlign);
    double getSubstituteCost(DATATYPE* a, DATATYPE* b);
    double getGapCost(DATATYPE* a);

    CSetOfSeq* m_pOSS1; //reference to the objects
    CSetOfSeq* m_pOSS2;
    CSetOfSeq* m_pSS1; // cloned local objects
    CSetOfSeq* m_pSS2;

    double m_fTotalScore;
    CModel* m_model;
    bool m_bTwoWay;
    int  m_iMaxStep;
    CAlignment* m_pAlign;

    //hash table manage: map pair <int,int> to an alignment
    std::map<pair<int,int>, CAlignment*> m_mapCachedMatching;// cached matching score 
    void ReleaseMatchingTable();
    void RemoveTable(int shapei, int shapej);
    void UpdateMatching();
protected:
    void showinfo(int bestP1, int bestP2, double fMaxScore);
    CAlignment* FindBestMatch(int& bestP1, int& bestP2, double& fMaxScore);
 
};

// dynamic matching where the shape context features are updated with matching proceeds
class CDynamicMatch : public CSWMatch
{
public:
    CDynamicMatch();
    ~CDynamicMatch();

    double DynaMatch(CSetOfSeq* pSS1, CSetOfSeq* pSS2, CAlignment* pASet, CModel* pModel, CModel* dynmodel, bool bTwoway = false);
    CModel* m_DynaModel;
    vector<CAlignment*> m_vCandidate;

protected:
    double FindAllCandidate();
    double InitCandidate(CAlignment*);
    double Verify();
    void   GetRef(vector<int>& ref1, vector<int>& ref2);
    void   GetMeanDist(float& dist1, float& dist2);
    
};
#endif
