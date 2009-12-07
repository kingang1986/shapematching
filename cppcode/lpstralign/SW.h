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


class CSWMatch
{
public:
    CSWMatch() 
    {
         m_pSS1 = m_pSS2 = NULL;
         m_model = NULL;
    }   
    ~CSWMatch(){};
    double Match(CSetOfSeq* pSS1, CSetOfSeq* pSS2, CAlignment* pASet, CModel* pModel, bool bTwoway = false);

protected:
    double InitMatch(CSetOfSeq* pSS1, CSetOfSeq* pSS2, CAlignment* pASet, CModel* model, bool bTwoWay);
    double CleanUpMatch();
    double MatchAStep();
    double MatchSequenceOneWay(CSequence* pSeqA, CSequence* pSeqB, CAlignment* pAlign, CModel* model);
    double MatchSequence(CSequence* pSeqA, CSequence* pSeqB, CAlignment* pAlign, CModel* model, bool bTwoway = false);
    double getSubstituteCost(DATATYPE* a, DATATYPE* b, CModel* model);
    double getGapCost(DATATYPE* a, CModel* model);
    std::map<pair<int,int>, CAlignment*> m_mapCachedMatching;// cached matching score 
    void ReleaseMatchingTable();
    void RemoveTable(int shapei, int shapej);
    CSetOfSeq* m_pOSS1;
    CSetOfSeq* m_pOSS2;
    CSetOfSeq* m_pSS1;
    CSetOfSeq* m_pSS2;
    void UpdateMatching();
    double m_fTotalScore;

    CModel* m_model;
    bool m_bTwoWay;
    CAlignment* m_pAlign;
    
 
};

#endif
