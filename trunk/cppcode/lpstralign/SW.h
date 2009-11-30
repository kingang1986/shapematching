#ifndef SW_H_
#define SW_H_
#include <vector>
#include "sequence.h"
#include "model.h"
#include "alignment.h"
#include <map>



using namespace std;

//double SmithWaterman(CSetOfSeq* pSS1, CSetOfSeq* pSS2, int iSSIndex1, int iSSIndex2, CSequence* pSeqA, CSequence* pSeqB, CAlignment* pAlign, CModel* model);
//double SmithWaterman_SetOfSeq(CSetOfSeq* pSS1, CSetOfSeq* pSS2, CAlignment* pASet, CModel* model);

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
    double Match(CSetOfSeq* pSS1, CSetOfSeq* pSS2, CAlignment* pASet, CModel* pModel);
    double Match(CSequence* pSeqA, CSequence* pSeqB, CAlignment* pAlign, CModel* model);

protected:
    double getSubstituteCost(DATATYPE* a, DATATYPE* b, CModel* model);
    double getGapCost(DATATYPE* a, CModel* model);
    std::map<pair<int,int>, CAlignment*> m_mapCachedMatching;// cached matching score 
    void ReleaseMatchingTable();
    void RemoveTable(int shapei, int shapej);
    CSetOfSeq* m_pOSS1;
    CSetOfSeq* m_pOSS2;
    CSetOfSeq* m_pSS1;
    CSetOfSeq* m_pSS2;
    CModel* m_model;
    void UpdateMatching();
};

#endif
