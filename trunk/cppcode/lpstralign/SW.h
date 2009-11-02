#ifndef SW_H_
#define SW_H_
#include <vector>
#include "Data.h"



using namespace std;

double SmithWaterman(CSetOfSeq* pSS1, CSetOfSeq* pSS2, int iSSIndex1, int iSSIndex2, CSequence* pSeqA, CSequence* pSeqB, CAlignment* pAlign, CModel* model);
double SmithWaterman_SetOfSeq(CSetOfSeq* pSS1, CSetOfSeq* pSS2, CAlignment* pASet, CModel* model);
//double GetPhi(CAlignment* pAlign, double* phi, int iParamDim, CModel* model);

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
#endif
