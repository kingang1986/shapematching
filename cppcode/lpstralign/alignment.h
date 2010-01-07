#ifndef ALIGN_HEADER
#define ALIGN_HEADER
#include <vector>
#include "math.h"
#include <string>
#include "sequence.h"

#define MAPTYPE_CONSTANT 0
#define MAPTYPE_UNCHANGE 1
#define MAPTYPE_ABS 2
#define MAPTYPE_EU 3
#define MAPTYPE_ABSEU 4
#define MAPTYPE_MATCH 101
#define MAPTYPE_GAP 102
#define MATCH       1
#define SUBST       2
#define DELET       3
#define INSRT       4

class CModel;
class CAlignment
{
public:
    CAlignment();
    ~CAlignment();
    double  AddAlignment(CAlignment& align);
    void    Clean(); //clear all info
    double  GetPhi(double* phi, int iParamDim, CModel* model);
    void    GetBound(int& start1, int& end1, int& start2, int& end2);
    int     GetOperNum();
    void    Print();
#ifndef SWIG
    int GetOper(int iIndex, int& oper, int& seq1, int& seq2, int&pt1, int& pt2, int& layer);
#else
    int GetOper(int iIndex, int& OUTPUT, int& OUTPUT, int& OUTPUT, int& OUTPUT, int& OUTPUT, int& OUTPUT);
#endif
    
    CSetOfSeq* m_pSS1;
    CSetOfSeq* m_pSS2;
    
    vector<int> m_SeqIndex1;
    vector<int> m_SeqIndex2;
    vector<int> m_PointIndex1;
    vector<int> m_PointIndex2;
    vector<int> m_operation;
    vector<int> m_layer;
    double m_fScore;
    int m_iStart1, m_iEnd1, m_iStart2, m_iEnd2;
    int m_nLength1, m_nLength2; // the lengths of two sequences in matching
protected:
    int m_iLayerCount;
};
#endif
