%module mss 
%include typemaps.i
%{
#include "alignment.h"
%}

class CAlignment
{
public:


    CAlignment();
    ~CAlignment();
    double  AddAlignment(CAlignment& align);
    double  GetPhi(double* phi, int iParamDim, CModel* model);
    void    GetBound(int& start1, int& end1, int& start2, int& end2);
    int GetOperNum();
    int GetOper(int iIndex, int& OUTPUT, int& OUTPUT, int& OUTPUT, int&OUTPUT, int& OUTPUT);
    CSetOfSeq* m_pSS1;
    CSetOfSeq* m_pSS2;
    
    double m_fScore;
};
