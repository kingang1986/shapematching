#ifndef SEQUENCE_H_
#define SEQUENCE_H_
#include <vector>
#include "math.h"
#include <string>
#include "stdlib.h"
#include "stdio.h"
using namespace std;

#define DATATYPE float
void chompstr(char* line);
int count_line(const char* file);
int count_column(char* line, char delim);
char* getnextstring(char* pstart, char* buffer, char delim);

class CMSSPoint
{
public:
    static int m_iFeatureDim;
public:
    CMSSPoint();
    virtual ~CMSSPoint();
    DATATYPE* m_pFeature;
    DATATYPE  m_fX;
    DATATYPE  m_fY;
    int Allocate(int iFeatureDim)
    {
        if (m_iFeatureDim != 0 && m_iFeatureDim != iFeatureDim)
        {
            fprintf(stderr, "Feature dimemsion inconsistant %d %d . ", m_iFeatureDim, iFeatureDim);
            exit(0);
        }
        if (m_pFeature != NULL)
        {
            delete [] m_pFeature;
        }
        m_pFeature  = new DATATYPE[iFeatureDim];
        m_iFeatureDim = iFeatureDim;
    }

    int     m_iOriginalSeqIdx;
    int     m_iOriginalPtIdx;
};

class CSequence
{
public:
    CSequence();
    CSequence(const CSequence& seq);
    ~CSequence();
    void Release();

    int AddPoint(CMSSPoint* pt);
    void Reverse();
    int GetPointNum() { return (int) m_vPoints.size();}
    DATATYPE* GetPointValue(int iIndex);
    vector<DATATYPE*> m_vFeature; // just reference, memory managed by CMSSPoint

    vector<CMSSPoint*> m_vPoints; // own the points;
    int m_iID;

    int m_iOriginalSeqId; 
    int m_iStartPos;
};

class CSetOfSeq
{
public:
    CSetOfSeq();
    CSetOfSeq(const CSetOfSeq& ss);
    CSetOfSeq(const char* strDataFile);
    ~CSetOfSeq();
    int AddSequence(CSequence* pSeq);
    int LoadSS(const char* strFile);
    int LoadSSBinary(const char* strFile);
    int SaveSSBinary(const char* strFile);

    int CheckSeq(const char* strFile, vector<int>& vSeqLength);
    void Print();
    void UpdateTotalPoints();
    void Release();
    int SplitSeq(int iSeqIndex, int iSplitPos1, int iSplitPos2);
    int SplitSeqByID(int iSeqID, int iSplitPos1, int iSplitPos2);
    int RemoveShortSeqs(int iMinLen);
    int GetXY(int iSeq, int iPt, float& x, float& y);
    int GetSeqNum() { return (int) m_vSeqs.size(); }
    int GetSeqLength(int iSeq) { if (iSeq < (int) m_vSeqs.size() ) return m_vSeqs[iSeq]->GetPointNum();}

//Data
public:    
//   int m_iSeqNum;
    int m_iTotalPoint;
    vector<CSequence*> m_vSeqs;
 //   vector<int> m_vSeqLength;
    string m_strFileName;
    int m_iShapeID;
    int m_iClassID;
    int m_iSeqIds;
};
#endif
