#ifndef SEQUENCE_H_
#define SEQUENCE_H_
#include <vector>
#include "math.h"
#include <string>
#include "stdlib.h"
#include "stdio.h"
#include <map>
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
    static const char* GetFeatureName(int iIdx);
    static int GetFeatureIdx(const char* ); 
    static int AddFeature(const char* newf);
    static vector<string>  m_vFeatureName;
   

public:
    CMSSPoint();
    virtual ~CMSSPoint();
    void Reverse();
    DATATYPE* m_pFeature;

   
    DATATYPE  m_fX;
    DATATYPE  m_fY;
    int Allocate()
    {
       
        if (m_iFeatureDim <= 0)
        {
            fprintf(stderr, "Feature dimemsion inconsistant %d . ", m_iFeatureDim);
            exit(0);
        }
        if (m_pLFeature != NULL)
        {
            delete [] m_pLFeature;
        }
        m_pLFeature  = new DATATYPE[m_iFeatureDim];
        if (m_pRFeature != NULL)
        {
            delete [] m_pRFeature;
        }
        m_pRFeature  = new DATATYPE[m_iFeatureDim];
        m_pFeature = m_pLFeature; //initially, left way features as the default feature 
        return m_iFeatureDim;
    }

    int     m_iOriginalSeqIdx;
    int     m_iOriginalPtIdx;

    DATATYPE* m_pRFeature; //one way feature , feature startswith 'f'
    DATATYPE* m_pLFeature; //reverse way feature, feature name startswith 'g'

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
    bool m_bOwner;

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
#ifdef SWIG
    int GetXY(int iSeq, int iPt, float& OUTPUT, float& OUTPUT);
#else
    int GetXY(int iSeq, int iPt, float& x, float& y);
#endif
    int GetX(int iSeq, int iPt);
    int GetY(int iSeq, int iPt);
    int GetSeqNum() { return (int) m_vSeqs.size(); }
    int GetSeqLength(int iSeq) { if (iSeq < (int) m_vSeqs.size() ) return m_vSeqs[iSeq]->GetPointNum();}
    CSequence* GetSeq(int iSeq) { return m_vSeqs[iSeq];}
    void SetFeatureValue(int seq, int pt, int idx, DATATYPE value, int bidirecti = 0); //bidrect : 0->both, 1->left, 2 -. right

public:    
    int m_iTotalPoint;
    vector<CSequence*> m_vSeqs;
    string m_strFileName;
    int m_iShapeID;
    int m_iClassID;
    int m_iSeqIds;
};
#endif
