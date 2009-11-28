#ifndef DATA_H_
#define DATA_H_
#include <vector>
#include "math.h"
#include <string>
using namespace std;

#define DATATYPE float

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




void chompstr(char* line);
int count_line(const char* file);
int count_column(char* line, char delim);
char* getnextstring(char* pstart, char* buffer, char delim);

class CSequence
{
public:
    CSequence();
    CSequence(const CSequence& seq);
    ~CSequence();
    int Allocate(int nPoint, int nFeatureDim);
    void Release();
    void SetPointValue(int iPoint, DATATYPE* vFeature);
    DATATYPE* GetPoint(int iIndex);

//data
public: 
    int m_iFeatureDim;
    int m_iPoint;
    DATATYPE** m_vFeature;
    DATATYPE* m_vX;
    DATATYPE* m_vY;
    int m_iID;

// for chopped sequence
    int m_iOriginalSeqId; 
    int m_iStartPos;
};

class CSetOfSeq
{

public:
    CSetOfSeq();
    CSetOfSeq(const CSetOfSeq& ss);
    ~CSetOfSeq();
    int AddSequence(CSequence* pSeq);
    int LoadSS(const char* strFile);
    int LoadSSBinary(const char* strFile);
    int SaveSSBinary(const char* strFile);

    int CheckSeq(const char* strFile);
    void Print();
    void Update();
    void Release();
    int SplitSeq(int iSeqIndex, int iSplitPos1, int iSplitPos2);
    int SplitSeqByID(int iSeqID, int iSplitPos1, int iSplitPos2);
    int RemoveShortSeqs(int iMinLen);


//Data
public:    
    int m_iSeqNum;
    int m_iFeatureDim;
    int m_iTotalPoint;
    vector<CSequence*> m_vSeqs;
    vector<int> m_vSeqLength;
    string m_strFileName;
    int m_iShapeID;
    int m_iClassID;
    int m_iSeqIds;
    
};


class CModel
{
protected:
    void Release();
public:
    CModel();
    ~CModel();

public: 
    int Read(const  char* strFile);
    int Write(const char* strFile);
    void Init(int iParamDim);
    void InitTheta(int iPatternNum);
    void Print();

    // info about training 
    int m_iPatternNum;

    //features
    int m_iFeatureDim; // the number of features
//    vector<string> m_vFeatureName;

    //weights i.e. theta
    int m_iParamDim; // the number of Params, i.e., the length of weights
    int* m_vFeatureIndex;
    int* m_vWeightIndex;
    int* m_vMapType; // the method to calculate similarity measurement 
    int* m_vMatchOrGap; //MAPTYPE_MATCH or MAPTYPE_GAP

    //
    int m_iGapCount;
    int* m_vGap;
    int m_iMatchCount;
    int* m_vMatch;
    double* m_vTheta;

    //model weights
    int m_iMapNum;
    double* m_vWeight;
    int* m_vSign; // the constraints for the weights

};

class CAlignment
{
public:
    CSetOfSeq* m_pSS1;
    CSetOfSeq* m_pSS2;
    vector<int> m_SeqIndex1;
    vector<int> m_SeqIndex2;
    vector<int> m_PointIndex1;
    vector<int> m_PointIndex2;
    vector<int> m_operation;
    double m_fScore;

    CAlignment();
    ~CAlignment();
    double AddAlignment(CAlignment& align);
    double GetPhi(double* phi, int iParamDim, CModel* model);
    void GetBound(int& start1, int& end1, int& start2, int& end2);
    bool m_bSameClass;
};


#endif /*DATA_H_*/
