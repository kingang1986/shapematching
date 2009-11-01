#ifndef DATA_H_
#define DATA_H_
#include <vector>
#include "math.h"
using namespace std;
#define DATATYPE float

#define MAPTYPE_CONSTANT 0
#define MAPTYPE_UNCHANGE 1
#define MAPTYPE_ABS 2
#define MAPTYPE_EU 3
#define MAPTYPE_ABSEU 4
#define MAPTYPE_MATCH 101
#define MAPTYPE_GAP 102



void chompstr(char* line);
int count_line(const char* file);
int count_column(char* line, char delim);
char* getnextstring(char* pstart, char* buffer, char delim);

class CSequence
{
public:
    CSequence();
    ~CSequence();
    int Allocate(int nPoint, int nFeatureDim);
    DATATYPE* GetPoint(int iIndex);
	int m_iFeatureDim;
	int m_iPoint;
	DATATYPE** m_vFeature;
    void Release();

};

class CSetOfSeq
{

public:
    CSetOfSeq();
    ~CSetOfSeq();
    int AddSequence(CSequence* pSeq);
    int LoadSS(const char* strFile);
    int LoadSSBinary(const char* strFile);
    int SaveSSBinary(const char* strFile);

    int CheckSeq(const char* strFile);
    void Print();
    void Release();
//Data
    int m_iSeqNum;
    int m_iFeatureDim;
    int m_iTotalPoint;
    vector<CSequence*> m_vSeqs;
    vector<int> m_vSeqLength;

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
};


class CModel
{
protected:
	void Release();
public:
	CModel();
	~CModel();

public: //data
	int Read(const  char* strFile);
	int Write(const char* strFile);
	void Init(int iParamDim);
	void InitTheta(int iPatternNum);
	void Print();

	//data
	int m_iPatternNum;
	int m_iFeatureDim;
	int m_iParamDim;
	int* m_vFeatureIndex;
	int* m_vWeightIndex;
	int* m_vMapType;
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
	int* m_vSign;

};

/*
class CModel
{
protected:
	void Release();

public:
	CModel();
	~CModel();
	int Read(const  char* strFile);
	int Write(const char* strFile);
	void Init(int iParamDim, int iFeatureDim);
	void InitTheta(int iPatternNum);
	void Print();

	//data
	int m_iPatternNum;
	int m_iFeatureDim;
	int m_iParamDim;
	double* m_vWeight;
	int* m_vSign;
	int* m_vGapMap;
	double* m_vTheta;
	int* m_vMatchMap;
	bool m_bGapExtension;

};
*/
#endif /*DATA_H_*/
