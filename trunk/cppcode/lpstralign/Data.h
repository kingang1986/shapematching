#ifndef DATA_H_
#define DATA_H_
#include <vector>
#include "math.h"
#include <string>
#include "sequence.h"
using namespace std;


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
