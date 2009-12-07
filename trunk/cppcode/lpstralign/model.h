#ifndef MODEL_HEADER
#define MODEL_HEADER
#include "alignment.h"

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
    void Default(int iFeatureNum);
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
#endif
