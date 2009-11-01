#ifndef LP_H_
#define LP_H_

#include <glpk.h>

//A pattern is an original shape, along with the homologs and decoys

class CPattern
{
public:
    int m_id;
    CSetOfSeq* m_original;
    vector<CSetOfSeq*> m_vHomolog;
    vector<CSetOfSeq*> m_vDecoy;
    int m_iHomolog;
    int m_iDecoy;
    int m_iBestHomolog;
    CAlignment m_homo_align;
    double m_fHomologScore;
    int m_iBestDecoy;
    CAlignment m_decoy_align;
    double m_fDecoyScore;
};

//////////////////////////////////////////////////////
//Linear constraints
//////////////////////////////////////////////////////
class CConstraints
{

public:
    CConstraints();
    ~CConstraints();
    int  Init(CModel* model, int iPatternNum, double fEpsilon);
    int Add(double* pWeight, int iPatternIndex);
    void Clear();
    int PrintMathProg(const char* filename);
    int Save(const char* file);
    int Load(const char* file);

    int GLPK_lp(CModel* pmodel);

public: //data
    int m_iWeightLength; //|w|
    int m_iPatternNum;// |epsilon|
    vector<double*> m_vWeights;
    vector<int> m_vPatternIndex;

    double m_fEpsilon, m_fC, m_fDistance;

};

class CSample
{
protected:
    CSetOfSeq* LoadSoS(const  char* strSoSFile);
    int Release();


public:
    int m_nShapes; //not uniq
    int m_nPattern;
    int m_nSetOfSeq;
    int m_iFeatureDim;
    bool m_bBinaryData;

    vector<char*> m_vFileNames;
    vector<CSetOfSeq*> m_vSS;
    vector<CPattern*> m_vPatterns;
    char m_szFolder[500];
    void UpdateHomologScore(CModel* pModel);

public: //members
    CSample();
    ~CSample();
    int LoadSample(const  char* strFile);
    //int AlignPatterns(CModel* pModel);
    int AlignSamples(CModel* pModel, const char* outputfile);

    double AlignDecoy(CConstraints* pCC, CModel* pModel, bool bActiveOnly, int iMinNewConstraint, int& iNumVio);
    int AlignHomolog(CModel* pModel);
    double AddConstraint(CConstraints* pCC, CModel* pModel, bool bActiveOnly);

    //active patterns
    vector<int> m_vActive;
    int SetAllActive();
    int GetActiveNum();

};


class CStructureLearning
{
public:
    CStructureLearning();
    ~CStructureLearning();
    int Learn(const char* );
    int Init(const char* datafile, const char* szpath, const char* modelfile, bool bBinaryData);

    int FindMostViolate(CPattern* p);

public:
    bool m_bLoadModelFromFile;
    char* m_strModelFile;
    CSample m_Sample;
    CModel m_Model;



    double m_fEpsilon;
    double m_fDistance;
    double m_fC;
    int m_iMaxIteration;
    int m_iMaxStep;
    int m_iMinNewConstraint;

};

#endif /*LP_H_*/
