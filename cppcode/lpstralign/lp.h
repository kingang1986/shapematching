#ifndef LP_H_
#define LP_H_

#include <glpk.h>

//A pattern is an original shape, along with the homologs and decoys

class CPattern
{
public:
    static int  m_iTopK;
    static int  m_iTotalClass;

public:
    CPattern();
    ~CPattern();

    void Align();
    double GetLabelLoss();

    CSetOfSeq* m_original;
    vector<vector<CSetOfSeq*> > m_vShapeClass;
    vector<vector<CAlignment> > m_vShapeAlign;
    vector<vector<CAlignment*> > m_vSortedShapeAlign;

    void   AlignDecoy(CModel* pModel);
    void   AlignHomolog(CModel* pModel);
    int    AlignClass(int iClassId, CModel* pModel);
    double GetLabelPhi(double*, int iParamDim, CModel*);
    int GetDecoyPhi(double*, int iParamDim, CModel*);
    double GetLabelLoss(CModel*);

    double Classify(CModel* pModel) ;
    int m_id;

protected:
    DATATYPE* m_pLabelPhi;
    double GetClassPhi(double* pw, int iParamDim, int iClassID, CModel* pModel) ;
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
    int Add(double* pWeight, int iPatternIndex, double loss);
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
    vector<double> m_vLoss;

    double m_fEpsilon, m_fC, m_fDistance;

};

class CSample
{
protected:
    CSetOfSeq* LoadSoS(const  char* strSoSFile);
    CSetOfSeq* FindSoS(const  char* strSoSFile);
    int Release();


public:
    int m_iFeatureDim;
    bool m_bBinaryData;

    vector<char*> m_vFileNames;
    vector<CSetOfSeq*> m_vSS;
    vector<CPattern*> m_vPatterns;
    char m_szFolder[500];

public: //members
    CSample();
    ~CSample();
    int AlignHomolog(CModel* pModel);
    int Classify(CModel* pModel, const char* outputfile);

    double UpdateConstraint(CConstraints* pCC, CModel* pModel, bool bActiveOnly, int iMinNewConstraint, int& iNumVio, double fEpsilon);

    //active patterns
    vector<int> m_vActive;
    int SetAllActive();
    int GetActiveNum();
    int m_iLastUpdated;

    int LoadShape(const char* strFile);
    int LoadPattern(const char* strFile);

};


class CStructureLearning
{
public:
    CStructureLearning();
    ~CStructureLearning();
    int Learn(const char* );
    int Init(const char* datafile, const char* patternfile, const char* szpath, const char* modelfile, bool bBinaryData);

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
