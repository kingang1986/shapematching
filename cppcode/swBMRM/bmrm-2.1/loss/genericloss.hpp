#ifndef _GENERICLOSS_HPP_
#define _GENERICLOSS_HPP_

#include "common.hpp"
#include "sml.hpp"
#include "info.hpp"
#include "data.hpp"
#include "genericdata.hpp"
#include "loss.hpp"
#include <string>
#include "timer.hpp"
#include "math.h"

#include <vector>

#include <blitz/array.h>
using namespace blitz;
using namespace std;

//#define SAVEMATCHES

/** 
 * Class for graphmatch classification loss.
 */
class CGenericLoss : public CLoss
{
  public:
    CGenericLoss(CModel* &model, CGenericData* &data);
    virtual ~CGenericLoss();

    // Interfaces
    void Usage();
    void ComputeLoss(Scalar& loss)
    {
      throw CBMRMException("ERROR: not implemented!\n", "CGraphMatchLoss::ComputeLoss()");
    }
    void ComputeLossAndGradient(Scalar& loss, TheMatrix& grad);
    void Predict(CModel *model);
    void Evaluate(CModel *model);

    Scalar LabelLoss(int search_ind, int* scene_inds);

    void LoadModel(std::string modelFilename="");
    void SaveModel(std::string modelFilename="");

    int* assignment(int n1, Scalar* weights, adjmatrix** res_adjs, int test);

  protected:
    void Phi(int n1, int* Knearest, adjmatrix** matches, Scalar* res);
    void Phi(int n1, int n2, adjmatrix* matches, Scalar* res);
    void Phi1(Point p1, Point p2, Scalar* a, Scalar* b, Scalar* res);
    int dimOfFeature;
    int dimOfWeight;
    int _Kneighbours;
    string _LossType; // choices: AVG_K, KNN
    void KNNLoss(Scalar& loss, TheMatrix& grad);
    void AvgKLoss(Scalar& loss, TheMatrix& grad);
    adjmatrix* GetDP(int n1, int n2, Scalar* weights, Scalar & cost);
    int* findTopK(int n1, int c, Scalar* weight, adjmatrix** res_adj, Scalar& cost);
    int* findTopKLabeled(int n1, int c, Scalar* weight, adjmatrix** res_adj, Scalar& cost);


    CGenericData* _data;
};

#endif
