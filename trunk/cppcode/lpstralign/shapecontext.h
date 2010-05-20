#ifndef SHAPECONTEXT_HEADER
#define SHAPECONTEXT_HEADER
#include "sequence.h"
#include "alignment.h"

    
class CShapeContext
{
public: 

//parameters;
    static int BINS, THETA; 
    static vector<int> FEATIDX;
//    static int ExtractFeature(CSetOfSeq& mss);
    static int ExtractFeature(CSetOfSeq& mss, vector<int>& vShapeRef, double fMeanDist);
 //   static int ExtractFeature(CSetOfSeq& mss, CAlignment& align);
    static int Clean(CSetOfSeq& mss);
    static int GetFeatureIdx();
};
#endif
