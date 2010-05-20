#ifndef APP_HEADER
#define APP_HEADER
#include "SW.h"
#include "sequence.h"
#include "model.h"
#include "alignment.h"
#include "lp.h"
#include "shapecontext.h"

double match2shapes(const char* shape1, const char* shape2, const char* modelfile, bool twoway = false)
{
    CSetOfSeq S1, S2;
    S1.Load(shape1);
    S2.Load(shape2);
    CModel m;
    if (strcmp(modelfile, "") == 0)
    {
       m.Default(CMSSPoint::m_iFeatureDim);
//       m.Write("tmp.model.txt");
    }
    else
       m.Read(modelfile);
    //CShapeContext::GetFeatureIdx();
    CDynamicMatch match;
    CAlignment align;
    match.DynaMatch(&S1,&S2,&align, &m, &m, twoway);  
    printf("matching score %f \n", align.m_fScore);
    return align.m_fScore; 
}



#endif
