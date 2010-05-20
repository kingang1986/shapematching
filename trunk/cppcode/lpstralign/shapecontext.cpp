#include "shapecontext.h"
#include "sequence.h"
#include "alignment.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "memory.h"
#define PI 3.1415923
int CShapeContext::BINS = 5;
int CShapeContext::THETA = 12;
vector<int> CShapeContext::FEATIDX;


int CShapeContext::GetFeatureIdx()
{
    FEATIDX.clear();
    for (int i = 0; i < BINS * THETA; i ++)
    {
        char tmpb[40];
        sprintf(tmpb, "sc%2d", i);
        int idx = CMSSPoint::GetFeatureIdx(tmpb);
        if (idx < 0)
            fprintf(stderr, "[ERROR] feature idx %s -> %d\n", tmpb, idx);
        FEATIDX.push_back(idx);
    }
    return BINS*THETA;
}

int CShapeContext::ExtractFeature(CSetOfSeq& mss, vector<int>& vShapeRef, double fMeanDist)
{
    fprintf(stderr, " mean dist %f\n", fMeanDist);
    int nSample = mss.m_iTotalPoint;
    int nRef = (int) vShapeRef.size();
    double r_inner = log(1.0/8);
    double r_outer = log(2.0);
    float* X = new float[nSample];
    float* Y = new float[nSample];
    float* logspace = new float[BINS];
    
    for (int i = 0; i < BINS; i ++)
    {
        logspace[i] = exp(r_inner + i * (r_outer - r_inner) / (BINS - 1));
        //fprintf(stderr, "%d: %f \n", i, logspace[i]);
    }

    int icount = 0;
    for (int s = 0; s < mss.GetSeqNum(); s ++)
    {
        for (int p = 0; p < mss.GetSeqLength(s); p ++)
        {
            mss.GetXY(s, p, X[icount], Y[icount]);
            icount ++;
        }
    }
    DATATYPE* buffer = new DATATYPE[BINS*THETA]; 
    float ang = 2 * PI / THETA; 
    for (int i = 0; i < nSample; i ++)
    {
        memset(buffer, 0, sizeof(DATATYPE) * (BINS*THETA));
        for (int j = 0;j < nRef; j ++) 
        {
            int idx = vShapeRef[j];
            float dy = Y[i] - Y[idx]; 
            float dx = X[i] - X[idx]; 
            float theta = atan2(dy, dx);
            if (theta < 0) theta += 2* PI;
            int nTheta = int(theta / ang);
            if (nTheta >= THETA) fprintf(stderr, "ntheta outof bound!!!\n");
            float dist = sqrt(dx * dx + dy * dy)/fMeanDist;
            int nDist = 0;
            for (int k = 0; k < BINS; k ++)
                if (dist < logspace[k]) nDist ++;
            if (nDist > 0)
            {
                int iFeatureIdx = nTheta * BINS + nDist - 1;  
                if (iFeatureIdx < 0 || iFeatureIdx >= 60) fprintf(stderr, "[ERROR] ! %d %d %d\n", iFeatureIdx, nTheta, nDist);
                buffer[iFeatureIdx] += 1.0;
            }
        } 
        CMSSPoint* pt = mss.GetPoint(i);
        for (int s = 0; s < BINS*THETA; s ++)
        {
//            if (s == 34) fprintf(stderr, "feat idx %d %d %d %.4f\n", i, s, FEATIDX[s], buffer[s]/nRef);
            pt->m_pLFeature[FEATIDX[s]] = pt->m_pRFeature[FEATIDX[s]] = buffer[s]/nRef;
        }
    }  
    return BINS*THETA;
    
}
int CShapeContext::Clean(CSetOfSeq& mss)
{
    for (int i = 0; i < mss.m_iTotalPoint; i ++)
    {
        CMSSPoint* pt = mss.GetPoint(i);
        for (int s = 0;  s < BINS*THETA; s ++ )
        {
//            fprintf(stderr, "pt %d, %d, feat %d \n", i, s,  FEATIDX[s]); 
            pt->m_pLFeature[FEATIDX[s]] = pt->m_pRFeature[FEATIDX[s]] = 0;
        }
    } 
    return BINS*THETA;
}

