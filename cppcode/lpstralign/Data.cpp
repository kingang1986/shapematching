# include "Data.h"
#include <cstdlib>
#include <string.h>
#include <stdio.h>


void chompstr(char* line)
{
    int len = strlen(line);
    if (line[len - 1] == '\n')
    {
        line[len - 1] = 0;
    }
    if (line[len - 2] == '\r')  line[len-2] = 0;
    return;
}


int count_line(const char* file) // by longbin
{
    FILE* f;
    char line[32768];
    int counter;
    if ((f = fopen (file, "r")) == NULL)
    { perror (file); exit (1); }
    counter = 0;
    while (!feof(f))
    {
        fgets(line,  32768, f);
        chompstr(line);
        if (strlen(line) > 2)
        {
            counter ++;
        }
        line[0] = 0;
    }
    fclose(f);
    return counter;

}


int count_column(char* line, char delim)
{
    int count = 0;

    char* pbuffer = line;
    while(pbuffer != NULL)
    {
        pbuffer = strchr(pbuffer, delim);
        if (pbuffer != NULL)
        {
            pbuffer ++;
        }
        count ++;
    }
    return count;
}

char* getnextstring(char* pstart, char* buffer, char delim)
{
    char* pnext  = NULL;
    if (pstart == NULL)
    {
        buffer[0] = 0;
        return NULL;
    }
    pnext = strchr(pstart, delim);
    if (pnext == NULL)
    {
        strcpy(buffer, pstart);
    }
    else
    {
        memcpy(buffer, pstart, pnext - pstart);
        buffer[pnext - pstart  ] = 0;
        pnext ++;
    }
    return pnext;
}


/////////////////////////////////////////////////////////////
//class CSequence
/////////////////////////////////////////////////////////////
CSequence::CSequence()
{
    m_vX = NULL;
    m_vY = NULL;
    m_vFeature = NULL;
}
CSequence::~CSequence()
{
    Release();
}

void CSequence::Release()
{
    if (m_vFeature != NULL)
    {
        for (int i = 0; i < m_iPoint; i ++)
        {
            DATATYPE* pPoint = m_vFeature[i];
            delete [] pPoint;
        }
        delete [] m_vFeature;
    }
    if (m_vX != NULL) delete [] m_vX;
    if (m_vY != NULL) delete [] m_vY;
    m_vFeature = NULL;
    m_vX = m_vY = NULL;
    m_iPoint = 0;
    m_iFeatureDim = 0;
    return ;
}

int CSequence::Allocate(int nPoint, int nFeatureDim)
{
    Release();
    m_iFeatureDim = nFeatureDim;
    m_iPoint = nPoint;
    m_vFeature = new DATATYPE*[nPoint];
    m_vX = new DATATYPE[nPoint];
    m_vY = new DATATYPE[nPoint];
    for (int i = 0; i < m_iPoint; i ++)
    {
        DATATYPE* p = new DATATYPE[m_iFeatureDim];
        m_vFeature[i] = p;
    }
    return nPoint * nFeatureDim;

}

DATATYPE* CSequence::GetPoint(int iIndex)
{
    return m_vFeature[iIndex];
}

/////////////////////////////////////////////////////////////
//class CSetOfSeq
/////////////////////////////////////////////////////////////
CSetOfSeq::CSetOfSeq()
{
    m_iSeqNum = 0;
    m_iTotalPoint = 0;
}
void CSetOfSeq::Release()
{
    for (int i = 0; i < (int) m_vSeqs.size(); i ++)
    {
        CSequence* pSeq = m_vSeqs[i];
        delete pSeq;
    }
    m_vSeqs.clear();
    m_iSeqNum = 0;
    m_iTotalPoint = 0;

}
CSetOfSeq::~CSetOfSeq()
{
    Release();
}

int CSetOfSeq::AddSequence(CSequence* pSeq)
{
    m_vSeqs.push_back(pSeq);
    m_iSeqNum ++;
    m_iTotalPoint += (int)(pSeq->m_iPoint);
    return m_iSeqNum;
}

int CSetOfSeq::CheckSeq(const char* file)
{
    FILE* f;
    char line[32768];
    char tmp[32768];
    int iFeatureDim = -1;
    //fprintf(stderr, "checking MSS file\n");
    if ((f = fopen(file, "r")) == NULL)
    {
        fprintf(stderr, "can't open file %s\n", file);
        exit(1);
    }
    int linecount = 0;
    fgets(line, 32768, f); // head line
    m_vSeqLength.clear();
    while(!feof(f))
    {
        line[0] = 0;
        fgets(line, 32768, f);
        if (strlen(line) > 2)
        {
             chompstr(line);
             //get first two columns
             memcpy(tmp, line, 32768);
             char tmpbuffer[512];
             char* pbuffer = tmp;
             pbuffer = getnextstring(pbuffer, tmpbuffer, '\t');
             int seq = atoi(tmpbuffer);
             if (seq >= (int)m_vSeqLength.size())
             m_vSeqLength.push_back(0);
             pbuffer = getnextstring(pbuffer, tmpbuffer, '\t');
             int ptn = atoi(tmpbuffer);
             //fprintf(stderr, "%d %d \n", seq, ptn);
             m_vSeqLength[seq] = ptn + 1;
             int iColumn = count_column(line, '\t');
             if (iFeatureDim != -1 && iColumn != iFeatureDim)
             {
                 fprintf(stderr, "inconsistant feature length (%d , %d) at line %d\n|%s| ", iColumn, iFeatureDim, linecount, line);
                 exit(1);
             }
             
             iFeatureDim = iColumn;
        }
   }
   fclose(f);
   m_iFeatureDim = iFeatureDim - 4;
   m_iSeqNum = (int)m_vSeqLength.size();
//   for (int k = 0; k < m_vSeqLength.size(); k ++)
//       fprintf(stderr, "%d of %d, feat len %d\n", k, m_vSeqLength.size(), m_iFeatureDim);
   return (int)m_vSeqLength.size();
}


int CSetOfSeq::LoadSSBinary(const char* strFile)
{
    Release();

    FILE* f = fopen(strFile, "r");
    if (!f)
    {
        fprintf(stderr, "Can't open file %s to read data", strFile);
    }
    fread(&m_iFeatureDim, sizeof(int), 1, f);
    fread(&m_iTotalPoint, sizeof(int), 1, f);
    fread(&m_iSeqNum, sizeof(int), 1, f);
    for (int i = 0; i < m_iSeqNum; i ++)
    {
        int t;
        fread(&t, sizeof(int), 1, f);
        m_vSeqLength.push_back(t);
    }
    for (int i = 0; i < m_iSeqNum; i ++)
    {
        CSequence* pSeq = new CSequence();
        pSeq->m_iFeatureDim = m_iFeatureDim;
        pSeq->m_iPoint = m_vSeqLength[i];
        pSeq->m_vFeature = new DATATYPE*[m_vSeqLength[i]];
        for (int j = 0;j < m_vSeqLength[i];j ++)
        {
            DATATYPE* pdata = new DATATYPE[m_iFeatureDim];
            fread(pdata, sizeof(DATATYPE), m_iFeatureDim, f);
            pSeq->m_vFeature[j] = pdata;
        }
        m_vSeqs.push_back(pSeq);

    }
      fclose(f);
     return 1;

}
int CSetOfSeq::SaveSSBinary(const char* strFile)
{

    FILE* f=fopen(strFile, "w");
    //write file header
    fwrite(&m_iFeatureDim, sizeof(int), 1, f);
    fwrite(&m_iTotalPoint, sizeof(int), 1, f);
    fwrite(&m_iSeqNum, sizeof(int), 1, f);
    for (int i = 0; i < m_iSeqNum; i ++)
    {
        int t = m_vSeqLength[i];
        fwrite(&t, sizeof(int), 1, f);
    }
    for (int i = 0; i < m_iSeqNum; i ++)
    {
        CSequence* pSeq = m_vSeqs[i];
        for (int j = 0; j < m_vSeqLength[i]; j ++)
        {
            DATATYPE* pdata = pSeq->GetPoint(j);
            fwrite(pdata, sizeof(DATATYPE), m_iFeatureDim, f);
        }
    }
    fclose(f);
    return 1;
}


int CSetOfSeq::LoadSS(const char* strFile)
{
    FILE* f;
    char line[32768];
    char tmp[32768];
    //start here
    //fprintf(stderr, "cehck file");
    CheckSeq(strFile) ; //count the number of sequence separators
    if ((f = fopen(strFile, "r")) == NULL)
    {
        fprintf(stderr, "cant' open file %s", strFile);
        return -1;
    }
    fgets(line, 32768, f); //head line
    m_iTotalPoint = 0;
    for (int i = 0; i < (int)m_vSeqLength.size(); i ++)
    {
        CSequence* pSeq = new CSequence();
        pSeq->Allocate(m_vSeqLength[i], m_iFeatureDim);
        for (int j = 0; j < m_vSeqLength[i]; j ++)
        {
            line[0] = 0;
            fgets(line, 32768, f);
            chompstr(line);
            if (strlen(line) > 1)
            {
                memcpy(tmp, line, 32768);
                char tmpbuffer[512];
                char* pbuffer = tmp;
                pbuffer = getnextstring(pbuffer, tmpbuffer, '\t');
                int seq = atoi(tmpbuffer);
                pbuffer = getnextstring(pbuffer, tmpbuffer, '\t');
                int ptn = atoi(tmpbuffer);
                pbuffer = getnextstring(pbuffer, tmpbuffer, '\t');
                double x = atof(tmpbuffer); 
                pbuffer = getnextstring(pbuffer, tmpbuffer, '\t');
                double y = atof(tmpbuffer); 
    //            fprintf(stderr, "\n%d %d %f %f", seq, ptn, x, y);
                pSeq->m_vX[j] = x;
                pSeq->m_vY[j] = y;
                for (int k = 0; k < m_iFeatureDim; k ++)
                {
                    pbuffer = getnextstring(pbuffer, tmpbuffer, '\t');
                    pSeq->GetPoint(j)[k] = (DATATYPE)atof(tmpbuffer);
   //                 fprintf(stderr, "%f ", atof(tmpbuffer));
                }
            }
        }
        m_vSeqs.push_back(pSeq);
        m_iTotalPoint += m_vSeqLength[i];
    }
    fclose(f);
    //Print();
    return m_iTotalPoint;

}

void CSetOfSeq::Print()
{
    fprintf(stderr, "%d sequence, %d points, %d feature dim\n", (int) m_vSeqs.size(), m_iTotalPoint, m_iFeatureDim);
    for (int ss = 0; ss < (int)m_vSeqs.size(); ss ++)
    {
        fprintf(stderr, "%d sequence is %d long\n", ss + 1, m_vSeqs[ss]->m_iPoint);
        for (int k = 0; k < m_vSeqs[ss]->m_iPoint; k ++)
        {
            DATATYPE* pData = m_vSeqs[ss]->GetPoint(k);
            for (int t = 0; t < m_iFeatureDim; t++)
            {
                fprintf(stderr, "%f\t", pData[t]);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "##########\n");
    }

}

/////////////////////////////////////////////////////////
//class CAlignment
/////////////////////////////////////////////////////////
CAlignment::CAlignment()
{
    m_pSS1 = m_pSS2 = NULL;
    m_fScore = 0.0;
}

CAlignment::~CAlignment()
{

}

double CAlignment::AddAlignment(CAlignment& align)
{
    if ((m_pSS1 != NULL && m_pSS1 != align.m_pSS1) || (m_pSS2 != NULL && m_pSS2 != align.m_pSS2))
    {
        fprintf(stderr, "can't merge two alignment because the shapes are different\n");
        return -1;
    }
    m_pSS1 = align.m_pSS1;
    m_pSS2 = align.m_pSS2;
    for (int i = 0; i < (int) align.m_operation.size(); i ++)
    {
        m_operation.push_back(align.m_operation[i]);
        m_SeqIndex1.push_back(align.m_SeqIndex1[i]);
        m_SeqIndex2.push_back(align.m_SeqIndex2[i]);
        m_PointIndex1.push_back(align.m_PointIndex1[i]);
        m_PointIndex2.push_back(align.m_PointIndex2[i]);
    }
    m_fScore += align.m_fScore;
    return m_fScore;

}


/////////////////////////////////////////////////////////////
//class CModel
/////////////////////////////////////////////////////////////

CModel::CModel()
{
    m_vWeight = NULL;
    m_vFeatureIndex = NULL;
    m_vWeightIndex = NULL;
    m_vMapType = NULL;
    m_vMatchOrGap = NULL;
    m_vMatch = NULL;
    m_vGap = NULL;
    m_vSign = NULL;
    m_vTheta = NULL;
}

CModel::~CModel()
{
    Release();
}

void CModel::Print()
{
    /* Writes structural model sm to file file. */
    fprintf(stderr, "\nweight\tfeaidx\tw_idx\tmaptype\tm_or_g\n");
    for (int i = 0; i < m_iMapNum; i ++)
    {
        fprintf(stderr, "%d\t%d\t%d\t%d\n",   m_vFeatureIndex[i], m_vWeightIndex[i],m_vMapType[i], m_vMatchOrGap[i]);
    }
    for (int i = 0; i < m_iParamDim; i ++)
    {
        fprintf(stderr, "%f ", m_vWeight[i]);
    }
    fprintf(stderr, "\n");
    for (int i = 0; i < m_iParamDim; i ++)
    {
        fprintf(stderr, "%d ", m_vSign[i]);
    }

    fprintf(stderr, "\nmax weight index %d, mapping # %d\n", m_iParamDim, m_iMapNum);
}


int CModel::Write(const char* strFile)
{
    /* Writes structural model sm to file file. */
    FILE* fp = fopen(strFile, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "Can't open file %s for writing model.\n", strFile);
        exit(1);
    }
    fprintf(fp, "feaidx\tw_idx\tmaptype\tm_or_g\n");
    for (int i = 0; i < m_iMapNum; i ++)
    {
        fprintf(fp, "%d\t%d\t%d\t%d\n",  m_vFeatureIndex[i], m_vWeightIndex[i],m_vMapType[i], m_vMatchOrGap[i]);
    }
    for (int i = 0; i < m_iParamDim; i ++)
    {
        fprintf(fp, "%f ", m_vWeight[i]);
    }
    fprintf(fp, "\n");
    for (int i = 0; i < m_iParamDim; i ++)
    {
        fprintf(fp, "%d ", m_vSign[i]);
    }
    fprintf(fp, "\n");

    fclose(fp);
    return m_iParamDim;
}

int CModel::Read(const char* strFile)
{

    Release();

    int itotaline = count_line(strFile);
    char line[32768];
    //fprintf(stderr, "Reading struct model ...");
    FILE* fp = fopen(strFile, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "can't open model file %s\n", strFile);
        exit(0);
    }
    fgets(line, 32768, fp); //head line
    Init(itotaline - 1);
    int matchcount = 0;
    int gapcount = 0;
    int maxWeightIndex = -1;
    for (int i = 0; i < itotaline - 3; i ++)
    {
        int idx, maptype, matchgap,  wid;
        fscanf(fp, "%d %d %d %d", &idx, &wid, &maptype, &matchgap);
        //fprintf(stderr, "%f %d %d %d %d",idx, wid, maptype, matchgap);
        if (wid > maxWeightIndex) maxWeightIndex = wid;
        m_vFeatureIndex[i] = idx;
        m_vWeightIndex[i] = wid;
        m_vMapType[i] = maptype;
        m_vMatchOrGap[i] = matchgap;
        if (matchgap == 101) matchcount ++;
        if (matchgap == 102) gapcount ++;
    }

    m_iParamDim = maxWeightIndex + 1;
    fprintf(stderr, "param dim = %d\n",m_iParamDim);
    m_iMapNum = itotaline - 3;
    m_vWeight = new double[m_iParamDim];
    for (int i = 0; i < m_iParamDim; i ++)
    {
        float w;
        fscanf(fp, "%f", &w);
        fprintf(stderr, "%f\t", w);
        m_vWeight[i] = w;
    }
    fprintf(stderr, "\nsign of each weight\n");
    for (int i = 0; i < m_iParamDim; i ++)
    {
        int w;
        fscanf(fp, "%d", &w);
        fprintf(stderr, "%d\t", w);
        m_vSign[i] = w;
    }
    fclose(fp);
    fprintf(stderr, "\ndone\n");

    //build
    m_vMatch = new int[matchcount];
    m_vGap = new int[gapcount];
    int idx1 = 0;
    int idx2 = 0;
    for (int i = 0; i < itotaline - 2; i ++)
    {
        if (m_vMatchOrGap[i] == MAPTYPE_MATCH)
        {
            m_vMatch[idx1] = i; idx1 ++ ;
        }
        else if (m_vMatchOrGap[i] == MAPTYPE_GAP)
        {
            m_vGap[idx2] = i; idx2 ++;
        }
    }
    m_iMatchCount = matchcount;
    m_iGapCount = gapcount;

    return m_iParamDim;

}

void CModel::Release()
{
    if (m_vWeight != NULL)
    {
        delete [] m_vWeight;
    }
    if (m_vFeatureIndex != NULL)
    {
        delete [] m_vFeatureIndex;
    }
    if (m_vWeightIndex != NULL)
    {
        delete [] m_vWeightIndex;
    }
    if (m_vMatchOrGap != NULL)
    {
        delete [] m_vMatchOrGap;
    }
    if (m_vTheta != NULL)
    {
        delete [] m_vTheta;
    }
    if (m_vMatch!= NULL)
    {
        delete [] m_vMatch;
    }
    if (m_vGap!= NULL)
    {
        delete [] m_vGap;
    }
    if (m_vSign != NULL)
    {
        delete [] m_vSign;
    }
    m_vWeight = NULL;
    m_vMatch = m_vGap = NULL;
    m_vFeatureIndex = NULL;
    m_vWeightIndex = NULL;
    m_vMatchOrGap = NULL;
    m_vSign  = NULL;
    m_vTheta = NULL;
    return ;
}

void CModel::InitTheta(int iPatternNum)
{
    m_vTheta = new double[iPatternNum];
    m_iPatternNum = iPatternNum;
}

void CModel::Init(int iParamDim)
{
    Release();
    m_vWeight = new double[iParamDim];
    m_vMapType = new int[iParamDim];
    m_vFeatureIndex = new int[iParamDim];
    m_vWeightIndex = new int[iParamDim];
    m_vSign = new int[iParamDim];
    m_vMatchOrGap = new int[iParamDim];
    m_iParamDim = iParamDim;
    return;
}

