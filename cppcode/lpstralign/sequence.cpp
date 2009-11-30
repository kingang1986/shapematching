#include "Data.h"
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


int count_line(const char* file) 
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

CSequence::CSequence(const CSequence& seq) // deep copy
{
    m_vX = NULL;
    m_vY = NULL;
    m_vFeature = NULL;
    
    m_iFeatureDim = seq.m_iFeatureDim;
    m_iPoint = seq.m_iPoint;
    Allocate(m_iPoint, m_iFeatureDim);
    memcpy(m_vX, seq.m_vX, sizeof(DATATYPE) * m_iPoint);
    memcpy(m_vY, seq.m_vY, sizeof(DATATYPE) * m_iPoint);
    for (int i = 0; i < m_iPoint; i ++)
    {
         memcpy(m_vFeature[i], seq.m_vFeature[i], sizeof(DATATYPE) * m_iFeatureDim);
    }
    m_iID = seq.m_iID;
    m_iOriginalSeqId = seq.m_iOriginalSeqId;
    m_iStartPos = seq.m_iStartPos;

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
    //printf("%d ", nPoint);
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

void CSequence::SetPointValue(int iPoint, DATATYPE* vFeature)
{
    memcpy(m_vFeature[iPoint], vFeature, sizeof(DATATYPE) * m_iFeatureDim);
    return ;
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
    m_iSeqIds = 10001;
    m_iClassID = -1;
    m_iShapeID = -1;
}
CSetOfSeq::CSetOfSeq(const CSetOfSeq& ss)
{
    m_iSeqIds = ss.m_iSeqIds;
    m_iSeqNum = ss.m_iSeqNum;
    m_iFeatureDim = ss.m_iFeatureDim;   
    m_iTotalPoint = ss.m_iTotalPoint;
    for (int i = 0; i < (int) ss.m_vSeqs.size(); i ++)
    {
        CSequence* pSeq = new CSequence(*(ss.m_vSeqs[i]));
        m_vSeqs.push_back(pSeq); 
    }
    m_iClassID = ss.m_iClassID;
    
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
    pSeq->m_iID = m_iSeqIds;
    m_iSeqIds ++;
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
   //for (int k = 0; k < m_vSeqLength.size(); k ++)
    //   fprintf(stderr, "%d of %d, feat len %d\n", k, m_vSeqLength.size(), m_iFeatureDim);
   Update();
   return (int)m_vSeqLength.size();
}

int CSetOfSeq::RemoveShortSeqs(int iMinLen)
{
    vector<CSequence*>::iterator iter = m_vSeqs.begin();
    while( iter != m_vSeqs.end() )
    {
      if ((*iter)->m_iPoint < iMinLen)
        iter = m_vSeqs.erase( iter );
      else
        ++iter;
    }
    Update();
    return m_iSeqNum;

}


int CSetOfSeq::SplitSeqByID(int iSeqID, int iSplitPos1, int iSplitPos2)
{
   
//    fprintf(stderr, "Split seq with ID %d \n", iSeqID);
    for (int i = 0; i < m_iSeqNum; i ++)
    {
        if (m_vSeqs[i]->m_iID == iSeqID)
            return SplitSeq(i, iSplitPos1 , iSplitPos2);
    }
    fprintf(stderr, "Can't find sequence %d!!!\n", iSeqID);
    return -1;
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
    Update();
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
    //strcpy(strFile, m_strFileName);
    //start here
    //fprintf(stderr, "cehck file");
    CheckSeq(strFile) ; //count the number of sequence separators
    if ((f = fopen(strFile, "r")) == NULL)
    {
        fprintf(stderr, "can't open file %s", strFile);
        return -1;
    }
    //    strcpy(m_strFileName, strFile);
    
    m_strFileName = strFile;
    //m_strFileName = "hello"; 
    fgets(line, 32768, f); //head line
    m_iTotalPoint = 0;
    for (int i = 0; i < (int)m_vSeqLength.size(); i ++)
    {
        CSequence* pSeq = new CSequence();
        pSeq->m_iID = i;
        pSeq->Allocate(m_vSeqLength[i], m_iFeatureDim);
        pSeq->m_iStartPos = 0;
        for (int j = 0; j < m_vSeqLength[i]; j ++)
        {
            line[0] = 0;
            fgets(line, 32760, f);
            chompstr(line);
            if (strlen(line) > 1)
            {
                memcpy(tmp, line, 32768);
                char tmpbuffer[512];
                char* pbuffer = tmp;
                pbuffer = getnextstring(pbuffer, tmpbuffer, '\t');
                //int seq = atoi(tmpbuffer);
                pbuffer = getnextstring(pbuffer, tmpbuffer, '\t');
                //int ptn = atoi(tmpbuffer);
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
        pSeq->m_iOriginalSeqId = (int)m_vSeqs.size();
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


void CSetOfSeq::Update()
{
    m_iSeqNum = (int) m_vSeqs.size();
    m_iTotalPoint = 0;
    for (int i = 0; i < (int)m_vSeqs.size(); i ++)
    {
       m_iTotalPoint += m_vSeqs[i]->m_iPoint; 
    }
    return; 
}

int CSetOfSeq::SplitSeq(int iSeqIndex, int iSplitPos1, int iSplitPos2)
{
    CSequence* pSeq = m_vSeqs[iSeqIndex];
    iSplitPos1 -= pSeq->m_iStartPos;
    iSplitPos2 -= pSeq->m_iStartPos;
//    fprintf(stderr, "split %d (len %d ) from %d to %d \n", iSeqIndex, pSeq->m_iPoint, iSplitPos1, iSplitPos2);
    int iLen = pSeq->m_iPoint;
    if (iSplitPos1 > 0)
    {
        CSequence* pSeq1 = new CSequence();
        pSeq1->Allocate(iSplitPos1, m_iFeatureDim);
        for (int i = 0; i < iSplitPos1; i ++)
        {
            pSeq1->SetPointValue(i, m_vSeqs[iSeqIndex]->GetPoint(i));
        }
        pSeq1->m_iOriginalSeqId = m_vSeqs[iSeqIndex]->m_iOriginalSeqId;
        pSeq1->m_iStartPos = m_vSeqs[iSeqIndex]-> m_iStartPos;
        AddSequence(pSeq1);
        
    } 
    if (iSplitPos2 < iLen - 1)
    {
        CSequence* pSeq2 = new CSequence();
    
        pSeq2->Allocate(iLen - iSplitPos2 - 1, m_iFeatureDim);
        for (int i = iSplitPos2 + 1; i < iLen; i ++)
        {
             pSeq2->SetPointValue(i - iSplitPos2 - 1, m_vSeqs[iSeqIndex]->GetPoint(i));
        }
        pSeq2->m_iOriginalSeqId = m_vSeqs[iSeqIndex]->m_iOriginalSeqId;
        pSeq2->m_iStartPos = m_vSeqs[iSeqIndex]-> m_iStartPos + iSplitPos2 + 1;
        AddSequence(pSeq2);
    }
    m_vSeqs.erase(m_vSeqs.begin() + iSeqIndex);
    delete  pSeq;
    Update();
 //   fprintf(stderr, " total %d seqs\n", m_vSeqs.size());
    return m_iSeqNum; 
 
}
