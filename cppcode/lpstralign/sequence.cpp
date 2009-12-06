#include <cstdlib>
#include <string.h>
#include <stdio.h>
#include "sequence.h"


void chompstr(char* line)
{
    int len = strlen(line);
    if (len <= 0) return;
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
//class CMSSPoint
/////////////////////////////////////////////////////////////

int CMSSPoint::m_iFeatureDim = 0;
CMSSPoint::CMSSPoint()
{
   m_pFeature = NULL;
   m_fX = m_fY = -1;
   m_iOriginalSeqIdx = m_iOriginalPtIdx = -1;
}
CMSSPoint::~CMSSPoint()
{
    if (m_pFeature != NULL)
    {
        delete [] m_pFeature;
    }
    m_pFeature = NULL;
}


/////////////////////////////////////////////////////////////
//class CSequence
/////////////////////////////////////////////////////////////
CSequence::CSequence()
{
   m_bOwner = true;
}

CSequence::CSequence(const CSequence& seq) // shallow copy
{
    for (int i = 0; i < (int)seq.m_vPoints.size(); i ++)
    {
        m_vPoints.push_back(seq.m_vPoints[i]);
        m_vFeature.push_back(seq.m_vFeature[i]);
    }
    m_iID = seq.m_iID;
    m_bOwner = false;
   

}
CSequence::~CSequence()
{
    Release();
}

void CSequence::Release()
{
    if (m_bOwner)
    {
        for (int i = 0; i < GetPointNum(); i ++)
        {
            CMSSPoint* pt = m_vPoints[i];
            delete pt;
        }
    }
    m_vPoints.clear();
    return ;
}

void CSequence::Reverse() // inplace reverse 
{
    int n = GetPointNum(); 
    for (int i = 0; i < n/2; i ++)
    {
        CMSSPoint* pt = m_vPoints[i];
        m_vPoints[i] = m_vPoints[n - i];
        m_vPoints[n - i] = pt;

        DATATYPE* dt = m_vFeature[i];
        m_vFeature[i] = m_vFeature[n - i];
        m_vFeature[n - i] = dt;
    }

}


DATATYPE* CSequence::GetPointValue(int iIndex)
{
    return m_vFeature[iIndex];
}


int CSequence::AddPoint(CMSSPoint* pt)
{
    m_vPoints.push_back(pt);
    m_vFeature.push_back(pt->m_pFeature);
    return (int) m_vPoints.size();
}

/////////////////////////////////////////////////////////////
//class CSetOfSeq
/////////////////////////////////////////////////////////////

CSetOfSeq::CSetOfSeq()
{
    m_iTotalPoint = 0;
    m_iSeqIds = 10001;
    m_iClassID = -1;
    m_iShapeID = -1;
}

CSetOfSeq::CSetOfSeq(const char* strDataFile)
{
    CSetOfSeq();
    LoadSS(strDataFile);

} 
CSetOfSeq::CSetOfSeq(const CSetOfSeq& ss)
{
    m_iSeqIds = ss.m_iSeqIds;
    m_iTotalPoint = ss.m_iTotalPoint;
    for (int i = 0; i < (int) ss.m_vSeqs.size(); i ++)
    {
        CSequence* pSeq = new CSequence(*(ss.m_vSeqs[i]));
        m_vSeqs.push_back(pSeq); 
    }
    m_iClassID = ss.m_iClassID;
    m_iShapeID = ss.m_iShapeID;
    
}

void CSetOfSeq::Release()
{
    for (int i = 0; i < (int) m_vSeqs.size(); i ++)
    {
        CSequence* pSeq = m_vSeqs[i];
        delete pSeq;
    }
    m_vSeqs.clear();
    m_iShapeID = m_iClassID = -1;
    m_iTotalPoint = 0;

}

CSetOfSeq::~CSetOfSeq()
{
    Release();
}

int CSetOfSeq::GetX(int iSeq, int iPt)
{
    double x =-1;
    if (iSeq < (int) m_vSeqs.size())
    {
        CSequence* pSeq = m_vSeqs[iSeq];
        if (iPt < pSeq->GetPointNum())
        {
            x = pSeq->m_vPoints[iPt]->m_fX;
            return x;
        }
    }
    return x;
}
int CSetOfSeq::GetY(int iSeq, int iPt)
{
    double y =-1;
    if (iSeq < (int) m_vSeqs.size())
    {
        CSequence* pSeq = m_vSeqs[iSeq];
        if (iPt < pSeq->GetPointNum())
        {
            y = pSeq->m_vPoints[iPt]->m_fY;
            return y;
        }
    }
    return y;
}
int CSetOfSeq::GetXY(int iSeq, int iPt, float& x, float& y)
{
    x = y = -1;
    if (iSeq < (int) m_vSeqs.size())
    {
        CSequence* pSeq = m_vSeqs[iSeq];
        if (iPt < pSeq->GetPointNum())
        {
            x = pSeq->m_vPoints[iPt]->m_fX;
            y = pSeq->m_vPoints[iPt]->m_fY;
            return 1;
        }
    }
    return -1;
}

int CSetOfSeq::AddSequence(CSequence* pSeq)
{
    pSeq->m_iID = m_iSeqIds;
    m_iSeqIds ++;
    m_vSeqs.push_back(pSeq);
    UpdateTotalPoints();
    return (int) m_vSeqs.size();
}

int CSetOfSeq::CheckSeq(const char* file, vector<int>& vSeqLength)
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
    vSeqLength.clear();
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
             if (seq >= (int)vSeqLength.size())
             vSeqLength.push_back(0);
             pbuffer = getnextstring(pbuffer, tmpbuffer, '\t');
             int ptn = atoi(tmpbuffer);
             //fprintf(stderr, "%d %d \n", seq, ptn);
             vSeqLength[seq] = ptn + 1;
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
   CMSSPoint::m_iFeatureDim = iFeatureDim - 4;
   //for (int k = 0; k < m_vSeqLength.size(); k ++)
    //   fprintf(stderr, "%d of %d, feat len %d\n", k, m_vSeqLength.size(), m_iFeatureDim);
   UpdateTotalPoints();
   return (int)vSeqLength.size();
}

int CSetOfSeq::RemoveShortSeqs(int iMinLen)
{
    vector<CSequence*>::iterator iter = m_vSeqs.begin();
    while( iter != m_vSeqs.end() )
    {
      if ((*iter)->GetPointNum() < iMinLen)
        iter = m_vSeqs.erase( iter );
      else
        ++iter;
    }
    UpdateTotalPoints();
    return (int) m_vSeqs.size();

}


int CSetOfSeq::SplitSeqByID(int iSeqID, int iSplitPos1, int iSplitPos2)
{
    for (int i = 0; i < (int) m_vSeqs.size() ; i ++)
    {
        if (m_vSeqs[i]->m_iID == iSeqID)
            return SplitSeq(i, iSplitPos1 , iSplitPos2);
    }
    fprintf(stderr, "Can't find sequence %d!!!\n", iSeqID);
    return -1;
}

int CSetOfSeq::LoadSSBinary(const char* strFile)
{
 /*
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
*/
    return 1;

}
int CSetOfSeq::SaveSSBinary(const char* strFile)
{
/*
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
            DATATYPE* pdata = pSeq->GetPointValue(j);
            fwrite(pdata, sizeof(DATATYPE), m_iFeatureDim, f);
        }
    }
    fclose(f);
*/
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
    vector<int> vSeqLength;
    CheckSeq(strFile, vSeqLength) ; //count the number of sequence separators
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
    for (int i = 0; i < (int)vSeqLength.size(); i ++)
    {
        CSequence* pSeq = new CSequence();
        AddSequence(pSeq);
        for (int j = 0; j < vSeqLength[i]; j ++)
        {
            line[0] = 0;
            fgets(line, 32760, f);
            chompstr(line);
            if (strlen(line) > 1)
            {
                memcpy(tmp, line, 32768);
                CMSSPoint* newpt = new CMSSPoint();
                newpt->m_iOriginalSeqIdx = pSeq->m_iID;
                newpt->m_iOriginalPtIdx = j; 
           
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
                newpt->m_fX = x;
                newpt->m_fY = y;
                newpt->Allocate();
                for (int k = 0; k < CMSSPoint::m_iFeatureDim; k ++)
                {
                    pbuffer = getnextstring(pbuffer, tmpbuffer, '\t');
                    newpt->m_pFeature[k] = (DATATYPE) atof(tmpbuffer); 
                }
                pSeq->AddPoint(newpt);
            }
        }
        m_iTotalPoint += vSeqLength[i];
    }
    fclose(f);
    UpdateTotalPoints();
    return m_iTotalPoint;

}

void CSetOfSeq::Print()
{
/*
    fprintf(stderr, "%d sequence, %d points, %d feature dim\n", (int) m_vSeqs.size(), m_iTotalPoint, m_iFeatureDim);
    for (int ss = 0; ss < (int)m_vSeqs.size(); ss ++)
    {
        fprintf(stderr, "%d sequence is %d long\n", ss + 1, m_vSeqs[ss]->m_iPoint);
        for (int k = 0; k < m_vSeqs[ss]->m_iPoint; k ++)
        {
            DATATYPE* pData = m_vSeqs[ss]->GetPointValue(k);
            for (int t = 0; t < m_iFeatureDim; t++)
            {
                fprintf(stderr, "%f\t", pData[t]);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "##########\n");
    }
*/
}


void CSetOfSeq::UpdateTotalPoints()
{
    m_iTotalPoint = 0;
    for (int i = 0; i < (int)m_vSeqs.size(); i ++)
    {
       m_iTotalPoint += m_vSeqs[i]->GetPointNum();
    }
    return; 
}
//sqq
int CSetOfSeq::SplitSeq(int iSeqIndex, int iSplitPos1, int iSplitPos2)
{
    CSequence* pSeq = m_vSeqs[iSeqIndex];
    fprintf(stderr, "split %d (len %d ) from %d to %d \n", iSeqIndex, pSeq->GetPointNum(), iSplitPos1, iSplitPos2);
    int iLen = pSeq->GetPointNum();
    if (iSplitPos1 > 0)
    {
        CSequence* pSeq1 = new CSequence();
        for (int i = 0; i < iSplitPos1; i ++)
        {
            pSeq1->m_vFeature.push_back(m_vSeqs[iSeqIndex]->m_vFeature[i]);
            pSeq1->m_vPoints.push_back(m_vSeqs[iSeqIndex]->m_vPoints[i]);
        }
        AddSequence(pSeq1);
        pSeq1->m_bOwner = m_vSeqs[iSeqIndex]->m_bOwner;
        
    } 
    if (iSplitPos2 < iLen - 1)
    {
        CSequence* pSeq2 = new CSequence();
    
        for (int i = iSplitPos2 + 1; i < iLen; i ++)
        {
            pSeq2->m_vFeature.push_back(m_vSeqs[iSeqIndex]->m_vFeature[i]);
            pSeq2->m_vPoints.push_back(m_vSeqs[iSeqIndex]->m_vPoints[i]);
        }
        AddSequence(pSeq2);
        pSeq2->m_bOwner = m_vSeqs[iSeqIndex]->m_bOwner;
    }
    
    m_vSeqs[iSeqIndex]->m_bOwner = false;
    m_vSeqs.erase(m_vSeqs.begin() + iSeqIndex);
    delete  pSeq;
    UpdateTotalPoints();
 //   fprintf(stderr, " total %d seqs\n", m_vSeqs.size());
    return (int) m_vSeqs.size();
 
}
