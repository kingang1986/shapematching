#ifndef _PARSER_H_
#define _PARSER_H_



#include <stdarg.h>
#include <vector>
#include <string>
#include <string.h>
#include <stdlib.h>
#define _MAX_PATH 1024

/**********************************************************
// class CCmdArg: basic cmd arg class
**********************************************************/
class CCmdArg
{
public:
	enum  CmdArgSyntax 
	{
		isOPT       = 0x00,  // argument is optional
		isREQ       = 0x01,  // argument is required
		isVALOPT    = 0x02,  // argument value is optional
		isVALREQ    = 0x04,  // argument value is required
		isLIST      = 0x20,  // argument is a list
		isPOS       = 0x40,  // argument is positional
		isHIDDEN    = 0x80,  // argument is not to be printed in usage
		isVALTAKEN  = (isVALREQ | isVALOPT),    // argument takes a value
		isOPTVALOPT = (isOPT | isVALOPT),
		isOPTVALREQ = (isOPT | isVALREQ),
		isPOSVALOPT = (isPOS | isVALOPT),
		isPOSVALREQ = (isPOS | isVALREQ),
	};

	CCmdArg(const char cOptChar, const char *szKeyword, const char *szValue, const char * szDescription, unsigned int syntax_flags)
	{
		m_cOptChar = cOptChar;
		m_szKeyword = szKeyword;
		m_szValue = szValue;
		m_szDescription = szDescription;
		m_dwSyntaxFlags = syntax_flags;
		m_bArgAvailable = false;
	}
	char GetCharName() { return m_cOptChar; }
	const char *GetKeywordName() { return m_szKeyword; }
	const char *GetValueName() { return m_szValue; }
	const char *GetDescription() { return m_szDescription; }
	unsigned int GetFlags() { return m_dwSyntaxFlags; }
	virtual int SetValue(char *argv[]) = 0;
	bool IsAvailable() { return m_bArgAvailable; }
	void SetAvailable() { m_bArgAvailable = true; }
	virtual int ValueNeeded() { return 1; }
private:
	char	m_cOptChar;
    const	char	*m_szKeyword;
    const	char	*m_szValue;
    const	char	*m_szDescription;
	unsigned int m_dwSyntaxFlags;
	bool	m_bArgAvailable;
};

class CCmdArgInt : public CCmdArg
{
public:
	CCmdArgInt(char cOptChar, 
	const	       char *szKeyword, 
	const		   char *szValue, 
    char iDefault,
	const		   char* szDescription, 
    unsigned int syntax_flags = CCmdArg::isOPTVALREQ)
	    : CCmdArg(cOptChar, szKeyword, szValue, szDescription, syntax_flags), m_iValue(iDefault) {}
	virtual ~CCmdArgInt() {}
	operator int(void)  const { return  m_iValue; }
	virtual int SetValue(char *argv[])
	{
		m_iValue = atoi(*argv++);
		return 0;
	}
private:
	int		m_iValue;
};

class CCmdArgIntList : public CCmdArg
{
public:
	CCmdArgIntList(char cOptChar, 
		       const char *szKeyword, 
		const 	   char *szValue, 
		const	   int iDefault,
		const	   int iIntNumber,
		const	   char* szDescription, 
			   unsigned int syntax_flags = CCmdArg::isOPTVALREQ)
		: CCmdArg(cOptChar, szKeyword, szValue, szDescription, syntax_flags)
	{
		m_arrValue.resize(iIntNumber);
		for (int i = 0; i < iIntNumber; i++)
			m_arrValue[i] = iDefault;
	}
	virtual ~CCmdArgIntList(){}
	int operator[] (int i) { return m_arrValue[i]; }
	virtual int SetValue(char *argv[])
	{
		for (unsigned int i = 0; i < m_arrValue.size(); i++)
		{
			m_arrValue[i] = atoi(*argv++);
		}
		return 0;
	}
	virtual int ValueNeeded() { return (int)m_arrValue.size(); }
private:
	std::vector<int>	m_arrValue;
};

class CCmdArgBool: public CCmdArg
{
public:
	CCmdArgBool(char cOptChar,
		const		char *szKeyword,
		const		bool bDefault,
		const		char* szDescription, 
		const		unsigned int syntax_flags = CCmdArg::isOPT)
		: CCmdArg(cOptChar, szKeyword, NULL, szDescription, syntax_flags), m_bValue(bDefault) {}
	virtual ~CCmdArgBool() {}		
	operator bool(void)  const { return  m_bValue; }
	virtual int SetValue(char *argv[])
	{
		m_bValue = !m_bValue;
		return 0;
	}
	virtual int ValueNeeded() { return 0; }
private:
	bool	m_bValue;
};

class CCmdArgFloat : public CCmdArg
{
public:
	CCmdArgFloat(char cOptChar, 
		  const       char *szKeyword, 
		const		 char *szValue, 
		const		 double lfDefault,
		const		 char* szDescription, 
				 unsigned int syntax_flags = CCmdArg::isOPTVALREQ)
		: CCmdArg(cOptChar, szKeyword, szValue, szDescription, syntax_flags), m_lfValue(lfDefault) {}
	operator double(void)  const { return  m_lfValue; }
	operator float(void)  const { return  (float)m_lfValue; }
	virtual int SetValue(char *argv[])
	{
		m_lfValue = atof(*argv++);
		return 0;
	}
private:
	double	m_lfValue;
};

class CCmdArgFloatList : public CCmdArg
{
public:
	CCmdArgFloatList(char cOptChar, 
		    const   char *szKeyword, 
			const   char *szValue, 
			double lfDefault,
			   int iFloatNumber,
			const   char* szDescription, 
			   unsigned int syntax_flags = CCmdArg::isOPTVALREQ)
		: CCmdArg(cOptChar, szKeyword, szValue, szDescription, syntax_flags)
	{
		m_arrValue.resize(iFloatNumber);
		for (int i = 0; i < iFloatNumber; i++)
			m_arrValue[i] = lfDefault;
	}
	virtual ~CCmdArgFloatList(){}
	double operator[] (int i) { return m_arrValue[i]; }
	virtual int SetValue(char *argv[])
	{
		for (unsigned int i = 0; i < m_arrValue.size(); i++)
		{
			m_arrValue[i] = atof(*argv++);
		}
		return 0;
	}
	virtual int ValueNeeded() { return (int)m_arrValue.size(); }
private:
	std::vector<double>	m_arrValue;
};

class CCmdArgString : public CCmdArg
{
public:
	CCmdArgString(char cOptChar, 
		         const  char *szKeyword, 
			const	  char *szValue, 
			const	  char *szDefault,
			const	  char* szDescription, 
				  unsigned int syntax_flags = CCmdArg::isOPTVALREQ)
		: CCmdArg(cOptChar, szKeyword, szValue, szDescription, syntax_flags), m_str(szDefault) {}
	operator const char *(void)  const { return  m_str.c_str(); }
	virtual ~CCmdArgString(){}
	virtual int SetValue(char *argv[])
	{
		m_str = *argv++;
		return 0;
	}
private:
	std::string m_str;
};

class CCmdArgStringList : public CCmdArg
{
public:
	CCmdArgStringList(char cOptChar, 
const		          char *szKeyword, 
const				  char *szValue, 
const				  char *szDefault,
				  int iStringNumber,
const				  char* szDescription, 
				  unsigned int syntax_flags = CCmdArg::isOPTVALREQ)
		: CCmdArg(cOptChar, szKeyword, szValue, szDescription, syntax_flags)
	{
		m_arrStr.resize(iStringNumber);
		for (int i = 0; i < iStringNumber; i++)
			m_arrStr[i] = szDefault;
	}
	virtual ~CCmdArgStringList(){}
	const char * operator[] (int i) { return m_arrStr[i].c_str(); }
	virtual int SetValue(char *argv[])
	{
		for (unsigned int i = 0; i < m_arrStr.size(); i++)
		{
			m_arrStr[i] = *argv++;
		}
		return 0;
	}
	virtual int ValueNeeded() { return (int)m_arrStr.size(); }
private:
#ifdef _UNICODE
	std::vector<std::wstring> m_arrStr;
#else
	std::vector<std::string> m_arrStr;
#endif
};

/**********************************************************
class CTaskItem: include a set of params
***********************************************************/
class CTaskItem
{
	std::vector<CCmdArg*> m_arrArgList;
public:
	const char*	m_szDescription;	//e.g. "detect folder: detect the face in the whole image folder
    const char*	m_szTaskName;		//e.g. "Folder"
	int		(*m_MethodFunction)();
	CTaskItem()
	{
		m_szDescription = m_szTaskName = NULL;
	}
	CTaskItem(const char* szDescription, const char* szTaskName, int (*pFun)(), int iCount, ...);
	int Size() { return (int)m_arrArgList.size(); }
	CCmdArg &GetAt(int i) { return *m_arrArgList[i]; }
	CCmdArg &operator[] (int i) { return *m_arrArgList[i]; }
};

CTaskItem::CTaskItem(const char* szDescription, const char* szTaskName, int (*pFun)(), int iCount, ...)
{
	m_szDescription = szDescription;
	m_szTaskName =  szTaskName;
	m_MethodFunction = pFun;
	
	va_list marker;
	va_start(marker, iCount);     /* Initialize variable arguments. */
	CCmdArg* pParam = NULL;
	for (int i = 0; i < iCount; i ++)
	{
		pParam = (CCmdArg*) va_arg(marker, void*);
		m_arrArgList.push_back(pParam);
	}
	va_end( marker );              /* Reset variable arguments. */
}

/**********************************************************
class CParser: Parse the arguments of .exe file
***********************************************************/
class CParser
{
public:
	CParser() { m_iTaskIndex = -1;};
	~CParser() {};
	CParser(const char* szProgramName, int iCount, ...);
	bool	Parse(int argc, char* argv[]);
	void	SetVersion(const char* szVersion) 
	{
		strcpy(m_szVersion, szVersion);
	}
	int		Run();
	const char * operator[] (int i) { return m_arrArgs[i]; }
    int   GetArgNum() { return (int) m_arrArgs.size();}
	void		DisplayAllTask();
	void		DisplayTaskDetails();
protected:
	bool		ParseArgs(int argc, char* argv[]);
	bool		ParseTask(int argc, char* argv[]);
	std::vector<CTaskItem*>	m_arrTaskItems;
	int			m_iTaskIndex;
	const char*		szAppName;
	char		m_szVersion[_MAX_PATH];
    vector<const char*> m_arrArgs;
};

CParser::CParser(const char* szProgramName, int iCount, ...)
{
	int i = 0;
	szAppName = szProgramName;
	va_list marker;
	va_start(marker, iCount);     /* Initialize variable arguments. */
	CTaskItem* pItem = NULL; 
	for (i = 0;i < iCount;	i ++)
	{
		pItem = (CTaskItem*) va_arg(marker, void*), 
		m_arrTaskItems.push_back(pItem);
	}
	va_end( marker );					/* Reset variable arguments.      */
	return ;
}

int CParser::Run()
{
	if (m_iTaskIndex < 0)
	{
		return false;
	}
	int (*fun)();
	fun = m_arrTaskItems[m_iTaskIndex]->m_MethodFunction;
	return fun();
}

void CParser::DisplayAllTask()
{
	printf("Usage: %s [tasktype]\n", szAppName);
	printf("\n");
	printf("Task type: %d\n", m_arrTaskItems.size());
	for (unsigned int i = 0; i < m_arrTaskItems.size(); i ++)
	{
		printf("        %-10s %s\n", m_arrTaskItems[i]->m_szTaskName, m_arrTaskItems[i]->m_szDescription);
	}
	printf("\n");
	printf("%s\n", m_szVersion);
}

void CParser::DisplayTaskDetails()
{
	char buf[1024];
	//printf("Task: %s\n", m_arrTaskItems[m_iTaskIndex]->m_szTaskName);
	printf(" %-10s %s\n", m_arrTaskItems[m_iTaskIndex]->m_szTaskName, m_arrTaskItems[m_iTaskIndex]->m_szDescription);
	printf("%s %s", szAppName, m_arrTaskItems[m_iTaskIndex]->m_szTaskName);
	for (int i = 0; i < m_arrTaskItems[m_iTaskIndex]->Size(); i++)
	{
		CCmdArg *pArg = &m_arrTaskItems[m_iTaskIndex]->GetAt(i);
		if (pArg->GetCharName() == 0)
			printf(" [-%s", pArg->GetKeywordName());
		else
			printf(" [-%c", pArg->GetCharName());
		if (pArg->GetValueName() == NULL)
			printf("]");
		else
			printf(" %s]", pArg->GetValueName());
	}
	printf("\n");
	printf("\n");
	printf("Options/Arguments:\n");
//	printf("=============================================\n");
	for (int i = 0 ;i < m_arrTaskItems[m_iTaskIndex]->Size(); i++)
	{
		CCmdArg *pArg = &m_arrTaskItems[m_iTaskIndex]->GetAt(i);
		if (pArg->GetKeywordName() == 0)
			sprintf(buf, "        -%c", pArg->GetCharName());
		else if (pArg->GetCharName() == 0)
			sprintf(buf, "        -%s", pArg->GetKeywordName());
		else
			sprintf(buf, "        -%c|%s", pArg->GetCharName(), pArg->GetKeywordName());
		if (pArg->GetValueName() != NULL)
		{
			strcat(buf, " ");
			strcat(buf, pArg->GetValueName());
		}
		printf("%-31s %s\n", buf, pArg->GetDescription());
	}
	//printf("---------------------------------------------\n");
	printf("\n");
	printf("%s\n", m_szVersion);
	return ;
}

bool CParser::ParseArgs(int argc, char* argv[])
{
	if (argc <= 2)
		return false;
	for( int i = 2; i < argc; i++)
    {
        if ( argv[i][0] != '-') 
		{
            //printf("Switch '-' is expected before %s\n\n", argv[i]);
			//return false;
            m_arrArgs.push_back(argv[i]);
            continue;
		}
		CCmdArg *pArg = NULL;
		int iItemIndex = 0;
		for (iItemIndex = 0; iItemIndex < m_arrTaskItems[m_iTaskIndex]->Size(); iItemIndex++)
		{
			pArg = &m_arrTaskItems[m_iTaskIndex]->GetAt(iItemIndex);
			if (argv[i][2] != 0 && pArg->GetKeywordName() != NULL)
			{
				//if (strcmp(argv[i] + 1, pArg->GetKeywordName(), strlen(argv[i]) - 1) == 0)
				if (strcmp(argv[i] + 1, pArg->GetKeywordName()) == 0)
				{
					break;
				}
			}
			else if (argv[i][2] == 0 && argv[i][1] == pArg->GetCharName())
			{
				break;
			}
		}
		if (iItemIndex < m_arrTaskItems[m_iTaskIndex]->Size())
		{
			if (i + pArg->ValueNeeded() >= argc)
			{
				char buf[50];

				if (pArg->GetKeywordName() == 0)
					sprintf(buf, "-%c", pArg->GetCharName());
				else if (pArg->GetCharName() == 0)
					sprintf(buf, "-%s", pArg->GetKeywordName());
				else
					sprintf(buf, "-%c|%s", pArg->GetCharName(), pArg->GetKeywordName());
				printf("Insufficient values for option %s, which requirs %d values\n\n", 
					buf, pArg->ValueNeeded());
				return false;
			}
			else
			{
				pArg->SetValue(argv + i + 1);
				i += pArg->ValueNeeded();
				pArg->SetAvailable();
			}
		}
		else
		{
			printf("Unknown option %s \n\n", argv[i-1]);
			return false;
		}
	}

	return true;
}

bool CParser::ParseTask(int argc, char* argv[])
{
	if (argc < 2)
	{
		return false;
	}
	for (unsigned int i = 0; i < m_arrTaskItems.size(); i ++)
	{
		if (strcmp(argv[1], m_arrTaskItems[i]->m_szTaskName) == 0)
		{
			m_iTaskIndex = i;
			break;
		}
	}
	return true;
}

bool CParser::Parse(int argc, char* argv[])
{
	if (ParseTask(argc, argv))
	{
		if (ParseArgs(argc, argv))
		{
			return true;
		}
		else
		{
			DisplayTaskDetails();
			return false;
		}
	}
	else 
	{
		printf("----");
		DisplayAllTask();
		printf("++++");
		return false;
	}
	
}

#endif //_PARSER_H_
