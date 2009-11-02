#include "stdio.h"
#include <iostream>
#include <fstream>
#include <string>
#include "Data.h"
#include "SW.h"
#include "lp.h"
#include <glpk.h>
#include "cmdline.h"

int Train();
int Test();
int Convert();
int Convert2();

CCmdArgInt arg_maxIteration = CCmdArgInt('i', "iteration", "20", 20, "maximal number of iteration");
CCmdArgInt arg_knn = CCmdArgInt('i', "knn", "3", 3, "K of the KNN, default 3");
CCmdArgInt arg_maxStep = CCmdArgInt('t', "step", "10", 10, "maximal number of step in each iteration");
CCmdArgInt arg_minNewConstraint = CCmdArgInt('n', "NewConstraint", "20", 20, "minimal number of new constraints");
CCmdArgFloat arg_C = CCmdArgFloat('c', "c", "1000", 1000, "regularization between W and slack factors, f = w + C * slack");
CCmdArgFloat arg_e = CCmdArgFloat('e', "epsilon", "0.001", 0.001, "epsilon");
CCmdArgBool arg_binary = CCmdArgBool('b', "binary", false, "binary file/default true");
CCmdArgString arg_datafile = CCmdArgString('d', "datafile", "trainsample", "trainsample", "train files");
CCmdArgString arg_modelfile = CCmdArgString('m', "modelfile", "model", "model", "initial model file");
CCmdArgString arg_outputfile = CCmdArgString('o', "outputfile", "outputfile", "outputfile", "output model file");
CCmdArgString arg_path = CCmdArgString('p', "path", "", "", "data path");
CCmdArgString arg_input = CCmdArgString('i', "input", "input", "input", "inputfile");

CTaskItem itm_train = CTaskItem( "build models using structural linear programming",
	"train", Train, 10,
	&arg_datafile,
	&arg_modelfile,
	&arg_path,
	&arg_maxIteration,
	&arg_maxStep,
	&arg_C,
	&arg_e,
	&arg_outputfile,
	&arg_minNewConstraint,
	&arg_binary);
CTaskItem itm_test = CTaskItem( "test models using structural linear programming",
        "test", Test, 5,
        &arg_datafile,
        &arg_modelfile,
        &arg_path,
        &arg_outputfile,
        &arg_binary);

CTaskItem itm_convert = CTaskItem( "Convert txt file to binary file",
	"convert", Convert, 2,
	&arg_input,
	&arg_outputfile
	);
CTaskItem itm_convert2 = CTaskItem( "Convert binary file to text file",
	"convert2", Convert2, 2,
	&arg_input,
	&arg_outputfile
	);

CParser parser = CParser("lpstralign",4, &itm_train, &itm_test, &itm_convert, &itm_convert2);


int Convert()
{
	CSetOfSeq ss;
	ss.LoadSS(arg_input);
	ss.SaveSSBinary(arg_outputfile);
	return 1;
}
int Convert2()
{
	CSetOfSeq ss;
	ss.LoadSSBinary(arg_input);
	ss.SaveSSBinary(arg_outputfile);
	return 1;
}
int Train()
{
    CPattern::m_iTopK = arg_knn;
	CStructureLearning s1;
	//int Init(char* datafile, char* szpath, char* modelfile);
	s1.Init(arg_datafile, arg_path, arg_modelfile, (bool)arg_binary);
	s1.m_iMaxIteration = arg_maxIteration;
	s1.m_iMaxStep = arg_maxStep;
	s1.m_fC = arg_C;
	s1.m_fEpsilon = arg_e;
	s1.m_iMinNewConstraint = arg_minNewConstraint;
	s1.Learn(arg_outputfile);
	return 0;
}
int Test()
{
    CPattern::m_iTopK = arg_knn;
       CModel m;
       CSample sample;
       sample.m_bBinaryData = arg_binary;
       sample.LoadSample(arg_datafile);
       m.Read(arg_modelfile);
       sample.AlignSamples(&m, arg_outputfile);
	return 0;
}
int main( int argc, char* argv[] )
{

	parser.SetVersion("---------------------------------------\nLongbin Chen, longbinc@yahoo.com, 10/29/2009");
	if (parser.Parse(argc,argv))
	{
		parser.Run();
	}
	return 1;

}

