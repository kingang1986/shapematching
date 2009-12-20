#include "stdio.h"
#include <iostream>
#include <fstream>
#include <string>
#include "SW.h"
#include "lp.h"
#include <glpk.h>
#include "cmdline.h"
#include "app.h"
#include "shapecontext.h"

int Train();
int Match();
int Test();
int Convert();
int Convert2();

CCmdArgInt arg_maxIteration = CCmdArgInt('i', "iteration", "20", 20, "maximal number of iteration");
CCmdArgInt arg_knn = CCmdArgInt('k', "knn", "3", 3, "K of the KNN, default 3");
CCmdArgInt arg_class = CCmdArgInt('C', "class", "9", 9, "number of class, default 9");
CCmdArgInt arg_maxStep = CCmdArgInt('t', "step", "10", 10, "maximal number of step in each iteration");
CCmdArgInt arg_minNewConstraint = CCmdArgInt('n', "NewConstraint", "20", 20, "minimal number of new constraints");
CCmdArgFloat arg_C = CCmdArgFloat('c', "c", "1000", 1000, "regularization between W and slack factors, f = w + C * slack");
CCmdArgFloat arg_e = CCmdArgFloat('e', "epsilon", "0.001", 0.001, "epsilon");
CCmdArgBool arg_binary = CCmdArgBool('b', "binary", false, "binary file/default true");
CCmdArgString arg_datafile = CCmdArgString('d', "datafile", "trainsample", "trainsample", "train files");
CCmdArgString arg_patternfile = CCmdArgString('t', "pattern", "train pattern", "train pattern", "pattern files");
CCmdArgString arg_modelfile = CCmdArgString('m', "modelfile", "model", "", "initial model file");
CCmdArgString arg_outputfile = CCmdArgString('o', "outputfile", "outputfile", "outputfile", "output model file");
CCmdArgString arg_path = CCmdArgString('p', "path", "", "", "data path, default empty");
CCmdArgString arg_input = CCmdArgString('i', "input", "input", "input", "inputfile");

CTaskItem itm_train = CTaskItem( "build models using structural linear programming",
	"train", Train, 11,
	&arg_datafile,
	&arg_patternfile,
	&arg_modelfile,
	&arg_path,
	&arg_maxIteration,
    &arg_knn,
//    &arg_class,
	&arg_maxStep,
	&arg_C,
	&arg_e,
	&arg_outputfile,
	&arg_minNewConstraint,
	&arg_binary);
CTaskItem itm_test = CTaskItem( "test models using structural linear programming",
        "test", Test, 6,
        &arg_datafile,
        &arg_modelfile,
        &arg_knn,
        &arg_path,
        &arg_outputfile,
        &arg_binary);
CTaskItem itm_match = CTaskItem( "match two shapes: usage: match [options] shape1.mss shape2.mss",
	"match", Match, 2,
	&arg_modelfile,
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

CParser parser = CParser("lpstralign",5, &itm_train, &itm_test, &itm_match, &itm_convert, &itm_convert2);


int Convert()
{
	CSetOfSeq ss;
	ss.Load(arg_input);
//	ss.SaveBinary(arg_outputfile);
	return 1;
}
int Convert2()
{
	CSetOfSeq ss;
	//ss.LoadSSBinary(arg_input);
	//ss.SaveSSBinary(arg_outputfile);
	return 1;
}
int Train()
{
    CPattern::m_iTopK = arg_knn;
    CPattern::m_iTotalClass = arg_class;
	CStructureLearning s1;
	//int Init(char* datafile, char* szpath, char* modelfile);
	s1.Init(arg_datafile, arg_patternfile, arg_path, arg_modelfile, (bool)arg_binary);
	s1.m_iMaxIteration = arg_maxIteration;
	s1.m_iMaxStep = arg_maxStep;
	s1.m_fC = arg_C;
	s1.m_fEpsilon = arg_e;
	s1.m_iMinNewConstraint = arg_minNewConstraint;
	s1.Learn(arg_outputfile);
	return 0;
}

int Match()
{
    if (parser.GetArgNum() != 2)
    {
       parser.DisplayTaskDetails();    
       return 0;
    }       
    match2shapes(parser[0], parser[1], arg_modelfile, arg_binary);
    return 0;
}
int Test()
{
/*
    CPattern::m_iTopK = arg_knn;
    CModel m;
    CSample sample;
    Sample.m_bBinaryData = arg_binary;
    sample.LoadSample(arg_datafile);
    m.Read(arg_modelfile);
    sample.AlignSamples(&m, arg_outputfile);
*/
	return 0;
}
int main( int argc, char* argv[] )
{

	parser.SetVersion("---------------------------------------\nLongbin Chen, longbinc@yahoo.com, 11/29/2009");
	if (parser.Parse(argc,argv))
	{
		parser.Run();
	}
	return 1;

}

