#ifndef _GENERICDATA_CPP_
#define _GENERICDATA_CPP_

#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "genericdata.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "math.h"

#include "CImg.h"
using namespace cimg_library;
using namespace std;

adjmatrix::adjmatrix(int n1, int n2) : n1(n1), n2(n2)
{
  mat = new int* [n1];
  for (int i = 0; i < n1; i ++)
  {
    mat[i] = new int [n2];
    for (int j = 0; j < n2; j ++) mat[i][j] = 0;
  }
}

adjmatrix::~adjmatrix()
{
  for (int i = 0; i < n1; i ++)
    delete [] mat[i];
  delete [] mat;
}

CGenericData::CGenericData()
{
  Configuration &config = Configuration::GetInstance();
  string name;
  string name2;

  if      (config.GetString("Program.mode") == "LEARNING")   
     { name = config.GetString("Data.trainFile"); name2 = config.GetString("Data.trainPairs"); }
  else if (config.GetString("Program.mode") == "PREDICTION") 
     { name = config.GetString("Data.testFile"); name2 = config.GetString("Data.testPairs"); }

  FILE* infile;
  infile = fopen(name.c_str(), "r");

  if (!infile)
  {
    cout << "Could not find " << name.c_str() << endl;
    exit(0);
  }

  char* corner_name  = new char [100];
  char* feature_name = new char [100];
  char* label_name   = new char [100];
  char* image_name   = new char [100];
  int class_label = 0;

  int ind = 0;
  _Kneighbours = config.GetInt("KNN.Kneighbour");
  category = config.GetInt("Data.category");
  dimOfFeature = config.GetInt("Data.dimOfFeature");
  dimOfWeight = 2 * dimOfFeature + 1;

  //map<int,vector<int>* > categories;
  double* featmax = new double [dimOfFeature];
  for (int i = 0; i < dimOfFeature; i ++) featmax[i] = -1;

  while (fscanf(infile, "%s %s %s %d", corner_name, feature_name, image_name, &class_label) == 4)
  {
    printf("%s %s %s %d\n", corner_name, feature_name, image_name, class_label);

    if (categories.find(class_label) == categories.end()) categories[class_label] = new vector<int>();
    categories[class_label]->push_back(ind);

    /// Store the class label...
    classes.push_back(class_label);

    /// Store the name of the image...
    char* iname = new char [100];
    memcpy(iname, image_name, 100*sizeof(char));
    imnames.push_back(iname);

    /// Store the corners...
    // Always store the dummy node as the first corner.
    vector<Point>* corn = new vector<Point>();
    Point dummynode(0,0);
    dummynode.dummy = 1;
    corn->push_back(dummynode);

    FILE* f = fopen(corner_name, "r");
    float x = 0;
    float y = 0;
    int seq = 0;
    int pnt = 0;
    char strbuff[32768];
    fgets(strbuff, 32768, f);
    while (fscanf(f, "%d %d %f %f", &seq, &pnt, &x, &y) == 4)
      corn->push_back(Point(x,y));
    corners.push_back(corn);
    fclose(f);
    shape_sizes.push_back((int) corn->size());

    /// Store the features...
    f = fopen(feature_name, "r");
    vector<Scalar*>* scf = new vector<Scalar*>();

    // The dummy node has no features...
    scf->push_back(NULL);
    // Not a very elegant way of reading the features, but it'll do for now :-)
    char* line = new char[10000];
    for (int i = 0; i < (int) corn->size() - 1; i ++)
    {
      Scalar* feat = new Scalar [dimOfFeature];
      fgets(line, 10000, f);
      char* s2 = line;

      for (int j = 0; j < dimOfFeature; j ++)
      {
        sscanf(s2, "%lf", &feat[j]);
        if (isnan(feat[j]))
        {
          printf("Failed to read feature.\n");
          exit(0);
        }
        if (feat[j] > featmax[j]) featmax[j] = feat[j];
        while (s2[0] != ' ' and s2[0] != '\n' and s2[0] != '\0' and s2[0] != '\t') s2 ++;
        if (s2[0] == ' ' or s2[0] == '\t') s2 ++;
      }

      scf->push_back(feat);
    }
    features.push_back(scf);
    fclose(f);
    delete [] line;


    if (class_label == category || category == 0) categoryinds.push_back(ind);
    ind ++;
  }
  for (int i = 0; i < (int) features.size(); i ++)
    for (int j = 1; j < (int) features[i]->size(); j ++)
      for (int k = 0; k < dimOfFeature; k ++)
        if (featmax[k] > 0) features[i]->at(j)[k] /= featmax[k];
  delete [] featmax;

  // Build the negative set. Assumes that each category has the same number of instances.
  negativeset = new vector<int> [(int) corners.size()];
  positiveset = new vector<int> [(int) corners.size()];

  for (map<int,vector<int>* >::iterator it = categories.begin(); it != categories.end(); it ++)
  {
    int cat = it->first;
    vector<int>* vec = it->second;

    for (int i = 0; i < (int) vec->size(); i ++)
      for (int j = 0; j < (int) vec->size(); j ++)
        if (i != j) positiveset[vec->at(i)].push_back(vec->at(j));

    for (int i = 0; i < (int) vec->size(); i ++)
    {
      for (map<int,vector<int>* >::iterator it2 = categories.begin(); it2 != categories.end(); it2 ++)
      {
        int cat2 = it2->first;
        vector<int>* vec2 = it2->second;

        if (cat != cat2) negativeset[vec->at(i)].push_back(vec2->at(i));
      }
    }
  }

  for (map<int,vector<int>* >::iterator it = categories.begin(); it != categories.end(); it ++)
    delete it->second;

  fclose(infile);
  infile = fopen(name2.c_str(), "r");

  // Now start reading the labels...
  int pair1;
  int pair2;
  fprintf(stderr, "start loading labeled matching from %s ", name2.c_str());
  while (fscanf(infile, "%s %d %d", label_name, &pair1, &pair2) == 3)
  {
    printf("%s %d %d\n", label_name, pair1, pair2);
    /// Store the correct labelling...
    FILE* f = fopen(label_name, "r");

    int n = (int) corners[pair1]->size();
    int m = (int) corners[pair2]->size();

    adjmatrix* adj = new adjmatrix(n, m);

    // Initially, everything matches to the dummy node.
    for (int i = 0; i < n; i ++) adj->mat[i][0] = 1;
    for (int i = 0; i < m; i ++) adj->mat[0][i] = 1;
    adj->mat[0][0] = 0;

    int c1 = 0;
    int c2 = 0;
    // Add one to the index, since zero actually corresponds to the dummy node.
    while (fscanf(f, "%d %d", &c1, &c2) == 2)
    {
   
      c1 ++;
      c2 ++;

      adj->mat[c1][0] = 0;
      adj->mat[0][c2] = 0;
      adj->mat[c1][c2] = 1;
    }

    correct[pair<int,int>(pair1, pair2)] = adj;

    adjmatrix* adj_trans = new adjmatrix(m, n);
    for (int i = 0; i < n; i ++)
      for (int j = 0; j < m; j ++)
        adj_trans->mat[j][i] = adj->mat[i][j];
    correct[pair<int,int>(pair2, pair1)] = adj_trans;
    fclose(f);
  }

  delete [] corner_name;
  delete [] feature_name;
  delete [] label_name;
  delete [] image_name;

  fclose(infile);

  _N = (int) categoryinds.size();
  _total = (int) corners.size();

  printf("Finished reading data., samples in training %d, total %d\n", _N, _total);
}

CGenericData::~CGenericData()
{
  for (vector<vector<Point>*>::iterator   it = corners.begin();  it != corners.end();  it ++) delete *it;
  for (vector<vector<Scalar*>*>::iterator it = features.begin(); it != features.end(); it ++)
  {
    for (vector<Scalar*>::iterator it2 = (*it)->begin(); it2 != (*it)->end(); it2 ++) delete [] *it2;
    delete *it;
  }
  for (map<pair<int,int>,adjmatrix*>::iterator it = correct.begin(); it != correct.end(); it ++)
    delete it->second;

  for (vector<char*>::iterator it = imnames.begin(); it != imnames.end(); it ++) delete [] *it;

  delete [] negativeset;
  delete [] positiveset;
}

#endif
