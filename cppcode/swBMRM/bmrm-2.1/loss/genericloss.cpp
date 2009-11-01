#ifndef _GENERICLOSS_CPP_
#define _GENERICLOSS_CPP_

#include "genericdata.hpp"
#include "genericloss.hpp"
#include "configuration.hpp"
#include <fstream>
#include <sstream>

#include "math.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <new>


//#include <omp.h>

#include "CImg.h"
using namespace cimg_library;

/**
 * Classification loss.
 */
Scalar loss1(int search_class, int scene_class)
{
  if (search_class == scene_class) 
       return 0;
  return 1;
}

/**
 * First-order feature vector.
 * During the first stage, this is just the linear features.
 * During the second stage, this is a scalar, corresponding to the inner product of the linear features with their weights from the first stage.
 */
void CGenericLoss::Phi1(Point p1, Point p2, Scalar* a, Scalar* b, Scalar* res)
{
  for (int i = 0; i < dimOfWeight ; i ++) 
      res[i] = 0;

  if (p1.dummy && p2.dummy)
     return;

  if (p1.dummy)
  {
    res[dimOfFeature] = 1;
    for (int i = 0; i < dimOfFeature; i ++) res[dimOfFeature + 1 + i] = fabs(b[i]);
  }
  else if (p2.dummy)
  {
    res[dimOfFeature] = 1;
    for (int i = 0; i < dimOfFeature; i ++) res[dimOfFeature+ 1 + i] = fabs(a[i]);
  }
  else
  {
    for (int i = 0; i < dimOfFeature; i ++)
      if (a[i] == 0 and b[i] == 0) res[i] = 0;
      else res[i] = (a[i] - b[i])*(a[i] - b[i]) / (fabs(a[i]) + fabs(b[i]));
  }
/*
  for (int i = 0; i < dimOfWeight; i ++)
    if (isnan(res[i]))
    {
      printf("Entry %d of feature vector is invalid.\n", i);
      exit(0);
    }
*/
}

void CGenericLoss::Phi(int n1, int n2, adjmatrix* matches, Scalar* res)
{ 
    for (int i = 0; i < dimOfWeight; i ++) 
        res[i] = 0;
    Scalar* phi = new Scalar[dimOfWeight];
    for (int i = 0; i < dimOfWeight; i ++) 
        res[i] = 0;
    for (int x = 0; x < matches->n1; x ++)
    {
      for (int y = 0; y < matches->n2; y ++)
      {
        if (matches->mat[x][y] == 1)
        {
          Phi1(_data->corners[n1]->at(x), _data->corners[n2]->at(y),
               _data->features[n1]->at(x), _data->features[n2]->at(y), phi);
          for (int i = 0; i < dimOfWeight; i ++) 
             res[i] -= phi[i];
        }
      }
    }
    delete[] phi;

}

void CGenericLoss::Phi(int n1, int* Knearest, adjmatrix** matches, Scalar* res)
{
  for (int i = 0; i < dimOfWeight; i ++) res[i] = 0;
  Scalar* phi = new Scalar[dimOfWeight];

  for (int k = 0; k < _Kneighbours; k ++)
  {
    int n2 = Knearest[k];
    for (int x = 0; x < matches[k]->n1; x ++)
      for (int y = 0; y < matches[k]->n2; y ++)
        if (matches[k]->mat[x][y] == 1)
        {
          Phi1(_data->corners[n1]->at(x), _data->corners[n2]->at(y),
               _data->features[n1]->at(x), _data->features[n2]->at(y), phi);
          for (int i = 0; i < dimOfWeight; i ++) 
             res[i] -= phi[i];
        }
  }
  delete [] phi;
}

/**
 * Dynamic programming solution to the matching problem.
 */
adjmatrix* dp(Scalar** matchcosts, int n, int m, Scalar* cost)
{
  Scalar** costmatrix = new Scalar*[n];
  int** dirmatrix = new int*[n];

  for (int i = 0; i < n; i ++)
  {
    costmatrix[i] = new Scalar[m];
    dirmatrix[i] = new int[m];
  }

  costmatrix[0][0] = 0;
  dirmatrix[0][0] = DIAGONAL;

  for (int i = 1; i < n; i ++)
  {
    costmatrix[i][0] = costmatrix[i-1][0] + matchcosts[i][0];
    dirmatrix[i][0] = UP;
  }

  for (int j = 1; j < m; j ++)
  {
    costmatrix[0][j] = costmatrix[0][j-1] + matchcosts[0][j];
    dirmatrix[0][j] = LEFT;
  }

  for (int i = 1; i < n; i ++)
  {
    for (int j = 1; j < m; j ++)
    {
      Scalar c_up = costmatrix[i-1][j] + matchcosts[i][0];
      Scalar c_left = costmatrix[i][j-1] + matchcosts[0][j];
      Scalar c_diag = costmatrix[i-1][j-1] + matchcosts[i][j];

      if (c_up < c_left and c_up < c_diag)
      {
        costmatrix[i][j] = c_up;
        dirmatrix[i][j] = UP;
      }
      else if (c_left < c_diag)
      {
        costmatrix[i][j] = c_left;
        dirmatrix[i][j] = LEFT;
      }
      else
      {
        costmatrix[i][j] = c_diag;
        dirmatrix[i][j] = DIAGONAL;
      }
    }
  }
  // Trace back through the arrows to get the adjacency matrix...
  adjmatrix* res = new adjmatrix(n,m);
  int cn = n - 1;
  int cm = m - 1;

  while (cn >= 0 and cm >= 0)
  {
    if (dirmatrix[cn][cm] == UP)
    {
      res->mat[cn][0] = 1;
      cn --;
    }
    else if (dirmatrix[cn][cm] == LEFT)
    {
      res->mat[0][cm] = 1;
      cm --;
    }
    else
    {
      res->mat[cn][cm] = 1;
      cn --;
      cm --;
    }
  }

  *cost = costmatrix[n-1][m-1];

  for (int i = 0; i < n; i ++)
  {
    delete [] costmatrix[i];
    delete [] dirmatrix[i];
  }
  delete [] costmatrix;
  delete [] dirmatrix;

  return res;
}

int compare(const void * a, const void * b)
{
  pair<Scalar,int> v1 = *(pair<Scalar,int>*) a;
  pair<Scalar,int> v2 = *(pair<Scalar,int>*) b;

  Scalar diff  = v1.first - v2.first;
  return diff < 0 ? -1 : (diff > 0 ? 1 : 0);
}

/**
 * Find the most similar shapes, and return their indices.
 */
int* CGenericLoss::assignment(int n1, Scalar* weights, adjmatrix** res_adjs, int test)
{
  int examps = _data->_total - 1;
  if (!test) examps = ((int) _data->positiveset[n1].size()) + ((int) _data->negativeset[n1].size());

  map<int,adjmatrix*> adjs;
  pair<Scalar,int>* costs = new pair<Scalar,int> [examps];
  Scalar* phi = new Scalar [dimOfWeight];

  int current = 0;

  // If test == 1, we will skip the example when (n1 == n2), so we must add 1 to the loop boundary.
  for (int id = 0; id < examps + test; id ++)
  {
    int n2 = id;
    if (!test)
    {
      if (id < (int) _data->positiveset[n1].size()) 
           n2 = _data->positiveset[n1][id];
      else 
           n2 = _data->negativeset[n1][id - (int) _data->positiveset[n1].size()];
    }
    if (n1 == n2) continue;

    int n = _data->shape_sizes[n1];
    int m = _data->shape_sizes[n2];

    Scalar totalcost = 0;
    Scalar** costmatrix = new Scalar*[n];
    for (int i = 0; i < n; i ++)
      costmatrix[i] = new Scalar[m];

    for (int i = 0; i < n; i ++)
      for (int j = 0; j < m; j ++)
      {
        Scalar cost = 0;
        Phi1(_data->corners[n1]->at(i), _data->corners[n2]->at(j),
             _data->features[n1]->at(i), _data->features[n2]->at(j), phi);
        if (weights != NULL)
          for (int k = 0; k < dimOfWeight; k ++) cost += weights[k]*phi[k];
        else
          for (int k = 0; k < dimOfWeight; k ++) cost += phi[k];
        costmatrix[i][j] = cost;
      }

    adjmatrix* a = dp(costmatrix,n,m,&totalcost);
    adjs[n2] = a;

    costs[current] = pair<Scalar,int>(totalcost,n2);

    if (!test and _data->classes[n1] != _data->classes[n2])
      costs[current].first -= 1.0 / _Kneighbours;

    for (int i = 0; i < n; i ++)
      delete [] costmatrix[i];
    delete [] costmatrix;
    current ++;
  }

  qsort(costs, examps, sizeof(pair<Scalar,int>), compare);

  int* res = new int[_Kneighbours];
  for (int k = 0; k < _Kneighbours; k ++)
  {
    res[k] = costs[k].second;
    res_adjs[k] = adjs[res[k]];
  }
  for (int k = _Kneighbours; k < examps; k ++)
    delete adjs[costs[k].second];

  delete [] costs;
  delete [] phi;

  return res;
}

/** Constructor. */
CGenericLoss::CGenericLoss(CModel* &model, CGenericData* &data)
  : CLoss(model, 1, data->dimOfWeight, data->bias()),
    _data(data)
{
  srand(0);
  Configuration &config = Configuration::GetInstance();

  verbosity = config.GetInt("Loss.verbosity");
  dimOfFeature = data->dimOfFeature;
  dimOfWeight = data->dimOfWeight;
  _Kneighbours = data->_Kneighbours;
  _LossType = config.GetString("Loss.Type");
}

/** Destructor. */
CGenericLoss::~CGenericLoss() {}

void CGenericLoss::Usage() {}
void CGenericLoss::ComputeLossAndGradient(Scalar& loss, TheMatrix& grad)
{
    if (_LossType == "KNN") 
        KNNLoss(loss, grad);
    if (_LossType == "AVG_K")
        AvgKLoss(loss, grad);
    return ;
}
/* Calculate the average matching score of top K samples from each class c, Cost_k
   Please note that the matching score to the same class uses the labeled match
    the matching score to other classes uses optimal match with current w */
//find top K results from class c that are most similar to sample n1

adjmatrix* CGenericLoss::GetDP(int n1, int n2, Scalar* weights, Scalar& cost)
{
    cost = 0;
    int n = _data->shape_sizes[n1];
    int m = _data->shape_sizes[n2];
    //create cost matrix
    Scalar* phi = new Scalar [dimOfWeight];
    Scalar** costmatrix = new Scalar*[n];
    for (int i = 0; i < n; i ++)
    {
        costmatrix[i] = new Scalar[m];
    }
    for (int i = 0; i < n; i ++)
    {
      for (int j = 0; j < m; j ++)
      {
        Scalar cost = 0;
        Phi1(_data->corners[n1]->at(i), _data->corners[n2]->at(j),
             _data->features[n1]->at(i), _data->features[n2]->at(j), phi);
        if (weights != NULL)
          for (int k = 0; k < dimOfWeight; k ++) cost += weights[k]*phi[k];
        else
          for (int k = 0; k < dimOfWeight; k ++) cost += phi[k];
        costmatrix[i][j] = cost;
      }
    }
    adjmatrix* a = dp(costmatrix, n, m, &cost);
    for (int i = 0; i < n; i ++)
      delete [] costmatrix[i];
    delete [] costmatrix;
    delete[] phi;
    return a;
}
int* CGenericLoss::findTopKLabeled(int n1, int c, Scalar* weights, adjmatrix** res_adjs, Scalar& avgcost)
{
  vector<int> testsamples;
  for (int i = 0; i < _data->classes.size(); i ++)
  {
      if (_data->classes[i] == c)
          testsamples.push_back(i);
  } 
  int examps = (int)testsamples.size();
  map<int,adjmatrix*> adjs;
  pair<Scalar,int>* costs = new pair<Scalar,int> [examps];
  Scalar* phi = new Scalar [dimOfWeight];

  int current = 0;
  //assert((int) testsamples.size() == 11);
  //printf ("%d ", testsamples.size());
  for (int i = 0; i < (int)testsamples.size(); i ++)
  {
      int n2 = testsamples[i];
      if (n1 == n2) continue;
      adjmatrix* a = _data->correct[pair<int,int>(n1, n2)];
   //   printf("%d,%d: %d %d\n", n1, n2, a->n1, a->n2);
      assert(a);
      Phi(n1, n2, a,  phi);
      Scalar cost = 0.0;
      for (int s = 0; s < dimOfWeight; s ++ ) cost += weights[s] * phi[s];
      adjs[n2] = a;
      costs[current] = pair<Scalar,int>(cost, n2);
    //  printf("add %d %f\n", n2, cost);
      current ++;
  }

  qsort(costs, current , sizeof(pair<Scalar,int>), compare);
  avgcost = 0;
  int* res = new int[_Kneighbours];
  for (int k = 0; k < _Kneighbours; k ++)
  {
      avgcost += costs[k].first;
      Scalar f = 0; 
 //      for (int ss = 0; ss < _Kneighbours; ss ++ ) f +=  
      res[k] = costs[k].second;
      res_adjs[k] = adjs[res[k]];
     // printf("sorted %d %d %f %d\n", k, res[k], costs[k].first, res_adjs[k]->n1);
  }
  delete [] phi;
  delete[] costs;
  return res;
}

int* CGenericLoss::findTopK(int n1, int c, Scalar* weights, adjmatrix** res_adjs, Scalar& avgcost)
{
  vector<int> testsamples;
  for (int i = 0; i < _data->classes.size(); i ++)
  {
      if (_data->classes[i] == c)
          testsamples.push_back(i);
  } 
  int examps = (int)testsamples.size();
  map<int,adjmatrix*> adjs;
  pair<Scalar,int>* costs = new pair<Scalar,int> [examps];
  Scalar* phi = new Scalar [dimOfWeight];

  int current = 0;

  for (int i = 0; i < (int)testsamples.size(); i ++)
  {
      int n2 = testsamples[i];
      Scalar cost;
      adjmatrix* a = GetDP(n1, n2, weights, cost);
      adjs[n2] = a;
      costs[current] = pair<Scalar,int>(cost, n2);
      current ++;
  }

  qsort(costs, current, sizeof(pair<Scalar,int>), compare);
  avgcost = 0;

  int* res = new int[_Kneighbours];
  for (int k = 0; k < _Kneighbours; k ++)
  {
      avgcost += costs[k].first;
      res[k] = costs[k].second;
      res_adjs[k] = adjs[res[k]];
  }
  for (int k = _Kneighbours; k < examps; k ++)
  {
      delete adjs[costs[k].second];
  }
  delete [] costs;
  delete [] phi;
  return res;
}


void CGenericLoss::AvgKLoss(Scalar& loss, TheMatrix& grad)
{
  TheMatrix &w = _model->GetW();
  loss = 0;
  grad.Zero();

  Scalar* dat = w.Data();
  Scalar* raw_g = grad.Data();
  for (int i = 0; i < dimOfWeight; i ++) 
      raw_g[i] = 0;  
  Scalar* lossi = new Scalar [_data->_N];
  Scalar** raw_gi = new Scalar* [_data->_N];

  int i;
  
  for (i = 0; i < _data->_N; i ++)
  {
    int n1 = _data->categoryinds[i]; // n1 is the class of sample i 
    int c1 = _data->classes[n1]; //to determn 
    Scalar* resy = new Scalar [dimOfWeight]; // phi
    Scalar* resybar = new Scalar [dimOfWeight]; // phi
    Scalar  mincost = 1e10; 
    Scalar  samecost = 1e10;
    //find the worst negative category
    for (int c = 1; c < (int)_data->categories.size() + 1; c ++)
    {
        if (c != c1) 
        {
            Scalar avgcost;
            adjmatrix** ybar_labelling = new adjmatrix* [_Kneighbours];
            int* ybar = findTopK(i, c, dat, ybar_labelling, avgcost); 
            if (mincost > avgcost)
            {
                Phi(i, ybar, ybar_labelling, resybar);
                mincost = avgcost;
            }
            //release temp memory
            for (int k = 0; k < _Kneighbours; k ++) 
               delete ybar_labelling[k];
            delete [] ybar;
            delete [] ybar_labelling;
        }
    }
    //find the best positive category
    adjmatrix** correct_adjs = new adjmatrix* [_Kneighbours];
    int* correct_neighbours = findTopKLabeled(n1, c1, dat, correct_adjs, samecost); 

    printf(".");
    fflush(stdout);

    Phi(n1, correct_neighbours, correct_adjs, resy);

    Scalar inp = 0;
    Scalar cst1 = 0;
    Scalar cst2 = 0;
    for (int j = 0; j < dimOfWeight; j ++)
    {
      inp += dat[j]*(resybar[j]-resy[j]);
      cst1 += dat[j] * resybar[j];
      cst2 += dat[j] * resy[j];
    }
    Scalar labloss = 1; 
    raw_gi[i] = new Scalar [dimOfWeight];
    
    if (inp + labloss > 0)
    {
        lossi[i] = inp + labloss;
        for (int j = 0; j < dimOfWeight; j ++)
        {
            raw_gi[i][j] = (1.0/_data->_N)*(resybar[j]-resy[j]);
        }
    }
    else for (int j = 0; j < dimOfWeight; j ++)
    {
        lossi[i] = 0; 
       // raw_gi[i][j] = (1.0/_data->_N)*(resybar[j]-resy[j]);
        raw_gi[i][j] = 0; 
    }

    printf("cost: same cls %f vs diff clas %f, loss %f = 1 + %f, %f, %f\n", samecost, mincost, lossi[i], inp, cst2, cst1);
    delete [] resy;
    delete [] resybar;

    delete [] correct_neighbours;
    delete [] correct_adjs;
  }
  printf("\r");

//   for (int i = 0; i < Weights; i ++) cout << dat[i] << endl;
//   cout << endl;

  for (int j = 0; j < _data->_N; j ++)
  {
    loss += lossi[j];
    for (int k = 0; k < dimOfWeight; k ++)
      raw_g[k] += raw_gi[j][k];

    delete [] raw_gi[j];
  }
  delete [] raw_gi;
  delete [] lossi;
  loss = loss/_data->_N;
  printf("total loss = %f \n", loss);
}
void CGenericLoss::KNNLoss(Scalar& loss, TheMatrix& grad)
{
  TheMatrix &w = _model->GetW();
  loss = 0;
  grad.Zero();

  Scalar* dat = w.Data();
  Scalar* raw_g = grad.Data();

  //Configuration &config = Configuration::GetInstance();

  //Scalar lambda = config.GetDouble("BMRM.lambda");

  //Evaluate(w);

  for (int i = 0; i < dimOfWeight; i ++) 
      raw_g[i] = 0;   //lambda*dat[i];

  Scalar* lossi = new Scalar [_data->_N];
  Scalar** raw_gi = new Scalar* [_data->_N];

  int i;
  //#pragma omp parallel for
  for (i = 0; i < _data->_N; i ++)
  {
    int n1 = _data->categoryinds[i];

    Scalar* resy = new Scalar [dimOfWeight];
    Scalar* resybar = new Scalar [dimOfWeight];

    adjmatrix** ybar_labelling = new adjmatrix* [_Kneighbours];

    int* correct_neighbours = new int [_Kneighbours];
    adjmatrix** correct_adjs = new adjmatrix* [_Kneighbours];

    int* ybar = assignment(n1, dat, ybar_labelling, 0);
/*
    if (_Kneighbours < _data->positiveset[n1].size())
    {
      printf("Too few neighbours... %d < %d\n", _data->positiveset[n1].size(), _Kneighbours);
      exit(0);
    }
*/
    for (int id = 0; id < _Kneighbours; id ++)
    {
      int n2 = _data->positiveset[n1][id];
      correct_neighbours[id] = n2;
      correct_adjs[id] = _data->correct[pair<int,int>(n1,n2)];
    }

    printf(".");
    fflush(stdout);

    Phi(n1, correct_neighbours, correct_adjs, resy);
    Phi(n1, ybar, ybar_labelling, resybar);

    Scalar inp = 0;
    for (int j = 0; j < dimOfWeight; j ++)
      inp += dat[j]*(resybar[j]-resy[j]);

    Scalar labloss = LabelLoss(n1, ybar);
    //if (inp + labloss > 0)
      lossi[i] = labloss;

    raw_gi[i] = new Scalar [dimOfWeight];
    for (int j = 0; j < dimOfWeight; j ++)
    {
      //if (inp + labloss > 0)
      {
        lossi[i] += dat[j]*(resybar[j]-resy[j]);
        raw_gi[i][j] = (1.0/_data->_N)*(resybar[j]-resy[j]);
      }
    }

    delete [] resy;
    delete [] resybar;

    for (int k = 0; k < _Kneighbours; k ++) delete ybar_labelling[k];

    delete [] ybar;
    delete [] ybar_labelling;
    delete [] correct_neighbours;
    delete [] correct_adjs;
  }
  printf("\r");

//   for (int i = 0; i < Weights; i ++) cout << dat[i] << endl;
//   cout << endl;

  for (int j = 0; j < _data->_N; j ++)
  {
    loss += lossi[j];
    for (int k = 0; k < dimOfWeight; k ++)
      raw_g[k] += raw_gi[j][k];

    delete [] raw_gi[j];
  }
  delete [] raw_gi;
  delete [] lossi;

  loss = loss/_data->_N;
}

void CGenericLoss::Predict(CModel *model)
{
  Evaluate(model);
}

Scalar mean(Scalar* data, int n)
{
  Scalar total = 0;
  for (int i = 0; i < n; i ++) total += data[i];
  return total/n;
}

Scalar sde(Scalar* data, int n)
{
  Scalar avg = 0;
  Scalar avg2 = 0;

  for (int i = 0; i < n; i ++)
  {
    avg += data[i];
    avg2 += data[i]*data[i];
  }
  avg /= n;
  avg2 /= n;

  Scalar sd = 0;
  if (avg2 - avg*avg > 0) sd = sqrt(avg2 - avg*avg);
  return sd / sqrt((float) n);
}

void saveimage(char* filename, char* im1, char* im2, vector<Point>* corn1, vector<Point>* corn2, adjmatrix* match, char* tname)
{
  FILE* txt = fopen(tname, "w");
  int nodes1 = (int) corn1->size();
  int nodes2 = (int) corn2->size();

  CImg<unsigned char> image_lbw(im1);
  CImg<unsigned char> image_rbw(im2);

  CImg<unsigned char> image_l(image_lbw.dimx(), image_lbw.dimy(), 1, 3, 0);
  CImg<unsigned char> image_r(image_rbw.dimx(), image_rbw.dimy(), 1, 3, 0);

  int xoffl = 0;
  int yoffl = 0;
  int xoffr = image_l.dimx() + 5;
  int yoffr = 0;

  int diff = image_l.dimy() - image_r.dimy();
  int dimy = 0;
  if (diff > 0) { yoffr = diff/2; dimy = image_l.dimy(); }
  else          { yoffl = -diff/2; dimy = image_r.dimy(); }

  CImg<unsigned char> image_out(xoffr + image_r.dimx(), dimy, 1, 3, 0);
  //const unsigned char red[]   = { 255,0,0 };
  //const unsigned char green[] = { 0,255,0 };
  //const unsigned char blue[]  = { 0,0,255 };

  for (int y = 0; y < image_l.dimy(); y ++)
    for (int x = 0; x < image_l.dimx(); x ++)
      for (int c = 0; c < 3; c ++)
        image_l(x,y,c) = image_lbw(x,y);

  for (int y = 0; y < image_r.dimy(); y ++)
    for (int x = 0; x < image_r.dimx(); x ++)
      for (int c = 0; c < 3; c ++)
        image_r(x,y,c) = image_rbw(x,y);

  fprintf(txt, "Left image:\n");
  for (int i = 1; i < nodes1; i ++)
  {
    int plx = (int) corn1->at(i).x;
    int ply = (int) corn1->at(i).y;
    //image_l.draw_triangle(plx, ply-2, plx+2, ply+2, plx-2, ply+2,blue);
    fprintf(txt, "%d %d\n", plx, ply);
  }

  fprintf(txt, "Right image:\n");
  for (int i = 1; i < nodes2; i ++)
  {
    int prx = (int) corn2->at(i).x;
    int pry = (int) corn2->at(i).y;
    //image_r.draw_triangle(prx, pry-2, prx+2, pry+2, prx-2, pry+2,blue);
    fprintf(txt, "%d %d\n", prx+xoffr, pry+yoffr);
  }

  for (int y = 0; y < image_l.dimy(); y ++)
    for (int x = 0; x < image_l.dimx(); x ++)
      for (int c = 0; c < 3; c ++)
        image_out(x+xoffl,y+yoffl,c) = image_l(x,y,c);

  for (int y = 0; y < image_r.dimy(); y ++)
    for (int x = 0; x < image_r.dimx(); x ++)
      for (int c = 0; c < 3; c ++)
        image_out(x+xoffr,y+yoffr,c) = image_r(x,y,c);


  fprintf(txt, "Lines:\n");
  for (int i = 1; i < nodes1; i ++)
  {
    int plx = (int) corn1->at(i).x + xoffl;
    int ply = (int) corn1->at(i).y + yoffl;
    for (int j = 1; j < nodes2; j ++)
      if (match->mat[i][j] == 1)
      {
        int prx = (int) corn2->at(j).x + xoffr;
        int pry = (int) corn2->at(j).y + yoffr;

        //image_out.draw_line(plx,ply,prx,pry,blue);
        fprintf(txt, "%d %d %d %d (%d %d)\n", plx, ply, prx, pry, i, j);
      }
  }
  image_out.save(filename);
  fclose(txt);
}

/**
 * Compare the performance of the model after learning to the model before learning.
 * The performance before learning is found using a constant weight vector.
 */
void CGenericLoss::Evaluate(CModel *model)
{
  //Configuration &config = Configuration::GetInstance();
  TheMatrix &w = _model->GetW();

  Scalar* dat = w.Data();

  Scalar* lossn = new Scalar [_data->_N];
  Scalar* lossw = new Scalar [_data->_N];

  //#pragma omp parallel for
  for (int i = 0; i < _data->_N; i ++)
  {
    Scalar ll = 0;

    int n1 = _data->categoryinds[i];

    adjmatrix** adjs1 = new adjmatrix* [_Kneighbours];
    int* ybar1 = assignment(n1, dat, adjs1, 1);
    ll = LabelLoss(n1, ybar1); lossw[i] = ll;
    cout << n1 << ": Label loss(with weight) = " << ll << endl;

    adjmatrix** adjs0 = new adjmatrix* [_Kneighbours];
    int* ybar0 = assignment(n1, NULL, adjs0, 1);
    ll = LabelLoss(n1, ybar0); lossn[i] = ll;
    cout << n1 << ": Label loss(w/o weight) = " << ll << endl;

    char* name0 = new char [50]; char* name0_txt = new char [50];
    char* name1 = new char [50]; char* name1_txt = new char [50];
    int n2_0 = ybar0[0];
    int n2_1 = ybar1[0];
    sprintf(name0, "matches/noweight_%d_%d.jpg", n1, n2_0); sprintf(name0_txt, "matches/noweight_%d_%d.txt", n1, n2_0);
    sprintf(name1, "matches/weight_%d_%d.jpg", n1, n2_1); sprintf(name1_txt, "matches/weight_%d_%d.txt", n1, n2_1);

    //saveimage(name0, _data->imnames[n1], _data->imnames[n2_0], _data->corners[n1], _data->corners[n2_0], adjs0[0], name0_txt);
    //saveimage(name1, _data->imnames[n1], _data->imnames[n2_1], _data->corners[n1], _data->corners[n2_1], adjs1[0], name1_txt);

    delete [] name0; delete [] name0_txt;
    delete [] name1; delete [] name1_txt;

    delete [] ybar1;
    delete [] ybar0;
    for (int k = 0; k < _Kneighbours; k ++)
    {
      delete adjs1[k];
      delete adjs0[k];
    }
    delete [] adjs1;
    delete [] adjs0;
  }

  cout << "Error with weights: " << mean(lossw, _data->_N) << endl;
  cout << "Error without weights: " << mean(lossn, _data->_N) << endl;

  FILE* res_temp = fopen("res_temp.txt", "w");
  fprintf(res_temp, "%f %f %f %f\n", mean(lossw, _data->_N), sde(lossw, _data->_N), mean(lossn, _data->_N), sde(lossn, _data->_N));
  fclose(res_temp);

  delete [] lossw;
  delete [] lossn;
}

/** Loss function. */
Scalar CGenericLoss::LabelLoss(int search_ind, int* scene_inds)
{
  Scalar loss = 0;

//   printf("%d (%d)\t%d (%d)\n", search_ind, _data->classes[search_ind], scene_inds[0], _data->classes[scene_inds[0]]);
//   for (int i = 1; i < Kneighbours; i ++)
//     printf("\t%d (%d)\n", scene_inds[i], _data->classes[scene_inds[i]]);

  for (int k = 0; k < _Kneighbours; k ++)
    loss += loss1(_data->classes[search_ind], _data->classes[scene_inds[k]]) / _Kneighbours;
  return loss;
}

#endif
