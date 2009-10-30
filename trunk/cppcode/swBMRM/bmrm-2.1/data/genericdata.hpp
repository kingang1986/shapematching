#ifndef _GENERICDATA_HPP_
#define _GENERICDATA_HPP_

#include "common.hpp"
#include "info.hpp"
#include "sml.hpp"
#include "data.hpp"

#include <vector>
#include <map>


//#define NODE_WEIGHTS_HALF 78
//#define NODE_WEIGHTS (2*NODE_WEIGHTS_HALF + 1)

//#define Kneighbours 4

using namespace std;

enum { UP, LEFT, DIAGONAL };

/**
 * An adjacency matrix. Just a two-dimensional array of ints.
 */
class adjmatrix
{
  public:
    adjmatrix(int n1, int n2);
    ~adjmatrix();

    int** mat;
    int n1;
    int n2;
};

/**
 * A two-dimensional point. Can be a `dummy' or `articulation' point if these variables are set to 1.
 */
class Point
{
  public:
    Point() : x(0), y(0) {};
    Point(Scalar x, Scalar y) : x(x), y(y) {dummy = 0; articulation = 0;};
    ~Point() {};

    Point(const Point& p)
    {
      x = p.x;
      y = p.y;
      dummy = p.dummy;
      articulation = p.articulation;
    }

    Scalar dist(Point p);

    Scalar x;
    Scalar y;

    // Features are stored in the data container...
    //int* features;

    int dummy;
    int articulation;
};

class CGenericData: public CData
{
  public:
    // constructor
    CGenericData();
    virtual ~CGenericData();

    // Present for *all* examples.
    vector<vector<Point>*> corners;                   // The points in the scene (plus a dummy point).
    vector<vector<Scalar*>*> features;                // The features for the points.
    vector<int> classes;                              // The class of each point set.

    int category;
    vector<int> categoryinds;
    int _total;
    unsigned int _Kneighbours;

    // Present only for some pairs of training examples.
    vector<int> shape_sizes; // The size of each point set.
    vector<int>* negativeset;
    vector<int>* positiveset;
    map<pair<int,int>,adjmatrix*> correct; // Labelled ordering of the points in this set.
    map<int,vector<int>*> categories;

    // Information about the images, for output.
    vector<char*> imnames;
    

    // Various other information, required by BMRM...
    int _N;

    unsigned int numOfExample;
    unsigned int numOfAllExample;
    unsigned int dimOfFeature;
    unsigned int dimOfWeight; 

    virtual bool bias(void) const { return false; }
    virtual unsigned int dim(void) const
    {
      return dimOfWeight; 
    }

    virtual unsigned int slice_size(void) const { return _N; }
    virtual unsigned int size(void) const { return _N; }

  //protected:
    
 

  private:
    CGenericData(const CGenericData&);
};

#endif
