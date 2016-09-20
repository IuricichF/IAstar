#ifndef VERTEX_H
#define VERTEX_H

#include <vector>
#include <map>

#include "topsimplex.h"

using namespace std;

class Vertex
{

private:
    map<int, vector<int> > partial_coboundaryTop;
    vector<float> coordinates;

public:
    Vertex();
    Vertex(vector<float>);

    void changeCoordinate(float,int);
    vector<float>& getCoordinates();
    float getCoordinate(int) const;

    int getCoboundaryTopNum();
    int getCoboundaryTopNum(int);
    int getCoboundaryMaxDim();

    void addPartialCoboundaryTop(int dim, int index);
    map<int, vector<int> >& getPartialCoboundaryTopRelations();
    vector<int>* getPartialCoboundaryTop(int);

    float euclideanDistance(Vertex& v);


};

#endif // VERTEX_H
