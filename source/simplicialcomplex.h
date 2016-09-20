#ifndef SIMPLICIALCOMPLEX_H
#define SIMPLICIALCOMPLEX_H

#include <vector>
#include <forward_list>
#include <queue>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <list>
#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "topsimplex.h"
#include "vertex.h"



class SimplicialComplex
{


protected:

    map<int,int> realIndex;
    vector<vector<TopSimplex> > topSimplexes;
    vector<Vertex> vertices;
    vector<forward_list<int> > adjRelations;

    //EmbeddingSpace == number of coordinates for each vertex
    //SimplicialComplexSpace == dimension of the biggest top simplex

public:
    SimplicialComplex();

    //------------------------------------------------------------------------------------------
    //functions for reding input files. Supported formats [OFF,GMV]
    void readOFF(char* file); //read OFF file format
    void readTS(char* file); //read TS  file format
    void readGMV(char* file); //read GMV file format
    void readSimplexes(char* file); //TODO //set of simplexes. Vertices have no coordinates
    void readPoints(char* name_in_file, float threshold, uint dimension); //set of points on which a Vietoris-Rips complex is computed
    void readGraph(char* name_in_file); //set of arcs of a graph based on which clicks are computed
    void readIA(char* file); //read IA file format

    void saveIA(char* file); // write IA file format

    int getVerticesNum();
    //return: the number of vertices

    int getTopSimplexesNum(int d);
    //return: the number of top-simplexes of dimension d

    set<int> getTopSimplexesSet();
    //return: return the number of top-simplexes of each dimension (dimensions having 0 top-simplexes are ignored)

    Vertex& getVertex(int i);
    //return: vertex of index i

    vector<TopSimplex>& getTopSimplices(int d);
    //return: set of top-simplexes of dimendion d

    TopSimplex& getTopSimplex(explicitS);
    //return: top-simplex corresponding to the input simplex encoded in the explicit way

    vector<int> getVertices(const implicitS& simpl, TopSimplex& top);
    //return: return the indexes of the vertices on the boundary of simpl

    //------------------------------------------------------------------------------------------
    //functions for simplexes explicitly encoded in the IA* data structure (explicitS)

    implicitS toImplicit(const explicitS&);
    //return: implicit representation of a simplex computed from its explicit representation

    vector<implicitS>* link(const explicitS&);
    //return: the set of simplexes in the link of s.
    //NOTE: only the simplexes of highets dimension are returned, i.e. if edge e is in the link of s, e is returned by the function but its vertices are not.

    vector<implicitS>* boundaryk(const explicitS& s,uint k);
    //return: simplexes of dimension k on the boundary of s

    vector<implicitS>* coboundaryk(const explicitS& s,uint k);
    //return: simplexes of dimension k on the coboundary of s

    vector<explicitS>* topAdjacent(const explicitS& s, uint i);
    //return: simplexes adjacent to s and sharing the i-th face of s

    vector<explicitS>* topStar(const explicitS& s);
    //return: top-simplexes on the coboundary of s (or, in other words, incident to s)
    //NOTE: valid only for vertices (a top-simplex cannot have simplexes on their coboundary)


    vector<explicitS>* topStar(const explicitS& s, int d);
    //return: top-simplexes of dimension d on the coboundary of s (or, in other words, incident to s)
    //NOTE: valid only for vertices (top-simplex cannot have simplexes on their coboundary)


    //------------------------------------------------------------------------------------------
    //functions for simplexes not encoded in the IA* and implicitly represented (implicitS) [see implicitS.h for details]

    explicitS toExplicit(const implicitS&);
    //return: explicit representation of a simplex computed from its implicit representation
    //NOTE: only a vertex or a top simplex have a valid explicit representation;

    implicitS toImplicitInsideTop(const implicitS& s, const explicitS& tops);
    //return: implicit representation of s inside the top-simplex tops
    //NOTE: returned implicit representation is valid if and only if s is incident in tops

    bool theSame(const implicitS& s1,const implicitS& s2);
    //return: true if s1 and s2 represent the same simplex
    //NOTE: use theSame() when you do not know if s1 and s2 are represented based on the same top-simplex. Otherwise operator== is faster. [see implicitS.h for details]

    vector<implicitS>* link(const implicitS& s);
    //return: the set of simplexes in the link of s.
    //NOTE: only the simplexes of highets dimension are returned, i.e. if edge e is in the link of s, e is returned by the function but its vertices are not.

    vector<implicitS>* boundaryk(const implicitS& s,uint k);
    //return: simplexes of dimension k on the boundary of s

    vector<implicitS>* coboundaryk(const implicitS&,uint);
    //return: simplexes of dimension k on the coboundary of s

    vector<implicitS>* adjacents(const implicitS& s);
    //return: simplexes adjacent to s

    vector<explicitS>* topStar(const implicitS& s);
    //return: top-simplexes on the coboundary of s (or, in other words, incident to s)

    //------------------------------------------------------------------------------------------



protected:
    void buildDataStructure();
    //function for initializing the IA*

    implicitS toImplicitInsideTop(explicitS vertex, explicitS simpl);
    //return: implicit representation of a vertex (inside top simplex simpl) computed from its explicit representation

    forward_list<explicitS>* incidentCluster(explicitS v,explicitS s);
    //(given a vertex v, and a top simplex of dimension k incident in v, retrieve all the k-1 connected simplexes still incident in v)

    void recursive_insert(implicitS simplex, const implicitS& original,uint pos, uint dim, forward_list<implicitS>* ret);
    //recursive function used for computing boundary and coboundary simplexes.

    void build_top_simplex(set<uint> setR, set<uint> setP, set<uint> setX, const vector<set<uint>* >& arcs, map<int, list<TopSimplex>* >* top_simplexes_local);
    //Function for computing the top simplexes given a graph

    //Function used to compare simplexes in their implicit representation inside a set
    bool cmp(const implicitS& lhs, const implicitS& rhs)
    {
        uint lhsv = lhs.simpl.find_first();
        uint rhsv = rhs.simpl.find_first();

        TopSimplex lhsTop;
        TopSimplex rhsTop;

        if(lhs.getDimension() != rhs.getDimension())
            return lhs.getDimension() < rhs.getDimension();

        if(lhs.topDim==0)
            lhsTop = TopSimplex(lhs.topIndex);
        else
            lhsTop = getTopSimplex(explicitS(lhs.topDim,lhs.topIndex));

        if(rhs.topDim==0)
            rhsTop = TopSimplex(rhs.topIndex);
        else
            rhsTop = getTopSimplex(explicitS(rhs.topDim,rhs.topIndex));

        while(lhsv < lhs.simpl.size() && rhsv < rhs.simpl.size()){
            if(lhsTop.getVertexIndex(lhsv) != rhsTop.getVertexIndex(rhsv)){
                return lhsTop.getVertexIndex(lhsv) < rhsTop.getVertexIndex(rhsv);
            }

            lhsv = lhs.simpl.find_next(lhsv);
            rhsv = rhs.simpl.find_next(rhsv);
        }
        return false;
    }


};


#endif // SIMPLICIALCOMPLEX_H
