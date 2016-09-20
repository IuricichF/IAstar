#ifndef TOPSIMPLEX_H
#define TOPSIMPLEX_H

#include <vector>
#include <iostream>
#include <map>
#include <climits>
#include <boost/dynamic_bitset.hpp>

using namespace std;

typedef unsigned int uint;

//representation for explicit simplexes, i.e. those simplexes that are explicitly encoded (vertices and top-simplexes)
struct explicitS{
    int first;
    int second;

    explicitS() : first(-1), second(-1){}
    explicitS(int f, int s) : first(f), second(s){}

    friend ostream& operator<<(ostream& out, const explicitS& simplex){
        out << "(" << simplex.first << ", " << simplex.second << ") ";
        out.flush();
        return out;
    }

    bool operator==(explicitS simpl){return this->first == simpl.first && this->second == simpl.second;}
    bool operator<(const explicitS simpl) const{
        if(this->first == simpl.first)
            return this->second < simpl.second;
        else{
            return this->first < simpl.first;
        }
    }
    bool isValid() const{return first!=-1;}
};

struct implicitS{
    int topDim;
    int topIndex;
    boost::dynamic_bitset<> simpl;

    implicitS() : topDim(-1), topIndex(-1), simpl(boost::dynamic_bitset<>(0)){}
    implicitS(int d, int i, int size) : topDim(d), topIndex(i) { simpl = boost::dynamic_bitset<>(size); simpl.flip();} // in initialization an implicitS represents the topSimplex from which it has been created
    implicitS(int d, int i, boost::dynamic_bitset<> s) : topDim(d), topIndex(i), simpl(s) { }
    implicitS(const implicitS& simplex) : topDim(simplex.topDim), topIndex(simplex.topIndex), simpl(simplex.simpl) { }
    implicitS(const explicitS& top, const vector<int>& simplexV) : topDim(top.first), topIndex(top.second), simpl(boost::dynamic_bitset<>(top.first+1)) {
        for(uint i=0; i<simplexV.size();i++)
            simpl[simplexV[i]]=1;
    }

    uint getDimension() const {return simpl.count()-1;}

    bool isValid(){return topDim!=-1;}

    bool operator==(implicitS simpl){
        return this->topDim == simpl.topDim &&
               this->topIndex == simpl.topIndex &&
               this->simpl == simpl.simpl;}

    bool operator<(const implicitS simpl) const{
        if(this->topDim == simpl.topDim){
            if(this->topIndex == simpl.topIndex)
                return this->simpl < simpl.simpl;
            else
                return this->topIndex < simpl.topIndex;
        }
        else
            return this->topDim < simpl.topDim;
    }

    friend ostream& operator<<(ostream& os, const implicitS& simplex)
    {
        os << "(" << simplex.topDim << " " << simplex.topIndex << " " << simplex.simpl << ") ";
        return os;
    }
};



class TopSimplex
{

private:
    vector<int> vertices;
    vector<int> adjacents;

    //dimension of the top simplex vertices.size() OR adjacents.size();

public:
    TopSimplex();
    TopSimplex(vector<int> vertices);
    TopSimplex(int vertex);

    void setAdjacent(int fIndex, int top);

    int getDimension();
    int get_nVertices();
    vector<int>& getVertices();

    int getAdjacent(int fIndex);
    int getVertexIndex(int vertex);

    void print_debug();

};

#endif // TOPSIMPLEX_H
