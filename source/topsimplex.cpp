#include "topsimplex.h"

TopSimplex::TopSimplex()
{

}

TopSimplex::TopSimplex(vector<int> vertices)
{
    sort(vertices.begin(), vertices.end());
    this->vertices = vertices;
    adjacents = vector<int>(vertices.size(),INT_MAX-1);
}

TopSimplex::TopSimplex(int vertex)
{
    vertices.push_back(vertex);
}

int TopSimplex::getDimension()
{
    return vertices.size()-1;
}

int TopSimplex::get_nVertices()
{
    return vertices.size();
}

vector<int>& TopSimplex::getVertices()
{
    return vertices;
}


int TopSimplex::getAdjacent(int faceIndex)
{
    return adjacents[faceIndex];
}

int TopSimplex::getVertexIndex(int vertex)
{
    return vertices[vertex];
}

void TopSimplex::setAdjacent(int faceIndex, int simpl)
{
    adjacents[faceIndex]=simpl;
}




void TopSimplex::print_debug(){

    for(uint i=0; i<vertices.size(); i++)
        cout << vertices[i] << " ";
    cout << endl;

//    for(uint i=0; i<adjacents.size(); i++){
//        if(adjacents[i] < 0){
//            cout << "One adjacent: " << -(adjacents[i]+1) << endl;
//        }
//        else if(adjacents[i] == INT_MAX-1){
//            cout << "No adjacent" << endl;
//        }
//        else{
//            cout << "Multiple adjacents" << endl;
//        }
//    }

    cout << endl;
}
