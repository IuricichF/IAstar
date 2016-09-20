

#include "simplicialcomplex.h"
#include "Timer.h"
#include "Usage.h"

using namespace std;

int main(int argc, char* argv[])
{

    SimplicialComplex* iastar = new SimplicialComplex();

    //iastar->readGMV(argv[1]);
    Timer time;
    time.start();
    iastar->readOFF(argv[1]);
    time.stop();
    printf("Structure built\n");
    time.start();

    cout << "Save file" << endl;
    iastar->saveIA("file_dopo.ia");
    delete iastar;
    iastar = new SimplicialComplex();
    cout << "Read file" << endl;
    iastar->readIA("file_dopo.ia");
    cout << "Save file" << endl;
    iastar->saveIA("file_dopo2.ia");

//debug time
//TEST STAR OF A VERTEX
//    for(int i=0; i<iastar->getVerticesNum(); i++){
//        cout << "VERTICE - " << i << endl;
//        Vertex v = iastar->getVertex(i);
//        v.print_debug();

//        cout << "Complete Cluster: ";
//        vector<explicitS>* cluster = iastar->topStar(explicitS(0,i));
//        for(vector<explicitS>::iterator it2 =cluster->begin(); it2 != cluster->end(); it2++){
//            cout << "(" << it2->first << "," << it2->second << ") ";
//        }
//        cout << endl << endl << endl;
//        delete cluster;
//    }

//    cout << endl << endl;

//TEST ADJACENCIES
//    set<int> topSets = iastar->getTopSimplexesSet();
//    for(set<int>::iterator it = topSets.begin(); it != topSets.end(); it++){
//        for(int i=0; i<iastar->getTopSimplexesNum(*it); i++){
//            TopSimplex simpl = iastar->getTopSimplex(explicitS(*it,i));
//            simpl.print_debug();

//            for(int f=0; f<(*it+1); f++){
//                cout << "Opposite to " << simpl.getVertexIndex(f) << endl;
//                vector<explicitS>* adjs = iastar->topAdjacent(explicitS(*it,i),f);
//                for(vector<explicitS>::iterator it2 =adjs->begin(); it2 != adjs->end(); it2++){
//                    cout << "(" << it2->first << "," << it2->second << ") ";
//                }
//                cout << endl;
//                delete adjs;
//            }


//TEST BOUNDARY
//    set<int> topSets = iastar->getTopSimplexesSet();
//    for(set<int>::iterator it = topSets.begin(); it != topSets.end(); it++){
//        for(int i=0; i<iastar->getTopSimplexesNum(*it); i++){
//            TopSimplex simpl = iastar->getTopSimplex(explicitS(*it,i));
//            simpl.print_debug();

//            for(int d=simpl.getDimension()-1; d>=0; d--){
//                cout << "Boundary-" << d << endl;
//                vector<implicitS>* boundary = iastar->boundaryk(explicitS(*it,i), d);
//                for(vector<implicitS>::iterator it2 = boundary.begin(); it2 != boundary.end(); it2++){
//                    cout << it2->simpl << endl;
//                }
//                cout << endl;
//            }
//            cout << endl;
//        }
//    }


//TEST adjacents OF A VERTEX
//    for(int i=0; i<iastar->getVerticesNum(); i++){
//        cout << "VERTICE - " << i << endl;
//        Vertex v = iastar->getVertex(i);
//        v.print_debug();

//        cout << "Compute link: " << endl;
//        vector<implicitS>* cluster = iastar->link(explicitS(0,i));
//        for(vector<implicitS>::iterator it2 =cluster->begin(); it2 != cluster->end(); it2++){
//            iastar->getTopSimplex(explicitS(it2->topDim,it2->topIndex) ).print_debug();
//            cout << it2->simpl << endl;
//        }
//        delete cluster;
//        cout << endl << endl << endl;
//    }



    time.stop();
    cout << "Everything computed in " << time.getElapsedTime() << " sec." << endl;

    return 0;
}

