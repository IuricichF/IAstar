#include "simplicialcomplex.h"

void SimplicialComplex::readOFF(char* file){
    //reader for the OFF file

    ifstream fStream (file);
      if (fStream.is_open())
      {
        string line;
        int nV, nT, x;

        getline(fStream,line);
        if(line.compare("OFF") != 0){
            cout << line << endl;
            printf("Wrong file: this is not an OFF file\n");
            exit(0);
        }

        //read number of vertices and number of top simplexes
        fStream >> nV;
        fStream >> nT;
        fStream >> x;

        vertices = vector<Vertex>(nV);
        getline(fStream,line);
        //read coordinates for each vertex
        for(int i=0; i<nV; i++){

            getline(fStream,line);
            istringstream iss(line);
            string coord;
            list<float> coords;
            while( getline(iss, coord, ' ')){
                coords.push_back(atof(coord.c_str()));
            }

            vertices[i] = Vertex(vector<float>(coords.begin(),coords.end()));

        }


        map<int, list<TopSimplex>* > topSimplexeslists;

        //read top simplexes
        int nVIndexes;
        for(int i=0; i<nT; i++){
            fStream >> nVIndexes;

            vector<int> topVertices(nVIndexes);

            for(int k=0; k<nVIndexes; k++){
                fStream >> topVertices[k];
            }

            TopSimplex topS(topVertices);
            if(topSimplexeslists.find(topS.getDimension()) == topSimplexeslists.end()){
                topSimplexeslists[topS.getDimension()]=new list<TopSimplex>();
            }
            topSimplexeslists[topS.getDimension()]->push_back(topS);
        }

        int dim=0;
        topSimplexes= vector<vector<TopSimplex> >(topSimplexeslists.size(), vector<TopSimplex>());
        for(map<int, list<TopSimplex>*>::iterator it=topSimplexeslists.begin(); it!=topSimplexeslists.end(); it++){

            realIndex[it->first]=dim;
            topSimplexes[dim]=vector<TopSimplex>(it->second->begin(), it->second->end());
            delete it->second;
            dim++;
        }

        fStream.close();
      }

      buildDataStructure();
}

void SimplicialComplex::readTS(char* file){
    //reader for the OFF file


    ifstream fStream (file);
      if (fStream.is_open())
      {
        string line;
        int nV, nT, x;

        //read number of vertices and number of top simplexes
        fStream >> nV;
        fStream >> nT;

        vertices = vector<Vertex>(nV);
        getline(fStream,line);
        //read coordinates for each vertex
        for(int i=0; i<nV; i++){

            getline(fStream,line);
            istringstream iss(line);
            string coord;
            list<float> coords;
            while( getline(iss, coord, ' ')){
                coords.push_back(atof(coord.c_str()));
            }

            vertices[i] = Vertex(vector<float>(coords.begin(),coords.end()));

        }

        //read top simplexes
        realIndex[3]=0;
        topSimplexes.push_back(vector<TopSimplex>(nT));

        for(int i=0; i<nT; i++){

            vector<int> topVertices(4);

            for(int k=0; k<4; k++){
                fStream >> topVertices[k];
            }

            TopSimplex topS(topVertices);
            topSimplexes[0][i]=topS;

        }

        fStream.close();
      }

      cout << "Building structure " << endl;
      buildDataStructure();
}


void SimplicialComplex::readGMV(char* file){
    //reader for the GMV file


    ifstream fStream (file);
      if (fStream.is_open())
      {
        string line;
        getline(fStream,line);
        getline(fStream,line);
        getline(fStream,line);
        getline(fStream,line);
        getline(fStream,line);

        int nV, nT;
        fStream >> line;
        fStream >> nV;

        vertices = vector<Vertex>(nV);
        float coord;
        //read x coordinates for each vertex
        for(int i=0; i<nV; i++){
            vector<float> coords(3);
            fStream >> coords[0];
            vertices[i] = Vertex(coords);
        }
        //read y coordinates for each vertex
        for(int i=0; i<nV; i++){
            fStream >> coord;
            getVertex(i).changeCoordinate(coord,1);
        }
        //read z coordinates for each vertex
        for(int i=0; i<nV; i++){
            fStream >> coord;
            getVertex(i).changeCoordinate(coord,2);
        }

        fStream >> line;
        fStream >> nT;

        map<int, list<TopSimplex>* > topSimplexeslists;

        //read top simplexes
        int nVIndexes;
        for(int i=0; i<nT; i++){
            fStream >> line;
            fStream >> nVIndexes;

            vector<int> topVertices(nVIndexes);

            for(int k=0; k<nVIndexes; k++){
                fStream >> topVertices[k];
                topVertices[k] = topVertices[k]-1;
            }

            TopSimplex topS(topVertices);
            if(topSimplexeslists.find(topS.getDimension()) == topSimplexeslists.end()){
                topSimplexeslists[topS.getDimension()]=new list<TopSimplex>();
            }
            topSimplexeslists[topS.getDimension()]->push_back(topS);
        }

        int dim=0;
        topSimplexes= vector<vector<TopSimplex> >(topSimplexeslists.size(), vector<TopSimplex>());
        for(map<int, list<TopSimplex>*>::iterator it=topSimplexeslists.begin(); it!=topSimplexeslists.end(); it++){

            realIndex[it->first]=dim;
            topSimplexes[dim]=vector<TopSimplex>(it->second->begin(), it->second->end());
            delete it->second;
            dim++;
        }

        fStream.close();
      }

      buildDataStructure();
}


void SimplicialComplex::readPoints(char* name_in_file, float threshold, uint dimension){

    fstream fStream(name_in_file, ios::in);

    if (fStream.is_open())
    {
        list<Vertex >* points=new list<Vertex >();

        string line;
        double x;
        while( getline ( fStream , line ) )
        {
          vector<float> point;
          std::istringstream iss( line );
          while(iss >> x) { point.push_back(x); }
          points->push_back(Vertex(point));
        }
        fStream.close();

        vertices = vector<Vertex>(points->begin(), points->end());

        delete points;


        vector<set<uint>* > arcs(vertices.size());
        set<uint> ps;

        for(uint i=0; i<vertices.size(); i++){
            ps.insert(i);
            arcs[i] = new set<uint>();
        }


        int arcs_count=0;

        #pragma omp parallel for
        for(uint i=0; i<vertices.size(); i++){
            for(uint j=i+1; j<vertices.size(); j++){
                if(threshold >= vertices[i].euclideanDistance(vertices[j])){

                    #pragma omp critical
                    arcs[i]->insert(j);

                    #pragma omp critical
                    arcs[j]->insert(i);

                    #pragma omp critical
                    arcs_count++;
                }
            }
        }


        map<int, list<TopSimplex>* > topSimplexeslists;
        build_top_simplex(set<uint>(), ps, set<uint>(), arcs, &topSimplexeslists);


        int dim=0;
        int totTop=0;
        topSimplexes= vector<vector<TopSimplex> >(topSimplexeslists.size(), vector<TopSimplex>());
        for(map<int, list<TopSimplex>*>::iterator it=topSimplexeslists.begin(); it!=topSimplexeslists.end(); it++){

            realIndex[it->first]=dim;
            topSimplexes[dim]=vector<TopSimplex>(it->second->begin(), it->second->end());
            delete it->second;
            totTop+=topSimplexes[dim].size();
            dim++;
        }

        buildDataStructure();
    }
}

void SimplicialComplex::readGraph(char* name_in_file){

    fstream fStream(name_in_file, ios::in);

    if (fStream.is_open())
    {
        string line;
        int vertN;
        int edges;
        getline ( fStream , line );

        std::istringstream iss( line );
        iss >> vertN;
        iss >> edges;

        vertices = vector<Vertex>(vertN);

        vector<set<uint>* > arcs(vertices.size());

        set<uint> points;
        for(uint i=0; i<vertices.size(); i++){
            arcs[i] = new set<uint>();
            points.insert(i);
        }


        map<int, int> vertMap;
        int v1,v2;
        int vReal=0;
        for(int i=0; i<edges; i++){
            getline ( fStream , line );

            std::istringstream iss( line );
            iss >> v1;
            iss >> v2;

            if(vertMap.find(v1) == vertMap.end()){
                vertMap[v1]=vReal;
                v1=vReal;
                vReal++;
            }
            else
                v1 = vertMap[v1];

            if(vertMap.find(v2) == vertMap.end()){
                vertMap[v2]=vReal;
                v2=vReal;
                vReal++;
            }
            else
                v2 = vertMap[v2];


            arcs[v1]->insert(v2);
            arcs[v2]->insert(v1);
        }


        map<int, list<TopSimplex>* > topSimplexeslists;
        build_top_simplex(set<uint>(), points, set<uint>(), arcs, &topSimplexeslists);

        int dim=0;
        int totTop=0;
        topSimplexes= vector<vector<TopSimplex> >(topSimplexeslists.size(), vector<TopSimplex>());
        for(map<int, list<TopSimplex>*>::iterator it=topSimplexeslists.begin(); it!=topSimplexeslists.end(); it++){

            realIndex[it->first]=dim;
            topSimplexes[dim]=vector<TopSimplex>(it->second->begin(), it->second->end());
            delete it->second;
            totTop+=topSimplexes[dim].size();
            dim++;
        }

        buildDataStructure();
    }
}

void SimplicialComplex::readIA(char* file){
    //reader for the OFF file

    ifstream fStream (file);
    if (fStream.is_open())
    {
        //number of vertices and number of different types of top simplexes
        int nV, nT;
        fStream >> nV;
        fStream >> nT;

        //for each type of top simplex, its real index in the vector of top simplexes
        int index, real;
        for(int i=0; i<nT; i++){
            fStream >> index;
            fStream >> real;
            realIndex[index]=real;
        }


        //for each vertex
        vertices = vector<Vertex>(nV);
        int coordDim;
        int cobDim, cobSize, cobInt, cobIndex;
        for(int i=0; i<nV; i++){

            //number of coordinates and values
            fStream >> coordDim;
            vector<float> coords(coordDim);
            for(int j=0; j<coordDim; j++)
                fStream >> coords[j];

            vertices[i] = Vertex(coords);
            fStream >> cobDim;
            for(int j=0; j<cobDim; j++){
                fStream >> cobInt;
                fStream >> cobSize;
                for(int k=0; k<cobSize; k++){
                    fStream >> cobIndex;
                    vertices[i].addPartialCoboundaryTop(cobInt,cobIndex);
                }
            }
        }

        //for each topsimplex dimension
        int topSize;
        int topDim;
        for(int i=0; i<nT; i++){
            //topSimplexes = vector<vector<TopSimplex> >(nT,vector<TopSimplex>());
            fStream >> topSize;
            vector<TopSimplex> tops(topSize);
            for(int j=0; j<topSize; j++){
                fStream >> topDim;
                vector<int> verts(topDim);
                for(int k=0; k<topDim; k++){
                    fStream >> verts[k];
                }
                tops[j]=TopSimplex(verts);
                int adjs;
                for(int k=0; k<topDim; k++){
                    fStream >> adjs;
                    tops[j].setAdjacent(k,adjs);
                }
            }
            topSimplexes.push_back(tops);
        }

        //for each non-manifold adjacent relation
        int adjRels, dim, indexes;
        fStream >> adjRels;
        adjRelations = vector<forward_list<int> >(adjRels);
        for(int i=0; i<adjRels; i++){
            forward_list<int> adj;
            fStream >> dim;
            for(int j=0; j<dim; j++){
                fStream >> indexes;
                adj.push_front(indexes);
            }
            adjRelations[i]=adj;
        }

        fStream.close();
    }
}


void SimplicialComplex::saveIA(char* file){

    ofstream fStream (file);
    if (fStream.is_open())
    {
        //number of vertices and number of different types of top simplexes
        fStream << vertices.size() << " " << realIndex.size() << endl;

        //for each type of top simplex, its real index in the vector of top simplexes
        for(map<int,int>::iterator it=realIndex.begin(); it!=realIndex.end(); it++){
            fStream << it->first << " " << it->second << endl;
        }
        fStream << endl;

        //for each vertex
        for(int i=0; i<vertices.size(); i++){
            vector<float> coords = vertices[i].getCoordinates();
            //number of coordinates and values
            fStream << coords.size() << " ";
            for(int j=0; j<coords.size(); j++){
                fStream << coords[j] << " ";
            }
            fStream << endl;

            map<int, vector<int> > cob = vertices[i].getPartialCoboundaryTopRelations();
            fStream << cob.size() << endl;
            for(map<int, vector<int> >::iterator it=cob.begin(); it!= cob.end(); it++){
                fStream << it->first << " " << it->second.size() << " ";
                vector<int> starC = it->second;
                for(int k=0; k<starC.size(); k++){
                    fStream << starC[k] << " ";
                }
                fStream << endl;
            }
        }

        //for each topsimplex dimension
        for(int i=0; i<topSimplexes.size(); i++){
            fStream << topSimplexes[i].size() << endl;

            //for each topsimplex
            for(int j=0; j<topSimplexes[i].size(); j++){

                //save vertex and adjacent simplex indices
                fStream << topSimplexes[i][j].getDimension()+1 << endl;
                for(int k=0; k<=topSimplexes[i][j].getDimension(); k++){
                    fStream << topSimplexes[i][j].getVertices()[k] << " ";
                }
                fStream << endl;

                for(int k=0; k<=topSimplexes[i][j].getDimension(); k++){
                    fStream << topSimplexes[i][j].getAdjacent(k) << " ";
                }
                fStream << endl;
            }
        }

        //for each non-manifold adjacent relation
        fStream << adjRelations.size() << endl;
        for(int i=0; i<adjRelations.size(); i++){
            vector<int> simplexes(adjRelations[i].begin(), adjRelations[i].end());
            //number of simplexes
            fStream << simplexes.size() << " ";

            //adjacent simplexes
            for(int j=0; j<simplexes.size(); j++){
                fStream << simplexes[j] << " ";
            }
            fStream << endl;
        }

        fStream.close();
    }
}



















