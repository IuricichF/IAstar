#include "simplicialcomplex.h"



implicitS SimplicialComplex::toImplicit(const explicitS& simpl){

    if(simpl.first == 0){
        map<int, vector<int> > partialCob = getVertex(simpl.second).getPartialCoboundaryTopRelations();
        for(map<int, vector<int> >::iterator it = partialCob.begin(); it != partialCob.end(); it++){
            for(uint k=0; k<it->second.size(); k++){
                return toImplicitInsideTop(simpl, explicitS(it->first, it->second[k]));
            }
        }
        return implicitS(0,simpl.second,1);
    }

    return implicitS(simpl.first,simpl.second,simpl.first+1);
}

implicitS SimplicialComplex::toImplicitInsideTop(const implicitS& simplex, const explicitS& top){

    if(simplex.topDim == top.first && simplex.topIndex == top.second)
        return simplex;

    vector<int> vertices_small = getVertices(simplex, getTopSimplex(explicitS(simplex.topDim, simplex.topIndex)));
    vector<int> vertices_big = getTopSimplex(top).getVertices();

    vector<int> internal_indices;
    for(int i=0; i<vertices_small.size(); i++){
        vector<int>::iterator it = lower_bound(vertices_big.begin(),vertices_big.end(),vertices_small[i]);
        internal_indices.push_back(it - vertices_big.begin());
    }

//    vector<implicitS>* bound = boundaryk(top,simplex.getDimension());
//    for(vector<implicitS>::iterator it=bound->begin(); it!=bound->end(); it++){

//        if(theSame(*it,simplex)){

//            implicitS impl = *it;
//            delete bound;

//            return impl;
//        }
//    }

    return implicitS(top,internal_indices);
}

implicitS SimplicialComplex::toImplicitInsideTop(explicitS vertex, explicitS top){

    implicitS implSimpl(top.first,top.second,top.first+1);
    implSimpl.simpl.flip();
    assert(implSimpl.simpl.none());
    for(int i=0; i<top.first+1; i++){
        if(getTopSimplex(top).getVertexIndex(i) == vertex.second){
            implSimpl.simpl[i]=1;
            assert(implSimpl.simpl.count() == 1);
            return implSimpl;
        }
    }

    return implicitS();
}

explicitS SimplicialComplex::toExplicit(const implicitS& simpl){

    int non_zero = simpl.simpl.count();
    if(non_zero == 1){
        if(simpl.topDim==0)
            return explicitS(0,simpl.topIndex);
        return explicitS(0, (topSimplexes[realIndex[simpl.topDim]][simpl.topIndex]).getVertices()[simpl.simpl.find_first()]);
    }
    else if(non_zero == simpl.topDim+1){
        return explicitS(simpl.topDim, simpl.topIndex);
    }
    else{
        printf("The encoded simplex is neither a vertex nor a top simplex\n");
    }

    return explicitS();
}

bool SimplicialComplex::theSame(const implicitS& s1,const implicitS& s2){

    uint lhsv = s1.simpl.find_first();
    uint rhsv = s2.simpl.find_first();
    TopSimplex& lhsTop = getTopSimplex(explicitS(s1.topDim,s1.topIndex));
    TopSimplex& rhsTop = getTopSimplex(explicitS(s2.topDim,s2.topIndex));

    while(lhsv < s1.simpl.size() && rhsv < s2.simpl.size()){
        if(lhsTop.getVertexIndex(lhsv) != rhsTop.getVertexIndex(rhsv))
            return false;

        lhsv = s1.simpl.find_next(lhsv);
        rhsv = s2.simpl.find_next(rhsv);
    }
    return true;
}


vector<explicitS>* SimplicialComplex::topStar(const implicitS& simpl){

    if(simpl.getDimension() == 0){
        return topStar(toExplicit(simpl));
    }
    else{
        int v=0;
        int maxDim=0;
        vector<int> vertices(simpl.simpl.count(),0);
        for(boost::dynamic_bitset<>::size_type i = simpl.simpl.find_first(); i < simpl.simpl.size(); i = simpl.simpl.find_next(i)){
            int vIndex = getTopSimplex(explicitS(simpl.topDim,simpl.topIndex)).getVertexIndex(i);
            if(maxDim < getVertex(vIndex).getCoboundaryMaxDim())
                maxDim = getVertex(vIndex).getCoboundaryMaxDim();
            vertices[v]=vIndex; v++;

        }
        maxDim++;


        map<explicitS, uint> top_found;

        #pragma omp parallel for shared(top_found)
        for(uint i=0; i<vertices.size(); i++){

            set<explicitS> top_per_vertex;

            for(int d=simpl.getDimension()+1; d<=getVertex(vertices[i]).getCoboundaryMaxDim(); d++){

                vector<explicitS>* top = topStar(explicitS(0,vertices[i]),d);
                top_per_vertex.insert(top->begin(), top->end());

                delete top;
            }

            #pragma omp critical
            {

                for(set<explicitS>::iterator it=top_per_vertex.begin(); it!=top_per_vertex.end(); it++){
                    map<explicitS, uint>::iterator it2 = top_found.find(*it);
                    if(it2 == top_found.end())
                        top_found[*it]=1;
                    else
                        top_found[*it]=top_found[*it]+1;
                }
            }
        }

        forward_list<explicitS> ret = forward_list<explicitS>();
        for(map<explicitS, uint>::iterator it = top_found.begin(); it!=top_found.end(); it++){
            if(it->second == vertices.size()){
                ret.push_front(it->first);
            }
        }

        return new vector<explicitS>(ret.begin(), ret.end());
    }
}

vector<implicitS>* SimplicialComplex::link(const implicitS& simplex){


    forward_list<implicitS> ret;
    vector<explicitS>* cob = topStar(simplex);

    #pragma omp parallel for shared(ret)
    for(int i=0; i<cob->size(); i++){
        implicitS inside = toImplicitInsideTop(simplex,(*cob)[i]);
        assert(inside.topDim!=-1);
        inside.simpl.flip();

        #pragma omp critical
        ret.push_front(inside);
    }
    delete cob;

    return new vector<implicitS>(ret.begin(), ret.end());
}



vector<implicitS>* SimplicialComplex::coboundaryk(const implicitS& simplex,uint dim){


    if(simplex.getDimension() >= dim){
        printf("No simplexes of dimension %d on the coboundary of a %d-simplex", dim, simplex.topDim);
        return new vector<implicitS>();
    }


    vector<explicitS>* star = topStar(simplex);

    auto foo = bind(&SimplicialComplex::cmp, this,_1,_2);
    set<implicitS, boost::function<bool(const implicitS &, const implicitS &)>> simplexesSet(foo);

    #pragma omp parallel for shared(simplexesSet)
    for(int i=0; i<star->size(); i++){


        forward_list<implicitS> coblocal;
        implicitS inTop(toImplicitInsideTop(simplex, (*star)[i]));

        implicitS newS(inTop);
        inTop.simpl = inTop.simpl.flip();
        uint first = inTop.simpl.find_first();
        recursive_insert(newS,inTop,first, dim, &coblocal);
        #pragma omp critical
        {
            for(forward_list<implicitS>::iterator it2 = coblocal.begin(); it2!=coblocal.end(); it2++){
                simplexesSet.insert(*it2);
            }
        }
    }

    delete star;

    return new vector<implicitS>(simplexesSet.begin(), simplexesSet.end());
}



vector<implicitS>* SimplicialComplex::boundaryk(const implicitS& simplex,uint dim){

    forward_list<implicitS> ret;
    if(simplex.getDimension() <= dim){
        printf("No simplexes of dimension %d on the boundary of a %d-simplex", dim, simplex.topDim);
        return new vector<implicitS>();
    }

    if(simplex.topDim == dim){
        implicitS newSimpl(simplex.topDim,simplex.topIndex,simplex.topDim+1);
        vector<implicitS>* vec= new vector<implicitS>(1);
        (*vec)[0] = newSimpl;
        return vec;
    }

    implicitS newS(simplex);
    newS.simpl = newS.simpl.reset();
    uint first = simplex.simpl.find_first();
    recursive_insert(newS,simplex,first, dim, &ret);

    return new vector<implicitS>(ret.begin(), ret.end());
}


void SimplicialComplex::recursive_insert(implicitS simplex, const implicitS& original,uint pos, uint dim, forward_list<implicitS>* ret){


    if(simplex.simpl.count() == dim){
        while(pos < original.simpl.size()){
            simplex.simpl[pos]=1;
            ret->push_front(simplex);
            simplex.simpl[pos]=0;
            pos = original.simpl.find_next(pos);
        }

    }

    else if(simplex.simpl.count() < dim){
        while(pos < original.simpl.size()){
            uint nextPos = original.simpl.find_next(pos);
            simplex.simpl[pos]=1;
            recursive_insert(simplex, original, nextPos, dim, ret);
            simplex.simpl[pos]=0;
            pos = nextPos;
        }
    }
}

vector<implicitS>* SimplicialComplex::adjacents(const implicitS& simpl){

    auto foo = bind(&SimplicialComplex::cmp, this,_1,_2);
    set<implicitS, boost::function<bool(const implicitS &, const implicitS &)>> simplexesSet(foo);

    if(simpl.getDimension() > 0){
        vector<implicitS>* boundary = boundaryk(simpl,simpl.getDimension()-1);

        for(vector<implicitS>::iterator it = boundary->begin(); it!=boundary->end(); it++){
            vector<implicitS>* cob = coboundaryk(*it,simpl.getDimension());
            simplexesSet.insert(cob->begin(), cob->end());
            delete cob;
        }
        delete boundary;
    }
    else{
        vector<implicitS>* coboundary = coboundaryk(simpl,1);
        for(vector<implicitS>::iterator it = coboundary->begin(); it!=coboundary->end(); it++){
            vector<implicitS>* cob = boundaryk(*it,0);
            simplexesSet.insert(cob->begin(), cob->end());
            delete cob;
        }

        delete coboundary;
    }

    simplexesSet.erase(simpl);
    return new vector<implicitS>(simplexesSet.begin(), simplexesSet.end());
}
