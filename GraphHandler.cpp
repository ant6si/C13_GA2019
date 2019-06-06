//
// Created by park on 19. 4. 16.
//

#include <fstream>
#include <iostream>
#include <list>
#include <bitset>
#include <math.h>

#include "GraphHandler.h"

using namespace std;

GraphHandler::GraphHandler(string filename){
    ifstream fs(filename.c_str());
    string line;
    getline(fs, line);


    sscanf(line.c_str(), "%d %d", &_V, &_E);
    _edges = new list<Edge*>();

    for (string line; getline(fs,line); ){
        int from, to, weight;

        if (line.length() <= 0){
            continue;
        }

        sscanf(line.c_str(), "%d %d %d", &from, &to, &weight);
        Edge* new_edge = new Edge(from-1,to-1,weight);
        _edges->push_back(new_edge);

        // For flipped vertex computation
        adj_mat[from-1].push_back(new_edge);
        adj_mat[to-1].push_back(new_edge);

    }
}



GraphHandler::~GraphHandler(){
    if(_edges != NULL){
        // To do
    }
}

int GraphHandler::get_V(){
    return _V;
}



void GraphHandler::print(){
    cout<<"V: "<<_V<<", E: "<<_E<<endl;
    list<Edge*>::iterator iter;
    int _weight;
    int _v_a = -1, _v_b=-1;
    for (iter = _edges->begin(); iter != _edges->end(); iter++){

        _v_a = (*iter)->_from;
        _v_b = (*iter)->_to;
        _weight = (*iter)->_weight;

        cout<<_v_a<<" "<<_v_b<<" "<<_weight<<endl;

    }

}

void GraphHandler::compute_score(Chromosome* chrom){
    list<Edge*>::iterator iter;
    int weight;
    int score=0;
    for(iter = _edges->begin(); iter != _edges->end(); iter++){
        weight = (*iter)->_weight;
        bitset<L> seq = chrom->_sequence;
        if( seq[(*iter)->_from] != seq[(*iter)->_to] ){
            score+=weight;
        }
    }
    chrom->_score = score;
}

int GraphHandler::compute_flipped_score(Chromosome *chrom, int score, int index) {
    // Assume chrom has un-flipped (origin) sequence
    // return  score
    if (adj_mat[index].empty()){
        // No adjacent edges.
        return score;
    }
    bitset<L> seq = chrom->_sequence;
    list<Edge *>::iterator iter;
    for (iter = adj_mat[index].begin(); iter != adj_mat[index].end(); iter++) {
        Edge *elem = (*iter);
        if( seq[ elem->_from] != seq[ elem->_to] ){
            score -= elem->_weight;
        }else{
            score += elem->_weight;
        }
    }
    return score;
}

int GraphHandler::compute_gain(Chromosome* chrom, int index){
    if (adj_mat[index].empty()){
        return 0;
    }
    bitset<L> seq = chrom->_sequence;
    list<Edge *>::iterator iter;
    int gain = 0;
    for (iter = adj_mat[index].begin(); iter != adj_mat[index].end(); iter++){
        Edge *elem = (*iter);
        if( seq[ elem->_from] != seq[ elem->_to]){
            gain -= elem->_weight;
        }else{
            gain += elem->_weight;
        }
    }
    return gain;
}



int GraphHandler::compute_locked_gain(Chromosome* chrom, int index, bool* isLocked){
    if (adj_mat[index].empty()){
        return 0;
    }
    bitset<L> seq = chrom->_sequence;
    list<Edge *>::iterator iter;
    int lg = 0;
    for (iter = adj_mat[index].begin(); iter != adj_mat[index].end(); iter++){
        Edge *elem = (*iter);
        int this_idx;
        int that_idx;
        if (elem->_from == index){
            this_idx = elem->_from;
            that_idx = elem->_to;
        }else{
            this_idx = elem->_to;
            that_idx = elem->_from;
        }
        if ( ! isLocked[that_idx]){
            continue;
        }
        if( seq[this_idx] != seq[that_idx]){
            lg -= elem->_weight;
        }else{
            lg += elem->_weight;
        }
    }
//    cout<<"Locked Gain: "<<lg<<endl;
    return lg;
}
