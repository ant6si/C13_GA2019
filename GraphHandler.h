//
// Created by park on 19. 5. 17.
//

#ifndef C13_GA2019_GRAPHHANDLER_H
#define C13_GA2019_GRAPHHANDLER_H

#include <bitset>
#include <list>
#include <vector>

#define L 1000

using namespace std;

struct Edge{
public:
    int _from, _to;
    int _weight;

public:
    Edge(int f, int t, int w):
            _from(f), _to(t), _weight(w)
    {
    }
};

struct Chromosome{
public:
    bitset<L> _sequence;
    long _score;
public:
    Chromosome(bitset<L> seq, long s):
            _sequence(seq),_score(s)
    {}
    Chromosome(){
        _sequence = bitset<L>(0);
        _score = 0;
    }
    /*
public:
    bool operator< (Chromosome& rval) {
        long s1 = this->_score;
        long s2 = rval._score;

        return (s1 < s2);
    }
     */
};

class GraphHandler {
private:
    int _V, _E;
    list<struct Edge*>* _edges = NULL;
    list<struct Edge*> adj_mat[L];
public:
    GraphHandler(string filename);
    ~GraphHandler();

public:
    int get_V();
    void print();
    void compute_score(Chromosome* chrom);
    int compute_flipped_score(Chromosome* chrom, int score, int index);
};

#endif //C13_GA2019_GRAPHHANDLER_H
