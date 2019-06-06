#include <cstdio>
#include <cstring>
#include <stdlib.h>
#include <limits.h>

#include <iostream>
#include <fstream>

#include <string>
#include <list>
#include <set>
#include <algorithm>
#include <vector>

#include <time.h>
#include <assert.h>

#include <bitset>
#include "GraphHandler.h"
#include "methods.h"


using namespace std;


void do_GA_1(string input_file, ofstream &file_out) {
//    srand(time(NULL));
    time_t st = time(NULL);
    GraphHandler gh = GraphHandler(input_file);

//    gh.print();
    vector<Chromosome *> *population = new vector<Chromosome *>();
    MAX_NUM = gh.get_V();
    gen_population_uniform(population, &gh);
//    gen_population_various(population, &gh);
    sort(population->begin(), population->end(), compare);
    time_t remain = TIME_LIMIT - (time(NULL) - st);
    int total_best = -999999999;
    int best_score = -999999999;
    int epoch = 0;
    float last_converge;
    int xover_per_generation = int(POPULATION_SIZE * XOVER_RATIO);
    ///GA Start
    while (remain > 0) {
        vector<Chromosome *> *offsprings = new vector<Chromosome *>();
        for (int count = 0; count < xover_per_generation; ++count) {
            Chromosome *offspring = new Chromosome();
            // Selection
            int p1 = select_random();
            int p2 = p1;
            while (p1 == p2) {
                p2 = select_random();
            }
//            cout<<p1<<", "<<p2<<endl;

            // Xover
            one_point_xover(offspring, population->at(p1), population->at(p2), &gh);
            // Mutation
            mutation(offspring);
            local_optimize_one_chrom(offspring, &gh);
//            max_locked_gain(offspring, &gh);
            //get score and regularize the new offspring
//            regularize(offspring, &gh);
            get_score(offspring, &gh);
            offsprings->push_back(offspring);
        }

        // Sort
        sort(population->begin(), population->end(), compare);
        // Replacement
//        replace_elitism(offsprings, population);
        replace_worse(offsprings, population);
        sort(population->begin(), population->end(), compare);
        //local optimization
//        do_local_optimize(population, &gh);
//        do_local_optimize_random(population, &gh);
        sort(population->begin(), population->end(), compare);

        // print best
        best_score = get_best_score(population);

        if (total_best < best_score) {
            total_best = best_score;
//            cout<<total_best<<endl;
        }

        epoch++;
        float converge = how_converge(population);
        //if (true){
        if (epoch % 10 == 0) {
            int ws = get_worst_score(population);
            int ms = get_median_score(population);
            cout << "time:" << (time(NULL) - st) << "/ epoch: " << epoch << "/ best_score: " << best_score
                 << "/ median_score: " << ms << "/ worst_score: " << ws << "/ converge: " << converge << endl;
//            file_out<< "time:"<<(time(NULL)-st)<<"/ epoch: " << epoch <<"/ best_score: "<< best_score<<"/ worst_score: "<< ws<<"/ median_score: "<< ms<<"/ converge: "<< converge <<endl;
            //cout<<best_score<<endl;
        }
        remain = TIME_LIMIT - (time(NULL) - st);
        last_converge = converge;

        // Added to free memory
        while (!offsprings->empty()) {
            delete offsprings->back();
            offsprings->pop_back();
        }

        ///CAUTION
        offsprings->clear();
        delete (offsprings);
    }
    sort(population->begin(), population->end(), compare);
    Chromosome *champ = population->back();
    file_out << (champ->_score) << endl;
//    file_out<< "time:"<<(time(NULL)-st)<<"/ epoch: " << epoch << "/ best_score: "<< best_score<<"/ converge: "<< last_converge <<endl;
/* For submission
    sort( population->begin(), population->end(), compare);
    Chromosome* champ = population->back();
    bitset<L> champ_seq = champ->_sequence;

    vector<int>* champ_vec = new vector<int>();

    for(int x=0; x<MAX_NUM; x++){
        if(champ_seq[x]==0){
            champ_vec->push_back((x+1));
        }
    }
    vector<int>::iterator iter;
    iter = champ_vec->begin();
    file_out<<(*iter);
    iter++;
    for(; iter != champ_vec->end(); iter++){
        file_out<<" "<< (*iter);
    }
    delete champ_vec;
*/
//printVec(population);
}

///Local_optimization
void do_MS_local_opt(string input_file, ofstream &file_out) {
//    srand((unsigned int)time(0));
    time_t st = time(NULL);
    GraphHandler gh = GraphHandler(input_file);

//    gh.print();
    vector<Chromosome *> *population = new vector<Chromosome *>();
    MAX_NUM = gh.get_V();
    gen_population_uniform(population, &gh);
//    printVec(population);

//    gen_population_various(population, &gh);
    sort(population->begin(), population->end(), compare);
    time_t remain = TIME_LIMIT - (time(NULL) - st);
    int total_best = -999999999;
    int best_score = -999999999;
    int epoch = 0;
    /// Local opt start

    while (remain > 0) {
        delete (population);
        vector<Chromosome *> *population = new vector<Chromosome *>();
        //local optimization
        gen_population_uniform(population, &gh);
//        printVec(population);
        sort(population->begin(), population->end(), compare);
        int before_best = get_best_score(population);
        do_local_optimize(population, &gh);
        sort(population->begin(), population->end(), compare);

        // print best
        best_score = get_best_score(population);

        if (total_best < best_score) {
            total_best = best_score;
//            cout<<total_best<<endl;
        }

        epoch++;
        //if (true){
        if (true) {
            cout << "time:" << (time(NULL) - st) << "/ epoch: " << epoch << "/ before best: " << before_best
                 << "/ best_score: " << best_score << endl;
            //cout<<best_score<<endl;
        }
        remain = TIME_LIMIT - (time(NULL) - st);
    }

    file_out << total_best << endl;
//    file_out<< "time:"<<(time(NULL)-st)<<"/ epoch: " << epoch << "/ best_score: "<< best_score<<"/ converge: "<< last_converge <<endl;

//printVec(population);
}

int main(int argc, char *argv[]) {
    srand((unsigned int)time(NULL));

    /// 제출시 인덱스 1 더하기!
    string input_file;
    string output_file;

    if (argc == 3) {
        input_file = string(argv[1]);
        output_file = string(argv[2]);
    } else {
        // chimera_946.txt cubic_1000.txt planar_800.txt toroidal_800.txt random_500.txt random_1000.txt
        input_file = "../input/chimera_946.txt";
        output_file = "hello.txt";
    }
    ofstream file_out;
    file_out.open(output_file.c_str());
/*
    /// BELOW CODES ARE FOR TESTING
    GraphHandler gh = GraphHandler(input_file);
    MAX_NUM = gh.get_V();
    cout<<input_file<<endl;
    int ccount = 0;
    for (int iter=0; iter<100; iter++){
        Chromosome* debug_chrom = gen_chromosome(0.5, &gh);
        Chromosome* chrom_one_bit = new Chromosome();
        Chromosome* chrom_lg = new Chromosome();
        chrom_one_bit->_sequence = debug_chrom->_sequence;
        chrom_lg ->_sequence = debug_chrom->_sequence;

//    cout<<(debug_chrom->_score)<<"/ "<<(debug_chrom->_sequence)<<endl;
        local_optimize_one_chrom(chrom_one_bit, &gh);
        gh.compute_score(chrom_one_bit);
        max_locked_gain(chrom_lg, &gh);
        gh.compute_score(chrom_lg);

        cout<< "Locked > Onebit: "<<(chrom_lg->_score > chrom_one_bit->_score)<< "/ locked: "<<chrom_lg->_score<<" / one bit: "<<chrom_one_bit->_score<<endl;
        if (chrom_lg->_score > chrom_one_bit->_score){
            ccount++;
        }
        delete debug_chrom;
        delete chrom_one_bit;
        delete chrom_lg;
    }
    cout<< "Good count: "<<ccount<<endl;
*/

    do_GA_1(input_file, file_out);
//for (int cc=0;cc<5;cc++){
//    do_MS_local_opt(input_file,file_out);
//
//}


/*
 * /// For parameter seach
    string input_files[] = {"../input/chimera_946.txt", "../input/cubic_1000.txt", "../input/planar_800.txt", "../input/toroidal_800.txt", "../input/random_500.txt", "../input/random_1000.txt"};
    input_file = input_files[4];
    output_file = "../output/2_selection_Result.txt";
    ofstream file_out;
    file_out.open(output_file.c_str());

    for (int GA_NUM = 2; GA_NUM<3; GA_NUM++){
        file_out<< "New method start: "<<GA_NUM<<endl;
        for(int file_num=0; file_num<5;file_num++){
            input_file = input_files[file_num];
            file_out<< "Input file: "<<input_file<<endl;
            for(int pi = 1; pi<4; pi++){
               for (float ci = 1; ci<4; ci++){
                   for(int oi = 2; oi<4; oi++){
                       POPULATION_SIZE = pi*80;
                       XOVER_RATIO = ci*0.05;
                       OPTIMIZE_RATIO = 0.5+oi*0.1;
                       file_out<< "POPULATION SIZE: "<<POPULATION_SIZE << "/ XOVER_RATIO: "<<XOVER_RATIO<<"/ OPT_RATIO: "<<OPTIMIZE_RATIO<<endl;
                       for(int it=0; it<5; it++){

                           if(GA_NUM==1){
                               do_GA_1(input_file, file_out);
                           }
                           if(GA_NUM==2){
                               do_GA_2(input_file, file_out);
                           }


                       }
                   }

               }
            }


        }
    }
    file_out.close();

*/

/*
    ///FOR SUBMISSION
    string input_file = "maxcut.in";
    string output_file = "maxcut.out";
    ofstream file_out;
    file_out.open(output_file.c_str());
    do_GA_1(input_file, file_out);

    //printVec(population);
    //maxout 에 인덱스 1더하기!!!
*/
    file_out.close();

    return 0;
}
