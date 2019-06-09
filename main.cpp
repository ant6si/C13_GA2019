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
    TIME_LIMIT = int( 500 * (MAX_NUM/3000.0)) -3;
//    cout << "TIME_LIMIT: " <<TIME_LIMIT <<endl;
    gen_population_uniform(population, &gh);
//    gen_population_various(population, &gh);
//    sort(population->begin(), population->end(), compare);
    time_t remain = TIME_LIMIT - (time(NULL) - st);
    int total_best = -999999999;
    int best_score = -999999999;
    int epoch = 0;
    int xover_per_generation = int(POPULATION_SIZE * XOVER_RATIO);
    ///GA Start
    while (remain > 0) {
        vector<Chromosome *> *offsprings = new vector<Chromosome *>();
        /// New idea

        if (epoch > 400){
            xover_per_generation = 6;
        }
        else {
            if (epoch % 8 ==0){
                xover_per_generation = 4;
            }
            xover_per_generation = 3;
        }

        for (int count = 0; count < xover_per_generation; ++count) {
            sort(population->begin(), population->end(), compare);
            Chromosome *offspring = new Chromosome();
            // Selection
            int p1 = select();//_random();
            int p2 = p1;
            while (p1 == p2) {
                p2 = select();//_random();
            }
//            cout<<p1<<", "<<p2<<endl;

            // Xover
//            one_point_xover(offspring, population->at(p1), population->at(p2), &gh);

            n_point_xover(int(MAX_NUM/100), offspring, population->at(p1), population->at(p2), &gh);
            // Mutation

//            MUTATION_RATE = (MAX_MUTATION_RATE - MIN_MUTATION_RATE) / (TIME_LIMIT) * (remain) + 0.001; // annealing
            mutation(offspring);
//            local_optimize_one_chrom(offspring, &gh);
            max_locked_gain(offspring, &gh);
//            get score and regularize the new offspring
              // alreay computed by max-lg
//            regularize(offspring, &gh);
//            get_score(offspring, &gh);
            replace_hybrid(offspring, population->at(p1), population->at(p2), population->front());
//            offsprings->push_back(offspring);
        }

        // Sort
//        sort(population->begin(), population->end(), compare);
        // Replacement
//        replace_elitism(offsprings, population);
//        replace_worse(offsprings, population);
//        sort(population->begin(), population->end(), compare);
        //local optimization
//        do_local_optimize(population, &gh);
//        do_local_optimize_random(population, &gh);
//        do_local_optimize_lg(population, &gh);
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
//            cout << "time:" << (time(NULL) - st) << "/ epoch: " << epoch << "/ best_score: " << best_score
//                 << "/ median_score: " << ms << "/ worst_score: " << ws << "/ converge: " << converge << endl;
//            file_out<< "time:"<<(time(NULL)-st)<<"/ epoch: " << epoch <<"/ best_score: "<< best_score<<"/ worst_score: "<< ws<<"/ median_score: "<< ms<<"/ converge: "<< converge <<endl;
            //cout<<best_score<<endl;
        }
        remain = TIME_LIMIT - (time(NULL) - st);

        // Added to free memory
        while (!offsprings->empty()) {
            delete offsprings->back();
            offsprings->pop_back();
        }

        ///CAUTION
        offsprings->clear();
        delete (offsprings);
    }
    int ws = get_worst_score(population);
    int ms = get_median_score(population);
    float converge = how_converge(population);
    cout << epoch << "\t\t" << best_score
         << "\t\t" << ms << "\t\t" << ws << "\t\t" << converge << endl;

    sort(population->begin(), population->end(), compare);
    Chromosome *champ = population->back();
//    file_out << (champ->_score) << endl;
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


void do_islancd_GA_1(string input_file, ofstream &file_out) {
//    srand(time(NULL));
    time_t st = time(NULL);
    GraphHandler gh = GraphHandler(input_file);

//    gh.print();
    vector<Chromosome *> *population1 = new vector<Chromosome *>();
    vector<Chromosome *> *population2 = new vector<Chromosome *>();
//    vector<Chromosome *> *population3 = new vector<Chromosome *>();
//    vector<Chromosome *> *population4 = new vector<Chromosome *>();
//    vector<Chromosome *> *population5 = new vector<Chromosome *>();

MAX_NUM = gh.get_V();
    TIME_LIMIT = int( 500 * (MAX_NUM/3000.0)) -3;
//    cout << "TIME_LIMIT: " <<TIME_LIMIT <<endl;
    gen_population_uniform(population1, &gh);
    gen_population_uniform(population2, &gh);
//    gen_population_uniform(population3, &gh);
//    gen_population_uniform(population4, &gh);
//    gen_population_uniform(population5, &gh);

    time_t remain = TIME_LIMIT - (time(NULL) - st);
    int total_best = -999999999;
    int best_score = -999999999;
    int epoch = 1;
    int xover_per_generation = int(POPULATION_SIZE * XOVER_RATIO);
    ///GA Start
    while (remain > 0) {
        do_one_generation(population1, &gh);
        do_one_generation(population2, &gh);
//        do_one_generation(population3, &gh);
//        do_one_generation(population4, &gh);
//        do_one_generation(population5, &gh);

        // print best
        Chromosome* best_chrom;
        best_chrom = get_best_in_all_island(population1, population2,&gh);
//        best_chrom = get_best_in_all_island(population1, population2,population3,population4,population5, &gh);
        best_score = best_chrom->_score;
        delete(best_chrom);
/*
        if (epoch %10 == 0){
            cout << (time(NULL)-st)<< "\t\t" <<epoch << "\t\t" <<best_score<<endl;
            cout <<"POP1: ";
            print_population_status(population1);
            cout <<"POP2: ";
            print_population_status(population2);
//            cout <<"POP3: ";
//            print_population_status(population3);
//            cout <<"POP4: ";
//            print_population_status(population4);
//            cout <<"POP5: ";
//            print_population_status(population5);
  //          cout<<endl;
        }
*/
        if (epoch % 30 == 0){
            /// Apply max-lg to best champ for each island
//            max_locked_gain(population1->back(), &gh);
//            max_locked_gain(population2->back(), &gh);
//            max_locked_gain(population3->back(), &gh);
//            max_locked_gain(population4->back(), &gh);
//            max_locked_gain(population5->back(), &gh);

            /// Interchange chromosome between island
 //           cout<<"move!!!"<<endl;
            move_best_to_neighbor(population1, population2);
            move_best_to_neighbor(population2, population1);
//            move_best_to_neighbor(population2, population3);

//            move_best_to_neighbor(population3, population4);
//            move_best_to_neighbor(population4, population5);
//            move_best_to_neighbor(population5, population1);
        }

        if (total_best < best_score) {
            total_best = best_score;
        }
        epoch++;
        /*
        float converge = how_converge(population);
        //if (true){
        if (epoch % 10 == 0) {
            int ws = get_worst_score(population);
            int ms = get_median_score(population);
//            cout << "time:" << (time(NULL) - st) << "/ epoch: " << epoch << "/ best_score: " << best_score
//                 << "/ median_score: " << ms << "/ worst_score: " << ws << "/ converge: " << converge << endl;
//            file_out<< "time:"<<(time(NULL)-st)<<"/ epoch: " << epoch <<"/ best_score: "<< best_score<<"/ worst_score: "<< ws<<"/ median_score: "<< ms<<"/ converge: "<< converge <<endl;
            //cout<<best_score<<endl;
        }
         */
        remain = TIME_LIMIT - (time(NULL) - st);

    }
    float conv1 = how_converge(population1);
    float conv2 = how_converge(population2);
    cout << epoch << "\t\t" << total_best << "\t\t"<< conv1 << "\t\t"<< conv2 <<endl;

    delete(population1);
    delete(population2);
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

      input_file = "maxcut.in";
 //       input_file = "../data/HW3/treecone_overlapped_3000.txt";
        output_file = "hello.txt";
    }
    ofstream file_out;
    file_out.open(output_file.c_str());

    /// BELOW CODES ARE FOR TESTING
    GraphHandler gh = GraphHandler(input_file);
    MAX_NUM = gh.get_V();



    for(int cycle=0; cycle < MAX_CYCLE; cycle++ ){
          do_GA_1(input_file, file_out);
    //    do_islancd_GA_1(input_file, file_out);

    }

    file_out.close();





    return 0;
}
