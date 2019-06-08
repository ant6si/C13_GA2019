//
// Created by park on 19. 5. 17.
//

#ifndef C13_GA2019_METHODS_H
#define C13_GA2019_METHODS_H
using namespace std;
/// Parameters
int TIME_LIMIT = 500; //500
int MAX_NUM; // valid gene length
int POPULATION_SIZE = 30; //100     //300;170
//for crossover
float XOVER_RATIO = 0.2; //0.02
// for selection
float MAX_FITNESS = 1.1;
float MIN_FITNESS = 1;
// for mutation
float MUTATION_RATE = 0.07 ;//0.01
// for annealing
float MAX_MUTATION_RATE = 0.07;
float MIN_MUTATION_RATE = 0.001;
// for replacement
float ELITISM_RATE = 0.1;
// for local optimization
float OPTIMIZE_RATIO = 0.7;
//Cycle_count
int MAX_CYCLE = 20;

bool compare(Chromosome* c1, Chromosome* c2){
    long sc1 = c1->_score;
    long sc2 = c2->_score;
    return sc1 < sc2;
}

void printVec(vector<Chromosome*>* vv){
    vector<Chromosome*>::iterator iter;
    for(iter = vv->begin(); iter!= vv->end() ; iter++ ){
//        cout<< (*iter)->_score<<endl;
        cout<< (*iter)->_score << "\t\t/"   << ((*iter)->_sequence)<<endl;
        //cout<<(*iter)<<endl;
    }
    cout<<endl;
}

bool is_equal(bitset<L> a, bitset<L> b, int max){
    bool flag = true;
    for(int i=0; i<max;i++){
        if(a[i] != b[i]) flag = false;
    }
    return flag;
}

float how_converge(vector<Chromosome*>* vv){
    /// 최적해와 같은 해의 비율
    Chromosome* best = vv->back();
    int c = 0; //count
    int s = vv->size(); // size
    vector<Chromosome*>::reverse_iterator riter(vv->rbegin());
    for(; riter != vv->rend(); ++riter){
        if ((*riter)->_score == (best->_score)) {
//        if ((*riter)->_score > (best->_score) * 0.9) {
            //cout<<"same"<<endl;
            c++;
        }
        else{
            break;
        }
    }
    float converge = float(c) / float(s);
    return converge;
}

void regularize(Chromosome* chrom, GraphHandler* gh){
    // Regularize chromosome and compute score
    // rule: last bit must be 1

    if (chrom->_sequence[0] == 0){
        //cout<<"before reg: "<<(chrom->_sequence)<<endl;
        for(int i=0;i<MAX_NUM;i++) {
            chrom->_sequence.flip(i);
        }

        //cout<<"after reg: "<< (chrom->_sequence) <<endl;
    }

    gh->compute_score(chrom);
}

void get_score(Chromosome* chrom, GraphHandler* gh){
    gh->compute_score(chrom);
}

int get_flipped_score (Chromosome* chrom, int index, GraphHandler* gh){
    int flipped_score = gh->compute_flipped_score(chrom, (chrom->_score), index);
    return flipped_score;
}

Chromosome* gen_chromosome(float threshold, GraphHandler* gh){
    bitset<L> new_seq(0);
//    double p = ((double) rand() / double(RAND_MAX));
    for(int gene_num=0; gene_num<MAX_NUM; gene_num++){
        double p = ( (double) rand() / double(RAND_MAX));
        if(p > threshold){
            new_seq.flip(gene_num);
        }
    }
    Chromosome* new_chrom = new Chromosome(new_seq, 0);
    regularize(new_chrom,gh);
    get_score(new_chrom, gh);
    return new_chrom;
}

void gen_population_uniform (vector<Chromosome*>* pop, GraphHandler* gh){
    // Generate new random population
    for(int chrom_num=0; chrom_num<POPULATION_SIZE; chrom_num++){
        Chromosome* new_chrom = gen_chromosome(0.5, gh);
        pop->push_back(new_chrom);
    }
}

void gen_population_various(vector<Chromosome*>* pop, GraphHandler* gh){
    for(int chrom_num=0; chrom_num<POPULATION_SIZE; chrom_num++){
        Chromosome* new_chrom;
        if(chrom_num < POPULATION_SIZE * 0.2){
            new_chrom = gen_chromosome(0.3, gh);
        }
        else if (chrom_num < POPULATION_SIZE * 0.5){
            new_chrom = gen_chromosome(0.4, gh);
        }else{
            new_chrom = gen_chromosome(0.5, gh);
        }
        pop->push_back(new_chrom);
    }
}

//void gen_population_uniform(vector<Chromosome*>* pop, GraphHandler* gh){
//    for(int chrom_num=0; chrom_num<POPULATION_SIZE; chrom_num++){
//        bitset<L> new_seq(0);
//        for(int gene_num=0; gene_num<MAX_NUM; gene_num++){
//            if(rand()%2==1){
//                new_seq[gene_num]=1;
//            }
//        }
//        Chromosome* new_chrom = new Chromosome(new_seq, 0);
//        get_score(new_chrom, gh);
//        pop->push_back(new_chrom);
//    }
//}


/// Selection
float compute_fitness(int rank){
    float fitness = MAX_FITNESS + (float(rank-1)) * (MIN_FITNESS-MAX_FITNESS)/(float(POPULATION_SIZE)-1.0);
    return fitness;
}
int select(){
    int num;
    float sum_of_fitness = float(POPULATION_SIZE) * (MAX_FITNESS+MIN_FITNESS)/2.0;
    double point = ((float) rand() / float(RAND_MAX)) * sum_of_fitness;
    float sum = 0;
    for (int i=POPULATION_SIZE-1; i>=0; i--){
        sum += compute_fitness(POPULATION_SIZE-i);
        if(point<sum){
            num = i;
            break;
        }
    }
    return num;
}

int select_random(){
    int num = rand()%POPULATION_SIZE;
    return num;
}

int select_random_up() {
    int num = rand() % (POPULATION_SIZE / 2);
    return num;
}
/// Crossover
void xover(Chromosome* offspring, Chromosome* p1, Chromosome* p2, GraphHandler* gh){
// p1, p2 : parents
// uniform crossover
    bitset<L> seq1 = p1->_sequence;
    bitset<L> seq2 = p2->_sequence;
    bitset<L> new_seq(0);
    for(int i=0; i<MAX_NUM; i++){
        int r = rand()%2;
        if(r == 0 ){
            new_seq.set( i, seq1[i] );
        }else{
            new_seq.set( i, seq2[i] );
        }
    }
    offspring->_sequence = new_seq;
}

void one_point_xover(Chromosome* offspring, Chromosome* p1, Chromosome* p2, GraphHandler* gh){
    /// Need to Check
    int point = rand()%MAX_NUM;
    for (int i=0; i<MAX_NUM; i++){
        if(i<point){
            offspring->_sequence.set(i,p1->_sequence[i]);
        }else{
            offspring->_sequence.set(i,p2->_sequence[i]);
        }
    }
}

void n_point_xover(int n, Chromosome* offspring, Chromosome* p1, Chromosome* p2, GraphHandler* gh){

    /// Generate cut points
    vector<int> cut_points;
    for(int idx = 0; idx < n; idx++){
        cut_points.push_back(rand()%MAX_NUM);
    }
    sort(cut_points.begin(), cut_points.end());
    for(int i=0; i<cut_points.size(); i++){
//        if (i>0 && cut_points[i] == cut_points[i-1]){
//            cut_points[i] +=1;
//        }
//        cout<< cut_points[i]<< " ";
    }
//    cout<<endl;

    bool from_p1 = true;
    int cut_count = 0;
    for (int ii = 0; ii<MAX_NUM; ii++){
        if(from_p1){
            offspring->_sequence.set(ii,p1->_sequence[ii]);
        }
        else{
            offspring->_sequence.set(ii,p2->_sequence[ii]);
        }
        if(cut_count < n && ii == cut_points[cut_count]){
            cut_count++;
            from_p1 = !(from_p1);
            while(cut_count< n && cut_points[cut_count] == cut_points[cut_count-1]){
                cut_count++;
//                cout<<"Skip!!"<<endl;
            }
//            cout<<"Find "<< cut_count<<"-th cut point: "<<cut_points[cut_count-1]<< " / from p1: "<< from_p1<<endl;
        }
    }
    /*
    /// Need to Check
    int point = rand()%MAX_NUM;
    for (int i=0; i<MAX_NUM; i++){
        if(i<point){
            offspring->_sequence.set(i,p1->_sequence[i]);
        }else{
            offspring->_sequence.set(i,p2->_sequence[i]);
        }
    }*/
}


/// Mutation
void mutation(Chromosome* offspring){
    bitset<L> seq = offspring->_sequence;
    double p = ((double) rand() / (RAND_MAX));
    for(int i=0; i<MAX_NUM; i++){
        if (p < MUTATION_RATE){
            seq.flip(i);
        }
    }
    offspring->_sequence = seq;
}

/// Replacement
void replace_worse(vector<Chromosome*>* offsprings, vector<Chromosome*>* population){
    //offsprings 의 사이즈 만큼 숫자를 골라서 (상위 2는 제외) 하위의 것들을 모두 교체
    int size = offsprings->size();
    for(int i=0; i<size; i++){
        Chromosome* old_c = population->at(i);
        Chromosome* new_c = offsprings->at(i);
//        cout<<"NEW: "<<new_c->_sequence<<endl;
        old_c->_sequence = new_c->_sequence;
        old_c->_score = new_c ->_score;
        //population->at(i) = offsprings->at(i);
    }

}

void replace_elitism(vector<Chromosome*>* offsprings, vector<Chromosome*>* population){
    int size = offsprings->size();
    int top_reserve = int(POPULATION_SIZE * ELITISM_RATE);
    for (int i=0; i<size; i++){
        int r = rand()%POPULATION_SIZE;
        while(r > POPULATION_SIZE-top_reserve){
            r = rand()%POPULATION_SIZE;
        }
        Chromosome* old_c = population->at(r);
        Chromosome* new_c = offsprings->at(i);
        old_c->_sequence = new_c->_sequence;
        old_c->_score = new_c -> _score;
    }
}

void replace_hybrid(Chromosome* offspring, Chromosome* p1, Chromosome* p2, Chromosome* worst){
//    cout<< "before/ p1: "<<p1->_score<< "/ p2: "<<p2->_score<< " / worst: "<<worst->_score<<endl;

    Chromosome* worse_p = p1;
    if(p1->_score > p2->_score){
        worse_p = p2;
    }

    if(offspring->_score >= worse_p->_score){
        worse_p->_sequence = offspring->_sequence;
        worse_p->_score = offspring -> _score;
    }
    else{
        worst->_sequence = offspring->_sequence;
        worst->_score = offspring -> _score;
    }
//    cout<< "after/ p1: "<<p1->_score<< "/ p2: "<<p2->_score<< " / worst: "<<worst->_score<<endl;

}

/// local optimization
void local_optimize_one_chrom(Chromosome* chrom, GraphHandler* gh){
    /// To Do
    //int init_score = chrom->_score;
    bitset<L> best_seq = chrom->_sequence;
    int best_score = chrom->_score;
    bool improved = true;
    while (improved){
        Chromosome* origin_chrom = new Chromosome();
        origin_chrom->_sequence =  best_seq;
        origin_chrom->_score = best_score;
        improved = false;
        for (int i=0; i<MAX_NUM; i++){
            int new_score = get_flipped_score(origin_chrom, i, gh);
            if(best_score < new_score) {
                best_seq = origin_chrom->_sequence.flip(i);
                origin_chrom->_sequence.flip(i);
                best_score = new_score;
//                cout<<"improved!!!"<<endl;
                improved = true;
//                i=0;
            }
        }
        delete(origin_chrom);
    }
//    cout<<"One local_optimize done"<<endl;
    chrom->_sequence = best_seq;
    chrom->_score = best_score;
    /*
    if(init_score<best_score){
        cout<<"LOCAL OPTIMIZATION WAS SUCCESSFUL"<<endl;
    }*/
}

void do_local_optimize(vector<Chromosome*>* population, GraphHandler* gh){
    int optimize_num = int(POPULATION_SIZE * OPTIMIZE_RATIO);
    vector<Chromosome*>::reverse_iterator riter(population->rbegin());
    riter = population->rbegin();
    for(int i=0; i<optimize_num; i++){
        local_optimize_one_chrom((*riter), gh);
        ++riter;
    }
}

void do_local_optimize_random(vector<Chromosome*>* population, GraphHandler* gh){
    int optimize_num = int(POPULATION_SIZE * OPTIMIZE_RATIO);
    vector<Chromosome*>::reverse_iterator riter(population->rbegin());
    riter = population->rbegin();
    for(int i=0; i<optimize_num; i++){
        if (i<optimize_num*0.8){
            local_optimize_one_chrom((*riter), gh);
            ++riter;
        }
        else{
            int r = rand()%POPULATION_SIZE;
            local_optimize_one_chrom(population->at(r), gh);
        }

    }
}

int find_maximum(int* arr, int size, bool* isLocked){
    int max=-9999999;
    int max_id=-1;
    for(int i=0; i<size; i++){
        if ( (!isLocked[i]) && max < arr[i]){
            max = arr[i];
            max_id = i;
        }
    }
    return max_id;
}


int lg_chain_analysis(int* lg_chain, int size, int* ptr_sum){
    int max = -999999999;
    int max_id = -1;
    int sum = 0;
    for(int idx = 0; idx<size; idx++){
        sum += lg_chain[idx];
        if (sum > max){
            max = sum;
            max_id = idx;
        }
    }
    ptr_sum[0] = max;
    assert(max_id != -1);
//    cout<<"computed sum is "<<max<<endl;
    return max_id;
}

void max_locked_gain(Chromosome* chrom, GraphHandler* gh){
    gh->compute_score(chrom);
    int origin_score =  chrom->_score;
    bool improved = true;
    while (improved){
        /// Initialize
        int locked_gain[MAX_NUM];
        bool isLocked[MAX_NUM];
        for(int idx = 0; idx<MAX_NUM; idx++){
            locked_gain[idx] = -999999; // 0?
            isLocked[idx] = false;
        }
        /// Find the initial vertex to start ( the vertex with the max gain)

        int max_gain = -999999;
        int init_vertex = -1;
        for(int v_id = 0; v_id <MAX_NUM; v_id++){
            int temp_gain = gh->compute_gain(chrom, v_id);
            if (temp_gain > max_gain){
                max_gain = temp_gain;
                init_vertex = v_id;
            }
        }
        isLocked[init_vertex] = true;
        chrom->_sequence.flip(init_vertex);
//        cout<<"Gain: "<<max_gain<<"/ vertex: "<<init_vertex<<endl;
        list<Edge *>::iterator iter;
        int target = init_vertex;

        improved = false;
        int locked_gain_chain[MAX_NUM];
        int vertex_chain[MAX_NUM];

        locked_gain_chain[0] = 0;
        vertex_chain[0] = init_vertex;
        /// Fill locked gain chain
        for (int ii =1; ii<MAX_NUM; ii ++){
            locked_gain_chain[ii] = -9999999;
            vertex_chain[ii] = -1;
        }

        for (int it = 1; it<MAX_NUM; it++){
            for (iter = gh->adj_mat[target].begin(); iter != gh->adj_mat[target].end(); iter++) {
                /// Update lock gain for adjacent vertices
                Edge *elem = (*iter);
                int that_idx;
                if (elem->_from == target){
                    that_idx = elem->_to;
                }else{
                    that_idx = elem->_from;
                }
                int tmp_locked_gain = gh->compute_locked_gain(chrom, that_idx, isLocked);
                locked_gain[that_idx] = tmp_locked_gain;
            }

            int max_locked_gain_idx = find_maximum(locked_gain, MAX_NUM, isLocked);
            int max_locked_gain = locked_gain[max_locked_gain_idx];
//            cout<< max_locked_gain_idx <<"/ max locked gain: "<< max_locked_gain<<endl;
//            assert(max_locked_gain>=0);
            locked_gain_chain[it] = max_locked_gain;
            vertex_chain[it] = max_locked_gain_idx;

            target = max_locked_gain_idx;
            isLocked[target] = true;
//            cout<<"target: "<<target<<endl;
            chrom->_sequence.flip(target);
        }
        int ptr_sum[1] ={-9999};
        int chain_max_id = lg_chain_analysis(locked_gain_chain, MAX_NUM, ptr_sum);

//        cout<<chain_max_id<< "/ "<<ptr_sum[0]<<endl;
        for(int n=MAX_NUM; n > chain_max_id; n--){
            int to_flip = vertex_chain[n];
            chrom->_sequence.flip(to_flip);
        }
        gh->compute_score(chrom);
        if(origin_score < chrom->_score){
            origin_score = chrom->_score;
            improved = true;
//            cout<< "improved, new score; "<<origin_score<<endl;
        }
    }
}

void do_local_optimize_lg(vector<Chromosome*>* population, GraphHandler* gh){
    int optimize_num = int(POPULATION_SIZE * OPTIMIZE_RATIO);
    vector<Chromosome*>::reverse_iterator riter(population->rbegin());
    riter = population->rbegin();
    for(int i=0; i<optimize_num; i++){
        max_locked_gain((*riter), gh);
        ++riter;
    }
}

int get_best_score(vector<Chromosome*>* population){
    Chromosome* best = population->back();
    return best->_score ;
}

int get_worst_score(vector<Chromosome *> *population) {
    Chromosome *worst = population->front();
    return worst->_score;
}

int get_median_score(vector<Chromosome *> *population) {
    Chromosome *median = population->at(POPULATION_SIZE / 2);
    return median->_score;
}
/// For Island GA
void do_one_generation(vector<Chromosome *> *population, GraphHandler* gh){
    /// return GA results with sorted version
    int xover_per_generation = int(POPULATION_SIZE * XOVER_RATIO);
    for (int count = 0; count < xover_per_generation; ++count) {
        sort(population->begin(), population->end(), compare);
        Chromosome *offspring = new Chromosome();
        // Selection
        int p1 = select_random();
        int p2 = p1;
        while (p1 == p2) {
            p2 = select_random();
        }
        // Xover
        n_point_xover(int(MAX_NUM/100), offspring, population->at(p1), population->at(p2), gh);
        // Mutation
//            MUTATION_RATE = (MAX_MUTATION_RATE - MIN_MUTATION_RATE) / (TIME_LIMIT) * (remain) + 0.001; // annealing
        mutation(offspring);
        local_optimize_one_chrom(offspring,gh);
//        max_locked_gain(offspring, gh);
//            get score and regularize the new offspring
        regularize(offspring, gh);
        get_score(offspring, gh);
        replace_hybrid(offspring, population->at(p1), population->at(p2), population->front());
        delete(offspring);
    }
    sort(population->begin(), population->end(), compare);
}


void move_best_to_neighbor(vector<Chromosome *> *from_pop, vector<Chromosome *> *to_pop){
    /// Assume from_pop and to_pop is sorted
    /// Replace to_pop's worst chromosome with form_pop's best chromosome.
    Chromosome* best = from_pop->back();
    Chromosome* worst = to_pop->front();
    worst->_sequence = best->_sequence;
    worst->_score = best->_score;
}

Chromosome* get_best_in_all_island(vector<Chromosome *> *population1, vector<Chromosome *> *population2,vector<Chromosome *> *population3,vector<Chromosome *> *population4,vector<Chromosome *> *population5,GraphHandler* gh){
    ///Find the best chromosome in all island!
    Chromosome* champ1 = population1->back();
    Chromosome* champ2 = population2->back();
    Chromosome* champ3 = population3->back();
    Chromosome* champ4 = population4->back();
    Chromosome* champ5 = population5->back();
    Chromosome* real_champ = new Chromosome();
    int max_score = -9999999;
    if (max_score < champ1->_score){
        real_champ->_sequence = champ1->_sequence;
        real_champ->_score = champ1->_score;
        max_score = champ1->_score;
    }
    if (max_score < champ2->_score){
        real_champ->_sequence = champ2->_sequence;
        real_champ->_score = champ2->_score;
        max_score = champ2->_score;
    }
    if (max_score < champ3->_score){
        real_champ->_sequence = champ3->_sequence;
        real_champ->_score = champ3->_score;
        max_score = champ3->_score;
    }
    if (max_score < champ4->_score){
        real_champ->_sequence = champ4->_sequence;
        real_champ->_score = champ4->_score;
        max_score = champ4->_score;
    }
    if (max_score < champ5->_score){
        real_champ->_sequence = champ5->_sequence;
        real_champ->_score = champ5->_score;
        max_score = champ5->_score;
    }
    return real_champ;

}



void print_population_status(vector<Chromosome *> *population){
    float converge = how_converge(population);
    int best_score = get_best_score(population);
    int ws = get_worst_score(population);
    int ms = get_median_score(population);
        cout << "/ best_score: " << best_score << "/ median_score: " << ms << "/ worst_score: "
        << ws << "/ converge: " << converge << endl;
}




#endif //C13_GA2019_METHODS_H
