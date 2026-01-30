#include <iostream>
#include <fstream> 
#include <vector>   
#include <string>   
#include <algorithm>  
#include <iomanip>    
#include <chrono>    
#include <random>
#include <ctime>
#include <cstdint>
using namespace std;

struct Element {
    int id = 0;
    double weight = 0.0;
};

struct Item {
    int id = 0;
    double profit = 0.0;
    vector<int> required_elements;
    double profit_per_element = 0.0;
};

struct Candidate {
    int item_id = 0;
    double score = 0.0;
};

static double tiempo_segundos(){
    using clock = std::chrono::steady_clock;
    static const auto t0 = clock::now();
    auto now = clock::now();
    std::chrono::duration<double> diff = now - t0;
    return diff.count();
}

// Leer instancia
static void leer_instancia(const string &filename,int &m, int &n, double &capacity, vector<double> &item_profits, vector<double> &element_weights, vector<vector<int>> &matrix,vector<int> &element_counts)
{
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error al abrir archivo: " << filename << "\n";
        exit(1);
    }

    if (!(file >> m)) { cerr << "Error leyendo m en " << filename << "\n"; exit(1); }
    if (!(file >> n)) { cerr << "Error leyendo n en " << filename << "\n"; exit(1); }
    if (!(file >> capacity)) { cerr << "Error leyendo capacidad en " << filename << "\n"; exit(1); }

    item_profits.assign(m, 0.0);
    element_weights.assign(n, 0.0);
    matrix.assign(m, vector<int>(n, 0));
    element_counts.assign(m, 0);

    for(int i = 0; i < m; i++){
        if(!(file >> item_profits[i])){
            cerr << "Error leyendo profits en " << filename << "\n";
            exit(1);
        }
    }

    for(int j = 0; j < n; j++){
        if(!(file >> element_weights[j])){
            cerr << "Error leyendo pesos en " << filename << "\n";
            exit(1);
        }
    }

    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            if(!(file >> matrix[i][j])){
                cerr << "Error leyendo matrix[" << i << "][" << j << "] en " << filename << "\n";
                exit(1);
            }
            if(matrix[i][j] == 1) element_counts[i]++;
        }
    }
}

// inicializar items y elementos 
static void initialize_items_and_elements(int m, int n, const vector<double> &profits, const vector<double> &weights, const vector<vector<int>> &matrix, const vector<int> &elem_count, vector<Element> &elements, vector<Item> &items)
{
    elements.assign(n, Element{});
    for(int j = 0; j < n; j++){
        elements[j].id = j;
        elements[j].weight = weights[j];
    }

    items.assign(m, Item{});
    for(int i = 0; i < m; i++){
        items[i].id = i;
        items[i].profit = profits[i];

        items[i].required_elements.clear();
        items[i].required_elements.reserve(elem_count[i]);

        for(int j = 0; j < n; j++){
            if(matrix[i][j] == 1) items[i].required_elements.push_back(j);
        }

        if(elem_count[i] > 0) items[i].profit_per_element = profits[i] / (double)elem_count[i];
        else items[i].profit_per_element = 0.0;
    }
}

// algoritmo principal 
static double probabilistic_algorithm(int m, int n, const vector<Element> &elements, const vector<Item> &items, double capacity, vector<int> &solution, uint32_t seed)
{
    mt19937 rng(seed);

    double current_weight = 0.0;
    solution.assign(n, 0);

    vector<int> blocked(m, 0);
    vector<Candidate> candidates;
    candidates.reserve(m);

    // PARAMETROS PARA EL ALGORITMO PROBABILIDAD (no se buscó los mejores)
    const double ALPHA = 0.7; // peso para ganancia marginal
    const int K_BEST = 5;     // elegir al azar entre los K mejores

    while(true){
        candidates.clear();

        for(int i = 0; i < m; i++){
            if(blocked[i]) continue;

            double extra_weight = 0.0;
            bool missing = false;

            for(int elem : items[i].required_elements){
                if(!solution[elem]){
                    missing = true;
                    extra_weight += elements[elem].weight;
                }
            }

            if(!missing){
                blocked[i] = 1;
                continue;
            }

            if(current_weight + extra_weight > capacity){
                blocked[i] = 1;
                continue;
            }

            double marginal = items[i].profit / extra_weight;
            double combined_score = ALPHA * marginal + (1.0 - ALPHA) * items[i].profit_per_element;

            candidates.push_back({i, combined_score});
        }

        if(candidates.empty()) break;

        sort(candidates.begin(), candidates.end(),[](const Candidate &a, const Candidate &b){ return a.score > b.score; });

        int limit = min((int)candidates.size(), K_BEST);
        uniform_int_distribution<int> dist(0, limit - 1);
        int chosen_item = candidates[dist(rng)].item_id;

        double added = 0.0;
        for(int elem : items[chosen_item].required_elements){
            if(!solution[elem]) {
                solution[elem] = 1;
                added += elements[elem].weight;
            }
        }

        current_weight += added;
        blocked[chosen_item] = 1;
    }

    double total_profit = 0.0;
    for(int i = 0; i < m; i++){
        bool full = true;
        for(int elem : items[i].required_elements){
            if(!solution[elem]) { full = false; break; }
        }
        if(full) total_profit += items[i].profit;
    }

    return total_profit;
}

int main(){
    const int GROUPS = 3;
    const int PER_GROUP = 10;
    const int RUNS = 100;       // corridas por instancia

    bool printed_header = false;

    for(int g = 1; g <= GROUPS; g++){
        for(int k = 1; k <= PER_GROUP; k++){

            string name = to_string(g) + "_" + to_string(k) + ".txt";

            int m = 0, n = 0;
            double cap = 0.0;
            vector<double> profits, weights;
            vector<vector<int>> matrix;
            vector<int> elem_count;

            leer_instancia(name, m, n, cap, profits, weights, matrix, elem_count);

            vector<Element> elements;
            vector<Item> items;
            vector<int> solution;

            initialize_items_and_elements(m, n, profits, weights, matrix, elem_count, elements, items);

            double best = 0.0, sum = 0.0;
            double t1 = tiempo_segundos();

            for(int r = 0; r < RUNS; r++){
                int inst_id = (g - 1) * PER_GROUP + k;

                uint32_t seed = (uint32_t)time(nullptr) + (uint32_t)(r * 123457) + (uint32_t)(inst_id * 99991);

                double val = probabilistic_algorithm(m, n, elements, items, cap, solution, seed);

                if(r == 0 || val > best) best = val;
                sum += val;
            }

            double t2 = tiempo_segundos();
            double avg = sum / RUNS;
            double total_time = t2 - t1;
            double per_run = total_time / RUNS;

            if(!printed_header){
                cout << "-------------------------------------------------------------------------------------------------------------\n";
                cout << "| " << left << setw(10) << "Instancia"
                     << " | " << setw(16) << "Mejor Ganancia"
                     << " | " << setw(12) << "Promedio"
                     << " | " << setw(17) << "Tiempo Total (s)"
                     << " | " << setw(20) << "Tiempo por corrida (s)"
                     << " |\n";
                cout << "-------------------------------------------------------------------------------------------------------------\n";
                printed_header = true;
            }

            cout << "| " << left << setw(10) << name
                 << " | " << right << setw(16) << fixed << setprecision(6) << best
                 << " | " << setw(12) << fixed << setprecision(6) << avg
                 << " | " << setw(17) << fixed << setprecision(6) << total_time
                 << " | " << setw(20) << fixed << setprecision(6) << per_run
                 << " |\n";
        }
    }

    return 0;
}