#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <cstdint>
using namespace std;

struct Element {
    int id = 0;
    double weight = 0.0;
    double efficiency = 0.0;
};

struct Item {
    int id = 0;
    double profit = 0.0;
    vector<int> required_elements;
    double profit_per_element = 0.0;
};

// Para el tiempo
static double tiempo_segundos() {
    using clock = std::chrono::steady_clock;
    static const auto t0 = clock::now();
    auto now = clock::now();
    std::chrono::duration<double> diff = now - t0;
    return diff.count();
}

//Leer instancia
static void leer_instancia(const string &filename, int &m, int &n, double &capacity, vector<double> &item_profits, vector<double> &element_weights, vector<vector<int>> &matrix, vector<int> &element_counts)
{
    ifstream file(filename);
    if(!file.is_open()){
        cerr << "Error al abrir archivo: " << filename << "\n";
        exit(1);
    }

    if(!(file >> m)){ cerr << "Error leyendo m en " << filename << "\n"; exit(1); }
    if(!(file >> n)){ cerr << "Error leyendo n en " << filename << "\n"; exit(1); }
    if(!(file >> capacity)){ cerr << "Error leyendo capacidad en " << filename << "\n"; exit(1); }

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

// inicializar items 
static void initialize_items_and_elements(int m, int n, const vector<double> &profits, const vector<double> &weights, const vector<vector<int>> &matrix, const vector<int> &elem_count, vector<Element> &elements, vector<Item> &items)
{
    elements.assign(n, Element{});
    for(int j = 0; j < n; j++){
        elements[j].id = j;
        elements[j].weight = weights[j];
        elements[j].efficiency = 0.0;
    }

    items.assign(m, Item{});
    for(int i = 0; i < m; i++){
        items[i].id = i;
        items[i].profit = profits[i];

        items[i].required_elements.clear();
        items[i].required_elements.reserve(elem_count[i]);

        double total_weight = 0.0;
        for(int j = 0; j < n; j++){
            if(matrix[i][j] == 1){
                items[i].required_elements.push_back(j);
                total_weight += weights[j];

                if(elem_count[i] > 0){
                    elements[j].efficiency += profits[i] / (double)elem_count[i];
                }
            }
        }

        items[i].profit_per_element = (total_weight > 0.0) ? (profits[i] / total_weight) : 0.0;
    }

    for(int j = 0; j < n; j++){
        if(elements[j].weight > 0.0) elements[j].efficiency /= elements[j].weight;
    }
}

static double greedy_knapsack(int m, int n,  const vector<Element> &elements, const vector<Item> &items, double capacity, vector<int> &solution)
{
    solution.assign(n, 0);

    double total_profit = 0.0;
    double current_weight = 0.0;

    // 0 = no usado, 1 = elegido, -1 = imposible
    vector<int> item_state(m, 0);

    while(true){
        int best_item = -1;
        double best_ratio = -1.0;
        bool best_is_free = false;

        for(int i = 0; i < m; i++){
            if(item_state[i] != 0) continue;
            if(items[i].profit <= 0.0) continue;

            double extra_weight = 0.0;

            for(int elem : items[i].required_elements){
                if(!solution[elem]){
                    extra_weight += elements[elem].weight;
                    if(current_weight + extra_weight > capacity) break;
                }
            }

            if(current_weight + extra_weight > capacity){
                item_state[i] = -1;
                continue;
            }

            if(extra_weight == 0.0){
                if(!best_is_free ||
                    (best_item == -1) ||
                    items[i].profit > items[best_item].profit){
                    best_item = i;
                    best_is_free = true;
                }
            } else if(!best_is_free){
                double ratio = items[i].profit / extra_weight;
                if(ratio > best_ratio){
                    best_ratio = ratio;
                    best_item = i;
                }
            }
        }

        if(best_item == -1) break;

        item_state[best_item] = 1;

        for(int elem : items[best_item].required_elements){
            if(!solution[elem]){
                solution[elem] = 1;
                current_weight += elements[elem].weight;
            }
        }

        total_profit += items[best_item].profit;
    }

    return total_profit;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    const int GROUPS = 3;
    const int PER_GROUP = 10;
    const int RUNS = 1;

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
                double val = greedy_knapsack(m, n, elements, items, cap, solution);

                if (r == 0 || val > best) best = val;
                sum += val;
            }

            double t2 = tiempo_segundos();

            double avg = sum / RUNS;
            double total_time = t2 - t1;
            double per_run = total_time / RUNS;

            if (!printed_header) {
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