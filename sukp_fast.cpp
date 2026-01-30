#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <chrono>

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
    int element_count = 0;
    double profit_per_element = 0.0;
};

struct ItemOrder {
    int idx = 0;
    double key = 0.0;
};

static double tiempo_segundos(){
    using clock = std::chrono::steady_clock;
    static const auto t0 = clock::now();
    auto now = clock::now();
    std::chrono::duration<double> diff = now - t0;
    return diff.count();
}

static void leer_instancia(const string &filename,int &m, int &n, double &capacity,vector<double> &item_profits, vector<double> &element_weights, vector<vector<int>> &matrix,  vector<int> &element_counts)
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

static void initialize_items_and_elements( int m, int n,const vector<double> &item_profits,const vector<double> &element_weights,const vector<vector<int>> &matrix, const vector<int> &element_counts, vector<Element> &elements,vector<Item> &items) {

    elements.assign(n, Element{});
    for(int j = 0; j < n; j++){
        elements[j].id = j;
        elements[j].weight = element_weights[j];
        elements[j].efficiency = 0.0;
    }

    items.assign(m, Item{});
    for (int i = 0; i < m; i++){
        items[i].id = i;
        items[i].profit = item_profits[i];
        items[i].element_count = element_counts[i];

        items[i].required_elements.clear();
        items[i].required_elements.reserve(element_counts[i]);

        double total_weight = 0.0;
        int count = 0;

        for(int j = 0; j < n; j++){
            if(matrix[i][j] == 1){
                items[i].required_elements.push_back(j);
                total_weight += elements[j].weight;
                count++;

                if(element_counts[i] > 0){
                    elements[j].efficiency += item_profits[i] / (double)element_counts[i];
                }
            }
        }

        if(count != items[i].element_count){
            items[i].element_count = count;
        }

        items[i].profit_per_element = (total_weight > 0.0)
                                        ? (item_profits[i] / total_weight)
                                        : 0.0;
    }
    for(int j = 0; j < n; j++){
        if(elements[j].weight > 0.0){
            elements[j].efficiency /= elements[j].weight;
        }
    }
}

static double deterministic_fast_fast( int m, int n,const vector<Element> &elements,const vector<Item> &items,double capacity, vector<unsigned char> &used_elements)
{
    double current_weight = 0.0;
    double total_profit  = 0.0;

    used_elements.assign(n, 0);

    // ordenar items por profit_per_element descendente
    vector<ItemOrder> ord(m);
    for(int i = 0; i < m; i++){
        ord[i].idx = i;
        ord[i].key = items[i].profit_per_element;
    }

    sort(ord.begin(), ord.end(), [](const ItemOrder &a, const ItemOrder &b){ return a.key > b.key; });

    // recorrer items greedy
    for(int t = 0; t < m; t++){
        int item_id = ord[t].idx;
        const Item &it = items[item_id];

        double extra_weight = 0.0;

        // calcular peso adicional si agrego este item
        for(int elem_id : it.required_elements){
            if(!used_elements[elem_id]){
                extra_weight += elements[elem_id].weight;
                if(current_weight + extra_weight > capacity){
                    extra_weight = -1.0;
                    break;
                }
            }
        }

        if(extra_weight < 0.0) continue; // no cabe

        // marcar elementos usados
        for(int elem_id : it.required_elements){
            if(!used_elements[elem_id]) used_elements[elem_id] = 1;
        }

        current_weight += extra_weight;
        total_profit   += it.profit;
    }

    return total_profit;
}

int main(){
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
            initialize_items_and_elements(m, n, profits, weights, matrix, elem_count, elements, items);

            vector<unsigned char> used_elements;

            double best = 0.0, sum = 0.0;

            double t1 = tiempo_segundos();

            for(int r = 0; r < RUNS; r++){
                double val = deterministic_fast_fast(m, n, elements, items, cap, used_elements);
                if (r == 0 || val > best) best = val;
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