#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <random>
#include <chrono>
#include <stdexcept>
#include <cstdint>
#include <cstdio>
#include <algorithm>
using namespace std;

struct Instance {
    int m = 0;
    int n = 0;
    double capacity = 0.0;
    vector<double> item_profits;      // size m
    vector<double> element_weights;   // size n
    vector<uint8_t> matrix_flat;      // size m*n, values 0/1

    inline uint8_t A(int i, int e) const{
        return matrix_flat[(size_t)i * (size_t)n + (size_t)e];
    }
};

struct Solution {
    vector<uint8_t> selected_items;     // size m
    vector<uint8_t> selected_elements;  // size n
    double used_capacity = 0.0;
    double total_profit = 0.0;
};

// Leer instancia.
static void read_instance(const string &filename, Instance &inst){
    ifstream file(filename);
    if(!file){
        throw runtime_error("Error al abrir el archivo: " + filename);
    }

    file >> inst.m >> inst.n >> inst.capacity;
    inst.item_profits.assign(inst.m, 0.0);
    inst.element_weights.assign(inst.n, 0.0);
    inst.matrix_flat.assign((size_t)inst.m * (size_t)inst.n, 0);

    for(int i = 0; i < inst.m; i++) file >> inst.item_profits[i];
    for(int e = 0; e < inst.n; e++) file >> inst.element_weights[e];

    for(int i = 0; i < inst.m; i++){
        for(int e = 0; e < inst.n; e++){
            int x; file >> x;
            inst.matrix_flat[(size_t)i * (size_t)inst.n + (size_t)e] = (uint8_t)(x != 0);
        }
    }
}

// Iniciar la solucion, unicamente al inicio.
static void init_solution(const Instance &inst, Solution &sol){
    sol.selected_items.assign(inst.m, 0);
    sol.selected_elements.assign(inst.n, 0);
    sol.used_capacity = 0.0;
    sol.total_profit = 0.0;
}

// Lo mismo que init_solution, solo que se asume que el tamaño es correcto y es mas eficiente.
static void reset_solution(const Instance &inst, Solution &sol) {
    fill(sol.selected_items.begin(), sol.selected_items.end(), 0);
    fill(sol.selected_elements.begin(), sol.selected_elements.end(), 0);
    sol.used_capacity = 0.0;
    sol.total_profit = 0.0;
}

// copia la solucion, separado de la original
static void copy_solution(const Solution &src, Solution &dst) {
    dst.selected_items = src.selected_items;
    dst.selected_elements = src.selected_elements;
    dst.used_capacity = src.used_capacity;
    dst.total_profit = src.total_profit;
}

// Recalcula consistencia: item esta seleccionado si TODOS sus elementos requeridos estan seleccionados
static void finalize_solution(const Instance &inst, Solution &sol) {
    sol.total_profit = 0.0;

    for(int i = 0; i < inst.m; i++){
        bool ok = true;
        for(int e = 0; e < inst.n; e++){
            if(inst.A(i, e) && !sol.selected_elements[e]){ ok = false; break;}
        }
        sol.selected_items[i] = ok ? 1 : 0;
        if(ok) sol.total_profit += inst.item_profits[i];
    }

    sol.used_capacity = 0.0;
    for(int e = 0; e < inst.n; e++){
        if(sol.selected_elements[e]) sol.used_capacity += inst.element_weights[e];
    }
}

// elem_items[e] = lista de items que usan el elemento e, para no tener que recorrer la matriz muchas veces
static vector<vector<int>> build_elem_items(const Instance &inst){
    vector<vector<int>> elem_items(inst.n);
    for(int e = 0; e < inst.n; e++){
        elem_items[e].reserve(64); // heurístico; no crítico
        for(int i = 0; i < inst.m; i++){
            if(inst.A(i, e)) elem_items[e].push_back(i);
        }
        elem_items[e].shrink_to_fit();
    }
    return elem_items;
}

/*
Construye una solucion inicial seleccionando elementos de forma grasp
En cada iteracion evalua elementos factibles segun ganancia inmediata y potencial por unida,
elige uno al azar desde una lista restringida de candidatos y actualiza
los ítems completados hasta que no sea posible agregar más elementos útiles.
*/
static void build_solution_randomized_elements(const Instance &inst,const vector<vector<int>> &elem_items,Solution &sol,double alpha,double rcl_alpha,mt19937 &rng){
    const int m = inst.m;
    const int n = inst.n;
    const double cap = inst.capacity;

    vector<int> missing_count(m, 0);
    vector<double> score(n, 0.0);
    vector<int> rcl;
    rcl.reserve(n);

    // missing_count[i] = cuantos elementos requiere el item i
    for(int i = 0; i < m; i++){
        int cnt = 0;
        for(int e = 0; e < n; e++) if(inst.A(i, e)) cnt++;
        missing_count[i] = cnt;
    }

    reset_solution(inst, sol);

    // items que no requieren elementos = ganancia gratis
    for(int i = 0; i < m; i++){
        if(missing_count[i] == 0){
            sol.selected_items[i] = 1;
            sol.total_profit += inst.item_profits[i];
        }
    }

    while (true){
        double best_score = 0.0;
        bool any = false;

        for(int e = 0; e < n; e++){
            if(sol.selected_elements[e]){ score[e] = 0.0; continue;}

            double w = inst.element_weights[e];
            if(sol.used_capacity + w > cap){ score[e] = 0.0; continue;}

            double immediate = 0.0;
            double potential = 0.0;

            // solo items que usan el elemento e
            for(int i : elem_items[e]){
                int miss = missing_count[i];
                if (miss <= 0) continue;

                if (miss == 1) immediate += inst.item_profits[i];
                else           potential += inst.item_profits[i] / (double)miss;
            }

            if(immediate == 0.0 && potential == 0.0){ score[e] = 0.0; continue;}

            double s = (immediate + alpha * potential) / w;
            score[e] = s;

            if(s > best_score) best_score = s;
            if(s > 0.0) any = true;
        }

        if(!any || best_score <= 0.0) break;

        double threshold = rcl_alpha * best_score;
        rcl.clear();
        for(int e = 0; e < n; e++){
            if(!sol.selected_elements[e] && score[e] >= threshold){
                rcl.push_back(e);
            }
        }
        if(rcl.empty()) break;

        uniform_int_distribution<int> dist(0, (int)rcl.size() - 1);
        int chosen_e = rcl[dist(rng)];

        sol.selected_elements[chosen_e] = 1;
        sol.used_capacity += inst.element_weights[chosen_e];

        // actualizar missing_count y sumar ganancias cuando un item queda completo
        for(int i : elem_items[chosen_e]){
            if(missing_count[i] > 0){
                missing_count[i]--;
                if(missing_count[i] == 0){
                    sol.selected_items[i] = 1;
                    sol.total_profit += inst.item_profits[i];
                }
            }
        }
    }
}

// Loop principal de GRASP: genera múltiples soluciones y conserva la mejor.
static void grasp_elements(const Instance &inst,const vector<vector<int>> &elem_items,int iterations,double alpha,double rcl_alpha,Solution &best_sol,mt19937 &rng){
    Solution current;
    init_solution(inst, current);

    reset_solution(inst, best_sol);
    best_sol.total_profit = 0.0;

    for(int it = 0; it < iterations; it++){
        reset_solution(inst, current);

        build_solution_randomized_elements(inst, elem_items, current, alpha, rcl_alpha, rng);
        finalize_solution(inst, current);

        if(current.total_profit > best_sol.total_profit){
            copy_solution(current, best_sol);
        }
    }
}

int main(){    

    // leer mis instancias mas rapido
    const int GROUPS = 3; 
    const int PER_GROUP = 10;
    
    // Parametros 
    const int GRASP_ITERS = 550;
    const double ALPHA = 0.5;
    const double RCL_ALPHA = 0.8;

    // RNG
    std::random_device rd;
    std::mt19937 rng(rd());

    bool printed_header = false;

    for(int g = 1; g <= GROUPS; g++){
        for(int k = 1; k <= PER_GROUP; k++){

            char buf[64];
            snprintf(buf, sizeof(buf), "%d_%d.txt", g, k);
            string filename(buf);

            try {
                Instance inst;
                read_instance(filename, inst);

                Solution best_sol;
                init_solution(inst, best_sol);

                auto elem_items = build_elem_items(inst);

                auto start = chrono::high_resolution_clock::now();

                grasp_elements(inst, elem_items, GRASP_ITERS, ALPHA, RCL_ALPHA, best_sol, rng);
                finalize_solution(inst, best_sol);

                auto end = chrono::high_resolution_clock::now();
                double tiempo = chrono::duration<double>(end - start).count();

                int count_items = 0, count_elems = 0;
                for(int i = 0; i < inst.m; i++) if (best_sol.selected_items[i]) count_items++;
                for(int e = 0; e < inst.n; e++) if (best_sol.selected_elements[e]) count_elems++;

                if(!printed_header){
                    cout << "---------------------------------------------------------------------------------------------\n";
                    cout << "| " << left << setw(10) << "Instancia"
                         << " | " << right << setw(10) << "Items"
                         << " | " << right << setw(10) << "Elems"
                         << " | " << right << setw(12) << "Ganancia"
                         << " | " << right << setw(12) << "Tiempo (s)"
                         << " |\n";
                    cout << "---------------------------------------------------------------------------------------------\n";
                    printed_header = true;
                }

                cout << "| " << left << setw(10) << filename
                     << " | " << right << setw(10) << count_items
                     << " | " << right << setw(10) << count_elems
                     << " | " << right << setw(12) << fixed << setprecision(2) << best_sol.total_profit
                     << " | " << right << setw(12) << fixed << setprecision(4) << tiempo
                     << " |\n";

            } catch(const exception &ex){
                cerr << "Error en " << filename << ": " << ex.what() << "\n";
                return 1;
            }
        }
    }

    return 0;
}