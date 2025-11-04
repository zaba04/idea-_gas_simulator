//
// Created by Ann Kod on 15/10/2025.
//

#include <algorithm>
#include <vector>
#include "atom.h"
#include "parameters.h"
#include "system_params.h"
#include <cmath>

#ifndef ALGORYTM_2_H
#define ALGORYTM_2_H



class algorytm_2 {

    public:
    algorytm_2() = default;
    ~algorytm_2() = default;
    void calculate_forces_for_atom(std::vector<atom>& atoms, parameters& params, int i);
    double calculate_total_energy(std::vector<atom>& atoms, parameters& params);
    double calculate_pressure(std::vector<atom>& atoms, parameters& params);


};



#endif //ALGORYTM_2_H
