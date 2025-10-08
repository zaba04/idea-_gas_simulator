//
// Created by Ann Kod on 08/10/2025.
//

#include "parameters.h"

parameters::parameters(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file " << filename << std::endl;
        exit(1);
    }
    std::vector<double> tmpparams;
    double tmp;

    while (file >> tmp) {
        tmpparams.push_back(tmp);
    }

    n = tmpparams[0];
    m = tmpparams[1];
    e = tmpparams[2];
    R = tmpparams[3];
    f = tmpparams[4];
    L = tmpparams[5];
    a = tmpparams[6];
    t_zero = tmpparams[7];
    tau = tmpparams[8];
    s_o = tmpparams[9];
    s_d = tmpparams[10];
    s_out = tmpparams[11];
    s_xyz = tmpparams[12];

}
