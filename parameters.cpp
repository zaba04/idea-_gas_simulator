//
// Created by Ann Kod on 08/10/2025.
//

#include "parameters.h"

#include <sstream>

parameters::parameters(const std::string& filename) {
    // Ustaw wartości domyślne
    n = 5;
    m = 39.948;
    e= 1.0;
    R = 0.38;
    f = 10000.0;
    L = 2.3;
    a = 0.38;
    t_zero = 100.0;
    tau = 0.002;
    s_o = 500;
    s_d = 8000;
    s_out = 10;
    s_xyz = 10;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open file " << filename
                  << ". Using default parameters." << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;

        if (iss >> key) {
            if (key == "n") iss >> n;
            else if (key == "m") iss >> m;
            else if (key == "epsilon" || key == "e") iss >> e;
            else if (key == "R") iss >> R;
            else if (key == "f") iss >> f;
            else if (key == "L") iss >> L;
            else if (key == "a") iss >> a;
            else if (key == "T0") iss >> t_zero;
            else if (key == "tau") iss >> tau;
            else if (key == "S_o") iss >> s_o;
            else if (key == "S_d") iss >> s_d;
            else if (key == "S_out") iss >> s_out;
            else if (key == "S_xyz") iss >> s_xyz;
        }
    }
    file.close();
}
