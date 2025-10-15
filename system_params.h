//
// Created by Ann Kod on 15/10/2025.
//

#ifndef SYSTEM_PARAMS_H
#define SYSTEM_PARAMS_H
#include <vector>


class system_params {
    double V;
    double P;
    double H;
    double T;
    public:
    system_params() : V(0.0), P(0.0) {} ;
    ~system_params() = default;
    void set_params(double V, double P, double H, double T);
    void set_params(double V, double P);
    double get_V() {return V;};
    double get_P() {return P;};
    double get_H() {return H;};
    double get_T() {return T;};
    void setV(double V) {this->V = V;};
    void setP(double P) {this->P = P;};
};



#endif //SYSTEM_PARAMS_H
