//
// Created by Ann Kod on 08/10/2025.
//

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <iostream>
#include <fstream>
#include <string>


class parameters {
    int n; //amount of atoms on 1/3
    double m; // mass i guess
    int e; //funno
    float R; //radius of cube
    double f; //dunno
    double L; //raius of vessel
    float a; //cube siza
    long t_zero; //temperatur start
    double tau; //time step
    int s_o;
    long s_d;
    int s_out;
    int s_xyz;
    public:
    parameters()
        : n(5), m(39.948), e(1), R(0.38), f(10000), L(2.3), a(0.38), t_zero(100), tau(0.002), s_o(500), s_d(8000), s_out(10), s_xyz(10) {}

    parameters(const std::string& filename);
    ~parameters(){}
    int get_n() { return n; }
    double get_m() { return m; }
    int get_ee() { return e; }
    float get_r() { return R; }
    double get_f() { return f; }
    double get_l() { return L; }
    float get_a() { return a; }
    long get_t_zero() { return t_zero; }
    double get_tau() { return tau; }
    int get_s_o() { return s_o; }
    long get_s_d() { return s_d; }
    int get_s_out() { return s_out; }
    int get_s_xyz() { return s_xyz; }

    void set_n(int value) { n = value; }
    void set_m(double value) { m = value; }
    void set_e(double value) { e = value; }
    void set_R(double value) { R = value; }
    void set_f(double value) { f = value; }
    void set_L(double value) { L = value; }
    void set_a(double value) { a = value; }
    void set_t_zero(double value) { t_zero = value; }
    void set_tau(double value) { tau = value; }
    void set_S_o(int value) { s_o = value; }
    void set_S_d(int value) { s_d = value; }
    void set_S_out(int value) { s_out = value; }
    void set_S_xyz(int value) { s_xyz = value; }

};



#endif //PARAMETERS_H
