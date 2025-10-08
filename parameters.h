//
// Created by Ann Kod on 08/10/2025.
//

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <iostream>
#include <fstream>
#include <string>


class parameters {
    int n; //amount of atoms
    int m; // i dunno know
    int e; //funno
    float R; //radius of cube
    double f; //dunno
    double L; //raius of vessel
    float a;
    long t_zero; //time start
    double tau; //time step
    int s_o;
    long s_d;
    int s_out;
    int s_xyz;
    public:
    parameters(int, int, int , float , double , double , float , long , double , int s_o, long s_d, int s_out, int s_xyz)
        : n(3), m(1), e(1), R(0.38), f(10000), L(1.2), a(0.38), t_zero(1000), tau(0.002), s_o(100), s_d(2000), s_out(10), s_xyz(10) {};
    parameters(const std::string& filename);
    ~parameters();
};



#endif //PARAMETERS_H
