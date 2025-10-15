#include <iostream>
#include <cmath>
#include <random>
#include "atom.h"
#include "parameters.h"
#include "algorytm_2.h"


#define k_b 0.00831
int main() {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(std::nextafter(0.0, 1.0), 1.0);
    std::uniform_int_distribution<int> sig(-1, 1);


    parameters parameters;

    std::vector<atom> atoms;
    system_params system;

    algorytm_2 algorytm_2;

    //part setting start xyz values
    for (int i = 0; i < parameters.get_n(); ++i) {
        for (int j = 0; j < parameters.get_n(); ++j) {
            for (int k = 0; k < parameters.get_n(); ++k) {
                float tmpx = 0;
                float tmpy = 0;
                float tmpz = 0;

                double b0[3] = {parameters.get_a(), 0, 0};
                double b1[3] = {parameters.get_a()/2., parameters.get_a()*sqrt(3)/2. , 0};
                double b2[3] = {parameters.get_a()/2., parameters.get_a()*sqrt(3)/6. , parameters.get_a()*sqrt(2/3.)};

                float tmpp = (parameters.get_n()-1);

                tmpx = (i - tmpp/2) * b0[0] + (j - tmpp/2) * b1[0] + (k - tmpp/2) * b2[0];
                tmpy = (i - tmpp/2) * b0[1] + (j - tmpp/2) * b1[1] + (k - tmpp/2) * b2[1];
                tmpz = (i - tmpp/2) * b0[2] + (j - tmpp/2) * b1[2] + (k - tmpp/2) * b2[2];


                atoms.push_back(atom(tmpx, tmpy, tmpz));
                // std::cout << tmpx << "  " << tmpy << "  " << tmpz << std::endl;
            }
        }
    }

    double sumMomX = 0.;
    double sumMomY = 0.;
    double sumMomZ = 0.;
    //generating momenta
    for (size_t i = 0; i <= atoms.size(); ++i) {
        double tmpEx = -0.5 * (k_b * parameters.get_t_zero()* std::log(dis(gen)));
        double tmpEy = -0.5 * (k_b * parameters.get_t_zero()* std::log(dis(gen)));
        double tmpEz = -0.5 * (k_b * parameters.get_t_zero()* std::log(dis(gen)));

        int sign = sig(gen);
        double tmpMomx = sign * sqrt(2* parameters.get_m() * tmpEx);
        sumMomX += tmpMomx;
        sign = sig(gen);
        double tmpMomy = sign * sqrt(2* parameters.get_m() * tmpEy);
        sumMomY += tmpMomy;
        sign = sig(gen);
        double tmpMomz = sign * sqrt(2* parameters.get_m() * tmpEz);
        sumMomZ += tmpMomz;

        //sprawdzić artefakty układające się po krzyżu

        atoms[i].set_p(tmpMomx, tmpMomy, tmpMomz);
    }

    for (size_t i = 0; i <= atoms.size(); ++i) {
        double tmpx = atoms[i].get_px() - 1./atoms.size()*sumMomX;
        double tmpy = atoms[i].get_py() - 1./atoms.size()*sumMomY;
        double tmpz = atoms[i].get_pz() - 1./atoms.size()*sumMomZ;
        atoms[i].set_p(tmpx, tmpy, tmpz);
        // std::cout << tmpx << " " << tmpy << " " << tmpz << std::endl;
    }

    algorytm_2.algo_2(atoms, parameters, system);

    //
    // std::cout << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;
    //
    // std::cout << "V:    " << system.get_V() << std::endl;


    return 0;
}