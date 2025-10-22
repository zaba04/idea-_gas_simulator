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

    for (size_t i = 0; i <= atoms.size(); ++i) {
        algorytm_2.algo_2(atoms, parameters, system, i);
    }

    //dynamika
    double tmppx = 0.0;
    double tmppy = 0.0;
    double tmppz = 0.0;
    double tmpFx = 0.0;
    double tmpFy = 0.0;
    double tmpFz = 0.0;
    double tmpRx = 0.0;
    double tmpRy = 0.0;
    double tmpRz = 0.0;
    double Tmean = 0.0;
    double Pmean = 0.0;
    double Hmean = 0.0;
    double tau = parameters.get_tau();
    double m = parameters.get_m();

    for (int s = 0; s < parameters.get_s_d()+parameters.get_s_o(); s++) {
        double tmpEkin = 0.0;
        for (size_t i = 0; i < atoms.size(); i++) {
            tmppx = atoms[i].get_px();
            tmppy = atoms[i].get_py();
            tmppz = atoms[i].get_pz();
            tmpFy = atoms[i].get_FwallY() + atoms[i].get_FatomY();
            tmpFx = atoms[i].get_FwallX() + atoms[i].get_FatomX();
            tmpFz = atoms[i].get_FwallZ() + atoms[i].get_FatomZ();
            tmpRx = atoms[i].get_x();
            tmpRy = atoms[i].get_y();
            tmpRz = atoms[i].get_z();

            //momentas
            double tmppx12 = tmppx + 0.5*tmpFx*tau;
            double tmppy12 = tmppy + 0.5*tmpFy*tau;
            double tmppz12 = tmppz + 0.5*tmpFz*tau;

            //positions
            double tmprx = tmpRx + 1./m * tmppx12 * tau;
            double tmpry = tmpRy + 1./m * tmppy12 * tau;
            double tmprz = tmpRz + 1./m * tmppz12 * tau;

            atoms[i].set_p(tmppx12, tmppy12, tmppz12);
            atoms[i].set_pos(tmprx, tmpry, tmprz);

            //calculating new
            algorytm_2.algo_2(atoms, parameters, system, i);

            tmpFy = atoms[i].get_FwallY() + atoms[i].get_FatomY();
            tmpFx = atoms[i].get_FwallX() + atoms[i].get_FatomX();
            tmpFz = atoms[i].get_FwallZ() + atoms[i].get_FatomZ();

            double tmppx1 = tmppx12 + 0.5*tmpFx*tau;
            double tmppy1 = tmppy12 + 0.5*tmpFy*tau;
            double tmppz1 = tmppz12 + 0.5*tmpFz*tau;

            atoms[i].set_p(tmppx1, tmppy1, tmppz1);

            tmpEkin += (tmppx1*tmppx1 + tmppy1*tmppy1 + tmppz1*tmppz1) /(2 * m);
        }

        double tmpT = 2./(3.*k_b*atoms.size()) * tmpEkin;
        double tmpH = tmpEkin + system.get_V();
        system.setH(tmpH);
        system.setT(tmpT);

        if (s%parameters.get_s_out()==0) {
            //zapis t,H,P,V
        }
        if (s%parameters.get_s_xyz()==0) {
            //zapis współrzędnych
        }
        if (s >= parameters.get_s_o()) {
            Tmean += tmpT / s;
            Pmean += system.get_P() / s;
            Hmean += tmpH / s;
        }
    }


    return 0;
}