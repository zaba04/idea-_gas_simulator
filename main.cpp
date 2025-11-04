#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include "atom.h"
#include "parameters.h"
#include "algorytm_2.h"


#define k_b 0.00831

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(std::nextafter(0.0, 1.0), 1.0);
    std::uniform_int_distribution<int> sign(0, 1);



    parameters parameters;
    std::vector<atom> atoms;
    system_params system;
    algorytm_2 algorytm_2;

    std::ofstream param_file("params");
    if (!param_file.is_open()) {
        std::cerr << "Nie mogę utworzyć pliku params" << std::endl;
        return 1;
    }

    // MUSI BYĆ DOKŁADNIE 13 WARTOŚCI:
    param_file << parameters.get_n() << " "               // n
               << parameters.get_m() << " "               // mm
               << parameters.get_ee() << " "             // epseps
               << parameters.get_r() << " "               // RR
               << parameters.get_f() << " "               // ff
               << parameters.get_l() << " "        // L_sphere
               << parameters.get_a() << " "               // aa
               << parameters.get_t_zero() << " "          // TT
               << parameters.get_tau() << " "             // tautau
               << parameters.get_s_o() << " "             // So
               << parameters.get_s_d() << " "             // Sd
               << parameters.get_s_out() << " "           // Sout
               << parameters.get_s_xyz() << std::endl;    // Sxyz
    param_file.close();

    //part setting start xyz values
    for (int i = 0; i < parameters.get_n(); ++i) {
        for (int j = 0; j < parameters.get_n(); ++j) {
            for (int k = 0; k < parameters.get_n(); ++k) {
                float tmpx = 0;
                float tmpy = 0;
                float tmpz = 0;

                double b0[3] = {parameters.get_a(), 0, 0};
                double b1[3] = {parameters.get_a()/2., parameters.get_a()*sqrt(3.)/2. , 0};
                double b2[3] = {parameters.get_a()/2., parameters.get_a()*sqrt(3.)/6. , parameters.get_a()*sqrt(2./3.)};

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
    for (size_t i = 0; i < atoms.size(); ++i) {
        double tmpEx = -0.5 * (k_b * parameters.get_t_zero()* std::log(dis(gen)));
        double tmpEy = -0.5 * (k_b * parameters.get_t_zero()* std::log(dis(gen)));
        double tmpEz = -0.5 * (k_b * parameters.get_t_zero()* std::log(dis(gen)));

        int signe = 0;
        signe = (sign(gen) == 0) ? -1 : 1;
        double tmpMomx = signe * sqrt(2* parameters.get_m() * tmpEx);
        sumMomX += tmpMomx;
        signe = (sign(gen) == 0) ? -1 : 1;
        double tmpMomy = signe * sqrt(2* parameters.get_m() * tmpEy);
        sumMomY += tmpMomy;
        signe = (sign(gen) == 0) ? -1 : 1;
        double tmpMomz = signe * sqrt(2* parameters.get_m() * tmpEz);
        sumMomZ += tmpMomz;


        atoms[i].set_p(tmpMomx, tmpMomy, tmpMomz);
    }

    for (size_t i = 0; i < atoms.size(); ++i) {
        double tmpx = atoms[i].get_px() - 1./atoms.size()*sumMomX;
        double tmpy = atoms[i].get_py() - 1./atoms.size()*sumMomY;
        double tmpz = atoms[i].get_pz() - 1./atoms.size()*sumMomZ;
        atoms[i].set_p(tmpx, tmpy, tmpz);
        // std::cout << tmpx << " " << tmpy << " " << tmpz << std::endl;
    }

    for (size_t i = 0; i < atoms.size(); ++i) {
        algorytm_2.calculate_forces_for_atom(atoms, parameters, i);
    }
    double total_V = algorytm_2.calculate_total_energy(atoms, parameters);
    double total_P = algorytm_2.calculate_pressure(atoms, parameters);
    system.setV(total_V);
    system.setP(total_P);

    //opening files
    std::ofstream file_xyz("xyz.dat");
    std::ofstream file_out("out.dat");


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

        // Krok 1: p(t + τ/2) dla wszystkich atomów
        for (size_t i = 0; i < atoms.size(); i++) {
            double tmppx = atoms[i].get_px();
            double tmppy = atoms[i].get_py();
            double tmppz = atoms[i].get_pz();
            double tmpFx = atoms[i].get_FwallX() + atoms[i].get_FatomX();
            double tmpFy = atoms[i].get_FwallY() + atoms[i].get_FatomY();
            double tmpFz = atoms[i].get_FwallZ() + atoms[i].get_FatomZ();

            atoms[i].set_p(tmppx + 0.5*tmpFx*tau,
                           tmppy + 0.5*tmpFy*tau,
                           tmppz + 0.5*tmpFz*tau);
        }

        // Krok 2: r(t + τ) dla wszystkich atomów
        for (size_t i = 0; i < atoms.size(); i++) {
            double tmpRx = atoms[i].get_x();
            double tmpRy = atoms[i].get_y();
            double tmpRz = atoms[i].get_z();
            double tmppx12 = atoms[i].get_px();
            double tmppy12 = atoms[i].get_py();
            double tmppz12 = atoms[i].get_pz();

            atoms[i].set_pos(tmpRx + (1./m) * tmppx12 * tau,
                            tmpRy + (1./m) * tmppy12 * tau,
                            tmpRz + (1./m) * tmppz12 * tau);
        }

        // Oblicz nowe siły dla wszystkich atomów
        for (size_t i = 0; i < atoms.size(); i++) {
            algorytm_2.calculate_forces_for_atom(atoms, parameters, i);
        }

        // Oblicz energię i ciśnienie
        double total_V = algorytm_2.calculate_total_energy(atoms, parameters);
        double total_P = algorytm_2.calculate_pressure(atoms, parameters);
        system.setV(total_V);
        system.setP(total_P);

        // Krok 3: p(t + τ) i energia kinetyczna
        double tmpEkin = 0.0;
        for (size_t i = 0; i < atoms.size(); i++) {
            double tmppx12 = atoms[i].get_px();
            double tmppy12 = atoms[i].get_py();
            double tmppz12 = atoms[i].get_pz();
            double tmpFx = atoms[i].get_FwallX() + atoms[i].get_FatomX();
            double tmpFy = atoms[i].get_FwallY() + atoms[i].get_FatomY();
            double tmpFz = atoms[i].get_FwallZ() + atoms[i].get_FatomZ();

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

        if (s % parameters.get_s_out() == 0) {
            file_out << std::fixed << std::setprecision(6);
            file_out << s * tau << " "
                     << tmpT << " "
                     << tmpEkin << " "
                     << system.get_V() << " "
                     << (tmpEkin + system.get_V())
                     << std::endl;
        }

        if (s % parameters.get_s_xyz() == 0) {
            for (size_t i = 0; i < atoms.size(); i++) {
                file_xyz << atoms[i].get_x() << " " << atoms[i].get_y() << " " << atoms[i].get_z() << " ";
            }
            file_xyz << std::endl;
            file_xyz.flush();
        }

        if (s >= parameters.get_s_o()) {
            Tmean += tmpT;
            Pmean += system.get_P();
            Hmean += tmpH;
        }
    }


    int count_steps = parameters.get_s_d();
    Tmean /= count_steps;
    Pmean /= count_steps;
    Hmean /= count_steps;

    file_xyz.close();
    file_out.close();
    return 0;
}