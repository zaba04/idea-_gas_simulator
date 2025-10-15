//
// Created by Ann Kod on 15/10/2025.
//

#include "algorytm_2.h"


void algorytm_2::algo_2(std::vector<atom>& atoms, parameters& params, system_params& system) {
    double tmpVS = 0.0;
    double tmpVP = 0.0;
    double tmpFs = 0.0;

    for (size_t i = 0; i < atoms.size(); i++) {
        //for that calculates Vs ad Vp and Fsi Fpi
        double tmpFsi = 0.0;
        double tmpFpi = 0.0;

        double tmpR = atoms[i].get_R();

        if (tmpR >= params.get_l()) {
            tmpVS += 0.5 * params.get_f() * (tmpR - params.get_l()) * (tmpR - params.get_l());
            //force from walls ov vessel
            tmpFsi = params.get_f() * (params.get_l() - tmpR);
            tmpFs += fabs(tmpFsi);
        }

        for (size_t j = 0; j < i; j++) {
            //Vp section

            double dx = atoms[i].get_x() - atoms[j].get_x();
            double dy = atoms[i].get_y() - atoms[j].get_y();
            double dz = atoms[i].get_z() - atoms[j].get_z();
            double tmpRjk = sqrt(dx*dx + dy*dy + dz*dz);

            if (tmpRjk > 1e-10) {  // Unikaj dzielenia przez 0
                tmpVP += params.get_ee() * (pow(params.get_r()/tmpRjk, 12) - 2.0 * pow(params.get_r()/tmpRjk, 6));

                double force_magnitude = 12.0 * params.get_ee() / tmpRjk *
                    (pow(params.get_r()/tmpRjk, 12) - pow(params.get_r()/tmpRjk, 6));
                tmpFpi += force_magnitude;
            }
        }

        atoms[i].set_F(tmpFsi, tmpFpi);
    }
    system.setV(tmpVS+tmpVP);
    system.setP(1./(4*M_PI*params.get_l()*params.get_l()) * tmpFs);
}
