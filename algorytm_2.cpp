//
// Created by Ann Kod on 15/10/2025.
//

#include "algorytm_2.h"


void algorytm_2::algo_2(std::vector<atom>& atoms, parameters& params, system_params& system) {
    double tmpVS = 0.0;
    double tmpVP = 0.0;
    double tmpFs = 0.0;

    const double f = params.get_f();
    const double l = params.get_l();
    const double ee = params.get_ee();
    const double r = params.get_r();
    const double r_inv = 1.0 / r;

    const size_t n = atoms.size();

    for (size_t i = 0; i < n; i++) {
        //for that calculates Vs ad Vp and Fsi Fpi
        double tmpFsi = 0.0;
        double tmpFpi = 0.0;

        double tmpR = atoms[i].get_R();

        if (tmpR >= l) {
            const double delta = tmpR - l;

            tmpVS += 0.5 * f * delta * delta;
            //force from walls ov vessel
            tmpFsi = f * (l - tmpR);
            tmpFs += fabs(tmpFsi);
        }

        const double xi = atoms[i].get_x();
        const double yi = atoms[i].get_y();
        const double zi = atoms[i].get_z();

        for (size_t j = 0; j < i; j++) {
            //Vp section
            double dx = xi - atoms[j].get_x();
            double dy = yi - atoms[j].get_y();
            double dz = zi - atoms[j].get_z();
            const double r_sq = dx*dx + dy*dy + dz*dz;



            if (r_sq > 1e-20) {  // Unikaj dzielenia przez 0
                double tmpRjk = sqrt(r_sq);
                const double ratio = r / tmpRjk;

                const double ratio2 = ratio * ratio;        // r²
                const double ratio6 = ratio2 * ratio2 * ratio2;  // r⁶
                const double ratio12 = ratio6 * ratio6; //r^12

                // Potencjał Lennarda-Jonesa
                tmpVP += ee * (ratio12 - 2.0 * ratio6);

                double force_magnitude = 12.0 * ee / tmpRjk *
                    (ratio12 - ratio6);
                tmpFpi += force_magnitude;
            }
        }

        atoms[i].set_F(tmpFsi, tmpFpi);
    }
    system.setV(tmpVS+tmpVP);
    system.setP(1./(4*M_PI*l*l) * tmpFs);
}
