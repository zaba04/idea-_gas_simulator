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
    const size_t n = atoms.size();


    for (size_t i = 0; i < n; i++) {
        //for that calculates Vs ad Vp and Fsi Fpi
        double Fsi_x = 0.0;
        double Fsi_y = 0.0;
        double Fsi_z = 0.0;

        double tmpFpi = 0.0;

        const double xi = atoms[i].get_x();
        const double yi = atoms[i].get_y();
        const double zi = atoms[i].get_z();
        const double ri = atoms[i].get_R();

        if (ri >= l) {
            const double delta = ri - l;
            tmpVS += 0.5 * f * delta * delta;

            if (ri > 1e-10) {
                const double force_magnitude = f * (l - ri);

                Fsi_x = force_magnitude * xi / ri;
                Fsi_y = force_magnitude * yi / ri;
                Fsi_z = force_magnitude * zi / ri;

                double F_total = sqrt(Fsi_x*Fsi_x + Fsi_y*Fsi_y + Fsi_z*Fsi_z);
                tmpFs += F_total;
            }
        }

        double Fpi_x = 0.0;
        double Fpi_y = 0.0;
        double Fpi_z = 0.0;

        for (size_t j = 0; j < n; j++) {
            if (i == j) continue;

            //Vp section
            double dx = xi - atoms[j].get_x();
            double dy = yi - atoms[j].get_y();
            double dz = zi - atoms[j].get_z();
            const double r_sq = dx*dx + dy*dy + dz*dz;



            if (r_sq > 1e-20) {
                double rij = sqrt(r_sq);
                const double ratio = r / rij;

                const double ratio2 = ratio * ratio;        // r²
                const double ratio6 = ratio2 * ratio2 * ratio2;  // r⁶
                const double ratio12 = ratio6 * ratio6; //r^12

                if (j < i) {
                    tmpVP += ee * (ratio12 - 2.0 * ratio6);
                }

                const double force_factor = 12.0 * ee * (ratio12 - ratio6) / r_sq;

                Fpi_x += force_factor * dx;
                Fpi_y += force_factor * dy;
                Fpi_z += force_factor * dz;
            }
        }

        atoms[i].set_Fatoms(Fpi_x , Fpi_y, Fpi_z);
        atoms[i].set_Fwall(Fsi_x , Fsi_y, Fsi_z);
        std::cout << atoms[i].get_FatomX() << " " << atoms[i].get_FatomY() << " " << atoms[i].get_FatomZ() <<  std::endl;
    }
    system.setV(tmpVS+tmpVP);
    system.setP(tmpFs/(4.0*M_PI*l*l));
}
