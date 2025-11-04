#include "algorytm_2.h"
// Funkcja 1: Obliczanie sił dla pojedynczego atomu
void algorytm_2::calculate_forces_for_atom(std::vector<atom>& atoms, parameters& params, int i) {
    const double f = params.get_f();
    const double l = params.get_l();
    const double ee = params.get_ee();
    const double r = params.get_r();
    const size_t n = atoms.size();

    double Fsi_x = 0.0;
    double Fsi_y = 0.0;
    double Fsi_z = 0.0;

    const double xi = atoms[i].get_x();
    const double yi = atoms[i].get_y();
    const double zi = atoms[i].get_z();
    const double ri = atoms[i].get_R();

    // Siły od ścianek
    if (ri >= l) {
        if (ri > 1e-10) {
            const double force_magnitude = f * (l - ri);
            Fsi_x = force_magnitude * xi / ri;
            Fsi_y = force_magnitude * yi / ri;
            Fsi_z = force_magnitude * zi / ri;
        }
    }

    // Siły van der Waalsa
    double Fpi_x = 0.0;
    double Fpi_y = 0.0;
    double Fpi_z = 0.0;

    for (size_t j = 0; j < n; j++) {
        if (i == j) continue;

        double dx = xi - atoms[j].get_x();
        double dy = yi - atoms[j].get_y();
        double dz = zi - atoms[j].get_z();
        const double r_sq = dx*dx + dy*dy + dz*dz;

        if (r_sq > 1e-20) {
            double rij = sqrt(r_sq);
            const double term1 = pow(r/rij, 12);
            const double term2 = pow(r/rij, 6);
            const double force_factor = 12.0 * ee * (term1 - term2) / r_sq;

            Fpi_x += force_factor * dx;
            Fpi_y += force_factor * dy;
            Fpi_z += force_factor * dz;
        }
    }

    atoms[i].set_Fatoms(Fpi_x, Fpi_y, Fpi_z);
    atoms[i].set_Fwall(Fsi_x, Fsi_y, Fsi_z);
}

// Funkcja 2: Obliczanie całkowitej energii
double algorytm_2::calculate_total_energy(std::vector<atom>& atoms, parameters& params) {
    const double f = params.get_f();
    const double l = params.get_l();
    const double ee = params.get_ee();
    const double r = params.get_r();
    const size_t n = atoms.size();

    double total_VS = 0.0;
    double total_VP = 0.0;

    // Energia od ścianek
    for (size_t i = 0; i < n; i++) {
        const double ri = atoms[i].get_R();
        if (ri >= l) {
            const double delta = ri - l;
            total_VS += 0.5 * f * delta * delta;
        }
    }

    // Energia van der Waalsa (każda para raz)
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i+1; j < n; j++) {
            double dx = atoms[i].get_x() - atoms[j].get_x();
            double dy = atoms[i].get_y() - atoms[j].get_y();
            double dz = atoms[i].get_z() - atoms[j].get_z();
            const double r_sq = dx*dx + dy*dy + dz*dz;

            if (r_sq > 1e-20) {
                double rij = sqrt(r_sq);
                const double ratio = r / rij;
                const double ratio6 = pow(ratio, 6);
                const double ratio12 = ratio6 * ratio6;
                total_VP += ee * (ratio12 - 2.0 * ratio6);
            }
        }
    }

    return total_VS + total_VP;
}

// Funkcja 3: Obliczanie ciśnienia (equation 15)
double algorytm_2::calculate_pressure(std::vector<atom>& atoms, parameters& params) {
    const double l = params.get_l();
    double total_wall_force = 0.0;

    for (size_t i = 0; i < atoms.size(); i++) {
        double Fx = atoms[i].get_FwallX();
        double Fy = atoms[i].get_FwallY();
        double Fz = atoms[i].get_FwallZ();
        total_wall_force += sqrt(Fx*Fx + Fy*Fy + Fz*Fz);
    }

    return total_wall_force / (4.0 * M_PI * l * l);
}