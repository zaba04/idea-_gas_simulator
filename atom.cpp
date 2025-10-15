//
// Created by Ann Kod on 08/10/2025.
//

#include "atom.h"

#include <cmath>

void atom::set_p(double momx, double momy, double momz) {
    p[0] = momx;
    p[1] = momy;
    p[2] = momz;
}

double atom::get_R() {
    return sqrt(x*x + y*y + z*z);
}
