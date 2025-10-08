#include <iostream>
#include <cmath>
#include "atom.h"
#include "parameters.h"



int main() {

    parameters parameters;

    std::vector<atom> atoms;

    for (int i = 0; i < parameters.get_n(); ++i) {
        for (int j = 0; j < parameters.get_n(); ++j) {
            for (int k = 0; k < parameters.get_n(); ++k) {
                float tmpx = 0;
                float tmpy = 0;
                float tmpz = 0;

                double b0[3] = {parameters.get_a(), 0, 0};
                double b1[3] = {parameters.get_a()/2., (parameters.get_a()*sqrt(3))/2. , 0};
                double b2[3] = {parameters.get_a()/2., (parameters.get_a()*sqrt(3))/6. , parameters.get_a()*sqrt(2/3.)};

                float tmpp = (parameters.get_n()-1);

                tmpx = (i - tmpp/2) * b0[0] + (j - tmpp/2) * b1[0] + (k - tmpp/2) * b2[0];
                tmpy = (i - tmpp/2) * b0[1] + (j - tmpp/2) * b1[1] + (k - tmpp/2) * b2[1];
                tmpz = (i - tmpp/2) * b0[2] + (j - tmpp/2) * b1[2] + (k - tmpp/2) * b2[2];

                atoms.push_back(atom(tmpx, tmpy, tmpz));
                std::cout << tmpx << "  " << tmpy << "  " << tmpz << std::endl;
            }
        }


    }
    return 0;
}