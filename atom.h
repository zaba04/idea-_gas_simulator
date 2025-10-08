//
// Created by Ann Kod on 08/10/2025.
//

#ifndef ATOM_H
#define ATOM_H



class atom {
    double x;
    double y;
    double z;
    double p[3];
    public:
    atom(){};
    atom(double x, double y, double z) : x(x), y(y), z(z) {};
    ~atom(){};
    void set_p(double, double, double);
    double get_px(){return p[0];}
    double get_py(){return p[1];}
    double get_pz(){return p[2];}
};



#endif //ATOM_H
