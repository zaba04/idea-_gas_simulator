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
    double F[2]; //F[0] = Fwall, F[1] = Ffromotheratoms
    public:
    atom(){};
    atom(double x, double y, double z) : x(x), y(y), z(z) {};
    ~atom(){};
    void set_p(double, double, double);
    void set_F(double s, double p){F[0] = s; F[1] = p;};
    double get_px(){return p[0];}
    double get_py(){return p[1];}
    double get_pz(){return p[2];}
    double get_x(){return x;};
    double get_y(){return y;}
    double get_z(){return z;}
    double get_R();

};



#endif //ATOM_H
