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
    double Fwall[3];
    double Fatoms[3];
    public:
    atom(){};
    atom(double x, double y, double z) : x(x), y(y), z(z) {};
    ~atom(){};
    void set_p(double, double, double);
    void set_Fwall(double x, double y , double z){Fwall[0] = x;Fwall[1] = y;Fwall[2] = z;}
    void set_Fatoms(double x, double y , double z){Fatoms[0] = x;Fatoms[1] = y;Fatoms[2] = z;}
    double get_px(){return p[0];}
    double get_py(){return p[1];}
    double get_pz(){return p[2];}
    double get_x(){return x;}
    double get_y(){return y;}
    double get_z(){return z;}
    double get_FwallX(){return Fwall[0];}
    double get_FwallY(){return Fwall[1];}
    double get_FwallZ(){return Fwall[2];}
    double get_FatomX(){return Fatoms[0];}
    double get_FatomY(){return Fatoms[1];}
    double get_FatomZ(){return Fatoms[2];}
    double get_R();

};



#endif //ATOM_H
