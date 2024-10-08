#include <vector>

#ifndef ELEMENT_PROP_H
#define ELEMENT_PROP_H

class material{
public:
    double E;       // Youngs modulus
    double G;       // sheer modulus
    double nu;      // Poission ratio

    material(double Ev = 3E11, double Gv = 3E11/2.6, double nuv = 0.3){
        E = Ev;
        G = Gv;
        nu = nuv;
    };

    ~material(){};
};

class section_TSbeam{
public:
    double k;       
    double I;
    // double Iy;
    // double Ixy;
    double A;

    section_TSbeam(double Iv = 1E-14 / 12, double Av = 1E-6, double kv = 1.2){
        I = Iv;
        // Iy = Iyv;
        // Ixy = Ixyv;
        A = Av;
        k = kv;
    };

    ~section_TSbeam(){};
};

class beam_bnd{
public:
    double q_const;
    std::vector<std::pair<double, double>> p_j;
    std::vector<std::pair<double, double>> M_k;

    beam_bnd(double q_c = 0.0){
        q_const = q_c;
    };

    ~beam_bnd(){};
};


#endif