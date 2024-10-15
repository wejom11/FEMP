#include <vector>

#ifndef ELEMENT_PROP_H
#define ELEMENT_PROP_H

class material{
public:
    double E;       // Youngs modulus
    double G;       // sheer modulus
    double nu;      // Poission ratio

    material(double Ev = 2.1E7, double nuv = 0.30){
        E = Ev;
        nu = nuv;
        G = E / (1 + nu) / 2;
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

class plate_prop{
public:
    double thickness;
    double k_eqv;

    plate_prop(double t = 0.01, double k = 1.2){
        thickness = t;
        k_eqv = k;
    };

    ~plate_prop(){};
};

#endif