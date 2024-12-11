#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <vector>
#include <math.h>
#include <stdio.h>
#include "solver.h"

class Gauss_int{
private:
    int order;
    int dim;
    std::vector<int> selected_dim;

public:
    int int_pt_num;
    double* int_pt_coord;
    double* weightness;

    Gauss_int(int od = 1, int d = 1, std::vector<int> s_dim = {1}){
        order = od;
        dim = d;

        int sds = s_dim.size();
        if(d < s_dim.size()){
            printf("ERROR: dimension %i should greater than selectde dimension's array s_dim's size %i\n", d, sds);
        }
        selected_dim = s_dim;
        int_pt_coord = nullptr;
        weightness = nullptr;
    }

    /// @brief get order's number of gauss point and weightness with one-dimension
    /// @param GP_1 gauss point
    /// @param weight_1 weightness
    void getGP_1(double* GP_1, double* weight_1);

    /// @brief get order's number of gauss point and weightness with ${dim}-dimension
    void getGP();

    void del();

    ~Gauss_int(){
        this->del();
    };

};

class Lagrange_intp{
private:
    int order;
    int dim;

public:
    Gauss_int gintp;
    double* Ni;
    double* Jacobien;
    double* Ndx;

    Lagrange_intp(int od = 1, int d = 1){
        order = od;
        dim = d;
        Ni = nullptr;
        Jacobien = nullptr;
        Ndx = nullptr;
    };

    void set_integration(int int_od, std::vector<int> sd = {});

    void getN__();

    void getNdl__(double* Ndl);

    void getNdx__(std::vector<double*> xyz_c);

    void del();

    ~Lagrange_intp(){
        this->del();
    };

};

class Serendipity_intp{
private:
    int order;
    int dim;

public:
    Gauss_int gintp;
    double* Ni;
    double* Jacobien;
    double* Ndx;

    Serendipity_intp(int od = 2, int d = 2){
        order = od;
        dim = d;
        Ni = nullptr;
        Jacobien = nullptr;
        Ndx = nullptr;
    };

    void set_integration(int int_od, std::vector<int> sd = {});

    void getN__();

    void getNdl__(double* Ndl);

    void getNdx__(std::vector<double*> xyz_c);

    void del();

    ~Serendipity_intp(){
        this->del();
    };

};

#endif