#ifndef MESH_H
#define MESH_H

#include <stdio.h>
#include "elements.h"
#include "boundaries.h"

class Linemesh_demo{
private:
    int ele_num;
    double length;
    bool is_M;
    double f_val;
    int type;
public:
    Linemesh_demo(int num = 1, double len = 1, double val = 1.0, bool M = true, int t = 0){
        ele_num = num;
        length = len;
        is_M = M;
        f_val = val;
        type = t;
    };

    void generate(elements &eles, std::vector<material> &matlib, std::vector<section_TSbeam> &seclib, 
                  boundaries &dbnd, double* &xyz);

    ~Linemesh_demo(){};
};

class Square_Platemesh_demo{
private:
    int line_ele_num;
    double length;
    double p_val;
    int type;

public:
    Square_Platemesh_demo(int num = 2, double len = 1, double pval = 1.0, int tp = 0){
        line_ele_num = num;
        length = len;
        p_val = pval;
        type = tp;
    }

    void generate(elements &eles, std::vector<material> &matlib, std::vector<plate_prop> &seclib, 
                  boundaries &dbnd, double* &xyz);

};

#endif