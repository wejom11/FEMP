#ifndef MESH_H
#define MESH_H

#include <stdio.h>
#include "elements.h"

class disp_bnd{
public:
    std::vector<int> Anchorage_nid;
};

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
                  disp_bnd &dbnd, double* xyz);

    ~Linemesh_demo(){};
};


#endif