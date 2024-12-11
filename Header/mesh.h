#ifndef MESH_H
#define MESH_H

#include <stdio.h>
#include "elements.h"
#include "boundaries.h"
#include "solution.h"

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

class Trussmesh_demo{
private:
    int line_ele_num;       // 宽度方向
    int dis_order;
    double length;
    double dis_mag;        
    double sec_a;
    int type;
    int bnd_type;           // 22 - anchorage in two side; 11 - Simply supported at both ends
                            // 02 - one side anchorage and one side free
                            // 12 - one side anchorage and one side simply supported

public:
    Trussmesh_demo(int bndtype = 22, int num = 2, double len = 1, double pval = 1.0, int tp = 0, double sa = 1.0, int dis_od = 1){
        bnd_type = bndtype;
        line_ele_num = num;
        length = len;
        dis_mag = pval;
        type = tp;
        sec_a = sa;
        dis_order = dis_od;
    }

    void generate(elements &eles, std::vector<material> &matlib, std::vector<plate_prop> &seclib, 
                  boundaries &dbnd, solution& sln, double* &xyz, short ad);

};

double disp_disturb(short bnd_type, short order, double y, double l);

#endif