#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <mkl.h>
#include "elements.h"
#include "mesh.h"

class asb_manager{
public:
    elements eles;
    double* xyz_coord;
    std::vector<material> mater_lib;
    std::vector<section_TSbeam> section_lib;
    disp_bnd dbd;

    DenseMatrix K;
    double* Fout;

    double* Var;

    asb_manager(){
        xyz_coord = nullptr;
        Fout = nullptr;
        Var = nullptr;
    };

    void init_mesh(int num = 1, bool M = true, int type = 0, double len = 1, double val = 2.5);

    void init_KF();

    void add_bnd();

    void solve();

    void post();

    ~asb_manager(){
        delete[] xyz_coord; xyz_coord = nullptr;
        delete[] Fout; Fout = nullptr;
        delete[] Var; Var = nullptr;
    };
};

#endif