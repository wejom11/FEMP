#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <mkl.h>
#include "elements.h"
#include "mesh.h"
#include "stiffness.h"
#include "boundaries.h"

class asb_manager{
public:
    elements eles;
    double* xyz_coord;
    std::vector<material> mater_lib;
    std::vector<section_TSbeam> section_lib;
    std::vector<plate_prop> plate_prop_lib;
    boundaries bnds;

    std::string method;
    stiffness KF;

    double* Var;

    asb_manager(){
        xyz_coord = nullptr;
        Var = nullptr;
        method = 'S';
    };

    /// @brief beam demo mesh
    /// @param num element number
    /// @param M bending moment or sheering force
    /// @param type element type
    /// @param len beam length
    /// @param val magnitude of M or P
    void init_mesh(int num = 1, bool M = true, int type = 0, double len = 1, double val = 2.5);
    void init_mesh(int num = 1, int type = 0, double len = 1, double val = 2.0);

    void init_KF();

    void add_bnd();

    void solve();

    void solveD();

    void solveS();

    void post();

    ~asb_manager(){
        delete[] xyz_coord; xyz_coord = nullptr;
        delete[] Var; Var = nullptr;
    };
};

#endif