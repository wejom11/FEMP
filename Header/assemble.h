#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <mkl.h>
#include <string.h>
#include "elements.h"
#include "mesh.h"
#include "stiffness.h"
#include "boundaries.h"
#include "solution.h"

class asb_manager{
public:
    short analysis_dim;
    elements eles;
    double* xyz_coord;
    std::vector<material> mater_lib;
    std::vector<section_TSbeam> section_lib;
    std::vector<plate_prop> plate_prop_lib;
    boundaries bnds;

    std::string method;
    std::string task;
    stiffness KF;

    solution sln;

    asb_manager(std::string met = "S", std::string ts = "MLSM", short ad = 3){
        analysis_dim = ad;
        xyz_coord = nullptr;
        method = met;
        this->task = ts;
        // if(!_strcmpi(method.data(), "SN")){
        //     task = "BAEA";        // buckling analysis with eigenvalue analysis
        // }
    };

    /// @brief beam demo mesh
    /// @param num element number
    /// @param M bending moment or sheering force
    /// @param type element type
    /// @param len beam length
    /// @param val magnitude of M or P
    void init_mesh(int num = 1, bool M = true, int type = 0, double len = 1, double val = 2.5);
    void init_mesh(int num = 1, int type = 0, double len = 1, double val = 2.0);
    void init_mesh(int w_num = 2, int type = 1, int b_type = 22, double pcr = 1., double seca = 2.0,
                   double len = 20, int order = 1);

    /// @param parms parameters for solver
    void solve(int** parms = nullptr);

    void post();

    ~asb_manager(){
        delete[] xyz_coord; xyz_coord = nullptr;
    };
};

#endif