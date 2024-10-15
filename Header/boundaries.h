#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include <stdio.h>
#include <vector>
#include "solver.h"
#include "elements.h"

struct disp_bnd{
    std::vector<int> Anchorage_nid = {};
};

struct plate_dload{
    std::vector<int> ele_set;
    double val = 1.0;
};

struct plate_cload{
    std::vector<int> node_set;
    double Mn_val = 0.;
    double Ms_val = 0.;
    double Qn_val = 0.;
};


struct plate_bnd{
    std::vector<plate_dload> pbd_d;
    std::vector<plate_cload> pbd_c;
};

class boundaries{
public:
    disp_bnd dbd;
    plate_bnd pbd;

    boundaries(disp_bnd& db, plate_bnd& pb){
        dbd = db;
        pb = pb;
    }

    boundaries(){};

    void add_bndD(DenseMatrix& Kd, double* Fo, elements& eles);

    void add_bndS(SparseMatrix& Ks, double* Fo, elements& eles);

    ~boundaries(){};
};

#endif