#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include <stdio.h>
#include <vector>
#include "solver.h"
#include "elements.h"

struct disp_bnd{
    std::vector<int> Anchorage_nid = {};
    std::vector<std::vector<bool>> fixed_dof = {};
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

struct dload_2d{
    std::vector<double> norm_vec = {1., 0};
    double val = 0.;
    int type = 0;
    std::vector<std::pair<int, int>> elset_line;
};

struct cload_2d{
    std::vector<double> val;
    std::vector<int> nset;
};

class boundaries{
public:
    disp_bnd dbd;
    plate_bnd pbd;
    std::vector<dload_2d> dbd2ds;
    std::vector<cload_2d> cbd2ds;

    boundaries(disp_bnd& db, plate_bnd& pb, std::vector<dload_2d>& dd, std::vector<cload_2d>& cd){
        dbd = db;
        pb = pb;
        dbd2ds = dd;
        cbd2ds = cd;
    }

    boundaries(){};

    void add_bndD(DenseMatrix& Kd, double* Fo, elements& eles);

    void add_bndS(SparseMatrix& Ks, double* Fo, elements& eles, bool is_NL = false);

    ~boundaries(){};
};

#endif