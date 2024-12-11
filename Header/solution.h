#ifndef SOLUTION_H
#define SOLUTION_H

#include <stdio.h>
#include <iomanip>
#include "stiffness.h"
#include "boundaries.h"
#include "read.h"

class solution{
public:
    double* Var;

    solution(){
        Var = nullptr;
    }

    void init_sln(int dof);

    void solve_D(stiffness& KF, boundaries& bnd, elements& eles);

    void solve_S(stiffness& KF, boundaries& bnd, elements& eles);

    void solve_SN_MLSM(stiffness& KF, boundaries& bnd, elements& eles, double* xyz, int** parms = nullptr);

    void solve_SN_GALM(stiffness& KF, boundaries& bnd, elements& eles, double* xyz, int** parms = nullptr);

    bool solveeqnS(const SparseMatrix& A, double* b, double* x);

    // double eig_IPM(SparseMatrix& A, SparseMatrix& B, double* Phi = nullptr);

    void del();

    ~solution(){
        delete[] Var; Var = nullptr;
    }
};

void write_out(elements& eles, double* xyz, double* Vars, double load = -1., std::string path = "./",
                std::string name = "Disp_ans");

void write_point(const double* x, const double* y, bool refresh = false, int len = 1, std::string path = "./",
                    std::string name = "Load_Disp_curve");

#endif