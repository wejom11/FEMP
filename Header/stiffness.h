#ifndef STIFFNESS_H
#define STIFFNESS_H

#include <stdio.h>
#include <iostream>
#include <string>
#include "solver.h"
#include "elements.h"

class stiffness{
public:
    SparseMatrix K_sparse;
    DenseMatrix K_dense;
    double* Fout;
    double* Fint;

    stiffness(){
        Fout = nullptr;
        Fint = nullptr;
    };

    void init_KFD(elements& eles);

    void init_KFS(elements& eles);

    void init_KFSN(elements& eles, std::vector<short> option = {1111});

    void init_K_symbolic(elements& eles);

    void asb_spa_sym(double* K, double* F, std::vector<int>& node_id,
                     const std::vector<std::pair<int,int>>& sont);

    ~stiffness(){
        K_sparse.del();
        K_dense.dealloc();
        delete[] Fout; Fout = nullptr;
        delete[] Fint; Fint = nullptr;
    };
};

#endif