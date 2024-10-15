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

    stiffness(){
        Fout = nullptr;
    };

    void init_KFD(elements& eles);

    void init_KFS(elements& eles);

    void init_K_symbolic(elements& eles);

    ~stiffness(){
        K_sparse.del();
        K_dense.dealloc();
        delete[] Fout; Fout = nullptr;
    };
};


#endif