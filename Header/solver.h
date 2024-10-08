#ifndef SOLVER_H
#define SOLVER_H
#include <stdio.h>
#include <math.h>
#include <vector>

class SparseMatrix{
public:
    int type;
    double* val;
    int* col;
    int* row_st;

    SparseMatrix(int mtype = 2){
        type = mtype;
        val = nullptr;
        col = nullptr;
        row_st = nullptr;
    }

    void del(){
        delete[] val; val = nullptr;
        delete[] col; col = nullptr;
        delete[] row_st; row_st = nullptr;
    }

    /// @brief show this SparseMatrix
    /// @param row total rows of this matrix
    void show(int row);

    ~SparseMatrix(){
        if(!(val == nullptr && col == nullptr && row_st == nullptr)){
            this -> del();
        }
    }
};

class DenseMatrix{
public:
    int type;
    int rows;
    int cols;
    double* val;

    /// @brief initialize a dense matrix
    /// @param t type = 1: unsymmetric; type = 2: symmetric;
    /// @param row The number of rows of the matrix
    /// @param col The number of columns of the matrix
    /// @param ini_val The initial value of the matrix element(default 0)
    DenseMatrix(int row = 0, int col = 0, int t = 1, double ini_val = 0.){
        val = nullptr;
        type = t;
        rows = row;
        cols = col;
        if(row != 0 && col != 0){
            val = new double[row * col]{ini_val};
        }
    };

    void resize(int rows, int cols);

    void dealloc();

    void alloc();

    void show();

    ~DenseMatrix(){
        if(val != nullptr){
            delete[] val; val = nullptr;
        }
    };
};

class pardiso_cfg{
public:
    int maxfct;
    int mnum;
    int perm;
    int iparm[64];
    int msglvl;
    int error;

    pardiso_cfg(int &mtype){
        initial(mtype);
    };

    /// @brief intial the parameter of Intel Pardiso solver
    void initial(int &type);

    ~pardiso_cfg(){};
};

/// @brief get the length of line AB
/// @param ptA point A
/// @param ptB point B
/// @return length (|AB|)
double getlen(double* ptA, double* ptB);

/// @brief Determine which category the node belongs to
/// @param n_id Node number
/// @param sont category of node
/// @return category number
int wherend(int n_id, std::vector<std::pair<int, int>> &sont);

/// @brief find element's position in the array, or the interval which it belongs to.
/// (give the Lower bound of interval's position)
/// @param ele element
/// @param array array, must be an increasing series, which means array[i+1] > array[i] for arbitrary i.
/// @return position
// template <class TYPE>
int find(int &ele, std::vector<int> &array);

#endif