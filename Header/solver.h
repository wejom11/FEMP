#ifndef SOLVER_H
#define SOLVER_H
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <mkl.h>

class SparseMatrix{
public:
    int type;
    int rows;
    double* val;
    int* col;
    int* row_st;

    SparseMatrix(int mtype = 2){
        type = mtype;
        val = nullptr;
        col = nullptr;
        row_st = nullptr;
    }

    SparseMatrix(SparseMatrix& A);

    void del(){
        delete[] val; val = nullptr;
        delete[] col; col = nullptr;
        delete[] row_st; row_st = nullptr;
    }

    /// @brief show this SparseMatrix
    /// @param row total rows of this matrix
    /// @param is_w wheather write matrix in files ./mat.txt or ./${file_name}.txt
    void show(bool is_w = false, std::string file_name = "mat");

    /// @brief compute y = (this).x
    void product(const double* x, double* y);

    /// @brief A = A - B 
    void sub(SparseMatrix& B);

    /// @brief A = alpha * A 
    void mul(double alpha);

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

    void set_zero();

    void dealloc();

    void alloc();

    void show();

    void product(const CBLAS_TRANSPOSE tra, const DenseMatrix& B, const CBLAS_TRANSPOSE trb, DenseMatrix& C, double alpha = 1.0);

    void add(const DenseMatrix& B);

    void mul(double alpha);

    ~DenseMatrix(){
        if(val){
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
double getlen(double* ptA, double* ptB, int dim);

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

/// @brief get matrix M's determinate
/// @param M input Matrix
/// @return determinate
double det(const DenseMatrix &M);

/// @brief find (i,j) Matrix element's position in Sparse Matrix storage format 'val' array
/// @param row the i-th row dimension
/// @param col the j-th col dimension
/// @attention the row/col number and SparseMatrix is one-based indexing;
/// @return the position where this element stored in SparseMatrix.val Array.
int match(const int row, const int col,  SparseMatrix &SPM);

int decode_sont(const std::vector<std::pair<int, int>> &sont, std::vector<int>& noet,
                 std::vector<int>& dofet);

#endif