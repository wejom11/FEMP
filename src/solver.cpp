#include <math.h>
#include "solver.h"

void SparseMatrix::show(int row){
    int i, j, k, tiscol, lstcol;

    for(i = 0; i < row; i++){
        printf("{");
        lstcol = 1;
        for(j = this->row_st[i]-1; j < this->row_st[i+1]-1; j++){
            tiscol = this->col[j];
            for(k = lstcol; k < tiscol; k++){
                //printf("%.2f, ", 0.);
                printf("... ");
            }
            if(tiscol == row){
                printf("%.8f},", this->val[j]);
            }
            else{
                printf("%.8f, ", this->val[j]);
            }
            lstcol = tiscol + 1;
        }

        for(k = lstcol - 1; k < row; k++){
            if(k == row - 1){
                //printf("%.2f}, ",0.);
                printf("...}");
            }
            else{
                //printf("%.2f, ",0.);
                printf("... ");
            }
        }
    }
    printf("\n");
};

void pardiso_cfg::initial(int &type){
    int i = 0;
    maxfct = 1;
    mnum = 1;
    msglvl = 0;
    for(i = 0; i < 64; i++){
        iparm[i] = 0;
    }
    if(type == 2){
        iparm[0] = 1;
        iparm[1] = 2;
        iparm[9] = 13;
    }
    else if(type == 11 || type == 1){
        iparm[0] = 1;
        iparm[1] = 2;
        iparm[9] = 13;
    }
    else{
        printf("matrix type %i haven't supported yet!", type);
    }
}

void DenseMatrix::dealloc(){
    delete[] val;
    val = nullptr;
}

void DenseMatrix::alloc(){
    dealloc();
    if(rows != 0 && cols != 0){
        val = new double[rows * cols]{0};
    }
}

void DenseMatrix::resize(int row, int col){
    rows = row;
    cols = col;
    alloc();
}

void DenseMatrix::show(){
    int ROW;

    for(int i = 0; i < rows; i++){
        ROW = i * cols;
        printf("row %i: ", i+1);

        for(int j = 0; j < cols; j++){
            printf("%.4f ", val[ROW + j]);
        }
        printf("\n");
    }
}

double getlen(double* ptA, double* ptB){
    double len2 = 0.;
    for(int i = 0; i < 3; i++){
        len2 += pow(ptA[i] - ptB[i], 2);
    }
    double len = sqrt(len2);

    return len;
};

int wherend(int n_id, std::vector<std::pair<int, int>> &sont){
    int nt_size = sont.size();

    std::vector<int> noet(nt_size + 1, 1);
    for(int i = 0; i < nt_size; i++){
        noet.at(i + 1) = noet.at(i) + sont.at(i).first;
    }

    return find(n_id, noet);
};

// template <class TYPE>
int find(int &ele, std::vector<int> &array){
    bool done = false;
    int front_ptr = 0;
    int back_ptr = array.size() - 1;
    int ptr = (front_ptr + back_ptr) / 2;

    if(array.back() == ele){
        return back_ptr;
    }
    else if(array.front() > ele || array.back() < ele){
        return -1;
    }
    
    while(!done){
        if(array.at(ptr) == ele){
            return ptr;
        }
        else if(array.at(ptr) > ele){
            back_ptr = ptr;
            ptr = (front_ptr + back_ptr) / 2;
        }
        else{
            front_ptr = ptr;
            ptr = (front_ptr + back_ptr) / 2;
            if(ptr == front_ptr){
                return front_ptr;
            }
        }
    }

    return -1;
};

double det(const DenseMatrix &M){
    double* mat  = M.val;
    if(M.rows == 2){
        return mat[0]*mat[3] - mat[1]*mat[2];
    }
    else if(M.rows == 3){
        double d = mat[0]*mat[4]*mat[8] + mat[3]*mat[7]*mat[2] + mat[6]*mat[1]*mat[5] -
                   mat[2]*mat[4]*mat[6] - mat[0]*mat[5]*mat[7] - mat[8]*mat[1]*mat[3];
        return d;
    }
    return 0;
};

int match(const int row, const int col,  SparseMatrix &SPM){
    int start = SPM.row_st[row-1] - 1;
    int end = SPM.row_st[row] - 1;
    int mid = (start + end) / 2;
    int ptr = mid;
    bool is_done = false;

    while(!is_done){
        if(SPM.col[ptr] < col){
            start = mid;
            mid = (start + end) / 2;
            if(mid == ptr){
                is_done = true;
            }else{
                ptr = mid;
            }
        }
        else if(SPM.col[ptr] > col){
            end = mid;
            mid = (start + end) / 2;
            if(mid == ptr){
                is_done = true;
            }
            else{
                ptr = mid;
            }
        }
        else{
            return ptr;
        }
    }
    
    return -1;
};