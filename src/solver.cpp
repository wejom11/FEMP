#include <math.h>
#include "solver.h"

SparseMatrix::SparseMatrix(SparseMatrix& A){
    this->rows = A.rows;
    this->type = A.type;
    this->val = new double[A.row_st[A.rows] - 1]{0.};
    this->col = new int[A.row_st[A.rows] - 1]{0};
    this->row_st = new int[A.rows + 1]{0};

    int i;
    for(i = 0; i < A.row_st[A.rows] - 1; i++){
        this->val[i] = A.val[i];
        this->col[i] = A.col[i];
    }
    for(i = 0; i < A.rows + 1; i++){
        this->row_st[i] = A.row_st[i];
    }
};

void SparseMatrix::show(bool is_w, std::string file_name){
    int i, j, k, l, tiscol, lstcol;

    std::string fn = file_name;
    std::ofstream out_f_io;
    if (is_w) {
        fn.append(".txt");
        out_f_io.open(fn, std::ios_base::out);
        out_f_io.setf(out_f_io.scientific);
        out_f_io.precision(10);

        out_f_io << this->rows << std::endl;
        l = 1;
        for (i = 0; i < rows; i++) {
            lstcol = 1;
            for (j = this->row_st[i] - 1; j < this->row_st[i + 1] - 1; j++) {
                tiscol = this->col[j];
                for (k = lstcol; k < tiscol; k++) {
                    out_f_io << 0. << "   ";
                    if (l == 10) {
                        out_f_io << '\n';
                        l = 0;
                    }
                    l++;
                }
                if (tiscol == rows) {
                    out_f_io << this->val[j] << "\n";
                    l = 10;
                    if (l == 10) {
                        out_f_io << '\n';
                        l = 0;
                    }
                    l++;
                }
                else {
                    out_f_io << this->val[j] << "   ";
                    if (l == 10) {
                        out_f_io << '\n';
                        l = 0;
                    }
                    l++;
                }
                lstcol = tiscol + 1;
            }

            for (k = lstcol - 1; k < rows; k++) {
                if (k == rows - 1) {
                    out_f_io << 0. << "\n";
                    l = 10;
                    if (l == 10) {
                        out_f_io << '\n';
                        l = 0;
                    }
                    l++;
                }
                else {
                    out_f_io << 0. << "   ";
                    if (l == 10) {
                        out_f_io << '\n';
                        l = 0;
                    }
                    l++;
                }
            }
        }
    }
    else {
        for (i = 0; i < rows; i++) {
            printf("{");
            lstcol = 1;
            for (j = this->row_st[i] - 1; j < this->row_st[i + 1] - 1; j++) {
                tiscol = this->col[j];
                for (k = lstcol; k < tiscol; k++) {
                    //printf("%.2f, ", 0.);
                    printf("... ");
                }
                if (tiscol == rows) {
                    printf("%.8f},", this->val[j]);
                }
                else {
                    printf("%.8f, ", this->val[j]);
                }
                lstcol = tiscol + 1;
            }

            for (k = lstcol - 1; k < rows; k++) {
                if (k == rows - 1) {
                    //printf("%.2f}, ",0.);
                    printf("...}");
                }
                else {
                    //printf("%.2f, ",0.);
                    printf("... ");
                }
            }
        }
        printf("\n");
    }
};

void SparseMatrix::product(const double* x, double* y){
    if(type == 2){
        sparse_matrix_t A_handle;
        sparse_status_t info;
        matrix_descr mdescr;
        mdescr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
        mdescr.mode = SPARSE_FILL_MODE_UPPER;
        mdescr.diag = SPARSE_DIAG_NON_UNIT;

        info = mkl_sparse_d_create_csr(&A_handle, SPARSE_INDEX_BASE_ONE, rows, rows,
            row_st, row_st + 1, col, val);

        info = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, A_handle, mdescr, x, 0, y);
    }
    else{
        printf("\033[1;31mERROR\033[0m: wrong matrix type\n");
    }
}

void SparseMatrix::sub(SparseMatrix& B){
    for(int i = 0; i < this->row_st[this->rows] - 1; i++){
        this->val[i] -= B.val[i];
    }
}

void SparseMatrix::mul(double alpha){
    for(int i = 0; i < this->row_st[this->rows] - 1; i++){
        this->val[i] *= alpha;
    }
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
    else if(type == -2){
        iparm[0] = 1;
        iparm[9] = 13;
    }
    else{
        printf("matrix type %i haven't supported yet!", type);
    }
}

void DenseMatrix::dealloc(){
    delete[] val;
    val = nullptr;
    rows = 0;
    cols = 0;
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

void DenseMatrix::set_zero(){
    for(int i = 0; i < rows * cols; i++){
        val[i] = 0.;
    }
};

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

void DenseMatrix::product(const CBLAS_TRANSPOSE tra, const DenseMatrix& B, const CBLAS_TRANSPOSE trb, DenseMatrix& C, double alpha){
    int Acs = tra == CblasNoTrans ? this->cols : this -> rows;
    int Brs = trb == CblasNoTrans ? B.rows : B.cols;
    if(Acs != Brs){
        printf("\033[1;31mERROR\033[0m: matrices' dimension do not match!\n");
        exit(0);
    }

    cblas_dgemm(CblasRowMajor, tra, trb, C.rows, C.cols, Acs, alpha, this->val, 
                this->cols, B.val, B.cols, 0, C.val, C.cols);
};

void DenseMatrix::add(const DenseMatrix& B){
    int i;
    for(i = 0; i < rows * cols; i++){
        val[i] += B.val[i];
    }
};

void DenseMatrix::mul(double alpha){
    int i;
    for(i = 0; i < rows * cols; i++){
        val[i] *= alpha;
    }
};

double getlen(double* ptA, double* ptB, int dim){
    double len2 = 0.;
    for(int i = 0; i < dim; i++){
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

int decode_sont(const std::vector<std::pair<int, int>> &sont, std::vector<int>& noet,
                 std::vector<int>& dofoet){
    int ntype_num = sont.size();
    noet.resize(ntype_num + 1, 1);
    dofoet.resize(ntype_num + 1, 0);

    for(int i = 0; i < ntype_num; i++){
        noet.at(i + 1) = noet.at(i) + sont.at(i).first;
        dofoet.at(i + 1) = dofoet.at(i) + sont.at(i).second * sont.at(i).first;
    }

    return ntype_num;
};