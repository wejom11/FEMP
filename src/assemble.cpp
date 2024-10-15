#include "assemble.h"

void asb_manager::init_mesh(int num, bool M, int type, double len, double val){
    Linemesh_demo lmd(num, len, val, M, type);
    lmd.generate(eles, mater_lib, section_lib, bnds, xyz_coord);
}

void asb_manager::init_mesh(int num, int type, double len, double val){
    Square_Platemesh_demo spmd(num, len, val, type);
    spmd.generate(eles, mater_lib, plate_prop_lib, bnds, xyz_coord);
};

void asb_manager::init_KF(){
    if(!_strcmpi(method.data(), "D")){
        KF.init_KFD(eles);
    }
    else if(!_strcmpi(method.data(), "S")){
        KF.init_KFS(eles);
    }
};

void asb_manager::add_bnd(){
    if(!_strcmpi(method.data(), "D")){
        bnds.add_bndD(KF.K_dense, KF.Fout, eles);
    }
    else if(!_strcmpi(method.data(), "S")){
        bnds.add_bndS(KF.K_sparse, KF.Fout, eles);
    }
}

void asb_manager::solve(){
    if(!_strcmpi(method.data(), "D")){
        this->solveD();
    }
    else if(!_strcmpi(method.data(), "S")){
        this->solveS();
    }
};

void asb_manager::solveD(){
    int* ipiv = new int[KF.K_dense.rows]{0};
    int info;

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, KF.K_dense.rows, 1, KF.K_dense.val, KF.K_dense.rows, ipiv, KF.Fout, 1);

    if(info != 0){
        printf("ERROR: solving process error!\n");
        exit(1);
    }

    delete[] ipiv; ipiv = nullptr;
}

void asb_manager::solveS(){
    const int dof = KF.K_sparse.rows;
    int mtype = KF.K_sparse.type;
    int nrhs = 1;
    Var = new double[dof] {0};
    void* ptsolver[64];
    pardiso_cfg para = pardiso_cfg(mtype);
    double ddum;

    for(int i = 0; i < 64; i++){
        ptsolver[i] = 0; 
    }

    int phase = 11;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(KF.K_sparse.type), &phase, &dof, KF.K_sparse.val, 
            KF.K_sparse.row_st, KF.K_sparse.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), &ddum, &ddum, &(para.error));
    if(para.error != 0 ){
        printf ("ERROR during symbolic factorization: %i \n", para.error);
    }
    
    phase = 22;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(KF.K_sparse.type), &phase, &dof, KF.K_sparse.val, 
            KF.K_sparse.row_st, KF.K_sparse.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), &ddum, &ddum, &(para.error));
    if(para.error != 0 ){
        printf ("ERROR during numerical factorization: %i \n", para.error);
    }

    phase = 33;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(KF.K_sparse.type), &phase, &dof, KF.K_sparse.val, 
            KF.K_sparse.row_st, KF.K_sparse.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), KF.Fout, Var, &(para.error));
    if(para.error != 0 ){
        printf ("ERROR during solution: %i \n", para.error);
    }
    printf("solve completed ...\n");

    phase = -11;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(KF.K_sparse.type), &phase, &dof, KF.K_sparse.val, 
            KF.K_sparse.row_st, KF.K_sparse.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), KF.Fout, &ddum, &(para.error));
};
