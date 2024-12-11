#include "solution.h"

void solution::init_sln(int dof){
    Var = new double[dof]{0.};
}

void solution::solve_D(stiffness& KF, boundaries& bnd, elements& eles){
    KF.init_KFD(eles);
    bnd.add_bndD(KF.K_dense, KF.Fout, eles);

    int* ipiv = new int[KF.K_dense.rows]{0};
    int info;

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, KF.K_dense.rows, 1, KF.K_dense.val, KF.K_dense.rows, ipiv, KF.Fout, 1);

    if(info != 0){
        printf("\033[1;31mERROR\033[0m: solving process error!\n");
        exit(1);
    }

    delete[] ipiv; ipiv = nullptr;
}

void solution::solve_S(stiffness& KF, boundaries& bnd, elements& eles){
    KF.init_KFS(eles);
    bnd.add_bndS(KF.K_sparse, KF.Fout, eles);

    const int dof = KF.K_sparse.rows;
    Var = new double[dof] {0};

    this->solveeqnS(KF.K_sparse, KF.Fout, Var);
}

void solution::solve_SN_MLSM(stiffness& KF, boundaries& bnd, elements& eles, double* xyz, int** parms){
    int r0, ri, max_itr, max_incr,
        i, j, k, itr, incr;
    double dl_is1, dl_i, dl_min, max_inc_t,
           tol, resid, p/* , resid_his */;
    double *Var_his = nullptr, *Q0 = nullptr,
           *Var_inc = nullptr;
    bool is_end = false, is_singular, is_conver/* , is_write */,
         is_default = true, write_dis;
    std::string path, name;

    if(parms){
        if(*parms){
            if(**parms != 0){
                is_default = false;
            }
        }
    }

    if(!is_default){
        double a, b;
        int exp_digit = pow(10, abs(**parms));
        int sign = **parms > 0 ? 1 : -1;

        a = *parms[2] / exp_digit;
        b = *parms[2] % exp_digit;
        dl_is1 = abs(a) * pow(10, sign * b);

        a = *parms[3] / exp_digit;
        b = *parms[3] % exp_digit;
        max_inc_t = abs(a) * pow(10,sign * b);

        a = *parms[4] / exp_digit;
        b = *parms[4] % exp_digit;
        dl_min = abs(a) * pow(10, sign * b);

        r0 = *parms[5];
        max_itr = *parms[6];
        max_incr = *parms[7];

        a = *parms[8] / exp_digit;
        b = *parms[8] % exp_digit;
        tol = abs(a) * pow(10, sign * b);
    }
    else{
        dl_is1 = 50;
        r0 = 14;
        max_itr = 16;
        max_incr = 50;
        max_inc_t = 1.2;
        dl_min = 0.01;
        tol = 1E-8;
        write_dis = true;
        path = "../Result/nl_demo/";
        name = "Disp_ans_MLSM";
    }

    KF.init_KFSN(eles);
    try{
        Var_his = new double[KF.K_sparse.rows]{0.};
        Var_inc = new double[KF.K_sparse.rows] {0.};
        Q0 = new double[KF.K_sparse.rows]{0.};
    }
    catch(const std::bad_alloc& e){
        std::cerr << e.what() << '\n';
        exit(1);
    }
    bnd.add_bndS(KF.K_sparse, Q0, eles);
    bnd.dbd2ds.clear();
    bnd.dbd2ds.shrink_to_fit();

    ri = r0;
    p = 0;
    incr = 0;
    dl_i = sqrt(ri * 1.0 / r0) * dl_is1;
    while(!is_end){
        is_singular = false;
        for(i = 0; i < KF.K_sparse.rows; i++){
            Var_his[i] = Var[i];
        }

        itr = 0;
        // resid_his = 1E12;
        is_conver = false;
        // is_write = false;
        printf("\nIncrement Step %i start: \nCurrent: %.6f; Increment: %.6f\n", incr+1, p, dl_i);
        while(itr < max_itr){
            itr++;
            //for (i = 0; i < KF.K_sparse.rows; i++) {
            //    printf("%.8f\n", KF.Fint[i]);
            //}
            for(i = 0; i < KF.K_sparse.rows; i++){
                KF.Fout[i] = (p + dl_i) * Q0[i] - KF.Fint[i];
            }
            bnd.add_bndS(KF.K_sparse, KF.Fout, eles);
            resid = 0;
            for(i = 0; i < KF.K_sparse.rows; i++){
                resid += KF.Fout[i] * KF.Fout[i];
            }

            // if (resid > resid_his) is_write = true;
            // if (is_write) {
            //     write_out(eles, xyz, Var);
            // }

            is_singular = !solveeqnS(KF.K_sparse, KF.Fout, Var_inc);
            if(is_singular){
                printf("  Current Load: %.6f Unit\nIncrement set to %.6f and retry\n", p, dl_i/2.0);
                break;
                // KF.K_sparse.type = -2; itr--;
                // continue;
            }

            if (resid <= tol) {
                printf("  Residual: %.8f\nIncrement Step %i done\n", resid, incr+1);
                p = p + dl_i;
                is_conver = true;
                break;
            }
            printf("  NR iteration %i, Residual: %.8f\n", itr + 1, resid);
            // resid_his = resid;

            for (i = 0; i < KF.K_sparse.rows; i++) {
                Var[i] += Var_inc[i];
            }
            KF.init_KFSN(eles);
        }

        ri = itr;
        if(is_conver){
            incr++;
            dl_is1 = dl_i;
            double times = sqrt(r0 * 1.0 / ri);
            dl_i = (max_inc_t > times ? times : max_inc_t) * dl_is1;
        }
        else{
            for (i = 0; i < KF.K_sparse.rows; i++) {
                Var[i] = Var_his[i];
            }
            dl_i /= 2.0;
            // dl_i = sqrt(r0 * 1.0 / ri) * dl_i;
            KF.init_KFSN(eles);
        }
        if(dl_i <= dl_min || incr >= max_incr){
            is_end = true;
            write_out(eles, xyz, Var, p);
        }
    }
};

void solution::solve_SN_GALM(stiffness& KF, boundaries& bnd, elements& eles, double* xyz, int** parms){
    // bnd: 02 para: r0:14; max_itr:16; dl0:30; tol:1E-8; distur_magni: 1E-4 order1, max_inc_t: 1.3, dl_min: 0.005, alpha:0
    // bnd: 02 para: r0:14; max_itr:16; dl0:50; tol:1E-8; distur_magni: 1E-4 order2, max_inc_t: 1.2, dl_min: 0.005, alpha:0
    // bnd: 02 para: r0:14; max_itr:16; dl0:50; tol:1E-8; distur_magni: 3E-4 order3, max_inc_t: 1.2, dl_min: 0.005, alpha:0
    // bnd: 02 para: r0:14; max_itr:16; dl0:50; tol:1E-8; distur_magni: 5E-4 order4, max_inc_t: 1.2, dl_min: 0.005, alpha:0
    // bnd: 02 para: r0:14; max_itr:16; dl0:50; tol:1E-8; distur_magni: 6E-4 order5, max_inc_t: 1.2, dl_min: 0.005, alpha:0

    // bnd: 11 para: r0:14; max_itr:16; dl0:30; tol:1E-8; distur_magni: 3E-4 order1, max_inc_t: 1.2, dl_min: 0.005, alpha:0
    // bnd: 11 para: r0:14; max_itr:16; dl0:30; tol:1E-8; distur_magni: 4E-4 order2, max_inc_t: 1.2, dl_min: 0.005, alpha:0

    // bnd: 12 para: r0:14; max_itr:16; dl0:30; tol:1E-8; distur_magni: 1E-4 order1, max_inc_t: 1.2, dl_min: 0.005, alpha:0
    // bnd: 12 para: r0:14; max_itr:16; dl0:30; tol:1E-8; distur_magni: 2E-4 order2, max_inc_t: 1.2, dl_min: 0.005, alpha:0
    // bnd: 12 para: r0:14; max_itr:16; dl0:30; tol:1E-8; distur_magni: 3E-4 order3, max_inc_t: 1.2, dl_min: 0.005, alpha:0

    // bnd: 22 para: r0:14; max_itr:16; dl0:30; tol:1E-8; distur_magni: 3E-4 order1, max_inc_t: 1.2, dl_min: 0.005, alpha:0
    // bnd: 22 para: r0:14; max_itr:16; dl0:30; tol:1E-8; distur_magni: 3E-4 order2, max_inc_t: 1.2, dl_min: 0.005, alpha:0
    // bnd: 22 para: r0:14; max_itr:16; dl0:30; tol:1E-8; distur_magni: 8E-4 order3, max_inc_t: 1.2, dl_min: 0.005, alpha:0
    int r0, ri, max_itr, mode,
        i, itr, incr,
        max_incr;
    double alpha, dl_is1, dl_i, u_max,
           tol, resid, p_t0, dp_n, str_len,
           p_t1st0, A, B, C, u02, u0u, u0du,
           u2, du2, udu, delta, psln1, psln2,
           /* resid_his,  */max_inc_t, dl_min;
    double *Var_t0 = nullptr, *Q0 = nullptr,
           *Var_u = nullptr, *Var_t1st0 = nullptr,
           *Var_0 = nullptr;
    bool is_end = false, is_singular, is_conver/* , is_write */,
         is_refresh, plot_curve, is_default = true,
         write_dis;
    std::string path_cur, name_cur,
                path_dis, name_dis;

    if(parms){
        if(*parms){
            if(**parms != 0){
                is_default = false;
            }
        }
    }    

    if(*parms){
        double a, b;
        int exp_digit = pow(10, abs(**parms));
        int sign = **parms > 0 ? 1 : -1;

        a = *parms[1] / exp_digit;
        b = *parms[1] % exp_digit;
        alpha = abs(a) * pow(10,sign * b);

        a = *parms[2] / exp_digit;
        b = *parms[2] % exp_digit;
        dl_is1 = abs(a) * pow(10,sign * b);

        a = *parms[3] / exp_digit;
        b = *parms[3] % exp_digit;
        max_inc_t = abs(a) * pow(10,sign * b);

        a = *parms[4] / exp_digit;
        b = *parms[4] % exp_digit;
        dl_min = abs(a) * pow(10, sign * b);

        r0 = *parms[5];
        max_itr = *parms[6];
        max_incr = *parms[7];

        a = *parms[8] / exp_digit;
        b = *parms[8] % exp_digit;
        tol = abs(a) * pow(10, sign * b);

        mode = *parms[9];
        str_len = *parms[10];
        plot_curve = *parms[11];

        if(plot_curve){
            if(parms[12]) int2str(parms[12], str_len, path_cur);
            else path_cur = "../Result/nl_demo/";
            if(parms[13]) int2str(parms[13], str_len, name_cur);
            else name_cur = "Load_Disp_curve";
        }

        write_dis = *parms[14];

        if(write_dis){
            if(parms[15]) int2str(parms[15], str_len, path_dis);
            else path_dis = plot_curve ? path_cur : "../Result/nl_demo/";
            if(parms[16]) int2str(parms[16], str_len, name_dis);
            else name_dis = "Disp_ans_GALM";
        }
    }
    else{
        alpha = 0.;
        dl_is1 = 50;
        max_inc_t = 1.2;
        dl_min = 0.005;
        r0 = 14;
        max_itr = 16;
        max_incr = 40;
        tol = 1E-8;
        mode = 1;
        plot_curve = true;
        path_cur = "../Result/nl_demo/";
        name_cur = "Load_Disp_curve";
        write_dis = true;
        path_dis = path_cur;
        name_dis = "Disp_ans_GALM";
    }

    KF.init_KFSN(eles);
    try{
        Var_t0 = new double[KF.K_sparse.rows]{0.};
        Var_u = new double[KF.K_sparse.rows] {0.};
        Var_t1st0 = new double[KF.K_sparse.rows]{0.};
        Var_0 = new double[KF.K_sparse.rows] {0.};
        Q0 = new double[KF.K_sparse.rows]{0.};
    }
    catch(const std::bad_alloc& e){
        std::cerr << e.what() << '\n';
        exit(1);
    }
    bnd.add_bndS(KF.K_sparse, Q0, eles);
    bnd.dbd2ds.clear();
    bnd.dbd2ds.shrink_to_fit();

    ri = r0;
    p_t0 = 0;
    incr = 0;
    dl_i = sqrt(ri * 1.0 / r0) * dl_is1;
    is_refresh = true;
    while(!is_end){
        is_singular = false;
        for(i = 0; i < KF.K_sparse.rows; i++){
            Var_t0[i] = Var[i];
        }

        itr = 0;
        p_t1st0 = 0.;
        is_conver = false;
        // resid_his = 1E12;
        // is_write = false;
        printf("\nIncrement Step %i start: \nCurrent: %.6f; Arc_length: %.6f\n", incr+1, p_t0, dl_i);
        for(i = 0; i < KF.K_sparse.rows; i++){
            KF.Fout[i] = (p_t0 + p_t1st0) * Q0[i] - KF.Fint[i];
        }
        bnd.add_bndS(KF.K_sparse, KF.Fout, eles);
        while(itr < max_itr){
            itr++;

            is_singular = !solveeqnS(KF.K_sparse, Q0, Var_0);

            if(is_singular){
                printf("  Current Load: %.6f Unit\nStiffness Matrix type set to -2 and retry\n", p_t0);
                if(mode > 0){
                    KF.K_sparse.type = -2;
                    itr --;
                    continue;
                }
                else{
                    break;
                }
            }

            for(i = 0; i < KF.K_sparse.rows; i++){
                Var_t1st0[i] = Var[i] - Var_t0[i];
            }
            this->solveeqnS(KF.K_sparse, KF.Fout, Var_u);
            u02 = 0.;
            u0du = 0.;
            u0u = 0.;
            u2 = 0.;
            du2 = 0.;
            udu = 0.;
            for(i = 0; i < KF.K_sparse.rows; i++){
                u02 += Var_0[i] * Var_0[i];
                u0u += Var_0[i] * Var_u[i];
                u0du += Var_0[i] * Var_t1st0[i];
                u2 += Var_u[i] * Var_u[i];
                udu += Var_u[i] * Var_t1st0[i];
                du2 += Var_t1st0[i] * Var_t1st0[i];
            }

            A = 1 + alpha;
            B = 2 * (u0u + u0du) / u02 + 2*alpha*p_t1st0;
            C = (2*udu + u2 + du2) / u02 + alpha*p_t1st0*p_t1st0 - dl_i*dl_i;
            delta = pow(B, 2) - 4 * A * C;

            if(delta < 0){
                printf("delta le 0\n");
                break;
            }

            psln1 = (-B + sqrt(delta)) / 2 / A;
            psln2 = -(B + sqrt(delta)) / 2 / A;
            if (abs(u0du) < 1E-10) u0du = 1.;
            dp_n = u0du > 0 ? psln1 : psln2;
            p_t1st0 += dp_n;

            for (i = 0; i < KF.K_sparse.rows; i++) {
                Var[i] += Var_0[i] * dp_n + Var_u[i];
            }
            KF.init_KFSN(eles);

            for(i = 0; i < KF.K_sparse.rows; i++){
                KF.Fout[i] = (p_t0 + p_t1st0) * Q0[i] - KF.Fint[i];
            }
            bnd.add_bndS(KF.K_sparse, KF.Fout, eles);

            resid = 0;
            for(i = 0; i < KF.K_sparse.rows; i++){
                resid += KF.Fout[i] * KF.Fout[i];
            }
            // if (resid > resid_his) is_write = true;
            // if (is_write && incr >= 22) {
            //    write_out(eles, xyz, Var);
            // }
            // resid_his = resid;

            printf("  NR iteration %i, Residual: %.8f\n", itr + 1, resid);
            if (resid <= tol) {
                printf("  Residual: %.8f\nIncrement Step %i done\n", resid, incr+1);
                p_t0 = p_t0 + p_t1st0;
                is_conver = true;
                break;
            }
        }

        ri = itr;
        if(is_conver){
            incr++;
            dl_is1 = dl_i;
            double times = sqrt(r0 * 1.0 / ri);
            dl_i = (max_inc_t > times ? times : max_inc_t) * dl_is1;
            
            u_max = 0.;
            for(i = 0; i < KF.K_sparse.rows / 2; i++){
                if(u_max < abs(Var[2*i])){
                    u_max = abs(Var[2*i]);
                }
            }
            if(plot_curve) write_point(&u_max, &p_t0, is_refresh, 1, path_cur, name_cur);
            is_refresh = false;
        }
        else{
            for (i = 0; i < KF.K_sparse.rows; i++) {
                Var[i] = Var_t0[i];
            }
            dl_i = dl_i / 2;
            KF.init_KFSN(eles);
        }
        if(dl_i <= dl_min || incr >= max_incr){
            is_end = true;
            if(write_dis) write_out(eles, xyz, Var, p_t0, path_dis, name_dis);
        }
    }

    delete[] Var_0;
    delete[] Var_t0;
    delete[] Var_t1st0;
    delete[] Q0;
    delete[] Var_u;
};

bool solution::solveeqnS(const SparseMatrix& A, double* b, double* x){
    int mtype = A.type;
    int nrhs = 1;
    int dof = A.rows;
    void* ptsolver[64];
    pardiso_cfg para = pardiso_cfg(mtype);
    double ddum;

    for(int i = 0; i < 64; i++){
        ptsolver[i] = 0; 
    }

    int phase = 11;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(A.type), &phase, &dof, A.val, 
            A.row_st, A.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), &ddum, &ddum, &(para.error));
    if(para.error != 0 ){
        printf ("\033[1;31mERROR\033[0m during symbolic factorization: %i \n", para.error);
    }
    
    phase = 22;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(A.type), &phase, &dof, A.val, 
            A.row_st, A.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), &ddum, &ddum, &(para.error));
    if(para.error != 0 ){
        printf ("    \033[1;31mERROR\033[0m during numerical factorization: %i \n", para.error);
        if(para.error == -4){
            return false;
        }
        else{
            exit(0);
        }
    }

    phase = 33;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(A.type), &phase, &dof, A.val, 
            A.row_st, A.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), b, x, &(para.error));
    if(para.error != 0 ){
        printf ("\033[1;31mERROR\033[0m during solution: %i \n", para.error);
    }
    printf("    solve completed ...\n");

    phase = -1;
    pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(A.type), &phase, &dof, A.val, 
            A.row_st, A.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), b, &ddum, &(para.error));
    if(para.error != 0 ){
        printf ("\033[1;31mERROR\033[0m during free memory: %i \n", para.error);
    }
    return true;
};

// double solution::eig_IPM(SparseMatrix& A, SparseMatrix& B, double* Phi){
//     double* phi_n = nullptr;
//     bool is_tmp = true;

//     a_int nev, ncv, ldv, ido, k,
//         lworkl, info, ierr,
//         iparam[11], ipntr[11];
//     double sigma, *pcr = nullptr;

//     a_int maxitr = 1000,
//         ishfts = 1,
//         mode = 3;
//     double tol = 1E-6;
//     char bmat = 'G',
//          which[2] = {'L','M'};

//     double *x = nullptr, *workl = nullptr, *workd = nullptr,
//            *resid = nullptr, *D = nullptr, *V = nullptr;
//     a_int *select = nullptr;

//     nev = 1;
//     ncv = 2 * nev;
//     ldv = A.rows;
//     lworkl = ncv * (ncv + 8);
    
//     //allocate
//     if(!pcr){
//         pcr = new double[nev]{0.};
//     }
//     if(Phi){
//         phi_n = Phi;
//         is_tmp = false;
//     }
//     else{
//         phi_n = new double[A.rows * nev]{0.};
//     }

//     x = new double[A.rows]{0.};
//     D = new double[ncv * 2 + 1]{0.};
//     V = new double[ldv * ncv + 1]{0.};
//     workd = new double[3 * A.rows + 1]{0.};
//     workl = new double[lworkl + 1]{0.};
//     resid = new double[A.rows + 1]{0.};
//     select = new int[ncv + 1]{0};

//     info = 0; ido = 0;
//     iparam[0] = ishfts;
//     iparam[2] = maxitr;
//     iparam[6] = mode;
//     sigma = 0.0;

//     k = 0;
//     void* ptsolver[64];
//     int nrhs = 1;
//     pardiso_cfg para = pardiso_cfg(A.type);
//     double ddum;

//     for (int i = 0; i < 64; i++) {
//         ptsolver[i] = 0;
//     }

//     int phase = 11;
//     pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(A.type), &phase, &A.rows, A.val,
//         A.row_st, A.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), &ddum, &ddum, &(para.error));
//     if (para.error != 0) {
//         printf("\033[1;31mERROR\033[0m during symbolic factorization: %i \n", para.error);
//     }

//     phase = 22;
//     pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(A.type), &phase, &A.rows, A.val,
//         A.row_st, A.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), &ddum, &ddum, &(para.error));
//     if (para.error != 0) {
//         printf("\033[1;31mERROR\033[0m during numerical factorization: %i \n", para.error);
//     }

//     while(true){
//         k++;
//         dsaupd_c(&ido, &bmat, A.rows, which, nev, tol, resid, ncv, V, ldv,
//                  iparam, ipntr, workd, workl, lworkl, &info);

//         if(ido == -1){
//             B.product(workd + ipntr[0] - 1, x);

//             phase = 33;
//             pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(A.type), &phase, &A.rows, A.val, 
//                     A.row_st, A.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), x, workd+ipntr[1] - 1, &(para.error));
//             if(para.error != 0 ){
//                 printf ("\033[1;31mERROR\033[0m during solution: %i \n", para.error);
//             }
//         }
//         else if(ido == 1){
//             phase = 33;
//             pardiso(ptsolver, &(para.maxfct), &(para.mnum), &(A.type), &phase, &A.rows, A.val, 
//                     A.row_st, A.col, &(para.perm), &nrhs, para.iparm, &(para.msglvl), workd+ipntr[2] - 1, workd+ipntr[1] - 1, &(para.error));
//             if(para.error != 0 ){
//                 printf ("\033[1;31mERROR\033[0m during solution: %i \n", para.error);
//             }
//         }
//         else if(ido == 2){
//             B.product(workd + ipntr[0] - 1, workd + ipntr[1] - 1);
//         }
//         else{
//             break;
//         }
//     }
    
//     // exception
//     if(info < 0){
//         if(info == -7) printf("\033[1;31mERROR\033[0m: Lenth of WORKL is not sufficient\n");
//         if(info == -9) printf("\033[1;31mERROR\033[0m: ARNOLDI starting vectot if zero, \n\
// please try to again\n");

//         if(info < -9 || info > -7 || info == -8){
//             printf("\033[1;31mERROR\033[0m: %i\n", info);
//         }
//     }
//     else{
//         char howmny[3] = {'A'};
//         dseupd_c(true, howmny, select, D, V, ldv, sigma, &bmat, A.rows,
//                  which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd,
//                  workl, lworkl, &ierr);

//         if(ierr == 0){
//             for(int i = 0; i < nev; i++){
//                 pcr[i] = D[i];
//                 for(int j = 0; j < A.rows; j++){
//                     phi_n[i * A.rows + j] = V[j * nev + i];
//                 }
//             }
//         }
//     }

//      int j = 1;
//      for (int i = 0; i < A.rows/2; i++) {
//          printf("% .8f ", phi_n[2*i + 0]);
//          if (j == (1 + 1) * 2 + 1) {
//              printf("\n");
//              j = 0;
//          }
//          j++;
//      }

//     if(is_tmp) delete[] phi_n;
//     phi_n = nullptr;

//     delete[] x;
//     delete[] workd;
//     delete[] workl;
//     delete[] resid;
//     delete[] select;
//     delete[] D;
//     delete[] V;

//     return *pcr;
// };

void write_out(elements& eles, double* xyz, double* Vars, double load, std::string path, std::string name){
    std::ofstream ans_fio;
    trim(path); trim(name);

    for(std::string::iterator itchar = path.begin(); itchar != path.end(); itchar++){
        if(*itchar == '\\') *itchar = '/';
    }
    if(path.back() != '/' || path.back() != '\\') path = path.append({'/'});

    name.append(".out");
    ans_fio.open(path + name);
    ans_fio.setf(ans_fio.scientific);
    ans_fio.precision(10);
    std::vector<std::pair<int, int>> sont = eles.sont();
    std::vector<int> sot = eles.sot();
    int dof = 0, nnum = 0, elenum = 0,
        i, j, etype;
    for(std::vector<std::pair<int, int>>::iterator itsont = sont.begin(); itsont != sont.end(); itsont++){
        dof += itsont->first * itsont->second;
        nnum += itsont->first;
    }
    i = 0;
    for(std::vector<int>::iterator itsot = sot.begin(); itsot != sot.end(); itsot++){
        elenum += *itsot;
        if(*itsot != 0){
            etype = i;
        }
        i++;
    }

    PSL_ ele_ins;
    std::vector<PSL_>* eles_ref = eles.eleset_ptr(ele_ins);

    ans_fio << "GEOMETRY\n" << std::endl;
    ans_fio << "CURRENT_LOAD: " << load << std::endl;
    ans_fio << " COORDINATES " << nnum << std::endl;
    ans_fio << "! NODE  X   Y   Z" << std::endl;
    for(i = 0; i < nnum; i++){
        ans_fio << "    " << i+1 << "    " << xyz[2*i] << "    " << xyz[2*i + 1] << "    " << 0. << std::endl;
    }

    ans_fio << "\n ELEMENTS " << elenum << std::endl;
    ans_fio << "    ELEMENT_NODES" << std::endl;
    for(i = 0; i < elenum; i++){
        ans_fio << "    " << i+1;
        for(j = 0; j < eles_ref->front().nnum; j++){
            ans_fio << "    " << eles_ref->at(i).node_id.at(j);
        }
        ans_fio << std::endl;
    }

    ans_fio << "\n    ELEMENT_MATERIAL\n" << std::endl;
    ans_fio << " *** NODE_VALUE ***" << std::endl;
    ans_fio << "NODE    U    V" << std::endl;
    for(i = 0; i < nnum; i++){
        ans_fio << "    " << i+1 << "    " << Vars[2*i] << "    " << Vars[2*i + 1] << std::endl;
    }

    ans_fio << "\n*END" << std::endl;
    ans_fio.close();
};

void write_point(const double* x, const double* y, bool refresh, int len, std::string path, std::string name){
    int i;
    trim(path);
    trim(name);

    for(std::string::reverse_iterator itchar = path.rbegin(); itchar != path.rend(); itchar++){
        if(*itchar == '\\'){
            *itchar = '/';
        }
    }
    if(path.back() != '/' || path.back() != '\\') path = path.append({'/'});


    std::ofstream out_file_io;
    if(!refresh){
        out_file_io.open(path + name.append(".txt"), std::ios_base::ate | std::ios_base::app);
    }
    else{
        out_file_io.open(path + name.append(".txt"), std::ios::out);
    }
    out_file_io.setf(out_file_io.scientific);
    out_file_io.precision(10);

    for(i = 0; i < len; i++){
        out_file_io << x[i] << "   " << y[i] << std::endl;
    }
    out_file_io.close();
};

void solution::del(){
    delete[] Var; Var = nullptr;
}