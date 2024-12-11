#include "PSL_.h"
                                                                                               
void PSL_::set_order(bool ord){
    order = ord ? 1 : 0;
    if(order){
        gauss_num = 5;
        node_id.resize(9, 0); node_id.shrink_to_fit();
        xy.resize(9, nullptr); xy.shrink_to_fit();
        uv.resize(9, nullptr); uv.shrink_to_fit();
    }
    else{
        gauss_num = 3;
        node_id.resize(4, 0); node_id.shrink_to_fit();
        xy.resize(4, nullptr); xy.shrink_to_fit();
        uv.resize(4, nullptr); uv.shrink_to_fit();
    }
};

void PSL_::set_gauss_order(short order){
    gauss_num = order;
};

void PSL_::set_xy_uv(double* global_xy, double* global_uv, const std::vector<int>& n_id){
    node_id = n_id;

    set_xy_uv(global_xy, global_uv);
};

void PSL_::set_xy_uv(double* global_xy, double* global_uv){
    int i, posi;
    for(i = 0; i < node_id.size(); i++){
        posi = 2*(node_id.at(i) - 1);
        xy.at(i) = global_xy + posi;
        uv.at(i) = global_uv + posi;
    }
}

void PSL_::ele_init(){
    shape_info = Lagrange_intp(order + 1, 2);
    shape_info.set_integration(gauss_num);
    shape_info.getNdx__(xy);
};

void PSL_::get_BL0T(DenseMatrix& BL0T, short intp){
    int i = 0, indice;
    int row = intp * nnum * 2;

    for(i = 0; i < nnum; i++){
        indice = 6*i;
        BL0T.val[indice] = shape_info.Ndx[row + 2*i];
        BL0T.val[indice + 2] = shape_info.Ndx[row + 2*i + 1];
        BL0T.val[indice + 4] = shape_info.Ndx[row + 2*i + 1];
        BL0T.val[indice + 5] = shape_info.Ndx[row + 2*i];
    }
};

void PSL_::get_BL1T(DenseMatrix& BL1T, short intp){
    int i = 0, indice;
    int row = intp * nnum * 2;
    double L_cons[4]{0.};

    for(i = 0; i < nnum; i++){
        L_cons[0] += shape_info.Ndx[row + 2*i] * uv.at(i)[0];
        L_cons[1] += shape_info.Ndx[row + 2*i + 1] * uv.at(i)[0];
        L_cons[2] += shape_info.Ndx[row + 2*i] * uv.at(i)[1];
        L_cons[3] += shape_info.Ndx[row + 2*i + 1] * uv.at(i)[1];
    }

    for(i = 0; i < nnum; i++){
        indice = 6*i;
        BL1T.val[indice] = L_cons[0] * shape_info.Ndx[row + 2*i];
        BL1T.val[indice + 1] = L_cons[1] * shape_info.Ndx[row + 2*i + 1];
        BL1T.val[indice + 2] = L_cons[0] * shape_info.Ndx[row + 2*i + 1] + L_cons[1] * 
                               shape_info.Ndx[row + 2*i];
        BL1T.val[indice + 3] = L_cons[2] * shape_info.Ndx[row + 2*i];
        BL1T.val[indice + 4] = L_cons[3] * shape_info.Ndx[row + 2*i + 1];
        BL1T.val[indice + 5] = L_cons[2] * shape_info.Ndx[row + 2*i + 1] + L_cons[3] * 
                               shape_info.Ndx[row + 2*i];
    }
};

void PSL_::get_BNLT(DenseMatrix& BNLT, short intp){
    int i = 0, indice;
    int row = intp * nnum * 2;

    for(i = 0; i < nnum; i++){
        indice = 8*i;
        BNLT.val[indice] = shape_info.Ndx[row + 2*i];
        BNLT.val[indice + 1] = shape_info.Ndx[row + 2*i + 1];
        BNLT.val[indice + 6] = shape_info.Ndx[row + 2*i];
        BNLT.val[indice + 7] = shape_info.Ndx[row + 2*i + 1];
    }
};

void PSL_::get_S(DenseMatrix& St, short intp){
    DenseMatrix D(3, 3);
    D.val[0] = 1.0;
    D.val[1] = mat->nu;
    D.val[3] = mat->nu;
    D.val[4] = 1.0;
    D.val[8] = (1.0 - mat->nu) / 2;
    double du[2][2]{0.};
    
    int row = 2 * nnum * intp;
    int i, j, k;
    for(i = 0; i < 2; i++){
        for(j = 0; j < 2; j++){
            du[i][j] = 0;
            for(k = 0; k < nnum; k++){
                du[i][j] += shape_info.Ndx[row + 2*k + j] * uv.at(k)[i];
            }
        }
    }

    double eps[3]{0.};
    double S_bar[3]{0.};
    eps[0] = du[0][0] + 0.5*(pow(du[0][0], 2) + pow(du[1][0], 2));
    eps[1] = du[1][1] + 0.5*(pow(du[0][1], 2) + pow(du[1][1], 2));
    eps[2] = du[0][1] + du[1][0] + du[0][0]*du[0][1] + du[1][0]*du[1][1];
    for(i = 0; i < 3; i++){
        S_bar[i] = 0.;
        k = 3 * i;
        for(j = 0; j < 3; j++){
            S_bar[i] += D.val[k + j] * eps[j];
        }
    }

    for(i = 0; i < 2; i++){
        j = 10 * i;

        St.val[0 + j] = S_bar[0];
        St.val[1 + j] = S_bar[2];
        St.val[4 + j] = S_bar[2];
        St.val[5 + j] = S_bar[1];
    }
};

void PSL_::get_S(double* Sbar, short intp){
    DenseMatrix D(3, 3);
    D.val[0] = 1.0;
    D.val[1] = mat->nu;
    D.val[3] = mat->nu;
    D.val[4] = 1.0;
    D.val[8] = (1.0 - mat->nu) / 2;
    double du[2][2]{0.};
    
    int row = 2 * nnum * intp;
    int i, j, k;
    for(i = 0; i < 2; i++){
        for(j = 0; j < 2; j++){
            du[i][j] = 0;
            for(k = 0; k < nnum; k++){
                du[i][j] += shape_info.Ndx[row + 2*k + j] * uv.at(k)[i];
            }
        }
    }

    double eps[3]{0.};
    double eta[3]{0.};
    //eta[0] = 0.5 * (pow(du[0][0], 2) + pow(du[1][0], 2));
    //eta[1] = 0.5 * (pow(du[0][1], 2) + pow(du[1][1], 2));
    //eta[2] = (du[0][0] * du[0][1] + du[1][0] * du[1][1]);
    eps[0] = du[0][0] + 0.5*(pow(du[0][0], 2) + pow(du[1][0], 2));
    eps[1] = du[1][1] + 0.5*(pow(du[0][1], 2) + pow(du[1][1], 2));
    eps[2] = du[0][1] + du[1][0] + du[0][0]*du[0][1] + du[1][0]*du[1][1];
    for(i = 0; i < 3; i++){
        Sbar[i] = 0.;
        k = 3 * i;
        for(j = 0; j < 3; j++){
            Sbar[i] += D.val[k + j] * eps[j];
        }
    }
    //for (i = 0; i < 3; i++) {
    //    std::cout << std::setprecision(10) << eps[i] << std::endl;
    //}
};

void PSL_::get_KL0(DenseMatrix& KL0){
    DenseMatrix BL0T(nnum * 2, 3);
    DenseMatrix D(3, 3);
    D.val[0] = 1.0;
    D.val[1] = mat->nu;
    D.val[3] = mat->nu;
    D.val[4] = 1.0;
    D.val[8] = (1.0 - mat->nu) / 2;

    DenseMatrix M_temp(nnum * 2, 3);
    DenseMatrix KL0itr(nnum * 2, nnum * 2);
    int intp;
    for(intp = 0; intp < shape_info.gintp.int_pt_num; intp++){
        M_temp.set_zero();
        KL0itr.set_zero();

        get_BL0T(BL0T, intp);

        BL0T.product(CblasNoTrans, D, CblasNoTrans, M_temp, 
                     shape_info.Jacobien[intp] * shape_info.gintp.weightness[intp]);
        M_temp.product(CblasNoTrans, BL0T, CblasTrans, KL0itr);

        KL0.add(KL0itr);
    }
};

void PSL_::get_KL1(DenseMatrix& KL1, short option){
    bool is_u1 = option % 10 > 0, is_u2 = (option / 10) % 10 > 0;

    DenseMatrix BL0T(nnum * 2, 3);
    DenseMatrix BL1T(nnum * 2, 3);
    DenseMatrix D(3, 3);
    D.val[0] = 1.0;
    D.val[1] = mat->nu;
    D.val[3] = mat->nu;
    D.val[4] = 1.0;
    D.val[8] = (1.0 - mat->nu) / 2;
    int intp;

    DenseMatrix M_temp(nnum * 2, 3);
    DenseMatrix KL1itr(nnum * 2, nnum * 2);
    DenseMatrix Kuitr(nnum * 2, nnum * 2);
    for(intp = 0; intp < shape_info.gintp.int_pt_num; intp++){
        get_BL0T(BL0T, intp);
        get_BL1T(BL1T, intp);
        KL1itr.set_zero();

        //for (int i = 0; i < 3; i++) {
        //    for (int j = 0; j < 18; j++) {
        //        printf("%.15f\n", BL1T.val[3*j + i]);
        //    }
        //}

        if(is_u1){
            M_temp.set_zero();
            Kuitr.set_zero();

            BL0T.product(CblasNoTrans, D, CblasNoTrans, M_temp);
            M_temp.product(CblasNoTrans, BL1T, CblasTrans, Kuitr);
            KL1itr.add(Kuitr);

            M_temp.set_zero();
            Kuitr.set_zero();

            BL1T.product(CblasNoTrans, D, CblasNoTrans, M_temp);
            M_temp.product(CblasNoTrans, BL0T, CblasTrans, Kuitr);
            KL1itr.add(Kuitr);
        }

        if(is_u2){
            M_temp.set_zero();
            Kuitr.set_zero();

            BL1T.product(CblasNoTrans, D, CblasNoTrans, M_temp);
            M_temp.product(CblasNoTrans, BL1T, CblasTrans, Kuitr);
            KL1itr.add(Kuitr);
        }

        KL1itr.mul(shape_info.gintp.weightness[intp] * shape_info.Jacobien[intp]);
        KL1.add(KL1itr);
    }
};

void PSL_::get_KNL(DenseMatrix& KNL){
    DenseMatrix BNLT(nnum * 2, 4);
    DenseMatrix S(4, 4);

    DenseMatrix M_temp(nnum * 2, 4);
    DenseMatrix KNLitr(nnum * 2, nnum * 2);
    int intp;
    for(intp = 0; intp < shape_info.gintp.int_pt_num; intp++){
        M_temp.set_zero();
        KNLitr.set_zero();
        
        this->get_BNLT(BNLT, intp);
        this->get_S(S, intp);

        //for (int i = 0; i < 4; i++) {
        //    for (int j = 0; j < 18; j++) {
        //        printf("%.15f\n", BNLT.val[4*j + i]);
        //    }
        //}

        BNLT.product(CblasNoTrans, S, CblasNoTrans, M_temp,
                     shape_info.Jacobien[intp] * shape_info.gintp.weightness[intp]);
        M_temp.product(CblasNoTrans, BNLT, CblasTrans, KNLitr);

        KNL.add(KNLitr);
    }
};

void PSL_::get_F(double* Ft){
    DenseMatrix BLT(nnum * 2, 3);
    DenseMatrix M_temp(nnum * 2, 3);
    double Sb[3]{0.};
    int intp, i;

    for(intp = 0; intp < shape_info.gintp.int_pt_num; intp++){
        this->get_S(Sb, intp);
        for(i = 0; i < 3; i++){
            //std::cout << std::setprecision(10) << Sb[i] << std::endl;
            Sb[i] *= shape_info.Jacobien[intp] * shape_info.gintp.weightness[intp];
            //printf("%.18f\n", Sb[i]);
        }

        BLT.set_zero();
        M_temp.set_zero();
        this->get_BL0T(M_temp, intp);
        BLT.add(M_temp);
        M_temp.set_zero();
        this->get_BL1T(M_temp, intp);
        BLT.add(M_temp);

        //for (int i = 0; i < 3; i++) {
        //    for (int j = 0; j < 18; j++) {
        //        printf("%.15f\n", BLT.val[3*j + i]);
        //    }
        //}
        
        for(i = 0; i < nnum * 2; i++){
            Ft[i] += BLT.val[3*i] * Sb[0] + BLT.val[3*i+1] * Sb[1]
                     + BLT.val[3*i+2] * Sb[2];
        }
    }
};

void PSL_::get_KF(DenseMatrix& Ke, double* Fe, short option){
    int i;
    bool is_KL0 = (option / 1000) % 10,
         is_KNL = (option / 100) % 10;

    Ke.set_zero();
    for(i = 0; i < nnum * 2; i++){
        Fe[i] = 0.;
    }
    double D0 = mat->E / (1 - pow(mat->nu, 2));

    if(is_KL0) this->get_KL0(Ke);
    this->get_KL1(Ke, option);
    if(is_KNL) this->get_KNL(Ke);
    this->get_F(Fe);
    Ke.mul(pp->thickness * D0);
    for(i = 0; i < nnum * 2; i++){
        Fe[i] *= pp->thickness * D0;
    }
};