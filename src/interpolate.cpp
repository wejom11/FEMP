#include "interpolate.h"

void Gauss_int::getGP_1(double* GP_1, double* weight_1){
    if(order == 1){
        GP_1[0] = 0;
        weight_1[0] = 2;
    }
    else if(order == 2){
        GP_1[0] = -1/sqrt(3); GP_1[1] = 1/sqrt(3);
        weight_1[0] = 1; weight_1[1] = 1;
    }
    else if(order == 3){
        GP_1[0] = -sqrt(0.6); GP_1[1] = 0;
        GP_1[2] = sqrt(0.6);
        weight_1[0] = 5.0/9.0; weight_1[1] = 8.0/9.0;
        weight_1[2] = 5.0/9.0;
    }
    else{
        printf("integration order %i is not supported!\n", order);
    }
};

void Gauss_int::getGP(){
    double* GP1 = new double[order];
    double* wei1 = new double[order];
    this->getGP_1(GP1, wei1);
    int sds = selected_dim.size();

    int_pt_num = pow(order, sds);
    int_pt_coord = new double[int_pt_num * dim]{0.0};
    weightness = new double[int_pt_num]{0.0};

    int repeat = 0, i = 0, j = 0, k = 0,
        row = 0, col = 0, loop = 0, row_w = 0;
    repeat = 1;
    loop = pow(order, sds - 1);
    double wei_ini_val = pow(2, dim - sds);
    for(i = 0; i < int_pt_num; i++){
        weightness[i] = wei_ini_val;
    }

    for(std::vector<int>::iterator itsd = selected_dim.begin(); itsd != selected_dim.end(); itsd++){
        row = 0;
        row_w = 0;
        for(k = 0; k < loop; k++){
            for(i = 0; i < order; i++){
                for(j = 0; j < repeat; j++){
                    int_pt_coord[row + *itsd - 1] = GP1[i];
                    weightness[row_w] *= wei1[i];
                    row += dim;
                    row_w ++;
                }
            }
        }

        repeat *= order;
        loop /= order;
    }

    delete[] GP1; GP1 = nullptr;
    delete[] wei1; wei1 = nullptr;
}

void Gauss_int::del(){
    if(int_pt_coord){
        delete[] int_pt_coord; int_pt_coord = nullptr;
    }
    if(weightness){
        delete[] weightness; weightness = nullptr;
    }
}

void Lagrange_intp::set_integration(int int_od, std::vector<int> sd){
    std::vector<int> ssd;
    if(sd.size() == 0){
        ssd.resize(dim);
        for(int i = 0; i < dim; i++){
            ssd.at(i) = i+1;
        }
    };
    gintp = Gauss_int(int_od, dim, ssd);
    gintp.getGP();
};

void Lagrange_intp::getN__(){
    if(dim == 2){
        if(order == 2){
            Ni = new double[gintp.int_pt_num * 9];
            int i, j, k, row;
            int nc[18] = {0,0,2,0,2,2,0,2,1,0,2,1,1,2,0,1,1,1};
            for(i = 0; i < gintp.int_pt_num; i++){
                row = i * 9;
                double xi = gintp.int_pt_coord[2*i];
                double xi2 = xi * xi;
                double eta = gintp.int_pt_coord[2*i + 1];
                double eta2 = eta * eta;
                double l2xi[3] = {(xi2 - xi)/2, 1-xi2, (xi2 + xi)/2};
                double l2eta[3] = {(eta2 - eta)/2, 1-eta2, (eta2 + eta)/2};
                for(j = 0; j < 9; j++){
                    Ni[row + j] = l2xi[nc[2*j]] * l2eta[nc[2*j+1]];
                }
            }
        }
    }
};

void Lagrange_intp::getNdl__(double* Ndl){
    if(dim == 2){
        if(order == 2){
            int i, j, k, row, J3;
            int nc[18] = {0,0,2,0,2,2,0,2,1,0,2,1,1,2,0,1,1,1};
            for(i = 0; i < gintp.int_pt_num; i++){
                row = i * 18;
                double xi = gintp.int_pt_coord[2*i];
                double xi2 = xi * xi;
                double eta = gintp.int_pt_coord[2*i + 1];
                double eta2 = eta * eta;
                double l2xi[3] = {(xi2 - xi)/2, 1-xi2, (xi2 + xi)/2};
                double l2eta[3] = {(eta2 - eta)/2, 1-eta2, (eta2 + eta)/2};
                double l2xidxi[3] = {xi - 0.5, -2 * xi, xi + 0.5};
                double l2etadeta[3] = {eta - 0.5, -2 * eta, eta + 0.5};
                for(j = 0; j < 9; j++){
                    Ndl[row + 2*j] = l2xidxi[nc[2*j]] * l2eta[nc[2*j+1]];
                    Ndl[row + 2*j + 1] = l2xi[nc[2*j]] * l2etadeta[nc[2*j+1]];
                }
            }
        }
    }
};

void Lagrange_intp::getNdx__(std::vector<double*> xyz_c){
    int pt_num = pow(order + 1, dim);
    Ndx = new double[gintp.int_pt_num * pt_num * dim];
    Jacobien = new double[gintp.int_pt_num];

    getNdl__(Ndx);

    int i, j, k, intp, Jdim, row, id = 0;
    int* ipiv = new int[dim]{0};
    int info = 0;
    double val = 0.0;
    DenseMatrix Jaco(dim, dim);
    for(intp = 0; intp < gintp.int_pt_num; intp++){
        row = pt_num * dim * intp;
        id = 0;
        for(i = 0; i < dim; i++){
            for(j = 0; j < dim; j++){
                val = 0.;
                for(k = 0; k < pt_num; k++){
                    val += Ndx[row + k * dim + i] * xyz_c.at(k)[j];
                }
                Jaco.val[id] = val;
                id++;
            }
        }

        Jacobien[intp] = abs(det(Jaco));

        info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, Jaco.rows, Jaco.cols, Jaco.val, Jaco.cols, ipiv);
        if(info != 0){
            printf("ERROR: calculating Ndx failed!\n");
            exit(1);
        }

        for(j = 0; j < pt_num; j++){
            info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', Jaco.rows, 1, Jaco.val, Jaco.cols, ipiv, &Ndx[row + dim*j], 1);

            if(info != 0){
                printf("ERROR: calculating Ndx failed!\n");
                exit(info);
            }
        }
    }

    delete[] ipiv; ipiv = nullptr;
}

void Lagrange_intp::del(){
    if(Ni){
        delete[] Ni; Ni = nullptr;
    }
    if(Jacobien){
        delete[] Jacobien; Jacobien = nullptr;
    }
    if(Ndx){
        delete[] Ndx; Ndx = nullptr;
    }
};

void Serendipity_intp::set_integration(int int_od, std::vector<int> sd){
    std::vector<int> ssd;
    if(sd.size() == 0){
        ssd.resize(dim);
        for(int i = 0; i < dim; i++){
            ssd.at(i) = i+1;
        }
    };
    gintp = Gauss_int(int_od, dim, ssd);
    gintp.getGP();
};

void Serendipity_intp::getN__(){
    if(dim == 2){
        if(order == 2){
            Ni = new double[gintp.int_pt_num * 8];
            int i, row;
            for(i = 0; i < gintp.int_pt_num; i++){
                row = i * 8;
                double xi = gintp.int_pt_coord[2*i];
                double xi2 = xi * xi;
                double eta = gintp.int_pt_coord[2*i + 1];
                double eta2 = eta * eta;
                double N1bar = (1 - xi) * (1 - eta) / 4;
                double N2bar = (1 + xi) * (1 - eta) / 4;
                double N3bar = (1 + xi) * (1 + eta) / 4;
                double N4bar = (1 - xi) * (1 + eta) / 4;

                Ni[row + 4] = (1 - xi2) * (1 - eta) / 2;
                Ni[row + 5] = (1 - eta2) * (1 + xi) / 2;
                Ni[row + 6] = (1 - xi2) * (1 + eta) / 2;
                Ni[row + 7] = (1 - eta2) * (1 - xi) / 2;
                Ni[row] = N1bar - (Ni[row + 4] + Ni[row + 7]) / 2;
                Ni[row + 1] = N2bar - (Ni[row + 4] + Ni[row + 5]) / 2;
                Ni[row + 2] = N3bar - (Ni[row + 6] + Ni[row + 5]) / 2;
                Ni[row + 3] = N4bar - (Ni[row + 6] + Ni[row + 7]) / 2;
            }
        }
    }
};

void Serendipity_intp::getNdl__(double* Ndl){
    if(dim == 2){
        if(order == 2){
            int i, j, k, row, J3;
            for(i = 0; i < gintp.int_pt_num; i++){
                row = i * 16;
                double xi = gintp.int_pt_coord[2*i];
                double xi2 = xi * xi;
                double eta = gintp.int_pt_coord[2*i + 1];
                double eta2 = eta * eta;
                double N1bardx = - (1 - eta) / 4;
                double N2bardx = (1 - eta) / 4;
                double N3bardx = (1 + eta) / 4;
                double N4bardx = -(1 + eta) / 4;
                double N1barde = -(1 - xi) / 4;
                double N2barde = -(1 + xi) / 4;
                double N3barde = (1 + xi) / 4;
                double N4barde = (1 - xi) / 4;

                Ndl[row + 8] = -xi * (1 - eta);
                Ndl[row + 10] = (1 - eta2) / 2;
                Ndl[row + 12] = -xi * (1 + eta);
                Ndl[row + 14] = -(1 - eta2) / 2;
                Ndl[row + 9] = -(1 - xi2) / 2;
                Ndl[row + 11] = -eta * (1 + xi);
                Ndl[row + 13] = (1 - xi2) / 2;
                Ndl[row + 15] = -eta * (1 - xi);
                Ndl[row] = N1bardx - (Ndl[row + 8] + Ndl[row + 14]) / 2;
                Ndl[row + 2] = N2bardx - (Ndl[row + 8] + Ndl[row + 10]) / 2;
                Ndl[row + 4] = N3bardx - (Ndl[row + 12] + Ndl[row + 10]) / 2;
                Ndl[row + 6] = N4bardx - (Ndl[row + 12] + Ndl[row + 14]) / 2;
                Ndl[row + 1] = N1barde - (Ndl[row + 9] + Ndl[row + 15]) / 2;
                Ndl[row + 3] = N2barde - (Ndl[row + 9] + Ndl[row + 11]) / 2;
                Ndl[row + 5] = N3barde - (Ndl[row + 13] + Ndl[row + 11]) / 2;
                Ndl[row + 7] = N4barde - (Ndl[row + 13] + Ndl[row + 15]) / 2;
            }
        }
    }
};

void Serendipity_intp::getNdx__(std::vector<double*> xyz_c){
    int pt_num = pow(order + 1, dim) - pow(order - 1, dim);
    Ndx = new double[gintp.int_pt_num * pt_num * dim];
    Jacobien = new double[gintp.int_pt_num];

    getNdl__(Ndx);

    int i, j, k, intp, Jdim, row, id = 0;
    int* ipiv = new int[dim]{0};
    int info = 0;
    double val = 0.0;
    DenseMatrix Jaco(dim, dim);
    for(intp = 0; intp < gintp.int_pt_num; intp++){
        row = pt_num * dim * intp;
        id = 0;
        for(i = 0; i < dim; i++){
            for(j = 0; j < dim; j++){
                val = 0.;
                for(k = 0; k < pt_num; k++){
                    val += Ndx[row + k * dim + i] * xyz_c.at(k)[j];
                }
                Jaco.val[id] = val;
                id++;
            }
        }

        Jacobien[intp] = abs(det(Jaco));

        info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, Jaco.rows, Jaco.cols, Jaco.val, Jaco.cols, ipiv);
        if(info != 0){
            printf("ERROR: calculating Ndx failed!\n");
            exit(1);
        }

        for(j = 0; j < pt_num; j++){
            info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', Jaco.rows, 1, Jaco.val, Jaco.cols, ipiv, &Ndx[row + dim*j], 1);

            if(info != 0){
                printf("ERROR: calculating Ndx failed!\n");
                exit(info);
            }
        }
    }

    delete[] ipiv; ipiv = nullptr;
}

void Serendipity_intp::del(){
    if(Ni){
        delete[] Ni; Ni = nullptr;
    }
    if(Jacobien){
        delete[] Jacobien; Jacobien = nullptr;
    }
    if(Ndx){
        delete[] Ndx; Ndx = nullptr;
    }
};