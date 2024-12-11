#include "PQL_.h"

void PQL_::set(std::vector<int> nid, std::vector<double*> xc){
    if(nid.size() != 9 || xc.size() != 9){
        printf("ERROR: input vector's size must be 9!\n");
        return;
    }
    node_id = nid;
    xyz = xc;
};

void PQL_::get_Kb(double Kb[27][27]){
    int inte_b;
    if(inte_sche == 0){
        inte_b = 3;
    }
    else if(inte_sche == 1){
        inte_b = 3;
    }
    else if(inte_sche == 2){
        inte_b = 2;
    }

    Lagrange_intp shape_info(2, 2);
    shape_info.set_integration(inte_b);
    shape_info.getNdx__(xyz);

    double D0 = mat->E * pow(pp->thickness, 3) / 12 / (1 - pow(mat->nu, 2));
    double temp_1 = 0.5 * D0 * (1 - mat->nu);
    double temp_2 = D0 * mat->nu;
    double val1, val2, val3, val4, val5;

    int i, j, I3, J3, k, row,
        I2, J2;
    for(i = 0; i < 27; i++){
        for(j = 0; j < 27; j++){
            Kb[i][j] = 0.;
        }
    }
    
    for(k = 0; k < shape_info.gintp.int_pt_num; k++){
        row = k * 9 * 2;
        val5 = shape_info.gintp.weightness[k] * shape_info.Jacobien[k];
        for(i = 0; i < 9; i++){
            I3 = 3*i;
            I2 = 2*i;
            for(j = i; j < 9; j++){
                J2 = 2*j;
                J3 = 3*j;
                val1 = shape_info.Ndx[row + I2] * shape_info.Ndx[row + J2];
                val2 = shape_info.Ndx[row + I2 + 1] * shape_info.Ndx[row + J2 + 1];
                val3 = shape_info.Ndx[row + I2] * shape_info.Ndx[row + J2 + 1];
                val4 = shape_info.Ndx[row + I2 + 1] * shape_info.Ndx[row + J2];

                Kb[I3][J3] += (D0 * val1 + temp_1 * val2) * val5;
                Kb[I3][J3 + 1] += (temp_1 * val4 + temp_2 * val3) * val5;
                Kb[I3 + 1][J3 + 1] += (temp_1 * val1 + D0 * val2) * val5;
                if(i != j){
                    Kb[I3 + 1][J3] += (temp_1 * val3 + temp_2 * val4) * val5;
                }
            }
        }
    }
};

void PQL_::get_Ks(double Ks[27][27]){
    int inte_s;
    if(inte_sche == 0){
        inte_s = 3;
    }
    else if(inte_sche == 1){
        inte_s = 2;
    }
    else if(inte_sche == 2){
        inte_s = 2;
    }

    Lagrange_intp shape_info(2, 2);
    shape_info.set_integration(inte_s);
    shape_info.getN__();
    shape_info.getNdx__(xyz);

    double alpha = mat->G * pp->thickness / 1.2;
    double val1, val2;

    int i, j, I3, J3, k, row,
        I2, J2, row_n;
    for(i = 0; i < 27; i++){
        for(j = 0; j < 27; j++){
            Ks[i][j] = 0.;
        }
    }
    for(k = 0; k < shape_info.gintp.int_pt_num; k++){
        row = k * 9 * 2;
        row_n = k * 9;

        val2 = shape_info.gintp.weightness[k] * shape_info.Jacobien[k];
        for(i = 0; i < 9; i++){
            I3 = 3*i;
            I2 = 2*i;
            for(j = i; j < 9; j++){
                J2 = 2*j;
                J3 = 3*j;
                val1 = shape_info.Ni[row_n + i] * shape_info.Ni[row_n + j];

                Ks[I3][J3] += val1 * val2;
                Ks[I3][J3 + 2] += (-shape_info.Ni[row_n + i] * shape_info.Ndx[row + J2]) * val2;
                Ks[I3 + 1][J3 + 1] += val1 * val2;
                Ks[I3 + 1][J3 + 2] += (-shape_info.Ni[row_n + i] * shape_info.Ndx[row + J2 + 1]) * val2;
                Ks[I3 + 2][J3 + 2] += (shape_info.Ndx[row + I2]*shape_info.Ndx[row + J2] + 
                                       shape_info.Ndx[row + I2 + 1]*shape_info.Ndx[row + J2 + 1]) * val2;
                if(i != j){
                    Ks[I3 + 2][J3] += (-shape_info.Ni[row_n + j] * shape_info.Ndx[row + I2]) * val2;
                    Ks[I3 + 2][J3 + 1] += (-shape_info.Ni[row_n + j] * shape_info.Ndx[row + I2 + 1]) * val2;
                }
            }
        }
    }
    for(i = 0; i < 27; i++){
        for(j = 0; j < 27; j++){
            Ks[i][j] *= alpha;
        }
    }
};

void PQL_::get_K(double* K){
    double Kb[27][27]{0.}, Ks[27][27]{0.};
    this->get_Kb(Kb);
    this->get_Ks(Ks);

    int i, j, row;

    row = 0;
    for(i = 0; i < 27; i++){
        for(j = 0; j < 27; j++){
            K[row] = Kb[i][j] + Ks[i][j];
            row++;
        }
    }
};

void PQL_::asb_F_dload(double* Fo, double val, const std::vector<std::pair<int, int>> &sont){
    Lagrange_intp shape_info(2, 2);
    shape_info.set_integration(3);
    shape_info.getN__();
    shape_info.getNdx__(xyz);

    std::vector<int> noet, dofoet;
    int ntype_num = decode_sont(sont, noet, dofoet);

    int i = 0, k, Ni, row_n, rows;
    int where = find(node_id.front(), noet);
    double val_k;
    for(k = 0; k < shape_info.gintp.int_pt_num; k++){
        val_k = shape_info.gintp.weightness[k] * shape_info.Jacobien[k];
        row_n = 9 * k;

        for(i = 0; i < 9; i++){
            Ni = node_id.at(i);
            rows = 3 * (Ni - noet.at(where)) + dofoet.at(where) + 2;

            Fo[rows] += val_k * shape_info.Ni[row_n + i] * val;
        }
    }
};