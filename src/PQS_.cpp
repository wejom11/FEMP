#include "PQS_.h"

void PQS_::set(std::vector<int> nid, std::vector<double*> xc){
    if(nid.size() != 8 || xc.size() != 8){
        printf("ERROR: input vector's size must be 8!\n");
        return;
    }
    node_id = nid;
    xyz = xc;
};

void PQS_::get_Kb(double Kb[24][24]){
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

    Serendipity_intp shape_info(2, 2);
    shape_info.set_integration(inte_b);
    shape_info.getNdx__(xyz);

    double D0 = mat->E * pow(pp->thickness, 3) / 12 / (1 - pow(mat->nu, 2));
    double temp_1 = 0.5 * D0 * (1 - mat->nu);
    double temp_2 = D0 * mat->nu;
    double val1, val2, val3, val4, val5;

    int i, j, I3, J3, k, row,
        I2, J2;
    for(i = 0; i < 24; i++){
        for(j = 0; j < 24; j++){
            Kb[i][j] = 0.;
        }
    }
    
    for(k = 0; k < shape_info.gintp.int_pt_num; k++){
        row = k * 8 * 2;
        val5 = shape_info.gintp.weightness[k] * shape_info.Jacobien[k];
        for(i = 0; i < 8; i++){
            I3 = 3*i;
            I2 = 2*i;
            for(j = i; j < 8; j++){
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

void PQS_::get_Ks(double Ks[24][24]){
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

    Serendipity_intp shape_info(2, 2);
    shape_info.set_integration(inte_s);
    shape_info.getN__();
    shape_info.getNdx__(xyz);

    double alpha = mat->G * pp->thickness / 1.2;
    double val1, val2;

    int i, j, I3, J3, k, row,
        I2, J2, row_n;
    for(i = 0; i < 24; i++){
        for(j = 0; j < 24; j++){
            Ks[i][j] = 0.;
        }
    }
    for(k = 0; k < shape_info.gintp.int_pt_num; k++){
        row = k * 8 * 2;
        row_n = k * 8;

        val2 = shape_info.gintp.weightness[k] * shape_info.Jacobien[k];
        for(i = 0; i < 8; i++){
            I3 = 3*i;
            I2 = 2*i;
            for(j = i; j < 8; j++){
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
    for(i = 0; i < 24; i++){
        for(j = 0; j < 24; j++){
            Ks[i][j] *= alpha;
        }
    }
};

void PQS_::asb_KS(SparseMatrix& Kasb, const std::vector<std::pair<int, int>> &sont){
    int ntype_num = sont.size();
    std::vector<int> noet(ntype_num + 1, 1);
    std::vector<int> dofoet(ntype_num + 1, 0);
    for(int i = 0; i < ntype_num; i++){
        noet.at(i + 1) = noet.at(i) + sont.at(i).first;
        dofoet.at(i + 1) = dofoet.at(i) + sont.at(i).second * sont.at(i).first;
    }

    double Kb[24][24]{0.}, Ks[24][24]{0.};
    this->get_Kb(Kb);
    this->get_Ks(Ks);

    int i, j, row, col, I3, J3,
        m ,n, where, posi,
        rows, cols, I3pm, J3pm;
    where = find(node_id.front(), noet);
    for(i = 0; i < 8; i++){
        I3 = 3*i;
        for(j = i; j < 8; j++){
            J3 = 3*j;
            if(node_id.at(i) < node_id.at(j)){
                row = node_id.at(i);
                col = node_id.at(j);
                rows = dofoet.at(where) + 3 * (row - noet.at(where));
                cols = dofoet.at(where) + 3 * (col - noet.at(where));

                for(m = 0; m < 3; m++){
                    posi = match(rows + m + 1, cols + 1, Kasb);
                    I3pm = I3 + m;
                    for(n = 0; n < 3; n++){
                        Kasb.val[posi + n] += Kb[I3pm][J3 + n] + Ks[I3pm][J3 + n];
                    }
                }
            }
            else if(node_id.at(i) > node_id.at(j)){
                row = node_id.at(j);
                col = node_id.at(i);
                rows = dofoet.at(where) + 3 * (row - noet.at(where));
                cols = dofoet.at(where) + 3 * (col - noet.at(where));

                for(m = 0; m < 3; m++){
                    J3pm = J3 + m;
                    posi = match(rows + m + 1, cols + 1, Kasb);
                    for(n = 0; n < 3; n++){
                        Kasb.val[posi + n] += Kb[I3 + n][J3pm] + Ks[I3 + n][J3pm];
                    }
                }
            }
            else{
                row = node_id.at(i);
                rows = dofoet.at(where) + 3 * (row - noet.at(where));

                for(m = 0; m < 3; m++){
                    I3pm = I3 + m;
                    posi = match(rows + m + 1, rows + m + 1, Kasb);
                    for(n = m; n < 3; n++){
                        Kasb.val[posi + n - m] += Kb[I3pm][I3 + n] + Ks[I3pm][I3 + n];
                    }
                }
            }
        }
    }
};

void PQS_::asb_F_dload(double* Fo, double val, const std::vector<std::pair<int, int>> &sont){
    Serendipity_intp shape_info(2, 2);
    shape_info.set_integration(3);
    shape_info.getN__();
    shape_info.getNdx__(xyz);

    int ntype_num = sont.size();
    std::vector<int> noet(ntype_num + 1, 1);
    std::vector<int> dofoet(ntype_num + 1, 0);
    for(int i = 0; i < ntype_num; i++){
        noet.at(i + 1) = noet.at(i) + sont.at(i).first;
        dofoet.at(i + 1) = dofoet.at(i) + sont.at(i).second * sont.at(i).first;
    }

    int i = 0, k, Ni, row_n, rows;
    int where = find(node_id.front(), noet);
    double val_k;
    for(k = 0; k < shape_info.gintp.int_pt_num; k++){
        val_k = shape_info.gintp.weightness[k] * shape_info.Jacobien[k];
        row_n = 8 * k;

        for(i = 0; i < 8; i++){
            Ni = node_id.at(i);
            rows = 3 * (Ni - noet.at(where)) + dofoet.at(where) + 2;

            Fo[rows] += val_k * shape_info.Ni[row_n + i] * val;
        }
    }
};
