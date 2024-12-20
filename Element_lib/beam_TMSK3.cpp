#include "beam_TMSK3.h"

void B3TS::set_bnd(const beam_bnd &bd){
    bnd = bd;
}

void B3TS::set_node(std::vector<int> &nid, std::vector<double*> &x_c){
    node_id = nid;
    xyz = x_c;
};

void B3TS::get_KF(double K[6][6], double F[6]){
    double len = getlen(xyz[0], xyz[2], 3);
    double len2 = len * len;
    double bend_cof = mat->E * sec_prop->I / len / 3;
    double sheer_cof =  mat->G * sec_prop->A / len / sec_prop->k / 30;

    double Ksf_val[6][6] = {{70,     15*len,   -80,      20*len,   10,       -5*len},
                            {15*len, 4*len2,   -20*len,  2*len2,   5*len,    -len2},
                            {-80,    -20*len,  160,      0,        -80,      20*len},
                            {20*len, 2*len2,   0,        16*len2,  -20*len,  2*len2},
                            {10,     5*len,    -80,      -20*len,  70,       -15*len},
                            {-5*len, -len2,    20*len,   2*len2,   -15*len,  4*len2}};
    
    double Ksr_val[6][6] = {{70,     15*len,     -80,      20*len,     10,       -5*len},
                           {15*len, 10*len2/3,  -20*len,  10*len2/3,  5*len,    -5*len2/3},
                           {-80,    -20*len,    160,      0,          -80,      20*len},
                           {20*len, 10*len2/3,  0,        40*len2/3,  -20*len,  10*len2/3},
                           {10,     5*len,      -80,      -20*len,    70,       -15*len},
                           {-5*len, -5*len2/3,  20*len,   10*len2/3,  -15*len,  10*len2/3}};

    double Kb_val[6][6] = {{0,  0 ,  0,  0,  0,  0},
                           {0,  7 ,  0,  -8, 0,  1},
                           {0,  0 ,  0,  0,  0,  0},
                           {0,  -8,  0,  16, 0,  -8},
                           {0,  0,   0,  0,  0,  0},
                           {0,  1,   0,  -8, 0,  7}};

    double Pq_val = bnd.q_const * len / 6;
    double Pq[6] = {Pq_val, 0, 4*Pq_val, 0, Pq_val, 0};

    for(std::vector<std::pair<double, double>>::iterator itp = bnd.p_j.begin(); itp != bnd.p_j.end(); itp++){
        double xi = itp->second;
        Pq[0] += itp->first * (xi - 1) * xi / 2;
        Pq[2] += itp->first * (1 - xi*xi);
        Pq[4] += itp->first * (1 + xi) * xi / 2;
    }
    for(std::vector<std::pair<double, double>>::iterator itM = bnd.M_k.begin(); itM != bnd.M_k.end(); itM++){
        double xi = itM->second;
        Pq[1] += itM->first * (xi - 1) * xi / 2;
        Pq[3] += itM->first * (1 - xi*xi) / 2;
        Pq[5] += itM->first * (xi + 1) * xi / 2;
    }

    if(int_method == 0){
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                K[i][j] = Kb_val[i][j] * bend_cof + Ksf_val[i][j] * sheer_cof;
            }
            F[i] = Pq[i];
        }
    }
    else{
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                K[i][j] = Kb_val[i][j] * bend_cof + Ksr_val[i][j] * sheer_cof;
            }
            F[i] = Pq[i];
        }
    }

};

void B3TS::asb_KF(double* K, double* F, const std::vector<std::pair<int, int>> &sont){
    std::vector<int> noet, dofoet;
    int ntype_num = decode_sont(sont, noet, dofoet);

    double Ke[6][6]{0}, Pe[6]{0};
    this->get_KF(Ke, Pe);

    int Ni, Nj, wherei, wherej, row, col;
    int ROWI;
    for(int i = 0; i < 3; i++){
        Ni = node_id.at(i);
        wherei = find(Ni, noet);
        row = dofoet.at(wherei) + 2 * (Ni - noet.at(wherei));

        for(int j = 0; j < 3; j++){
            Nj = node_id.at(j);
            wherej = find(Nj, noet);
            col = dofoet.at(wherej) + 2 * (Nj - noet.at(wherej));

            for(int m = 0; m < 2; m++){
                ROWI = (row + m) * dofoet.back();
                for(int n = 0; n < 2; n++){
                    K[ROWI + col + n] += Ke[2*i+m][2*j+n];
                }
            }
        }

        F[row] += Pe[2*i];
        F[row + 1] += Pe[2*i + 1];
    }
};