#include "beam_TMSK2.h"

void B2TS::set_bnd(const beam_bnd &bd){
    bnd = bd;
}

void B2TS::set_node(std::vector<int> &nid, std::vector<double*> &x_c){
    node_id = nid;
    xyz = x_c;
};

void B2TS::get_KF(double K[4][4], double F[4]){
    double len = getlen(xyz[0], xyz[1]);
    double len2 = len * len;
    double bend_cof = mat->E * sec_prop->I / len;
    double sheer_cof =  mat->G * sec_prop->A / len / sec_prop->k;

    double Ks_val[4][4] = {{1,          0.5*len,    -1,         0.5*len},
                          {0.5*len,     len2/3,     -0.5*len,   len2/6},
                          {-1,          -0.5*len,   1,          -0.5*len},
                          {0.5*len,     len2/6,     -0.5*len,   len2/3}};

    double Kb_val[4][4] = {{0,  0 ,  0,  0},
                           {0,  1 ,  0,  -1},
                           {0,  0 ,  0,  0},
                           {0,  -1,  0,  1}};

    double Pq_val = bnd.q_const * len / 2;
    double Pq[4] = {Pq_val, 0, Pq_val, 0};

    for(std::vector<std::pair<double, double>>::iterator itp = bnd.p_j.begin(); itp != bnd.p_j.end(); itp++){
        Pq[0] += itp->first * (1 - itp->second) / 2;
        Pq[2] += itp->first * (1 + itp->second) / 2;
    }
    for(std::vector<std::pair<double, double>>::iterator itM = bnd.M_k.begin(); itM != bnd.M_k.end(); itM++){
        Pq[1] += itM->first * (1 - itM->second) / 2;
        Pq[3] += itM->first * (1 + itM->second) / 2;
    }

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            K[i][j] = Kb_val[i][j] * bend_cof + Ks_val[i][j] * sheer_cof;
        }
        F[i] = Pq[i];
    }
};

void B2TS::asb_KF(double* K, double* F, const std::vector<std::pair<int, int>> &sont){
    int ntype_num = sont.size();
    std::vector<int> noet(ntype_num + 1, 1);
    std::vector<int> dofoet(ntype_num + 1, 0);
    for(int i = 0; i < ntype_num; i++){
        noet.at(i + 1) = noet.at(i) + sont.at(i).first;
        dofoet.at(i + 1) = dofoet.at(i) + sont.at(i).second * sont.at(i).first;
    }

    double Ke[4][4]{0}, Pe[4]{0};
    this->get_KF(Ke, Pe);

    int Ni, Nj, wherei, wherej, row, col;
    int ROWI;
    for(int i = 0; i < 2; i++){
        Ni = node_id.at(i);
        wherei = find(Ni, noet);
        row = dofoet.at(wherei) + 2 * (Ni - noet.at(wherei));

        for(int j = 0; j < 2; j++){
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