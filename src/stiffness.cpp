#include "stiffness.h"

void stiffness::init_KFD(elements& eles){
    const std::vector<int> s_o_et = eles.sot();
    const std::vector<std::pair<int, int>> s_o_nt = eles.sont();

    int sot_size = s_o_et.size();
    int sont_size = s_o_nt.size();
    int i = 0, j = 0;
    int dof = 0;
    for(i = 0; i < sont_size; i++){
        dof += s_o_nt.at(i).first * s_o_nt.at(i).second;
    }

    K_dense.resize(dof, dof);
    Fout = new double[dof]{0.};

    for(i = 0; i < sot_size; i++){
        if(s_o_et.at(i) == 0){
            continue;
        }

        if(i == 0){
            B2TS ele_empty;
            std::vector<B2TS>* eles_ref = eles.eleset_ptr(ele_empty);

            for(std::vector<B2TS>::iterator itele = eles_ref->begin(); itele != eles_ref->end(); itele++){
                itele->asb_KF(K_dense.val, Fout, s_o_nt);
            }
        }
        else if(i == 1){
            B3TS ele_empty;
            std::vector<B3TS>* eles_ref = eles.eleset_ptr(ele_empty);

            for(std::vector<B3TS>::iterator itele = eles_ref->begin(); itele != eles_ref->end(); itele++){
                itele->asb_KF(K_dense.val, Fout, s_o_nt);
            }
        }
        else{
            printf("no such element type %i\n", i);
        }
    }
};

void stiffness::init_KFS(elements& eles){
    this->init_K_symbolic(eles);

    const std::vector<std::pair<int,int>> s_o_nt = eles.sont();
    const std::vector<int> s_o_t = eles.sot();

    Fout = new double[K_sparse.rows]{0.};

    int i;
    for(i = 0; i < s_o_t.size(); i++){
        if(s_o_t.at(i) == 0){
            continue;
        }

        if(i == 0){
            printf("WARNING: Not yet functional\n");
            exit(2);
        }
        else if(i == 1){
            printf("WARNING: Not yet functional\n");
            exit(2);
        }
        else if(i == 2){
            PQL_ ele_empty;
            std::vector<PQL_>* eles_ref = eles.eleset_ptr(ele_empty);

            for(std::vector<PQL_>::iterator itele = eles_ref->begin(); itele != eles_ref->end(); itele++){
                itele->asb_KS(K_sparse, s_o_nt);
            }
        }
        else if(i == 3){
            PQS_ ele_empty;
            std::vector<PQS_>* eles_ref = eles.eleset_ptr(ele_empty);

            for(std::vector<PQS_>::iterator itele = eles_ref->begin(); itele != eles_ref->end(); itele++){
                itele->asb_KS(K_sparse, s_o_nt);
            }
        }
        else{
            printf("no such element type %i\n", i);
        }
    }

};

void stiffness::init_K_symbolic(elements& eles){
    const std::vector<std::pair<int,int>> s_o_nt = eles.sont();
    const std::vector<int> s_o_t = eles.sot();
    int ntype_num = s_o_nt.size();
    std::vector<int> noet(ntype_num + 1, 1);
    std::vector<int> dofoet(ntype_num + 1, 0);
    for(int i = 0; i < ntype_num; i++){
        noet.at(i + 1) = noet.at(i) + s_o_nt.at(i).first;
        dofoet.at(i + 1) = dofoet.at(i) + s_o_nt.at(i).second * s_o_nt.at(i).first;
    }
    int pt_num = noet.back() - 1;
    int dof = dofoet.back();

    bool* K_symbol = new bool[pt_num * pt_num]{false};
    int i = 0, j ,k, l, row;
    for(std::vector<int>::const_iterator itsot = s_o_t.begin(); itsot != s_o_t.end(); itsot++){
        if(*itsot){
            if(i == 0){
                printf("WARNING: Not yet functional\n");
                exit(2);
            }
            else if(i == 1){
                printf("WARNING: Not yet functional\n");
                exit(2);
            }
            else if(i == 2){
                PQL_ empty_ele;
                std::vector<PQL_>* PQL_eles = eles.eleset_ptr(empty_ele);

                for(std::vector<PQL_>::iterator itele = PQL_eles->begin(); itele != PQL_eles->end(); itele++){
                    for(j = 0; j < 9; j++){
                        row = pt_num * (itele->node_id.at(j) - 1);
                        for(k = 0; k < 9; k++){
                            K_symbol[row + itele->node_id.at(k) - 1] = true;
                        }
                    }
                }
            }
            else if(i == 3){
                PQS_ empty_ele;
                std::vector<PQS_>* PQS_eles = eles.eleset_ptr(empty_ele);

                for(std::vector<PQS_>::iterator itele = PQS_eles->begin(); itele != PQS_eles->end(); itele++){
                    for(j = 0; j < 8; j++){
                        row = pt_num * (itele->node_id.at(j) - 1);
                        for(k = 0; k < 8; k++){
                            K_symbol[row + itele->node_id.at(k) - 1] = true;
                        }
                    }
                }
            }
            else{
                printf("ERROR: \n");
            }
        }
        i++;
    }

    int* nzn__ = new int[ntype_num * ntype_num]{0};
    int* nznd__ = new int[ntype_num]{0};
    int row_k, fst_id;
    for(k = 0; k < ntype_num; k++){
        row_k = ntype_num * k;
        for(i = noet[k] - 1; i < noet[k+1] - 1; i++){
            row = i * pt_num;
            if(K_symbol[row + i]){
                nznd__[k]++;
            }
    
            for(l = k; l < ntype_num; l++){
                fst_id = i+1 <= noet[l]-1 ? noet[l]-1 : i+1;
                for(j = fst_id; j < noet[l+1] - 1; j++){
                    if(K_symbol[row + j]){
                        nzn__[row_k + l] ++;
                    }
                }
            }
        }
    }

    int nzn = 0;
    for(i = 0; i < ntype_num; i++){
        row = i * ntype_num;
        l = s_o_nt.at(i).second;
        nzn += nznd__[i] * (l*l + l)/2;
        for(j = i; j < ntype_num; j++){
            nzn += nzn__[row + j] * l * s_o_nt.at(j).second;
        }
    }

    delete[] nzn__; nzn__ = nullptr;
    delete[] nznd__; nznd__ = nullptr;

    try{
        K_sparse.val = new double[nzn]{0};
        K_sparse.row_st = new int[dof + 1]{0};
        K_sparse.col = new int[nzn]{0};
        K_sparse.type = 2;
        K_sparse.rows = dof;
    }
    catch(const std::bad_alloc& e){
        std::cerr << e.what() << '\n';
        exit(1);
    }
    
    int icol = 0, IN, JN, m, n;
    K_sparse.row_st[0] = 1;
    for(k = 0; k < ntype_num; k++){
        for(i = noet[k]; i < noet[k+1]; i++){
            row = (i - 1) * pt_num;
            IN = s_o_nt.at(k).second * (i - noet[k]) + dofoet.at(k);

            for(l = k; l < ntype_num; l++){
                fst_id = i <= noet[l] ? noet[l] : i;
                for(j = fst_id; j < noet[l+1]; j++){
                    if(K_symbol[row + j - 1]){
                        JN = s_o_nt.at(l).second * (j - noet[l]) + dofoet.at(l);

                        for(m = 1; m < s_o_nt.at(l).second + 1; m++){
                            K_sparse.col[icol] = JN + m;
                            icol++;
                        }
                    }
                }
            }
            K_sparse.row_st[IN + 1] = icol + 1;

            for(m = 1; m < s_o_nt.at(k).second; m++){
                for(n = K_sparse.row_st[IN] + m - 1; n < K_sparse.row_st[IN + 1] - 1; n++){
                    K_sparse.col[icol] = K_sparse.col[n];
                    icol++;
                }
                K_sparse.row_st[IN + m + 1] = icol + 1;
            }
        }
    }

    delete[] K_symbol; K_symbol = nullptr;
    
};
