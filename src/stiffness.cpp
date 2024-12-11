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
    if(!K_sparse.val){
        this->init_K_symbolic(eles);
        Fout = new double[K_sparse.rows]{0.};
    }
    else{
        for(int i = 0; i < K_sparse.row_st[K_sparse.rows] - 1; i++){
            K_sparse.val[i] = 0.;
        }
        
        for(int i = 0; i < K_sparse.rows; i++){
            Fout[i] = 0.;
        }
    }  

    const std::vector<std::pair<int,int>> s_o_nt = eles.sont();
    const std::vector<int> s_o_t = eles.sot();
    int i;
    for(i = 0; i < s_o_t.size(); i++){
        if(s_o_t.at(i) == 0){
            continue;
        }

        if(i == 0){
            printf("\033[33mWARNING\033[0m: Not yet functional\n");
            exit(2);
        }
        else if(i == 1){
            printf("\033[33mWARNING\033[0m: Not yet functional\n");
            exit(2);
        }
        else if(i == 2){
            PQL_ ele_empty;
            std::vector<PQL_>* eles_ref = eles.eleset_ptr(ele_empty);
            double* Ke = new double[27 * 27]{0.};

            for(std::vector<PQL_>::iterator itele = eles_ref->begin(); itele != eles_ref->end(); itele++){
                itele->get_K(Ke);
                asb_spa_sym(Ke, nullptr, itele->node_id, s_o_nt);
            }
            delete[] Ke; Ke = nullptr;
        }
        else if(i == 3){
            PQS_ ele_empty;
            std::vector<PQS_>* eles_ref = eles.eleset_ptr(ele_empty);
            double* Ke = new double[24 * 24]{0.};

            for(std::vector<PQS_>::iterator itele = eles_ref->begin(); itele != eles_ref->end(); itele++){
                itele->get_K(Ke);
                asb_spa_sym(Ke, nullptr, itele->node_id, s_o_nt);
            }

            delete[] Ke; Ke = nullptr;
        }
        else if(i == 4){
            printf("\033[33mWARNING\033[0m: Not yet functional\n");
            exit(2);
        }
        else{
            printf("no such element type %i\n", i);
        }
    }
};

void stiffness::init_KFSN(elements& eles, std::vector<short> option){
    if(!K_sparse.val){
        this->init_K_symbolic(eles);
        Fout = new double[K_sparse.rows]{0.};
        Fint = new double[K_sparse.rows]{0.};
    }
    else{
        for(int i = 0; i < K_sparse.row_st[K_sparse.rows] - 1; i++){
            K_sparse.val[i] = 0.;
        }
        
        for(int i = 0; i < K_sparse.rows; i++){
            Fint[i] = 0.;
            Fout[i] = 0.;
        }
    }  

    const std::vector<std::pair<int,int>> s_o_nt = eles.sont();
    const std::vector<int> s_o_t = eles.sot();
    int i;
    for(i = 0; i < s_o_t.size(); i++){
        if(s_o_t.at(i) == 0){
            continue;
        }

        if(i == 4){
            PSL_ ele_empty;
            std::vector<PSL_>* eles_ref = eles.eleset_ptr(ele_empty);
            int eledof = eles_ref->front().nnum * 2;
            DenseMatrix Ke(eledof, eledof);
            double* Fe = new double[eledof]{0.};

            for(std::vector<PSL_>::iterator itele = eles_ref->begin(); itele != eles_ref->end(); itele++){
                itele->get_KF(Ke, Fe, option.front());
                asb_spa_sym(Ke.val, Fe, itele->node_id, s_o_nt);
            }

            Ke.dealloc();
            delete[] Fe; Fe = nullptr;
        }
        else if(i == 1){
            printf("\033[33mWARNING\033[0m: Not yet functional\n");
            exit(2);
        }
        else if(i == 2){
            printf("\033[33mWARNING\033[0m: Not yet functional\n");
            exit(2);
        }
        else if(i == 3){
            printf("\033[33mWARNING\033[0m: Not yet functional\n");
            exit(2);
        }
        else if(i == 0){
            printf("\033[33mWARNING\033[0m: Not yet functional\n");
            exit(2);            
        }
        else{
            printf("no such element type %i\n", i);
        }
    }
};

void stiffness::init_K_symbolic(elements& eles){
    const std::vector<std::pair<int,int>> s_o_nt = eles.sont();
    const std::vector<int> s_o_t = eles.sot();
    std::vector<int> noet, dofoet;
    int ntype_num = decode_sont(s_o_nt, noet, dofoet);
    int pt_num = noet.back() - 1;
    int dof = dofoet.back();

    bool* K_symbol = new bool[pt_num * pt_num]{false};
    int i = 0, j ,k, l, row;
    for(std::vector<int>::const_iterator itsot = s_o_t.begin(); itsot != s_o_t.end(); itsot++){
        if(*itsot){
            if(i == 0){
                printf("\033[33mWARNING\033[0m: Not yet functional\n");
                exit(2);
            }
            else if(i == 1){
                printf("\033[33mWARNING\033[0m: Not yet functional\n");
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
            else if(i == 4){
                PSL_ empty_ele;
                std::vector<PSL_>* PSL_eles = eles.eleset_ptr(empty_ele);

                for(std::vector<PSL_>::iterator itele = PSL_eles->begin(); itele != PSL_eles->end(); itele++){
                    for(j = 0; j < itele->nnum; j++){
                        row = pt_num * (itele->node_id.at(j) - 1);
                        for(k = 0; k < itele->nnum; k++){
                            K_symbol[row + itele->node_id.at(k) - 1] = true;
                        }
                    }
                }                
            }
            else{
                printf("\033[1;31mERROR\033[0m: \n");
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

void stiffness::asb_spa_sym(double* K, double* F, std::vector<int>& node_id,
                            const std::vector<std::pair<int,int>>& sont){
    std::vector<int> noet, dofoet;
    int ntype_num = decode_sont(sont, noet, dofoet);

    int i, j, row, col, ISD, JSD,
        m ,n, where, posi, sgdof,
        rows, cols, ISDpm, JSDpm;
    int dim = node_id.size();

    where = find(node_id.front(), noet);
    sgdof = sont.at(where).second;
    int dof = dim * sgdof;
    if(K){
        for(i = 0; i < dim; i++){
            ISD = sgdof*i;
            for(j = i; j < dim; j++){
                JSD = sgdof*j;
                if(node_id.at(i) < node_id.at(j)){
                    row = node_id.at(i);
                    col = node_id.at(j);
                    rows = dofoet.at(where) + sgdof * (row - noet.at(where));
                    cols = dofoet.at(where) + sgdof * (col - noet.at(where));

                    for(m = 0; m < sgdof; m++){
                        posi = match(rows + m + 1, cols + 1, K_sparse);
                        ISDpm = ISD + m;
                        for(n = 0; n < sgdof; n++){
                            K_sparse.val[posi + n] += K[ISDpm*dof + JSD + n];
                        }
                    }
                }
                else if(node_id.at(i) > node_id.at(j)){
                    row = node_id.at(j);
                    col = node_id.at(i);
                    rows = dofoet.at(where) + sgdof * (row - noet.at(where));
                    cols = dofoet.at(where) + sgdof * (col - noet.at(where));

                    for(m = 0; m < sgdof; m++){
                        JSDpm = JSD + m;
                        posi = match(rows + m + 1, cols + 1, K_sparse);
                        for(n = 0; n < sgdof; n++){
                            K_sparse.val[posi + n] += K[(ISD + n)*dof + JSDpm];
                        }
                    }
                }
                else{
                    row = node_id.at(i);
                    rows = dofoet.at(where) + sgdof * (row - noet.at(where));

                    for(m = 0; m < sgdof; m++){
                        ISDpm = ISD + m;
                        posi = match(rows + m + 1, rows + m + 1, K_sparse);
                        for(n = m; n < sgdof; n++){
                            K_sparse.val[posi + n - m] += K[ISDpm*dof + ISD + n];
                        }
                    }
                }
            }
        }
    }

    if(F){
        if(!Fint) printf("\033[1;31mERROR\033[0m: Fint haven't initialized!\n");
        for(i = 0; i < dim; i++){
            ISD = sgdof * i;
            row = node_id.at(i);
            rows = dofoet.at(where) + sgdof * (row - noet.at(where));

            for(j = 0; j < sgdof; j++){
                Fint[rows + j] += F[ISD + j];
            }
        }
    }
};
