#include "mesh.h"

void Linemesh_demo::generate(elements& eles, std::vector<material> &matlib, std::vector<section_TSbeam> &seclib, 
                             boundaries &dbnd, double* &xyz){
    if(type == 0 || type == 1){
        section_TSbeam sec_1;
        material mat_1;

        matlib.resize(1); matlib.front() = mat_1;
        seclib.resize(1); seclib.front() = sec_1;

        xyz = new double[3 * (ele_num + 1)];
        double dl = length / ele_num;
        for(int i = 0; i < ele_num + 1; i++){
            xyz[3*i] = i * dl;
            xyz[3*i + 1] = 0.;
            xyz[3*i + 2] = 0.;
        }
        eles.set_sont({{ele_num+1, 2}});

        beam_bnd bd;
        bd.p_j = {{is_M ? 0 : f_val, 1.0}};
        bd.M_k = {{is_M ? f_val : 0, 1.0}};
        dbnd.dbd.Anchorage_nid = {1};

        eles.resize({ele_num, 0, 0, 0});
        B2TS ele_empty;
        std::vector<B2TS>* b2ts_set = eles.eleset_ptr(ele_empty);

        if(type == 0){
            B2TS b2ts_ele(0,&matlib.front(), &seclib.front());
            for(int i = 0; i < ele_num; i++){
                std::vector<int> nnid = {i+1, i+2};
                std::vector<double*> xx_c = {xyz+3*i, xyz+3*i+3};
                b2ts_ele.set_node(nnid, xx_c);

                if(i == ele_num - 1){
                    b2ts_ele.set_bnd(bd);
                }
                
                b2ts_set->at(i) = b2ts_ele;
            }
        }
        else{
            B2TS b2ts_ele(1,&matlib.front(), &seclib.front());
            for(int i = 0; i < ele_num; i++){
                std::vector<int> nnid = {i+1, i+2};
                std::vector<double*> xx_c = {xyz+3*i, xyz+3*i+3};
                b2ts_ele.set_node(nnid, xx_c);

                if(i == ele_num - 1){
                    b2ts_ele.set_bnd(bd);
                }
                
                b2ts_set->at(i) = b2ts_ele;
            }
        }
    }
    else if(type == 2 || type == 3){
        section_TSbeam sec_1;
        material mat_1;

        matlib.resize(1); matlib.front() = mat_1;
        seclib.resize(1); seclib.front() = sec_1;

        xyz = new double[3 * (2*ele_num + 1)];
        double dl = length / ele_num / 2;
        for(int i = 0; i < 2*ele_num + 1; i++){
            xyz[3*i] = i * dl;
            xyz[3*i + 1] = 0.;
            xyz[3*i + 2] = 0.;
        }
        eles.set_sont({{2*ele_num+1, 2}});

        beam_bnd bd;
        bd.p_j = {{is_M ? 0 : f_val, 1.0}};
        bd.M_k = {{is_M ? f_val : 0, 1.0}};
        dbnd.dbd.Anchorage_nid = {1};

        eles.resize({0,ele_num,0,0});
        B3TS ele_empty;
        std::vector<B3TS>* b3ts_set = eles.eleset_ptr(ele_empty);
        if(type == 2){
            B3TS b3ts_ele(0, &matlib.front(), &seclib.front());
            for(int i = 0; i < ele_num; i++){
                std::vector<int> nnid = {2*i+1, 2*i+2, 2*i+3};
                std::vector<double*> xx_c = {xyz+6*i, xyz+6*i+3, xyz+6*i+6};
                b3ts_ele.set_node(nnid, xx_c);

                if(i == ele_num - 1){
                    b3ts_ele.set_bnd(bd);
                }
                
                b3ts_set->at(i) = b3ts_ele;
            }
        }
        else{
            B3TS b3ts_ele(1, &matlib.front(), &seclib.front());
            for(int i = 0; i < ele_num; i++){
                std::vector<int> nnid = {2*i+1, 2*i+2, 2*i+3};
                std::vector<double*> xx_c = {xyz+6*i, xyz+6*i+3, xyz+6*i+6};
                b3ts_ele.set_node(nnid, xx_c);

                if(i == ele_num - 1){
                    b3ts_ele.set_bnd(bd);
                }
                
                b3ts_set->at(i) = b3ts_ele;
            }
        }
    }
}

void Square_Platemesh_demo::generate(elements &eles, std::vector<material> &matlib, std::vector<plate_prop> &pplib, 
                                     boundaries &dbnd, double* &xyz){
    if(type == 0 || type == 1 || type == 2){
        material mat_1;
        plate_prop pp_1;

        matlib.resize(1);
        matlib.front() = mat_1;
        pplib.resize(1);
        pplib.front() = pp_1;

        double l_step = length / line_ele_num / 2;
        int line_node_num = 2 * line_ele_num + 1;
        dbnd.dbd.Anchorage_nid.resize(4 * line_node_num - 4);
        xyz = new double[3 * line_node_num * line_node_num]{0.};
        int i, j, indice, k, indice_a;
        indice = 0;
        indice_a = 0;
        for(i = 0; i < line_node_num; i++){
            for(j = 0; j < line_node_num; j++){
                xyz[indice] = j * l_step;
                xyz[indice + 1] = i * l_step;
                xyz[indice + 2] = 0;
                indice += 3;

                if(i == 0 || i == line_node_num - 1 || j == 0 || j == line_node_num - 1){
                    dbnd.dbd.Anchorage_nid.at(indice_a) = i*line_node_num + j + 1;
                    indice_a++;
                }
            }

        }

        eles.resize({0,0,line_ele_num * line_ele_num,0});
        eles.set_sont({{line_node_num * line_node_num, 3}});
        PQL_ ele_ins(type, &matlib.front(), &pplib.front());
        std::vector<PQL_>* eles_ref = eles.eleset_ptr(ele_ins);

        std::vector<int> n_id(9);
        std::vector<double*> xc(9, nullptr);
        indice = 0;
        int lnn2 = line_node_num * 2;
        for(i = 0; i < line_ele_num; i++){
            for(j = 0; j < line_ele_num; j++){
                n_id = {lnn2 * i + 2*j + 1, lnn2 * i + 2*j + 3, lnn2 * (i + 1) + 2*j + 3, lnn2 * (i + 1) + 2*j + 1,
                        lnn2 * i + 2 * (j+1), lnn2 * i + 2*j + 3 + line_node_num, lnn2 * (i+1) + 2 * (j+1), lnn2 * i + 2*j + 1 + line_node_num,
                        lnn2 * i + 2*j + 2 + line_node_num};
                for(k = 0; k < 9; k++){
                    xc.at(k) = xyz + 3*(n_id.at(k) - 1);
                }

                ele_ins.set(n_id, xc);

                eles_ref->at(indice) = ele_ins;
                indice++;
            }
        }
        dbnd.pbd.pbd_d.resize(1);
        dbnd.pbd.pbd_d.front().ele_set.resize(line_ele_num * line_ele_num + 1);
        dbnd.pbd.pbd_d.front().ele_set.front() = 2;
        dbnd.pbd.pbd_d.front().val = p_val;
        for(i = 0; i < line_ele_num * line_ele_num; i++){
            dbnd.pbd.pbd_d.front().ele_set.at(i+1) = i+1;
        }
    }
    else if(type == 3 || type == 4 || type == 5){
        material mat_1;
        plate_prop pp_1;

        matlib.resize(1);
        matlib.front() = mat_1;
        pplib.resize(1);
        pplib.front() = pp_1;

        double lx_step = 0,
               ly_step = length / line_ele_num / 2;
        int line_node_num = 2 * line_ele_num + 1;
        int pts_num = line_node_num * line_node_num - line_ele_num * line_ele_num;
        dbnd.dbd.Anchorage_nid.resize(4 * line_node_num - 4);
        xyz = new double[3 * pts_num]{0.};
        int i, j, indice, k, indice_a;
        indice = 0;
        indice_a = 0;
        for(i = 0; i < line_node_num; i++){
            lx_step = length / line_ele_num / (2 - i%2);
            for(j = 0; j < (2 - i%2) * line_ele_num + 1; j++){
                xyz[indice] = j * lx_step;
                xyz[indice + 1] = i * ly_step;
                xyz[indice + 2] = 0;
                indice += 3;

                if(i == 0 || i == line_node_num - 1 || j == 0 || 
                   j == (2 - i%2) * line_ele_num){
                    dbnd.dbd.Anchorage_nid.at(indice_a) = indice / 3;
                    indice_a++;
                }
            }

        }

        eles.resize({0,0,0,line_ele_num * line_ele_num});
        eles.set_sont({{pts_num, 3}});
        PQS_ ele_ins(type-3, &matlib.front(), &pplib.front());
        std::vector<PQS_>* eles_ref = eles.eleset_ptr(ele_ins);

        std::vector<int> n_id(8);
        std::vector<double*> xc(8, nullptr);
        indice = 0;
        int lnn2 = line_node_num * 2;
        for(i = 0; i < line_ele_num; i++){
            for(j = 0; j < line_ele_num; j++){
                n_id = {(3*line_ele_num + 2) * i + 2*j + 1, (3*line_ele_num + 2) * i + 2*j + 3, 
                        (3*line_ele_num + 2) * (i + 1) + 2*j + 3, (3*line_ele_num + 2) * (i + 1) + 2*j + 1,
                        (3*line_ele_num + 2) * i + 2 * (j+1), (3*line_ele_num + 2) * i + j + 2 + line_node_num, 
                        (3*line_ele_num + 2) * (i+1) + 2 * (j+1), (3*line_ele_num + 2) * i + j + 1 + line_node_num};
                for(k = 0; k < 8; k++){
                    xc.at(k) = xyz + 3*(n_id.at(k) - 1);
                }

                ele_ins.set(n_id, xc);

                eles_ref->at(indice) = ele_ins;
                indice++;
            }
        }
        dbnd.pbd.pbd_d.resize(1);
        dbnd.pbd.pbd_d.front().ele_set.resize(line_ele_num * line_ele_num + 1);
        dbnd.pbd.pbd_d.front().ele_set.front() = 3;
        dbnd.pbd.pbd_d.front().val = p_val;
        for(i = 0; i < line_ele_num * line_ele_num; i++){
            dbnd.pbd.pbd_d.front().ele_set.at(i+1) = i+1;
        }
    }
}
