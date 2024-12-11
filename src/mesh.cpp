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

void Trussmesh_demo::generate(elements &eles, std::vector<material> &matlib, std::vector<plate_prop> &pplib, 
                                     boundaries &dbnd, solution& sln, double* &xyz, short ad){
    material mat_1;
    double t = sec_a;
    plate_prop pp_1(t);

    matlib.resize(2);
    matlib.front() = mat_1;
    matlib.at(1) = material(2.1E11);
    pplib.resize(1);
    pplib.front() = pp_1;

    int nnum = type ? 9 : 4;
    int h_len = line_ele_num * round(length / sec_a);
    //int h_len = 10;
    double wide_step = sec_a / line_ele_num / (type + 1);
    double height_step = length / h_len / (type + 1);
    int line_node_num = (type + 1) * line_ele_num + 1;
    int h_lnn = (type + 1) * h_len + 1;

    int sidebnd1 = bnd_type % 10, sidebnd2 = (bnd_type / 10) % 10;
    if (type == 0 && (sidebnd1 == 1 || sidebnd2 == 1)) {
        if (line_ele_num % 2 != 0) {
            printf("\033[1;31mERROR\033[0m: please use even number element in width for element PSL4\n");
            exit(0);
        }
    }
    int anum1 = sidebnd1 == 2 ? line_node_num : sidebnd1,
        anum2 = sidebnd2 == 2 ? (2 + type) * line_node_num : sidebnd2;
    dbnd.dbd.Anchorage_nid.resize(anum1 + anum2);
    dbnd.dbd.fixed_dof.resize(anum1 + anum2,{});
    int dof = ad * line_node_num * (h_lnn + (sidebnd2 == 2 ? (type + 1) : 0));
    xyz = new double[dof]{0.};
    sln.Var = new double[dof]{0.};

    int i, j, indice, k, indice_a;
    indice = 0;
    indice_a = 0;
    for(i = 0; i < (h_lnn + (sidebnd2 == 2 ? (type + 1) : 0)); i++){
        for(j = 0; j < line_node_num; j++){
            if ( i >= 1 && i < h_lnn) {
                xyz[indice] = j * wide_step + dis_mag * disp_disturb(bnd_type, dis_order, i * height_step, length);
                xyz[indice + 1] = i * height_step;
            }
            else {
                xyz[indice] = j * wide_step;
                xyz[indice + 1] = i * height_step;
            }
            indice += 2;

            if(i == 0 && sidebnd1 == 2){
                dbnd.dbd.Anchorage_nid.at(indice_a) =  j + 1;
                indice_a++;
            }
            if(i >= h_lnn - 1 && sidebnd2 == 2){
                dbnd.dbd.Anchorage_nid.at(indice_a) = i*line_node_num + j + 1;
                dbnd.dbd.fixed_dof.at(indice_a) = {1, 0};
                indice_a++;
            }
        }
    }

    if(sidebnd1 == 1){
        dbnd.dbd.Anchorage_nid.at(indice_a) =  line_ele_num / (2 - type) + 1;
        indice_a++;
    }
    if(sidebnd2 == 1){
        dbnd.dbd.Anchorage_nid.at(indice_a) = (h_lnn-1)*line_node_num + line_ele_num / (2 - type) + 1;
        dbnd.dbd.fixed_dof.at(indice_a) = {1, 0};
        indice_a++;
    }

    eles.resize({ 0,0,0,0,line_ele_num * h_len + (sidebnd2 == 2 ? line_ele_num : 0) });
    eles.set_sont({{line_node_num * (h_lnn + (sidebnd2 == 2 ? (type + 1) : 0)), 2}});
    PSL_ ele_ins(&matlib.front(), &pplib.front(), type);
    std::vector<PSL_>* eles_ref = eles.eleset_ptr(ele_ins);

    std::vector<int> n_id(nnum);
    indice = 0;
    int lnn2 = line_node_num * (type + 1);
    for(i = 0; i < h_len + (sidebnd2 == 2 ? 1 : 0); i++){
        for(j = 0; j < line_ele_num; j++){
            if(type == 1){
                n_id = {lnn2 * i + 2*j + 1, lnn2 * i + 2*j + 3, lnn2 * (i + 1) + 2*j + 3, lnn2 * (i + 1) + 2*j + 1,
                        lnn2 * i + 2 * (j+1), lnn2 * i + 2*j + 3 + line_node_num, lnn2 * (i+1) + 2 * (j+1), lnn2 * i + 2*j + 1 + line_node_num,
                        lnn2 * i + 2*j + 2 + line_node_num};
                //n_id = {1,7,9,3,4,8,6,2,5};
            }
            else{
                n_id = {line_node_num * i + j + 1, line_node_num * i + j + 2, 
                        line_node_num * (i+1) + j + 2, line_node_num * (i+1) + j + 1};
            }

            ele_ins.set_xy_uv(xyz, sln.Var, n_id);
            if(i == h_len){
                ele_ins.mat = &matlib.at(1);
            }

            eles_ref->at(indice) = ele_ins;
            //eles_ref->at(indice).set_gauss_order(3);
            eles_ref->at(indice).ele_init();
            indice++;
        }
    }
    
    double Ix = pow(sec_a, 2) * t / 12;
    dbnd.dbd2ds.resize(1);
    dbnd.dbd2ds.front().elset_line.resize(line_ele_num);
    dbnd.dbd2ds.front().norm_vec = {0., -1.};
    dbnd.dbd2ds.front().type = type;
    //dbnd.dbd2ds.front().val = p_factor * mat_1.E * Ix / pow(length, 2);
    dbnd.dbd2ds.front().val = 1 * t;
    for(i = 0; i < line_ele_num; i++){
        dbnd.dbd2ds.front().elset_line.at(i) = {line_ele_num * (h_len - 1 + (sidebnd2 == 2 ? 1 : 0)) + i + 1, 3};
    }

    // dbnd.cbd2ds.resize(1);
    // dbnd.cbd2ds.front().val = {0, -p_factor * mat_1.E * Ix / pow(length, 2)};
    // dbnd.cbd2ds.front().nset = {line_node_num * (h_lnn - 1) + line_ele_num / (2 - type) + 1};
}


double disp_disturb(short bnd_type, short order, double y, double l){
    double x_d = 0.;
    double k = 0;
    double pi = 3.14159265358979323846;
    double phi = y / l;
    
    if(bnd_type == 02){
        k = (2*order - 1) * pi / 2.0;
        x_d = 1 - cos(k * phi);
        return x_d;
    }
    else if(bnd_type == 11){
        k = order * pi;
        x_d = sin(k * phi);
        return x_d;
    }
    else if(bnd_type == 12){
        if(order == 1){
            k = 4.4934094579090641753078809272803;
        }
        else if(order == 2){
            k = 7.725251836937707164195068933063;
        }
        else if(order == 3){
            k = 10.904121659428899827148702790189;
        }
        else if(order == 4){
            k = 14.066193912831473479978965600602;
        }
        x_d = 1 - phi - cos(k*phi) + sin(k*phi) / k;
        return x_d;
    }
    else if(bnd_type == 22){
        k = 2 * order * pi;
        x_d = 1 - cos(k * phi);
        return x_d;
    }
    else{
        printf("Waring: such boundary type %i havn't supported yet![No disturbation will be added]\n", bnd_type);
        return 0;
    }
};