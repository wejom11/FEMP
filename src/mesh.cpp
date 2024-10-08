#include "mesh.h"

void Linemesh_demo::generate(elements& eles, std::vector<material> &matlib, std::vector<section_TSbeam> &seclib, 
                             disp_bnd &dbnd, double* xyz){
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
        dbnd.Anchorage_nid = {1};

        if(type == 0){
            eles.resize({ele_num,0,0,0});
            B2TS ele_empty;
            std::vector<B2TS>* b2ts_set = eles.eleset_ptr(ele_empty);

            B2TS b2ts_ele(&matlib.front(), &seclib.front());
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
            eles.resize({0, ele_num,0,0});
            B2TSR ele_empty;
            std::vector<B2TSR>* b2tsr_set = eles.eleset_ptr(ele_empty);

            B2TSR b2tsr_ele(&matlib.front(), &seclib.front());
            for(int i = 0; i < ele_num; i++){
                std::vector<int> nnid = {i+1, i+2};
                std::vector<double*> xx_c = {xyz+3*i, xyz+3*i+3};
                b2tsr_ele.set_node(nnid, xx_c);

                if(i == ele_num - 1){
                    b2tsr_ele.set_bnd(bd);
                }
                
                b2tsr_set->at(i) = b2tsr_ele;
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
        dbnd.Anchorage_nid = {1};

        if(type == 2){
            eles.resize({0,0,ele_num,0});
            B3TS ele_empty;
            std::vector<B3TS>* b3ts_set = eles.eleset_ptr(ele_empty);

            B3TS b3ts_ele(&matlib.front(), &seclib.front());
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
            eles.resize({0,0,0,ele_num});
            B3TSR ele_empty;
            std::vector<B3TSR>* b3tsr_set = eles.eleset_ptr(ele_empty);

            B3TSR b3tsr_ele(&matlib.front(), &seclib.front());
            for(int i = 0; i < ele_num; i++){
                std::vector<int> nnid = {2*i+1, 2*i+2, 2*i+3};
                std::vector<double*> xx_c = {xyz+6*i, xyz+6*i+3, xyz+6*i+6};
                b3tsr_ele.set_node(nnid, xx_c);

                if(i == ele_num - 1){
                    b3tsr_ele.set_bnd(bd);
                }
                
                b3tsr_set->at(i) = b3tsr_ele;
            }
        }
    }
}
