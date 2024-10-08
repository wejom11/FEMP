#include "assemble.h"

void asb_manager::init_mesh(int num, bool M, int type, double len, double val){
    Linemesh_demo lmd(num, len, val, M, type);
    lmd.generate(eles, mater_lib, section_lib, dbd, xyz_coord);
}

void asb_manager::init_KF(){
    const std::vector<int> s_o_et = eles.sot();
    const std::vector<std::pair<int, int>> s_o_nt = eles.sont();

    int sot_size = s_o_et.size();
    int sont_size = s_o_nt.size();
    int i = 0, j = 0;
    int dof = 0;
    for(i = 0; i < sont_size; i++){
        dof += s_o_nt.at(i).first * s_o_nt.at(i).second;
    }

    K.resize(dof, dof);
    Fout = new double[dof]{0.};

    for(i = 0; i < sot_size; i++){
        if(s_o_et.at(i) == 0){
            continue;
        }

        if(i == 0){
            B2TS ele_empty;
            std::vector<B2TS>* eles_ref = eles.eleset_ptr(ele_empty);

            for(std::vector<B2TS>::iterator itele = eles_ref->begin(); itele != eles_ref->end(); itele++){
                itele->asb_KF(K.val, Fout, s_o_nt);
            }
        }
        else if(i == 1){
            B2TSR ele_empty;
            std::vector<B2TSR>* eles_ref = eles.eleset_ptr(ele_empty);

            for(std::vector<B2TSR>::iterator itele = eles_ref->begin(); itele != eles_ref->end(); itele++){
                itele->asb_KF(K.val, Fout, s_o_nt);
            }
        }
        else if(i == 2){
            B3TS ele_empty;
            std::vector<B3TS>* eles_ref = eles.eleset_ptr(ele_empty);

            for(std::vector<B3TS>::iterator itele = eles_ref->begin(); itele != eles_ref->end(); itele++){
                itele->asb_KF(K.val, Fout, s_o_nt);
            }
        }
        else if(i == 3){
            B3TSR ele_empty;
            std::vector<B3TSR>* eles_ref = eles.eleset_ptr(ele_empty);

            for(std::vector<B3TSR>::iterator itele = eles_ref->begin(); itele != eles_ref->end(); itele++){
                itele->asb_KF(K.val, Fout, s_o_nt);
            }
        }
        else{
            printf("no such element type %i\n", i);
        }
    }
};

void asb_manager::add_bnd(){
    const std::vector<std::pair<int,int>> sont = eles.sont();
    int ntype_num = sont.size();
    std::vector<int> noet(ntype_num + 1, 1);
    std::vector<int> dofoet(ntype_num + 1, 0);
    for(int i = 0; i < ntype_num; i++){
        noet.at(i + 1) = noet.at(i) + sont.at(i).first;
        dofoet.at(i + 1) = dofoet.at(i) + sont.at(i).second * sont.at(i).first;
    }

    int where_nd, row;
    int dof = dofoet.back();
    for(std::vector<int>::iterator itbd = dbd.Anchorage_nid.begin(); itbd != dbd.Anchorage_nid.end(); itbd++){
        where_nd = find(*itbd, noet);
        row = dofoet.at(where_nd) + sont.at(where_nd).second * (*itbd - noet.at(where_nd));

        for(int i = 0; i < sont.at(where_nd).second; i++){
            K.val[(row + i)*dof + row + i] = 1E10;
            Fout[row + i] = 0;
        }
    }
}

void asb_manager::solve(){
    int* ipiv = new int[K.rows]{0};
    int info;

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, K.rows, 1, K.val, K.rows, ipiv, Fout, 1);

    if(info != 0){
        printf("ERROR: solving process error!\n");
        exit(1);
    }

    delete[] ipiv; ipiv = nullptr;
}
