#include "assemble.h"

void asb_manager::init_mesh(int num, bool M, int type, double len, double val){
    Linemesh_demo lmd(num, len, val, M, type);
    lmd.generate(eles, mater_lib, section_lib, bnds, xyz_coord);
}

void asb_manager::init_mesh(int num, int type, double len, double val){
    Square_Platemesh_demo spmd(num, len, val, type);
    spmd.generate(eles, mater_lib, plate_prop_lib, bnds, xyz_coord);
};

void asb_manager::init_mesh(int w_num, int type, int b_type, double pcr, double seca,
                double len, int order){
    Trussmesh_demo tsmd(b_type, w_num, len, pcr, type, seca, order);
    tsmd.generate(eles, mater_lib, plate_prop_lib, bnds, sln, xyz_coord, analysis_dim);
};

void asb_manager::solve(int** parms){
    if(!_strcmpi(method.data(), "D")){
        sln.solve_D(KF, bnds, eles);
    }
    else if(!_strcmpi(method.data(), "S")){
        sln.solve_S(KF, bnds, eles);
    }
    else if(!_strcmpi(method.data(), "SN")){
        if(!_strcmpi(task.data(), "MLSM")){
            sln.solve_SN_MLSM(KF, bnds, eles, xyz_coord, parms);
        }
        else if(!_strcmpi(task.data(), "GALM")){
            sln.solve_SN_GALM(KF, bnds, eles, xyz_coord, parms);
        }
    }
};