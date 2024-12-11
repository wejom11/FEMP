#ifndef PSL__H
#define PSL__H

#include <vector>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include "element_prop.h"
#include "interpolate.h"

class PSL_{
public:
    short order; // 0 -- 1(PSL4); 1 -- 2(PSL9)
    int nnum;
    int gauss_num;

    plate_prop* pp;
    material* mat;
    std::vector<int> node_id;
    std::vector<double*> xy;
    std::vector<double*> uv;

    Lagrange_intp shape_info;

    PSL_(material* mat_ptr = nullptr, plate_prop* pp_ptr = nullptr, bool ord = true){
        pp = pp_ptr;
        mat = mat_ptr;
        order = ord ? 1 : 0;
        if(order){
            nnum = 9;
            gauss_num = 5;
            node_id.resize(9, 0);
            xy.resize(9, nullptr);
            uv.resize(9, nullptr);
        }
        else{
            nnum = 4;
            gauss_num = 3;
            node_id.resize(4, 0);
            xy.resize(4, nullptr);
            uv.resize(4, nullptr);
        }
    }

    /// @attention only set the element's order before set coordinate and displacement
    /// @param ord order: 1 --- PSL9; 0 -- PSL4
    void set_order(bool ord);

    /// @brief select gauss integration order
    /// @param order default: PSL4 -> 3; PSL9 -> 5
    void set_gauss_order(short order);

    /// @brief set coordinate, displacement and node_id
    /// @param global_xy golbal xy coordinate
    /// @param global_uv global displacement
    /// @param n_id node id
    void set_xy_uv(double* global_xy, double* global_uv, const std::vector<int>& n_id);
    /// @brief set coordinate, displacement(node_id is already setted)
    /// @param global_xy golbal xy coordinate
    /// @param global_uv global displacement
    void set_xy_uv(double* global_xy, double* global_uv);

    void ele_init();

    void get_BL0T(DenseMatrix& BL0T, short intp);

    void get_BL1T(DenseMatrix& BL1T, short intp);

    void get_BNLT(DenseMatrix& BNLT, short intp);

    void get_S(DenseMatrix& St, short intp);
    void get_S(double* Sbar, short intp);

    void get_KL0(DenseMatrix& KL0);
 
    /// @param option if value is 01, ignore the term B_{L1}^{T} D B_{L1}
    ///
    ///               if value is 10, ignore the term B_{L0}^{T} D B_{L1} + B_{L1}^{T} D B_{L0}
    ///
    ///               if value is 11, do not ignore any term
    void get_KL1(DenseMatrix& KL1, short option = 11);

    void get_KNL(DenseMatrix& KNL);

    void get_F(double* Ft);

    /// @brief get the element's stiffness and internal load
    /// @param Ke stiffness matrix
    /// @param Fe internal load vector
    /// @param option four digit int; representation wheather KL0, KNL, Ku2, Ku1 to be calculated
    void get_KF(DenseMatrix& Ke, double* Fe, short option = 1111);

    ~PSL_(){

    }
};

#endif