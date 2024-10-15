#include <stdio.h>
#include <vector>
#include "element_prop.h"
#include "solver.h"

#ifndef BEAM_TMSK2_H
#define BEAM_TMSK2_H

/// @brief Temoshenko beam element with 2 nodes
class B2TS{
public:
    int int_method;             // 0-Full integration; 1-reduced integration
    std::vector<int> node_id;
    material* mat;
    section_TSbeam* sec_prop;
    std::vector<double*> xyz;
    beam_bnd bnd;

// public:

    B2TS(int method = 0, material* mater = nullptr, section_TSbeam* sec = nullptr){
        mat = mater;
        sec_prop = sec;
        int_method = method;
    };

    /// @brief set the boundary information
    /// @param bd beam boundary
    void set_bnd(const beam_bnd &bd);

    /// @brief set the node information
    void set_node(std::vector<int> &nid, std::vector<double*> &x_c);

    /// @brief get the stiffness matrix and load vector of element
    /// @param K stiffness matrix
    /// @param F load vector
    void get_KF(double K[4][4], double F[4]);

    /// @brief assemble element stiffness matrix into global stiffness matrix
    /// @param K global stiffness matrix
    /// @param F global load vector
    /// @param sont number of each node type
    void asb_KF(double* K, double* F, const std::vector<std::pair<int, int>> &sont);

    ~B2TS(){};
};

#endif