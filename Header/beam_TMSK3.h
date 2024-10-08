#ifndef BEAM_TMSK3_H
#define BEAM_TMSK3_H

#include <stdio.h>
#include <vector>
#include "element_prop.h"
#include "solver.h"

class B3TS{
public:
    std::vector<int> node_id;
    material* mat;
    section_TSbeam* sec_prop;
    std::vector<double*> xyz;
    beam_bnd bnd;

    B3TS(material* mater = nullptr, section_TSbeam* sec = nullptr){
        mat = mater;
        sec_prop = sec;
    };

    /// @brief set the boundary information
    /// @param bd beam boundary
    void set_bnd(const beam_bnd &bd);

    /// @brief set the node information
    void set_node(std::vector<int> &nid, std::vector<double*> &x_c);

    /// @brief get the stiffness matrix and load vector of element
    /// @param K stiffness matrix
    /// @param F load vector
    void get_KF(double K[6][6], double F[6]);

    /// @brief assemble element stiffness matrix into global stiffness matrix
    /// @param K global stiffness matrix
    /// @param F global load vector
    /// @param sont number of each node type
    void asb_KF(double* K, double* F, const std::vector<std::pair<int, int>> &sont);

    ~B3TS(){};
};

#endif