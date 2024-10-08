#ifndef BEAM_TMSK3R_H
#define BEAM_TMSK3R_H

#include "beam_TMSK3.h"

class B3TSR: public B3TS
{
public:
    B3TSR(material* mater = nullptr, section_TSbeam* sec = nullptr):B3TS(mater, sec){};

    /// @brief get the stiffness matrix and load vector of element
    /// @param K stiffness matrix
    /// @param F load vector
    void get_KF(double K[6][6], double F[6]);

    /// @brief assemble element stiffness matrix into global stiffness matrix
    /// @param K global stiffness matrix
    /// @param F global load vector
    /// @param sont number of each node type
    void asb_KF(double* K, double* F, const std::vector<std::pair<int, int>> &sont);
};

#endif