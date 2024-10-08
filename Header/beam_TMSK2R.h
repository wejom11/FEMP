#ifndef BEAM_TMSK2R_H
#define BEAM_TMSK2R_H

#include "beam_TMSK2.h"

class B2TSR: public B2TS
{
public:
    B2TSR(material* mater = nullptr, section_TSbeam* sec = nullptr):B2TS(mater, sec){};

    /// @brief get the stiffness matrix and load vector of element
    /// @param K stiffness matrix
    /// @param F load vector
    void get_KF(double K[4][4], double F[4]);

    /// @brief assemble element stiffness matrix into global stiffness matrix
    /// @param K global stiffness matrix
    /// @param F global load vector
    /// @param sont number of each node type
    void asb_KF(double* K, double* F, const std::vector<std::pair<int, int>> &sont);
};

#endif