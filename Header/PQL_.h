#ifndef PQL__H
#define PQL__H

#include <vector>
#include <stdio.h>
#include "element_prop.h"
#include "interpolate.h"

/// order of node:
/// 4--7--3
/// |  |  |
/// 8--9--6
/// |  |  |
/// 1--5--2
class PQL_
{
public:
    short int inte_sche;      //0 - Full(PQLF); 1 - select(PQLS); 2 - reduced(PQLR)
    std::vector<int> node_id;
    std::vector<double*> xyz;
    material* mat;
    plate_prop* pp;

    PQL_(int itsh = 0, material* mater = nullptr, plate_prop* p_p = nullptr){
        inte_sche = itsh;
        mat = mater;
        pp = p_p;
    };

    /// @brief set elements properties
    /// @param nid node number
    /// @param xc xyz coordinate
    /// @param pbnd boundary conditions
    void set(std::vector<int> nid, std::vector<double*> xc);

    void get_Kb(double Kb[27][27]);

    void get_Ks(double Ks[27][27]);

    void asb_KS(SparseMatrix& Kasb, const std::vector<std::pair<int, int>> &sont);

    void asb_F_dload(double* Fo, double val, const std::vector<std::pair<int, int>> &sont);

    ~PQL_(){};
};

#endif