#include "boundaries.h"

void boundaries::add_bndD(DenseMatrix& Kd, double* Fo, elements& eles){
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
            Kd.val[(row + i)*dof + row + i] = 1E10;
            Fo[row + i] = 0;
        }
    }
};

void boundaries::add_bndS(SparseMatrix& Ks, double* Fo, elements& eles){
    const std::vector<std::pair<int,int>> sont = eles.sont();
    int ntype_num = sont.size();
    std::vector<int> noet(ntype_num + 1, 1);
    std::vector<int> dofoet(ntype_num + 1, 0);
    for(int i = 0; i < ntype_num; i++){
        noet.at(i + 1) = noet.at(i) + sont.at(i).first;
        dofoet.at(i + 1) = dofoet.at(i) + sont.at(i).second * sont.at(i).first;
    }

    int where_nd, row, posi;

    for(std::vector<plate_dload>::iterator itpd = pbd.pbd_d.begin(); itpd != pbd.pbd_d.end(); itpd++){
        if(itpd->ele_set.front() == 2){
            PQL_ ele_empty;
            std::vector<PQL_>* eles_ref = eles.eleset_ptr(ele_empty);

            for(std::vector<int>::iterator ites = itpd->ele_set.begin() + 1; ites != itpd->ele_set.end(); ites++){
                eles_ref->at(*ites - 1).asb_F_dload(Fo, itpd->val, sont);
            }
        }
        else if(itpd->ele_set.front() == 3){
            PQS_ ele_empty;
            std::vector<PQS_>* eles_ref = eles.eleset_ptr(ele_empty);

            for(std::vector<int>::iterator ites = itpd->ele_set.begin() + 1; ites != itpd->ele_set.end(); ites++){
                eles_ref->at(*ites - 1).asb_F_dload(Fo, itpd->val, sont);
            }
        }
    }

    for(std::vector<int>::iterator itbd = dbd.Anchorage_nid.begin(); itbd != dbd.Anchorage_nid.end(); itbd++){
        where_nd = find(*itbd, noet);
        row = dofoet.at(where_nd) + sont.at(where_nd).second * (*itbd - noet.at(where_nd));
        for(int i = 0; i < sont.at(where_nd).second; i++){
            posi = match(row + i + 1, row + i + 1, Ks);
            Ks.val[posi] *= 1E8;
            Fo[row + i] = 0;
        }
    }

    // for(std::vector<plate_cload>::iterator itpc = pbd.pbd_c.begin(); itpc != pbd.pbd_c.end(); itpc++){
    //     where_nd = find(itpc->node_set.front(), noet);
    //     for(std::vector<int>::iterator itns = itpc->node_set.begin(); itns != itpc->node_set.end(); itns++){
    //         row = dofoet.at(where_nd) + sont.at(where_nd).second * (*itns - noet.at(where_nd));
    //     }
    // }
};