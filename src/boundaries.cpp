#include "boundaries.h"

void boundaries::add_bndD(DenseMatrix& Kd, double* Fo, elements& eles){
    const std::vector<std::pair<int,int>> sont = eles.sont();
    std::vector<int> noet, dofoet;
    int ntype_num = decode_sont(sont, noet, dofoet);

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

void boundaries::add_bndS(SparseMatrix& Ks, double* Fo, elements& eles, bool is_NL){
    const std::vector<std::pair<int,int>> sont = eles.sont();
    std::vector<int> noet, dofoet;
    int ntype_num = decode_sont(sont, noet, dofoet);

    int where_nd, row, posi_s, ifix;

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

    for(std::vector<dload_2d>::iterator itdbd = dbd2ds.begin(); itdbd != dbd2ds.end(); itdbd++){
        double norm = pow(itdbd->norm_vec.at(0), 2) + pow(itdbd->norm_vec.at(1), 2);
        double val_xy[2]{0.}, length = 0.;
        int i, elenumber, ltag, face[3], posi_v[3],
            where;
        double wei[3];
        norm = sqrt(norm);
        for(i = 0; i < 2; i++){
            val_xy[i] = itdbd->val * itdbd->norm_vec.at(i) / norm;
        }

        PSL_ ele_empty;
        std::vector<PSL_>* eles_ref = eles.eleset_ptr(ele_empty);

        where = find(eles_ref->front().node_id.front(), noet);
        for(std::vector<std::pair<int, int>>::iterator itel = itdbd->elset_line.begin(); itel != itdbd->elset_line.end(); itel++){
            elenumber = itel->first;
            ltag = itel->second;

            if(itdbd->type == 0){
                wei[0] = 0.5;
                wei[1] = 0.5;
                face[0] = ltag - 1;
                face[1] = ltag % 4;
                length = getlen(eles_ref->at(elenumber-1).xy.at(face[0]), eles_ref->at(elenumber-1).xy.at(face[1]), 2);
                for(i = 0; i < 2; i++){
                    posi_v[i] = dofoet.at(where) + 2 * (eles_ref->at(elenumber-1).node_id.at(face[i]) - noet.at(where));
                    Fo[posi_v[i]] += val_xy[0] * length * wei[i];
                    Fo[posi_v[i] + 1] += val_xy[1] * length * wei[i];
                }
            }
            else{
                face[0] = ltag - 1; wei[0] = 1.0 / 6.0;
                face[1] = ltag + 3; wei[1] = 2.0 / 3.0;
                face[2] = ltag % 4; wei[2] = wei[0];
                length = getlen(eles_ref->at(elenumber-1).xy.at(face[0]), eles_ref->at(elenumber-1).xy.at(face[2]), 2);
                for(i = 0; i < 3; i++){
                    posi_v[i] = dofoet.at(where) + 2 * (eles_ref->at(elenumber-1).node_id.at(face[i]) - noet.at(where));
                    Fo[posi_v[i]] += val_xy[0] * wei[i] * length;
                    Fo[posi_v[i] + 1] += val_xy[1] * wei[i] * length;
                }
            }
        }
    }

    for(std::vector<cload_2d>::iterator itcd = cbd2ds.begin(); itcd != cbd2ds.end(); itcd++){
        where_nd = find(itcd->nset.front(), noet);
        for(std::vector<int>::iterator itnd = itcd->nset.begin(); itnd != itcd->nset.end(); itnd++){
            posi_s = dofoet.at(where_nd) + sont.at(where_nd).second * (*itnd - noet.at(where_nd));
            Fo[posi_s] += itcd->val[0];
            Fo[posi_s + 1] += itcd->val[1];
        }
    }

    ifix = 0;
    for(std::vector<int>::iterator itbd = dbd.Anchorage_nid.begin(); itbd != dbd.Anchorage_nid.end(); itbd++){
        where_nd = find(*itbd, noet);
        row = dofoet.at(where_nd) + sont.at(where_nd).second * (*itbd - noet.at(where_nd));
        for(int i = 0; i < sont.at(where_nd).second; i++){
            if(dbd.fixed_dof.at(ifix).size() > i){
                if(dbd.fixed_dof.at(ifix).at(i)){
                    posi_s = match(row + i + 1, row + i + 1, Ks);
                    Ks.val[posi_s] *= is_NL ? 1 : 1E12;
                    Fo[row + i] = 0;
                }
            }
            else{
                posi_s = match(row + i + 1, row + i + 1, Ks);
                Ks.val[posi_s] *= is_NL ? 1 : 1E12;
                Fo[row + i] = 0;
            }
        }
        ifix++;
    }

    // for(std::vector<plate_cload>::iterator itpc = pbd.pbd_c.begin(); itpc != pbd.pbd_c.end(); itpc++){
    //     where_nd = find(itpc->node_set.front(), noet);
    //     for(std::vector<int>::iterator itns = itpc->node_set.begin(); itns != itpc->node_set.end(); itns++){
    //         row = dofoet.at(where_nd) + sont.at(where_nd).second * (*itns - noet.at(where_nd));
    //     }
    // }
};