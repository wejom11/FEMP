#include <stdio.h>
#include <iostream>
#include "assemble.h"
#include "interpolate.h"

int main(){
    asb_manager asb;
    int elenum = 0;
    int eletype = 0;
    double len = 1;

    printf("Please select element type: \n\
(0 - PQLF; 1 - PQLS; 2 - PQLR;\n\
 3 - PQSF; 4 - PQSS; 5 - PQSR)\n");
    std::cin >> eletype;

    printf("Please line element number ln (total element number: ln^2): \n");
    std::cin >> elenum;

    printf("Please input length of square plate: \n");
    std::cin >> len;

    if(eletype < 0 || eletype > 5){
        std::cerr << "ERROR: element type " << eletype << " is not supported!";
        exit(0);
    }
    asb.init_mesh(elenum, eletype, len, 10.0);
    double start = clock();
    // for(int i = 0; i < 2 * (elenum + 1); i++){
    //     printf("%.4f\n", asb.Fout[i]);
    // }
    // asb.K.show();
    // DenseMatrix Kw = asb.K;
    // double zero = Kw.val[0];
    // for(int i = 0; i < Kw.rows * Kw.cols; i++){
    //     Kw.val[i] /= zero;
    // }
    // Kw.show();
    // asb.KF.K_sparse.show(asb.KF.K_sparse.rows);
    // for(int i = 0; i < asb.KF.K_sparse.rows + 1; i++){
    //     printf("%i\n", asb.KF.K_sparse.row_st[i]);
    // }
    // printf("\n");
    // asb.KF.K_sparse.show(asb.KF.K_sparse.rows);
    asb.solve();

    double D = asb.mater_lib.front().E * pow(asb.plate_prop_lib.front().thickness, 3) 
               / 12 /(1 - pow(asb.mater_lib.front().nu, 2));
    if(eletype == 1 || eletype == 0 || eletype == 2){
        int center = ((2*elenum+1)*(2*elenum+1) + 1)/2 - 1;
        printf("%.15f\n", asb.sln.Var[3*center + 2]*D/10.0/pow(len, 4));
    }
    else{
        int center = ((2 * elenum + 1) * (2 * elenum + 1) - elenum*elenum + 1) / 2 - 1;
        printf("%.15f\n", asb.sln.Var[3 * center + 2] * D / 10.0 / pow(len, 4));
    }

    double end = clock();
    printf("time: %.4f\n", end - start);
}