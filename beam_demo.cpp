#include <stdio.h>
#include <iostream>
#include "assemble.h"

int main(){
    asb_manager asb;
    int elenum = 0;
    int eletype = 0;
    int isM = 1;

    printf("Please select element type: \n\
(0 - 2-node Timoshenko beam element; \n\
 1 - 2-node Timoshenko beam element with reduced integration; \n\
 2 - 3-ndoe Timoshenko beam element; \n\
 3 - 2-node Timoshenko beam element with reduced integration)\n");
    std::cin >> eletype;

    printf("Please input element number: \n");
    std::cin >> elenum;

    printf("Please select load conditions: \n\
(1 - bending moment; \n\
 0 - sheering force)\n");
    std::cin >> isM;

    if(eletype < 0 || eletype > 3){
        std::cerr << "ERROR: element type " << eletype << " is not supported!";
        exit(0);
    }
    asb.init_mesh(elenum, isM, eletype, 1.0, 2.5);
    double start = clock();
    asb.init_KF();

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
    asb.add_bnd();
    asb.solveD();

    if(eletype == 1 || eletype == 0){
        for(int i = 0; i < 2*(elenum+1); i++){
            printf("%.8f\n", asb.KF.Fout[i]);
        }
    }
    else{
        for(int i = 0; i < 2*(2*elenum+1); i++){
            printf("%.8f\n", asb.KF.Fout[i]);
        }
    }

    printf("%.4f\n", asb.xyz_coord[3]);

    double end = clock();
    printf("time: %.4f\n", end - start);
}