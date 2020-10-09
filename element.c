/*=================================================================
 * element.c - example used to illustrate calculation of
 * elment matrices of beam elements.
 *
 * Input: E,rho,L,A,I,R,Ax,Ix
 * Output: Kg,Mg
 *
 * Creates element matrix for beam elements
 *
 * author : Suguang Dou, DTU Wind Energy
 * date   : October 9th, 2020
 * contact: sudou@dtu.dk; dousuguang@gmail.com
 *=================================================================*/

#include <math.h>
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize i, j, k, l; // int // size_t
    
#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble E, rho, L, A, I, *R, Ax, Ix, c, s;
    mxDouble AEL, EIL1, EIL2, EIL3, EIL4, rhoA;
    mxDouble *Mg, *Kg, Ke[6][6]={0}, Me[6][6]={0};
    mxDouble *y, *z, x;
#else
    double E, rho, L, A, I, *R, Ax, Ix, c, s;
    double AEL, EIL1, EIL2, EIL3, EIL4, rhoA;
    double *Mg, *Kg, Ke[6][6]={0}, Me[6][6]={0};
#endif
    
    // for (i=0; i<6; i++) {
    //     for (j=0; j<6; j++) {
    //         mexPrintf("%f", Me[i][j]);
    //     }
    //     mexPrintf("\n");
    // }
    
    
    // Get the value of the scalar input
    // (E,rho,L,A,I,R,Ax,Ix)
    E   = 	mxGetScalar(prhs[0]);
    rho = 	mxGetScalar(prhs[1]);
    L 	= 	mxGetScalar(prhs[2]);
    A 	= 	mxGetScalar(prhs[3]);
    I 	= 	mxGetScalar(prhs[4]);
    // Ax, Ix
    if (nrhs>6) {
        Ax 	= 	mxGetScalar(prhs[6]);
        Ix 	= 	mxGetScalar(prhs[7]);
        A   =   Ax;
        I   =   Ix;
    }
    // R
    R   =   mxGetPr(prhs[5]);
    // mexPrintf("%f\n", R[1]);
    
    
    c = R[0];
    s = R[1];
    double T[6][6]={{c,s,0,0,0,0},{-s,c,0,0,0,0},{0,0,1,0,0,0},{0,0,0,c,s,0},{0,0,0,-s,c,0},{0,0,0,0,0,1}};
    // for (i=0; i<6; i++) {
    //     for (j=0; j<6; j++) {
    //         mexPrintf("%f\n", T[i][j]);
    //     }
    //     mexPrintf("\n");
    // }
    
    
    AEL = A*E/L;
    EIL1 = 2*E*I/L;
    EIL2 = 6*E*I/pow(L,2);
    EIL3 = 12*E*I/pow(L,3);
    EIL4 = 4*E*I/L;
    
    double K[6][6]={{AEL,0,0,-AEL,0,0},{0,EIL3,EIL2,0,-EIL3,EIL2},{0,EIL2,EIL4,0,-EIL2,EIL1},{-AEL,0,0,AEL,0,0},{0,-EIL3,-EIL2,0,EIL3,-EIL2},{0,EIL2,EIL1,0,-EIL2,EIL4}};
    // for (i=0; i<6; i++) {
    //     for (j=0; j<6; j++) {
    //         mexPrintf("%10.1f", K[i][j]);
    //     }
    //     mexPrintf("\n");
    // }
    
    double M[6][6]={{140*L,0,0,70*L,0,0},{0,156*L,22*pow(L,2),0,54*L,-13*pow(L,2)},{0,22*pow(L,2),4*pow(L,3),0,13*pow(L,2),-3*pow(L,3)},{70*L,0,0,140*L,0,0},{0,54*L,13*pow(L,2),0,156*L,-22*pow(L,2)},{0,-13*pow(L,2),-3*pow(L,3),0,-22*pow(L,2),4*pow(L,3)}};
    rhoA = rho*A/420;
    for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
            M[i][j] = M[i][j]*rhoA;
        }
    }
    // for (i=0; i<6; i++) {
    //     for (j=0; j<6; j++) {
    //         mexPrintf("%10.1f", M[i][j]);
    //     }
    //     mexPrintf("\n");
    // }
    
    for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
            for (k=0; k<6; k++) {
                for (l=0; l<6; l++) {
                    Ke[i][j] = Ke[i][j] + T[k][i]*K[k][l]*T[l][j];
                    Me[i][j] = Me[i][j] + T[k][i]*M[k][l]*T[l][j];
                }
            }
        }
    }
    // for (i=0; i<6; i++) {
    //     for (j=0; j<6; j++) {
    //         mexPrintf("%10.1f", Ke[i][j]);
    //     }
    //     mexPrintf("\n");
    // }
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(6,6,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(6,6,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
#if MX_HAS_INTERLEAVED_COMPLEX
    Kg = mxGetDoubles(plhs[0]);
    Mg = mxGetDoubles(plhs[1]);
#else
    Kg = mxGetPr(plhs[0]);
    Mg = mxGetPr(plhs[1]);
#endif
    
    for (i=0; i<6; i++) {
        for (j=0; j<6; j++) {
            k = j*6+i;
            Kg[k]= Ke[i][j]; 
            Mg[k]= Me[i][j];
        }
    }
}