
#include <assert.h>
#include "CG.hpp"
#include <iostream>
double* comblineaire(int n, double a, double *A,
                     double b, double *B,double *R )
{
    for(int i=0;i< n;++i)
        R[i] = a*A[i] + b*B[i];
    return R;
}
double sdot(int n,double *A, double *B)
{
    double s=0;
    for(int i=0;i< n;++i)
        s += A[i]*B[i];
    return s;
}

int  CG( const MatVirt &A,const  MatVirt &C,
                    double *b, // second membre
                    double * x, // solution qui contient une initialisation
                    int nbitermax,double eps)
{
    int n = A.n;
    double * G  = new double[n];
    double * H  = new double[n];
    double * CG = new double[n];
    double * AH = CG;
    double * Ax = CG;
    double Cgg, Cggp, rho, gamma ;
    assert( n == A.m && n == C.n && n == C.m);
    comblineaire(n,1.,A.MatMul(x,Ax) , -1, b,G);
    C.MatMul(G,CG);
    
    for(int i=0;i< n;++i)
        H[i] = - CG[i];
    Cgg= sdot(n,G,CG);
    if( Cgg > eps)
    {
        for(int iter=0; iter <nbitermax ; ++ iter)
        {
            A.MatMul(H,AH);
            rho = - sdot(n,G,H)/sdot(n,H,AH);
            // Error comblineaire(n, 1.,x, rho, G,x);
            comblineaire(n, 1.,x, rho, H,x);
            comblineaire(n, 1.,G, rho, AH,G);
            Cggp = Cgg;
            C.MatMul(G,CG);
            Cgg =sdot(n,G,CG);
            gamma = Cgg / Cggp;
            comblineaire(n, -1., CG, gamma, H, H);
            //std::cout << iter<< " / " << nbitermax << " Cgg " << Cgg
                      //<< " , " << gamma  <<" , " << rho<<std::endl;
            if( Cgg < eps) break;
        }
    }
    
    delete[] G;
    delete[] H;
    delete[] CG;
    
    return Cgg < eps;
}


