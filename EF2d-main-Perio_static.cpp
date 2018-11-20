#include "EF2d-base.hpp"
#include <fstream>
#include "CG.hpp"
#include <cmath>
#include <typeinfo>
#include <vector>

using namespace std;


typedef double R;

class MatrixLap:public MatVirt {
public:
    Maillage2d &Th;
    R alpha,*beta;
    R *coef;
    int *dl;
    MatrixLap( Maillage2d &T,int *d,int np,R*cc, R a, R*b  ):MatVirt(np), Th(T), alpha(a),beta(b), coef(cc),dl(d)  {}
    R * AddMatMul(const R *x, R *Ax) const;
};


R* MatrixLap::AddMatMul(const R *u, R *Au) const
{ 
   //   += int_ alpha u*v + beta grad u . grad v
    const int  d = Triangle::d,nbv=Triangle::nbv;
    double coeff= alpha/(d+1) ; //modif
    R2 G[3];
    double ak[nbv][nbv];
    
    for(int k=0; k<Th.nt;++k)
    {   
        const Triangle & K=Th[k];      
        const R mes =K.mes, betames= beta[k]*mes , coefmes=coeff*mes;
        
        if(beta[k]) K.GradLambdaK(G); 
        for(int ip=0;ip<nbv;++ip)
            for(int jp=0; jp< nbv; ++jp)
                ak[ip][jp] = betames*(G[ip],G[jp]) + coefmes ;   //modif

        for(int ip=0; ip< nbv; ++ip)
        {
          int i= Th(K[ip]);
            int in =dl[i];
            R ci = coef[i];
          if( !K[ip].OnGamma())
           for(int jp=0;jp<nbv;++jp)
          {
              int j = Th(K[jp]);
              int jn = dl[j];
              R cj = coef[j] ;
              Au[in] += ci*cj*ak[ip][jp] * u[jn] ; 
              
        }}
    }
    return Au;
}



int main(int argc, const char ** argv)
{   
    cout << " Usage: " << argv[0] << " file.msh perio\n";
    assert(argc>2);
    Maillage2d Th(argv[1]); 
    int typeperio=atoi(argv[2]); 
    assert(abs(typeperio)==1);
    int nv  = Th.nv;
    int nt = Th.nt;
    int * dl = new int[nv];
    R * coef = new R[nv];  
    vector<int> t1; 
    vector<int> t2; 
    for(int i=0; i< nv; ++i)
    { 
	if(Th.v[i].lab==3)  // bord inferieur
           t1.push_back(i); 

        if(Th.v[i].lab==2)  // bord superieur
           t2.push_back(i); 
    }

    // Tri du tableau t2 en fonction de t1
    unsigned int tmp=0, i=0, j=0;
    assert(t1.size()==t2.size());
    for(i=0; i< t1.size(); ++i)
    { 
        int i1=0, i2=0;
        for( j=0; j< t2.size() ; ++j)
        { 
            i1=t1[j];
            i2=t2[i]; 
            double R1= sqrt( (Th.v[i1].x) * (Th.v[i1].x) + (Th.v[i1].y) * (Th.v[i1].y));
            double R2= sqrt( (Th.v[i2].x) * (Th.v[i2].x) + (Th.v[i2].y) * (Th.v[i2].y)); 
            if ( R1 == R2) 
                  { 
                    tmp = t2[j];
                    t2[j] = t2[i];
                    t2[i] = tmp;
                    j=i; 
                  }
        }
    }

    // Pour tout i appartenant au bord 3 on associe son image j du bord 2 par la rotation pi/6
    unsigned int s=t1.size();
    int Corr[s][2];
    for(i=0; i< s; ++i)
    {
        Corr[i][0] = t1[i];
        Corr[i][1] = t2[i];
    }

    /*for(unsigned int i=0; i< s; ++i)
       cout << Corr[i][0] << " est l image par R de " << Corr[i][1] << endl;*/

    //Renumerotation
    int count = 0;
    for( int i=0; i< nv; ++i)
    { 
        dl[i]=i-count;
        coef[i]= 1;
        for(unsigned j=0; j< s; ++j)
           if (i== Corr[j][0])
           { 
                dl[i]=Corr[j][1];
                coef[i]= typeperio;
                count++;
           }
    //cout << "dl" << i << " = " << dl[i] << " ; " << "c" << i << " = " << coef[i] << endl;  
    }

    //construction du tableau beta
    R epsilon = 5.1636e-04;
    R * beta = new R[nt]; 
    for(int l =0; l< nt; ++l)
    {
         const Triangle & K=Th[l];
         if ( (K.lab == 1) || (K.lab == 2)) // regions contenant du fer
             beta[l] = epsilon; 
         else 
             beta[l] = 1; 

         beta[l] /= (4*M_PI*1e-07);      
    }

    int n = dl[nv-1]+1; // nombre d'inconnues

    MatrixLap A(Th,dl,n,coef,0,beta);  

    // beta =0 pour le second membre
    R * beta_nul = new R[nt];
    for(int k=0 ; k < nt; ++k)
        beta_nul[k] =0.;

    MatrixLap M(Th,dl,n,coef,1,beta);  
    
    R * fh = new R[n];
    R * b  = new R[n];
    R * u  = new R[n]; // periodique
    R * un = new R[nv]; // non periodique

    // calcul des aires des regions
    int nb_regions= 8;
    int regions[8] = {0, 1, 2, 3, 19, 20, 21, 22};
    R  aires [8] = {0., 0., 0., 0., 0., 0., 0., 0.}, aire_totale=0.;
    for (int i=0; i < nb_regions; ++i)
    {
        for(int k=0; k < nt; ++k)
        {   
            const Triangle & K=Th[k];
            if (K.lab == regions[i])
            {
                const R mes =K.mes;
                aires[i]+= mes;
            }
        }
        aire_totale += aires[i] ;
    }

    double kappa = 80. ;//section transversale
    for (int i=0; i< n;++i)
       fh[i]=0;

    for(int k=0; k < nt; ++k)
    { 
        const Triangle & K=Th[k]; 
        for(int i=0; i< K.nbv; ++i)     
        {
            int in = dl[Th(K[i])]; 
            if (K.lab == 19 || K .lab == 21)  // domaine bobine d'alimentation
	       fh[in]= kappa*100/(aires[4]+aires[6]) ; 
        }
    }

    M.MatMul(fh,b); // assemblage second membre
    
    // initialisation
    for (int i=0; i< n;++i)
         u[i]=0.;
         
    R eps = 1e-7;
    int ok = CG(A,MatId(n),b,u,n,eps); // resolution avec le gradient conjugue
    if (ok)
       cout << "CG a converge pour epsilon " << eps << endl;
    else 
       cout << "CG n'a pas converge pour epsilon " << eps << endl;    

    // reconstuction du champ non periodic un, a partir du champ u.
    for(int i=0; i< nv;++i)
         un[i]=coef[i]*u[dl[i]]; 

    // sauvegarde de la solution periodique
    ofstream ff("sol.txt");
    for(int i=0; i< n; ++i)
       ff << u[i] << endl;   

    // pour ploter la solution   
    ofstream f("m.gp");
    for(int k=0; k < nt; ++k)
    { 
        for (int ip=0;ip< 4 ;++ ip)
        {
            int i= Th(k,ip%3);
             f << Th.v[i].x << " " <<Th.v[i].y << " " << un[i] << endl ;
        }
             f << endl << endl;
    }

    // clean memory ...
    delete [] b;
    delete [] u;
    delete [] un;
    delete [] fh;
    delete [] dl;
    delete [] coef ;
    delete [] beta;
    delete [] beta_nul;
}
