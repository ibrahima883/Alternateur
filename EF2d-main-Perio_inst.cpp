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
    R alpha,**beta;
    R *coef;
    int *dl;
    MatrixLap( Maillage2d &T,int *NN,int n,R*cc, R a, R**b  ):MatVirt(n), Th(T), alpha(a),beta(b), coef(cc),dl(NN)  {}
    R * AddMatMul(const R *x, R *Ax) const;
};


R* MatrixLap::AddMatMul(const R *u, R *Au) const
{ 
   //   += int_ alpha u*v + beta grad u . grad v
    const int  d = Triangle::d,nbv=Triangle::nbv;
    R2 G[3];
    double ak[nbv][nbv]; 
    double coeff= alpha/(d+1); 
    for(int k=0; k<Th.nt;++k)
    {   
        const Triangle & K=Th[k];
        const R coef0= beta[k][1]/(d+2)/(d+1);
        const R mes =K.mes, betames= beta[k][0]*mes , coefmes = coeff*mes, coef0mes= coef0*mes;  
        const R coef1 = beta[k][2]/(d+1), coef2 = beta[k][3]/(d+1);
        const R coef1mes=coef1*mes*mes, coef2mes=coef2*mes*mes;
        
        if(beta[k][0] )  K.GradLambdaK(G); 
        for(int ip=0;ip<nbv;++ip)
            for(int jp=0; jp< nbv; ++jp)
                 ak[ip][jp] = betames*(G[ip],G[jp]) + coef0mes*(1+ (ip==jp)) + coef1mes + coef2mes + coefmes;

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

double scalar(R *u, R *v, int n)
{ 
   R res=0;
   for(int i=0; i<n; ++i)
      res += u[i]*v[i]; 
   return res;
}

int main(int argc, const char ** argv)
{   
    cout << " Usage: " << argv[0] << " file.msh perio\n";
    assert(argc>2);
    Maillage2d Th(argv[1]);    
    int typeperio=atoi(argv[2]); 
    assert(abs(typeperio)==1);   
    int nv  = Th.nv;
    int * dl = new int[nv];
    R * coef = new R[nv];
    vector<int> t1; 
    vector<int> t2; 
    for(int i=0; i< nv; ++i)
    { 
	    if(Th.v[i].lab==3) //bord inferieur
           t1.push_back(i); 

        if(Th.v[i].lab==2) //bord superieur  
           t2.push_back(i);
     }
  
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
    unsigned int s=t1.size();
    int Corr[s][2];
    for(i=0; i< s; ++i)
    {
        Corr[i][0] = t1[i];
        Corr[i][1] = t2[i];
    }

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
    }

    // calcul des aires des regions
    int nb_regions= 8;
    int regions[8] = {0, 1, 2, 3, 19, 20, 21, 22};
    R  aires [8] = {0., 0., 0., 0., 0., 0., 0., 0.}, aire_totale=0.;
    for (int i=0; i <nb_regions; ++i)
    {
        for(int k=0; k<Th.nt;++k)
        {  
             const Triangle & K=Th[k];
            if (K.lab == regions[i])
            {const R mes =K.mes;
             aires[i]+= mes;}
        }
        aire_totale += aires[i] ;
    }
  
    int n = dl[nv-1]+1;
    int nt = Th.nt ;
    R epsilon = 5.1636e-04, kappa = 80.; 
    R sigma= 5*pow(10,7), Y11= 1, Y22= 1e-03, I01= 15, f= 50.;// Ri=1 ou 1000 ohm?, Ra=1 ohm et Va=15 volts
    R dt= 0.002, T= 5./f; // cinq periodes 

    // construction du tableau beta
	R **beta;
	beta = new R* [nt];
	for (int i=0; i < nt; ++i)
	   beta[i] = new R[4];

	for (int i=0; i < nt; ++i)
	  for (int j=0; j < 4; ++j)
	     beta[i][j] = 0;

   for(int k =0; k< nt; ++k)
   {
      const Triangle & K=Th[k];
      if ( (K.lab == 1) || (K.lab == 2)) 
         beta[k][0] = epsilon/(4*M_PI*1e-07); 
      else 
         beta[k][0] = 1./(4*M_PI*1e-07);
          
      if (K.lab == 3)
         beta[k][1]  = sigma/dt;    
    
      if (K.lab == 19 || K.lab == 21 )
         beta[k][2] = kappa*Y11*100*100/(aires[4]+aires[6])/(aires[4]+aires[6])/dt ;        

      if (K.lab == 20 || K.lab == 22 )
         beta[k][3] = kappa*Y22*60*60/(aires[5]+aires[7])/(aires[5]+aires[7])/dt ;

   }

   MatrixLap A(Th,dl,n,coef,0,beta); 

   // beta nul pour le second membre
   R **beta_nul;
	 beta_nul = new R* [nt];
   for (int i=0; i < nt; ++i)
	 beta_nul[i] = new R[4];
   for (int i=0; i < nt; ++i)
   for (int j=0; j < 4; ++j)
	   beta_nul[i][j] = 0;
   
   MatrixLap M(Th,dl,n,coef,1,beta_nul); 

   R * fh  = new R[n];
   R * u   = new R[n];
   R * un  = new R[nv];
   R * b   = new R[n];
   R * b0  = new R[n];
   R * b1  = new R[n];
   R * u0  = new R[n]; // old solution periodic
   R * u0n = new R[nv]; // old solution non periodic

   for (int i=0; i< n;++i)
      fh[i]=0;
   
   for(int k=0; k<Th.nt; ++k)
   {  
      const Triangle & K=Th[k]; 
      for(int i=0; i<K.nbv; ++i)     
      {
          int in = dl[Th(K[i])]; 
          if (K.lab == 19 || K .lab ==21)
                fh[in] = I01*100/(aires[4]+aires[6]); 
      }
   }

   M.MatMul(fh,b0); // calcul du terme constant du second membre

   R **new_beta; //pour le calcul du deuxeme terme du second membre
   new_beta = new R* [nt];
   for(int i=0; i < nt; ++i)
	  new_beta[i] = new R[4];
   for(int i=0; i < nt; ++i)
   for(int j=0; j < 4; ++j)
   {   
        if (j == 0)
	         new_beta[i][j] = 0;
        else
           new_beta[i][j] = beta[i][j];
    }

   MatrixLap M1(Th,dl,n,coef,1,new_beta);

    //pour le calcul du flux
    R * f1 = new R[nv];
    R * f2 = new R[nv]; 
    for (int i=0; i< nv; ++i)
    {   f1[i]=0; f2[i]=0; }
 
    for(int k=0; k<Th.nt; ++k)
    {  
          const Triangle & K=Th[k]; 
          const int  d = Triangle::d;
          const R mes =K.mes ;
          for(int i=0; i<K.nbv; ++i)     
          {  
             int num = Th(K[i]) ;
             if (K.lab == 19 || K .lab ==21)
                f1[num] = 100*mes/(aires[4]+aires[6])/(d+1); 

             if (K.lab == 20 || K .lab ==22)
                f2[num] = 60*mes/(aires[5]+aires[7])/(d+1); 

          }
   }

   /*ifstream fi("sol.txt");
   assert( fi); 
   for(int i=0; i<n; ++i)
   { 
       fi >> u0[i] ; 
       assert( fi.good());
   }*/

   for(int i=0; i<n; ++i)
        u0[i] =0. ;
   for(int i=0; i<nv;++i)
        u0n[i] = coef[i]*u0n[dl[i]];

   int n_temps = (int) (T/dt) ;

   R * i1 = new R [n_temps];
   R * i2 = new R [n_temps];
   R eps = 1e-7;

   for(R t=0.; t<= T ; t+=dt)
   {
      cout << "t= " << t << " s ; " ;
      M1.MatMul(u0,b1);
      for(int i=0; i<n;++i)
         b[i] = b0[i]+b1[i];  //second membre pour le pb instationnaire  
        
      int ok = CG(A,MatId(n),b,u,n,eps);
      if (ok)
         cout << "CG a converge pour epsilon " << eps << endl;
      else 
         cout << "CG n'a pas converge pour epsilon " << eps << endl;
    
     for(int i=0; i<nv;++i)
         un[i] = coef[i]*u[dl[i]]; 

     R flux1_old= kappa*scalar(f1,u0n,nv);
     R flux2_old= kappa*scalar(f2,u0n,nv);
     R flux1= kappa*scalar(f1,un,nv);
     R flux2= kappa*scalar(f2,un,nv);
   
     int idx = (int) (t/dt);
     i1[idx] = -Y11*( flux1 - flux1_old)/dt + I01 ; 
     i2[idx] = -Y22*( flux2 - flux2_old)/dt ;
   
     for(int i=0; i<nv;++i)
     {
       u0[i]  = u[i];
       u0n[i] = un[i];
     }


   } // fin boucle temps

   // pour ploter les intensites
   ofstream f0("intensites.gp");
   for(int i=0; i< n_temps; ++i)
   {
      R t = (i+1)*dt ;
      f0 <<  t <<  " " << i1[i] << " " << i2[i] << endl ;
   }
         
   // pour ploter le vecteur potentiel final
   ofstream ff1("m_inst.gp");
   for(int k=0;k<Th.nt;++k)
   { 
      for (int ip=0;ip< 4 ;++ ip)
      {
         int i= Th(k,ip%3);
         ff1 << Th.v[i].x << " " <<Th.v[i].y << " " << un[i] << endl ;
      }
         ff1 << endl<< endl;
   }

   // clean memory ...
   delete [] b;
   delete [] b0;
   delete [] u;
   delete [] un;
   delete [] fh;
   delete [] dl;
   delete [] coef;
   delete [] b1;
   delete [] i1;
   delete [] i2;
   delete [] u0;
   delete [] u0n;

   for (int i=0; i < nt; i++)
	  delete [] beta[i];
   delete [] beta; 
 
   for (int i=0; i < nt; i++)
	  delete [] beta_nul[i];
   delete [] beta_nul; 

   for (int i=0; i < nt; i++)
	  delete [] new_beta[i];
   delete [] new_beta; 

}
