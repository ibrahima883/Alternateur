#include "EF2d-base.hpp"
#include <fstream>
#include "CG.hpp"
#include <cmath>
#include <typeinfo>
#include <vector>

using namespace std;


typedef double R;


int main(int argc, const char ** argv)
{   
    cout << " Utilisation du code: " << argv[0] << " rotor.msh stator.msh\n";
    assert(argc>2);   
    Maillage2d Thm(argv[1],argv[2]);
    int nv  = Thm.nv;
    Maillage2d Thr(argv[1]);    
    Maillage2d Ths(argv[2]); 
    int nvr  = Thr.nv;
    int * dlr = new int[nvr];
    R * coefr = new R[nvr];
    int nvs  = Ths.nv;
    int * dls = new int[nvs];
    R * coefs = new R[nvs];

    vector <int> ligne; 
    for(int i=0; i< nv; ++i)
	if(Thm.v[i].lab==1)  //ligne de glissement
           ligne.push_back(i); 

    //for(unsigned i=0; i< ligne.size(); ++i)	
    // cout << ligne[i] << endl;

    unsigned int s=ligne.size()/2;

    for(int i=0; i< nv; ++i)
    {  
        if (i < nvr )
        {   dlr[i]=i; coefr[i]=1; }
        else 
        {   dls[i-nvr]=i; coefs[i-nvr]=1; }
    }

    R dt=1./19200. , T=16*dt ; // v=50 tours/s et 1 tour= 360° --> dt = (15/16)/(50*360) s
 
    R Rx[nvr] ;
    R Ry[nvr] ;

    for(R t=0.; t<= T ; t+=dt)
    {
       int l = (int) (t/dt) ;
       //cout << "Sommets du rotor dans la zone negative apres " << l << " pas de temps: " << endl;
       cout << "Sommets de la ligne de glissement nons communs apres " << l << " pas de temps: " << endl;
       int count=0;
       R phi = l*15*M_PI/16./180.; // angle en radian

       //rotation du rotor
       for(int i=0; i< nvr; ++i)
       {     
          double pente = 0.;
          R x = Thm.v[i].x, y = Thm.v[i].y ;
          Rx[i] = x*cos(phi) - y*sin(phi) ;
          Ry[i] = x*sin(phi) + y*cos(phi) ;
          if (Rx[i]) { pente = Ry[i]/Rx[i]; }
          if ((Thm.v[i].lab==1) & (abs(atan(pente)) - (M_PI/12.) > 1e-06))
          {  
              coefr[i] = -1;
              cout << i << " (rotor) ; coef= " << coefr[i] << endl;
          }
       }

       cout << endl ;
       cout << "Les points communs apres " << l << " pas de temps: " << endl;
       for(int i=nvr; i< nv; ++i)
       {   
          dls[i-nvr]= i-count;
          if ( Thm.v[i].lab == 1 )          
             for(unsigned j=0; j<s ; ++j)
             {
                int q = ligne[j] ; 
                if ( (abs(Rx[q] -Thm.v[i].x) < 1e-04 ) &  (abs(Ry[q] - Thm.v[i].y) < 1e-04)  )
                {      
                   dls[i-nvr] = dlr[q]; 
                   cout << i << " (global) <--> " << i-nvr << " (stator) <--> " << dls[i-nvr]  << " (rotor)" << endl;
                   count++;
                }
             }  
       } 

    cout << endl;    
 
    } // fin boucle temps 

    //for (int i=0; i< nvs; ++i)
    // cout << i << " ; " << dls[i] << endl;

   // clean memory ...
    delete [] dlr;
    delete [] coefr;
    delete [] dls;
    delete [] coefs;

}
