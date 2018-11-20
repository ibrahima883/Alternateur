
 struct  MatVirt
 {
   int n,m;
  virtual double *  AddMatMul(const double *x, double *Ax) const =0;
  double * MatMul(const double *x, double *Ax) const
     {//std::cout << "bla" <<std::endl;
         for (int i=0; i<n; ++i)
             Ax[i]=0.;
         return AddMatMul(x,Ax);
     }
 
     MatVirt(int nn): n(nn), m(nn) {}
     MatVirt(int nn,int mm): n(nn), m(mm) {}
     virtual ~MatVirt() {}
};


int CG( const MatVirt &A,
		     const MatVirt &C,
                     double *b, // second membre
		     double * x, // solution qui contient une initialisation
		     int nbitermax,double eps)
;

struct  MatId : public MatVirt {
    MatId(int n) :MatVirt(n) {}
    double *  AddMatMul(const double *x, double *Ax) const
    {
        for (int i=0; i<n; ++i)
            Ax[i]+=x[i];
        return Ax;
    }
    
};
