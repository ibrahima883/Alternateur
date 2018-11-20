#include "R2.hpp"
#include "clock.hpp"
#include <cassert>

class Label {
public:
  int lab; 
  int  OnGamma() const { return lab;}
  Label(int l=0) : lab(l) {} 
};

class Sommet :public R2,public  Label
{
 public:
  Sommet() {}   
};

inline std::ostream& operator <<(std::ostream& f, const  Sommet & P )
{ return  f << P.x << ' ' << P.y  << ' ' << P.lab << ' ' ;}
inline  std::istream& operator >>(std::istream& f,  Sommet & P)
{ return f >>  P.x >>  P.y >> P.lab ;  }

class Triangle :public  Label {
public:
    
  static const int nbv =3, d=2;
  Sommet * v[nbv]; 
  double mes;

  void set(Sommet *v0, int i0, int i1, int i2, int r){
    v[0] = v0+i0; v[1] = v0+i1; v[2] = v0+i2;
    lab =r;
    mes = det(*v[0], *v[1], *v[2]) * 0.5; 
    assert(mes>0) ;
   }

  Triangle(){ (v[0]=(v[1]=(v[2]=0)));}
  void build(Sommet *v0,int * I,int offset=0)
  {// I array of vertex number  
    for(int i=0; i < nbv; ++i) 
      v[i] =  v0 + I[i]+offset; 
    mes = det(*v[0], *v[1], *v[2]) * 0.5; 
    assert(mes>0) ;
  } 
  void GradLambdaK(R2 *G) const
  {
    double K2 = mes*2; 
    G[0] = R2(*v[1],*v[2]).perp()/K2;
    G[1] = R2(*v[2],*v[0]).perp()/K2;
    G[2] = R2(*v[0],*v[1]).perp()/K2;
  }
  
  Sommet & operator[](int i) { assert(i>=0 && i < nbv); return *(v[i]); }
  const Sommet & operator[](int i) const { assert(i>=0 && i < nbv); return *(v[i]);}
   
  const R2 &P(int i) const { return *v[i];}
    
  R2 operator()(R2 Phat) const { return (1-Phat.x-Phat.y)* P(0) + Phat.x* P(1) +  Phat.y* P(2) ; }

}; 

inline std::ostream& operator <<(std::ostream& f, const   Triangle & K )
{ return  f << K[0] << ' ' << K[1]  << ' ' << K[2]  << ' '<< K.lab << ' ' ;}
inline  std::istream& operator >>(std::istream& f,  Triangle & K)
{ return f >>  K[0] >>  K[1] >> K[2] >> K.lab ;  } 

class Maillage2d 
{
public:
  int nv,nt; 
  Sommet * v; 
  Triangle *t;
  Maillage2d(const char *  filename); 
  Maillage2d(const char *,const char *  );
  ~Maillage2d() { delete [] v; delete [] t; }
  // destuctor => careful with copie operator  
  // no copy operator
  // chech index number
  int CheckV(int i) const { assert( i>=0 && i < nv); return i; } 
  int CheckT(int i) const { assert( i>=0 && i < nt); return i; } 
  int operator()(const Sommet & vv) const { return CheckV(&vv-v);}
  int operator()(const  Triangle & tt) const  { return CheckT(&tt-t);}
  int operator()(const Sommet * vv)const  { return CheckV(vv-v);}  // (1)
  int operator()(const  Triangle * tt) const { return CheckT(tt-t);}
  Triangle & operator[](int k) { return t[CheckT(k)]; }
  const Triangle & operator[](int k) const { return t[CheckT(k)]; }
  int  operator()(int k,int i) const { return  operator()(t[k].v[i]); }// call (1)

private:
  Maillage2d(const Maillage2d &);
  Maillage2d & operator=(const Maillage2d &);
 
};
