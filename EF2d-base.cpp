#include "EF2d-base.hpp"
#include <fstream>

Maillage2d::Maillage2d(const char * filename)
{
  int i0, i1, i2, ir;
  std::ifstream  f(filename); 
  assert( f); 
  int unused; //I[4] ; 
  f >> nv >> nt >> unused ;
  assert( f.good());
  t = new Triangle[nt];
  v = new Sommet[nv];
  assert( t && v); 
  double mes =0; 
  for(int i=0;i<nv;++i)
    { 
      f >> v[i] ; 
      assert( f.good());
    }

  for(int k=0;k<nt;++k)
    { 
      //for(int i=0;i< 4; ++i)
	//f >> I[i] ; 
        f >> i0 >> i1 >> i2 >> ir ;
      //assert( f.good());
      assert( f.good() && i0 >0 && i0 <=nv && i1 >0 && i1 <=nv && i2 >0 && i2 <=nv);
      t[k].set(v,i0-1,i1-1,i2-1,ir);
      //t[k].build(v,I,-1);
      mes += t[k].mes; 
    }
  std::cout<< " End read " << nv << " " << nt << " mes =" << mes << std::endl; 
}
Maillage2d::Maillage2d(const char * filename1,const char * filename2)
{
    std::ifstream  f1(filename1);
    std::ifstream  f2(filename2);
    assert( f2 && f1);
    int unused, I[4] ;
    int nv1,nt1,nv2,nt2;
    f1 >> nv1 >> nt1 >> unused ;
    f2 >> nv2 >> nt2 >> unused ;
    nt = nt1+nt2;
    nv =nv1+nv2;
    assert( f1.good() && f2.good());
    t = new Triangle[nt];
    v = new Sommet[nv];
    assert( t && v);
    double mes =0;
    for(int i=0;i<nv1;++i)
    {
        f1 >> v[i] ;
        assert( f1.good());
    }
    for(int i=0;i<nv2;++i)
    {
        f2 >> v[nv1+i] ;
        assert( f2.good());
    }
    
    for(int k=0;k<nt1;++k)
    {
        for(int i=0;i< 4; ++i)
            f1 >> I[i] ;
        assert( f1.good());
        t[k].build(v,I,-1);
        mes += t[k].mes;
    }
    for(int k=0;k<nt2;++k)
    {
        for(int i=0;i< 4; ++i)
            f2 >> I[i] ;
        assert( f2.good());
        t[k+nt1].build(v,I,-1+nv1);
        mes += t[k].mes;
    }
    std::cout<< " End read " << nv << " " << nt << " mes =" << mes << std::endl;
}
