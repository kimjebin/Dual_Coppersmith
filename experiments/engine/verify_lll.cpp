// Streaming exact validation of Bareiss LLL data and coefficient bounds.
#include <gmpxx.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using Z=mpz_class;
static Z ceil_div(const Z &x,const Z &y) {
  Z r; mpz_cdiv_q(r.get_mpz_t(),x.get_mpz_t(),y.get_mpz_t()); return r;
}
static Z ceil_sqrt(const Z &x) {
  Z r; mpz_sqrt(r.get_mpz_t(),x.get_mpz_t()); if(r*r<x) ++r; return r;
}
int main(int argc,char **argv) {
  try {
    if(argc!=3) throw std::runtime_error("Usage: verify_lll PREPARATION INITIAL_BEST");
    std::ifstream in(argv[1]);
    std::string magic; int d; Z a,c;
    in>>magic>>d>>a>>c;
    if(magic!="EXACT-GMP-PREP-1" || d<1 || a<=0 || c<=0)
      throw std::runtime_error("Bad preparation header");
    Z ignored;
    for(int i=0;i<d;++i) in>>ignored;
    std::vector<Z> delta(d+1);
    for(auto &x:delta) { in>>x; if(x<=0) throw std::runtime_error("Nonpositive minor"); }
    if(delta[0]!=1) throw std::runtime_error("Bad initial minor");
    const Z denominator=Z(1)<<64;
    std::vector<std::vector<unsigned long>> mu(d);
    for(int i=1;i<d;++i) {
      mu[i].resize(i);
      Z x;
      for(int j=0;j<i;++j) {
        in>>x;
        Z absolute=abs(x);
        if(100*absolute>51*delta[j+1]) throw std::runtime_error("Size-reduction precondition fails");
        Z upper=ceil_div(absolute*denominator,delta[j+1]);
        if(!upper.fits_ulong_p()) throw std::runtime_error("mu bound overflow");
        mu[i][j]=upper.get_ui();
      }
      if(100*(delta[i+1]*delta[i-1]+x*x)<99*delta[i]*delta[i])
        throw std::runtime_error("Lovasz precondition fails");
    }
    if(!in) throw std::runtime_error("Truncated preparation");
    const Z initial(argv[2]);
    if(initial<=0) throw std::runtime_error("Bad incumbent");
    Z radius; Z numerator=a*(initial-1);
    mpz_fdiv_q(radius.get_mpz_t(),numerator.get_mpz_t(),c.get_mpz_t());
    std::vector<Z> bound(d);
    Z largest=0;
    for(int i=d-1;i>=0;--i) {
      Z weighted=0;
      for(int j=i+1;j<d;++j) weighted+=Z(mu[j][i])*bound[j];
      bound[i]=ceil_sqrt(ceil_div(radius*delta[i],delta[i+1]))+ceil_div(weighted,denominator);
      if(bound[i]>largest) largest=bound[i];
    }
    if(largest>=(Z(1)<<50)) throw std::runtime_error("Integer coefficient bound exceeds 2^50");
    std::cout<<"LLL_OK "<<d<<' '<<largest<<'\n';
  } catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 2; }
}
