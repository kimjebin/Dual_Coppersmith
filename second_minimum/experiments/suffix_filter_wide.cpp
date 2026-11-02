// Full-dimensional exact dyadic intervals; the older24-coordinate source is
// retained unchanged so existing certificates remain independently replayable.
#include <gmpxx.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using Z=mpz_class;
static Z floordiv(const Z &a,const Z &b) {
  Z out;mpz_fdiv_q(out.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return out;
}
int main(int argc,char **argv) {
  try {
    if(argc!=3) throw std::runtime_error("Usage: suffix_filter_wide PREPARATION INITIAL_BEST");
    std::ifstream in(argv[1]);std::string magic;int d;Z a,c,x;
    in>>magic>>d>>a>>c;
    if(magic!="EXACT-GMP-PREP-1" || d<1 || a<=0 || c<=0) throw std::runtime_error("Bad header");
    for(int i=0;i<d;++i) in>>x;
    std::vector<Z> delta(d+1);for(auto &z:delta) {in>>z;if(z<=0) throw std::runtime_error("Bad minor");}
    const Z denominator=Z(1)<<64;
    std::vector<std::vector<Z>> lower(d);
    for(int i=1;i<d;++i) {
      lower[i].resize(i);
      for(int j=0;j<i;++j) {in>>x;lower[i][j]=floordiv(x*denominator,delta[j+1]);}
    }
    if(!in) throw std::runtime_error("Truncated preparation");
    const Z initial(argv[2]);if(initial<=0) throw std::runtime_error("Bad incumbent");
    const Z radius=floordiv(a*(initial-1),c);
    std::cout<<"SUFFIX-BOUNDS-1 "<<d<<" 0 64 128\n";
    for(int i=0;i<d;++i) {
      Z r=radius==0 ? Z(0) : floordiv(delta[i+1]*(Z(1)<<128),delta[i]*radius);
      std::cout<<r<<' ';
    }
    std::cout<<'\n';
    for(const auto &row:lower) {for(const auto &z:row) std::cout<<z<<' ';std::cout<<'\n';}
  } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 2;}
}
