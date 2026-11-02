// Recompute exact fraction-free Schur complements from the independent Gram
// matrix, comparing every stored Bareiss pivot and off-diagonal coefficient.
#include <gmpxx.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using Z=mpz_class;
int main(int argc,char **argv) {
  try {
    if(argc!=3) throw std::runtime_error("Usage: verify_gram PREPARATION GRAM_INPUT");
    std::ifstream prep(argv[1]),input(argv[2]);
    std::string magic1,magic2;int d,n;Z a,c,aa,cc,initial,x;
    prep>>magic1>>d>>a>>c;input>>magic2>>n>>aa>>cc>>initial;
    if(magic1!="EXACT-GMP-PREP-1" || magic2!="FAST-LIFT-1" || d<1 || d!=n || a!=aa || c!=cc)
      throw std::runtime_error("Header mismatch");
    for(int i=0;i<d;++i) {Z p,q;prep>>p;input>>q;if(p!=q) throw std::runtime_error("Inner mismatch");}
    for(int i=0;i<d;++i) input>>x;
    std::vector<Z> delta(d+1);for(auto &z:delta) prep>>z;
    if(delta[0]!=1) throw std::runtime_error("Bad first minor");
    std::vector<std::vector<Z>> lower(d),H(d,std::vector<Z>(d));
    for(int i=1;i<d;++i) {lower[i].resize(i);for(auto &z:lower[i]) prep>>z;}
    for(auto &row:H) for(auto &z:row) input>>z;
    if(!prep || !input) throw std::runtime_error("Truncated input");
    for(int i=0;i<d;++i) for(int j=0;j<i;++j)
      if(H[i][j]!=H[j][i]) throw std::runtime_error("Nonsymmetric Gram matrix");
    for(int k=0;k<d;++k) {
      const Z pivot=H[k][k];
      if(pivot<=0 || pivot!=delta[k+1]) throw std::runtime_error("Bareiss pivot mismatch");
      for(int i=k+1;i<d;++i)
        if(H[i][k]!=lower[i][k]) throw std::runtime_error("Bareiss lower coefficient mismatch");
      for(int i=k+1;i<d;++i) for(int j=i;j<d;++j) {
        Z numerator=pivot*H[i][j]-H[i][k]*H[k][j];
        if(!mpz_divisible_p(numerator.get_mpz_t(),delta[k].get_mpz_t()))
          throw std::runtime_error("Nonintegral Schur complement");
        mpz_divexact(H[i][j].get_mpz_t(),numerator.get_mpz_t(),delta[k].get_mpz_t());
        H[j][i]=H[i][j];
      }
    }
    std::cout<<"GRAM_OK "<<d<<'\n';
  } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 2;}
}
