// Reuse the exact streaming validator, printing every certified box bound.
// This file includes its implementation in a separate namespace: the original
// verifier and its source hash remain unchanged.
#include <sstream>
#define main original_verifier_main
#include "verify_lll.cpp"
#undef main

int main(int argc,char **argv) {
  try {
    if(argc!=3) throw std::runtime_error("Usage: partition_bounds PREPARATION INITIAL_BEST");
    std::ifstream in(argv[1]);
    std::string magic; int d; Z a,c,ignored;
    in>>magic>>d>>a>>c;
    if(magic!="EXACT-GMP-PREP-1" || d<1 || a<=0 || c<=0)
      throw std::runtime_error("Bad preparation header");
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
        if(100*abs(x)>51*delta[j+1]) throw std::runtime_error("Size-reduction fails");
        Z upper=ceil_div(abs(x)*denominator,delta[j+1]);
        if(!upper.fits_ulong_p()) throw std::runtime_error("mu overflow");
        mu[i][j]=upper.get_ui();
      }
      if(100*(delta[i+1]*delta[i-1]+x*x)<99*delta[i]*delta[i])
        throw std::runtime_error("Lovasz fails");
    }
    if(!in) throw std::runtime_error("Truncated preparation");
    const Z initial(argv[2]); Z radius,numerator=a*(initial-1);
    if(initial<=0) throw std::runtime_error("Bad incumbent");
    mpz_fdiv_q(radius.get_mpz_t(),numerator.get_mpz_t(),c.get_mpz_t());
    std::vector<Z> bound(d);
    for(int i=d-1;i>=0;--i) {
      Z weighted=0;
      for(int j=i+1;j<d;++j) weighted+=Z(mu[j][i])*bound[j];
      bound[i]=ceil_sqrt(ceil_div(radius*delta[i],delta[i+1]))+ceil_div(weighted,denominator);
    }
    Z largest=0;
    for(const auto &b:bound) if(b>largest) largest=b;
    if(largest>=(Z(1)<<50)) throw std::runtime_error("Coefficient bound exceeds 2^50");
    std::cout<<"LLL_OK "<<d<<' '<<largest<<'\n'<<"BOUNDS";
    for(const auto &b:bound) std::cout<<' '<<b;
    std::cout<<'\n';
  } catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 2; }
}
