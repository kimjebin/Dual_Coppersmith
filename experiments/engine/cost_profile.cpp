// Heuristic full-sphere node cost from exact saved GSO minors. This program
// produces NO minimum or lower-bound certificate and never prunes a search.
#include <gmpxx.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using Z=mpz_class;
static long double logz(const Z &value) {
  if(value<=0) throw std::runtime_error("Nonpositive logarithm");
  long bits=mpz_sizeinbase(value.get_mpz_t(),2);
  Z top=value;if(bits>53) mpz_fdiv_q_2exp(top.get_mpz_t(),value.get_mpz_t(),bits-53);
  return std::log2(static_cast<long double>(top.get_d()))+std::max(0L,bits-53);
}
int main(int argc,char **argv) {
  try {
    if(argc!=3) throw std::runtime_error("Usage: cost_profile PREPARATION INITIAL_BEST");
    std::ifstream input(argv[1]);std::string magic;int d;Z a,c,x;
    input>>magic>>d>>a>>c;
    if(magic!="EXACT-GMP-PREP-1" || d<1) throw std::runtime_error("Bad header");
    for(int i=0;i<d;++i) input>>x;
    std::vector<Z> delta(d+1);for(auto &z:delta) input>>z;
    Z initial(argv[2]);Z radius=a*(initial-1)/c;
    if(radius<=0) {std::cout<<"TRIVIAL "<<d<<'\n';return 0;}
    long double logradius=logz(radius),sum=0,maximum=-1e100L;
    std::vector<long double> costs,diagonal;
    for(int i=0;i<d;++i) diagonal.push_back(logz(delta[i+1])-logz(delta[i]));
    for(int k=1;k<=d;++k) {
      sum+=diagonal[d-k];
      long double logball=(k/2.0L)*std::log2(std::acos(-1.0L))-
          std::lgamma(k/2.0L+1)/std::log(2.0L);
      long double cost=logball+(k/2.0L)*logradius-sum/2-1;
      costs.push_back(cost);maximum=std::max(maximum,cost);
    }
    long double scaled=0;for(auto cost:costs) scaled+=std::exp2(cost-maximum);
    long double estimate=maximum+std::log2(scaled);
    std::cout<<std::setprecision(16)<<"HEURISTIC "<<d<<' '<<estimate<<' '<<maximum<<'\n';
  } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 2;}
}
