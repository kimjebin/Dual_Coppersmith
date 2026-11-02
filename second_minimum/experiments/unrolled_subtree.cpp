// Experimental fixed-suffix adapter to MIT-licensed fplll 5.5.0 enumlib.
// Not used as certificate evidence until independent validation succeeds.
#define main unused_lift_main
#include "fast_lift.cpp"
#undef main
#include "generated/fixed_enumlib.h"
#include <cstdint>
#include <cmath>
#ifndef PROOF_DIM
#error PROOF_DIM is required
#endif
constexpr int dimension=PROOF_DIM;
using Kernel=fplll::fixed_enumlib::lattice_enum_t<dimension,1,1,1,false>;
using Globals=fplll::fixed_enumlib::globals_t<dimension>;
template<int level> void enter(Kernel &kernel,int free_dimension) {
  if(free_dimension==level+1) kernel.enumerate_recur(typename Kernel::template i_tag<level,true,-2,0>());
  else if constexpr(level>0) enter<level-1>(kernel,free_dimension);
  else throw std::runtime_error("Bad free dimension");
}
int main(int argc,char **argv) {
  try {
    if(argc<2 || argc-2>dimension) throw std::runtime_error("Usage: unrolled_subtree INPUT [suffix]");
    std::ifstream input(argv[1]);std::string magic;int d;Z a,c,best;
    input>>magic>>d>>a>>c>>best;
    if(magic!="FAST-LIFT-1" || d!=dimension || a<=0 || c<=0 || best<=0) throw std::runtime_error("Bad input");
    std::vector<Z> inner(d),initial(d);for(auto &z:inner) input>>z;for(auto &z:initial) input>>z;
    double rho;int precision=std::max(128,gso_min_prec(rho,d,LLL_DEF_DELTA,LLL_DEF_ETA)+32);
    FF::set_prec(precision);Matrix<FZ> H;H.resize(d,d);
    for(int i=0;i<d;++i) for(int j=0;j<d;++j) {Z z;input>>z;mpz_set(H[i][j].get_data(),z.get_mpz_t());}
    if(!input) throw std::runtime_error("Truncated input");
    Matrix<FZ> empty1,empty2;MatGSOGram<FZ,FF> gso(H,empty1,empty2,GSO_INT_GRAM);
    if(!gso.update_gso()) throw std::runtime_error("GSO failed");
    LiftEvaluator evaluator(d,gso,H,a,c,best,inner,initial);evaluator.init_delta_def(precision,rho,true);
    FF radius=evaluator.expanded_radius(),value;long normalizer=-1,exponent;
    for(int i=0;i<d;++i) {value=gso.get_r_exp(i,i,exponent);normalizer=std::max(normalizer,exponent+value.exponent());}
    evaluator.set_normexp(normalizer);radius.mul_2si(radius,-normalizer);
    Globals global;global.A=radius.get_d(GMP_RNDU);
    for(auto &signal:global.signal) signal=0;
    global.process_sol=[&](double distance,double *coordinates) {
      std::vector<FF> exact_coordinates(d);for(int i=0;i<d;++i) exact_coordinates[i]=coordinates[i];
      double bound=global.A.load();evaluator.eval_sol(exact_coordinates,distance,bound);return bound;
    };
    auto kernel=std::make_unique<Kernel>(global);kernel->threadid=0;kernel->_A=global.A;
    for(int i=0;i<d;++i) {
      kernel->pr[i]=kernel->pr2[i]=1.;kernel->_x[i]=kernel->_Dx[i]=0;kernel->_D2x[i]=-1;
      kernel->_r[i]=d-1;kernel->_l[i]=0;kernel->_counts[i]=0;
      value=gso.get_r_exp(i,i,exponent);value.mul_2si(value,exponent-normalizer);kernel->risq[i]=value.get_d();
      if(!(kernel->risq[i]>0) || !std::isfinite(kernel->risq[i])) throw std::runtime_error("Unsupported diagonal range");
      for(int j=0;j<d;++j) {kernel->_sigT[i][j]=0;kernel->muT[i][j]=0;}
      for(int j=i+1;j<d;++j) {gso.get_mu(value,j,i);kernel->muT[i][j]=value.get_d();}
    }
    kernel->_l[d]=0;kernel->_counts[d]=0;kernel->_update_AA();
    const int free_dimension=d-(argc-2);
    for(int j=free_dimension;j<d;++j) {
      Z z(argv[j-free_dimension+2]);if(abs(z)>=(Z(1)<<50)) throw std::runtime_error("Suffix range");
      kernel->_x[j]=z.get_si();
    }
    for(int i=d-1;i>=free_dimension;--i) {
      double center=0;for(int j=d-1;j>i;--j) center-=kernel->_x[j]*kernel->muT[i][j];
      const double residue=center-kernel->_x[i];
      kernel->_l[i]=kernel->_l[i+1]+residue*residue*kernel->risq[i];
    }
    auto began=std::chrono::steady_clock::now();
    if(free_dimension>0) {
      const int i=free_dimension-1;double center=0;
      for(int j=d-1;j>i;--j) center-=kernel->_x[j]*kernel->muT[i][j];
      kernel->_sigT[i][i]=center;
      enter<dimension-1>(*kernel,free_dimension);
    } else {
      bool nonzero=false;std::vector<FF> coordinates(d);
      for(int i=0;i<d;++i) {coordinates[i]=static_cast<double>(kernel->_x[i]);nonzero|=kernel->_x[i]!=0;}
      if(nonzero) {double dummy=global.A;evaluator.eval_sol(coordinates,kernel->_l[0],dummy);}
    }
    Z nodes=0;for(int i=0;i<d;++i) nodes+=Z(std::to_string(kernel->_counts[i]));
    auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
    std::cout<<"COMPLETE "<<evaluator.best<<' '<<nodes<<' '<<evaluator.leaves<<' '<<seconds<<' '<<precision<<'\n';
    for(const auto &z:evaluator.best_y) std::cout<<z<<' ';std::cout<<'\n';
  } catch(const std::exception &error) {std::cerr<<error.what()<<'\n';return 2;}
}
