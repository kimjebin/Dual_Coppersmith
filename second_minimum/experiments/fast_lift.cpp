// Error-bounded exhaustive enumeration with exact nearest-integer lifts.
// Uses fplll 5.5.0's public proved-enumeration error-bound API.
// Reference: fplll/svpcvp.cpp and enum/evaluator.h (LGPL-2.1-or-later).
#include <fplll/fplll.h>
#include <fplll/gso_gram.h>
#include <fplll/enum/enumerate.h>
#include <fplll/enum/evaluator.h>
#include <gmpxx.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace fplll;
using Z = mpz_class;
using FZ = Z_NR<mpz_t>;
using FF = FP_NR<mpfr_t>;

static Z floor_div(const Z &a,const Z &b) {
  Z q; mpz_fdiv_q(q.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t()); return q;
}

class LiftEvaluator : public ExactErrorBoundedEvaluator {
 public:
  const Matrix<FZ> &H;
  Z a,c,best;
  std::vector<Z> inner, best_y;
  unsigned long long leaves=0;
  LiftEvaluator(int d, MatGSOInterface<FZ,FF> &gso, const Matrix<FZ> &gram,
                const Z &aa,const Z &cc,const Z &bb,
                std::vector<Z> dd,std::vector<Z> initial)
      : ExactErrorBoundedEvaluator(d,gso,EVALMODE_COUNT),H(gram),a(aa),c(cc),best(bb),
        inner(std::move(dd)),best_y(std::move(initial)) { update_radius(); }
  void update_radius() {
    Z radius=floor_div(a*(best-1),c);
    mpz_set(int_max_dist.get_data(),radius.get_mpz_t());
  }
  FF expanded_radius() {
    FF r,error;
    r.set_z(int_max_dist,GMP_RNDU);
    if (!get_max_error_aux(r,true,error))
      throw std::runtime_error("fplll could not certify sphere error bound");
    r.add(r,error,GMP_RNDU);
    return r;
  }
  void eval_sol(const std::vector<FF> &coordinates,const enumf &,enumf &bound) override {
    ++leaves;
    std::vector<Z> y(d);
    FZ integer;
    bool nonzero=false;
    for(int i=0;i<d;++i) {
      integer.set_f(coordinates[i]);
      mpz_set(y[i].get_mpz_t(),integer.get_data());
      nonzero=nonzero || y[i]!=0;
    }
    if(!nonzero) throw std::runtime_error("Zero quotient candidate");
    Z q=0,dot=0;
    for(int i=0;i<d;++i) {
      Z entry;
      mpz_set(entry.get_mpz_t(),H[i][i].get_data());
      q+=y[i]*y[i]*entry;
      dot+=y[i]*inner[i];
      for(int j=0;j<i;++j) {
        mpz_set(entry.get_mpz_t(),H[i][j].get_data());
        q+=2*y[i]*y[j]*entry;
      }
    }
    const Z k=floor_div(2*dot+a,2*a);
    const Z residue=dot-k*a;
    Z numerator=c*q+residue*residue;
    if(!mpz_divisible_p(numerator.get_mpz_t(),a.get_mpz_t()))
      throw std::runtime_error("Nonintegral original-lattice lift norm");
    Z value=numerator/a;
    if(value<best) {
      best=value; best_y=std::move(y); update_radius();
      FF updated=expanded_radius();
      updated.mul_2si(updated,-normExp);
      // Both old and new bounds contain every remaining improvement.
      bound=std::min(bound,updated.get_d(GMP_RNDU));
    }
  }
};

int main(int argc,char **argv) {
  try {
    if(argc!=2) throw std::runtime_error("Usage: fast_lift INPUT");
    std::ifstream in(argv[1]);
    std::string magic; int d; Z a,c,best;
    in>>magic>>d>>a>>c>>best;
    if(magic!="FAST-LIFT-1" || d<1 || a<=0 || c<=0 || best<=0)
      throw std::runtime_error("Invalid input");
    std::vector<Z> inner(d),initial(d);
    for(auto &z:inner) in>>z;
    for(auto &z:initial) in>>z;
    double rho;
    int precision=std::max(128,gso_min_prec(rho,d,LLL_DEF_DELTA,LLL_DEF_ETA)+32);
    FF::set_prec(precision);
    Matrix<FZ> H; H.resize(d,d);
    for(int i=0;i<d;++i) for(int j=0;j<d;++j) {
      Z z; in>>z; mpz_set(H[i][j].get_data(),z.get_mpz_t());
    }
    if(!in) throw std::runtime_error("Truncated input");
    Matrix<FZ> empty1,empty2;
    MatGSOGram<FZ,FF> gso(H,empty1,empty2,GSO_INT_GRAM);
    if(!gso.update_gso()) throw std::runtime_error("GSO failure");
    LiftEvaluator evaluator(d,gso,H,a,c,best,inner,initial);
    evaluator.init_delta_def(precision,rho,true);
    FF radius=evaluator.expanded_radius();
    Enumeration<FZ,FF> enumeration(gso,evaluator);
    auto began=std::chrono::steady_clock::now();
    // Empty pruning and target vectors: complete SVP enumeration, sign classes.
    enumeration.enumerate(0,d,radius,0);
    auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
    Z nodes=0;
    const auto per_level=enumeration.get_nodes_array();
    for(int i=0;i<d;++i) nodes+=Z(std::to_string(per_level[i]));
    std::cout<<"COMPLETE "<<evaluator.best<<' '<<nodes<<' '
             <<evaluator.leaves<<' '<<seconds<<' '<<precision<<'\n';
    for(const auto &z:evaluator.best_y) std::cout<<z<<' ';
    std::cout<<'\n';
  } catch(const std::exception &e) {
    std::cerr<<e.what()<<'\n'; return 2;
  }
}
