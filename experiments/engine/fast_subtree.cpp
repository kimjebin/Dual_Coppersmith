// Subtree interface to the unchanged proved full-sphere lift evaluator.
#define main full_sphere_main
#include "fast_lift.cpp"
#undef main

int main(int argc,char **argv) {
  try {
    if(argc<2) throw std::runtime_error("Usage: fast_subtree INPUT [suffix entries]");
    std::ifstream in(argv[1]);
    std::string magic; int d; Z a,c,best;
    in>>magic>>d>>a>>c>>best;
    if(magic!="FAST-LIFT-1" || d<1 || a<=0 || c<=0 || best<=0 || argc-2>d)
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
    std::vector<enumxt> suffix;
    for(int i=2;i<argc;++i) {
      Z z(argv[i]);
      if(abs(z)>=(Z(1)<<50)) throw std::runtime_error("Suffix too large");
      suffix.push_back(z.get_d());
    }
    Matrix<FZ> empty1,empty2;
    MatGSOGram<FZ,FF> gso(H,empty1,empty2,GSO_INT_GRAM);
    if(!gso.update_gso()) throw std::runtime_error("GSO failure");
    LiftEvaluator evaluator(d,gso,H,a,c,best,inner,initial);
    evaluator.init_delta_def(precision,rho,true);
    FF radius=evaluator.expanded_radius();
    Enumeration<FZ,FF> enumeration(gso,evaluator);
    auto began=std::chrono::steady_clock::now();
    if(suffix.size()<static_cast<size_t>(d)) {
      // fplll disables sign symmetry for nonzero fixed suffixes itself.
      // subtree_reset=false retains the fixed suffix's full partial norm.
      enumeration.enumerate(0,d,radius,0,{},suffix,{},false,false);
    } else {
      std::vector<FF> coordinates(d);
      bool nonzero=false;
      for(int i=0;i<d;++i) {
        FZ z; z=static_cast<long>(suffix[i]); coordinates[i].set_z(z);
        nonzero=nonzero || suffix[i]!=0;
      }
      if(nonzero) { enumf dummy=radius.get_d(GMP_RNDU); evaluator.eval_sol(coordinates,0,dummy); }
    }
    auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
    Z nodes=0;
    if(suffix.size()<static_cast<size_t>(d)) {
      const auto per_level=enumeration.get_nodes_array();
      for(int i=0;i<d;++i) nodes+=Z(std::to_string(per_level[i]));
    }
    std::cout<<"COMPLETE "<<evaluator.best<<' '<<nodes<<' '
             <<evaluator.leaves<<' '<<seconds<<' '<<precision<<'\n';
    for(const auto &z:evaluator.best_y) std::cout<<z<<' ';
    std::cout<<'\n';
  } catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 2; }
}
