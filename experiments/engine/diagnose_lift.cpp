#define main unused_whole_sphere_main
#include "fast_lift.cpp"
#undef main
#include <sys/resource.h>
#include <execinfo.h>
#include <cstdlib>
#include <cstdio>
void *operator new(std::size_t size) {
  if(size>128UL*1024*1024) {
    std::fprintf(stderr,"LARGE_CPP_ALLOCATION %zu\n",size);
    void *trace[25];int count=backtrace(trace,25);backtrace_symbols_fd(trace,count,2);
  }
  void *out=std::malloc(size?size:1);
  if(!out) throw std::bad_alloc();
  return out;
}
void operator delete(void *pointer) noexcept {std::free(pointer);}
void operator delete(void *pointer,std::size_t) noexcept {std::free(pointer);}
class TraceEvaluator : public LiftEvaluator {
public:
  using LiftEvaluator::LiftEvaluator;
  void eval_sol(const std::vector<FF> &coordinates,const enumf &distance,enumf &bound) override {
    if(leaves<2) std::cerr<<"CALLBACK coordinates="<<coordinates.size()<<" evaluator_d="<<d<<" normExp="<<normExp<<std::endl;
    for(const auto &coordinate:coordinates)
      if(!mpfr_number_p(coordinate.get_data())) throw std::runtime_error("Nonfinite enumerated coordinate");
    LiftEvaluator::eval_sol(coordinates,distance,bound);
  }
};
static void checkpoint(const char *label) {
  struct rusage usage;getrusage(RUSAGE_SELF,&usage);
  std::cerr<<label<<" maxRSS_kB="<<usage.ru_maxrss<<std::endl;
}
int main(int argc,char **argv) {
  try {
    if(argc<2) return 1;
    std::ifstream in(argv[1]);std::string magic;int d;Z a,c,best;
    in>>magic>>d>>a>>c>>best;
    std::vector<Z> inner(d),initial(d);
    for(auto &z:inner) in>>z;for(auto &z:initial) in>>z;
    double rho;int precision=std::max(128,gso_min_prec(rho,d,LLL_DEF_DELTA,LLL_DEF_ETA)+32);
    std::cerr<<"dimension="<<d<<" precision="<<precision<<" rho="<<rho<<std::endl;
    FF::set_prec(precision);
    Matrix<FZ> H;H.resize(d,d);size_t max_bits=0;
    for(int i=0;i<d;++i) for(int j=0;j<d;++j) {
      Z z;in>>z;mpz_set(H[i][j].get_data(),z.get_mpz_t());
      max_bits=std::max(max_bits,mpz_sizeinbase(z.get_mpz_t(),2));
    }
    if(!in) throw std::runtime_error("Truncated input");
    std::cerr<<"maximum_Gram_entry_bits="<<max_bits<<std::endl;checkpoint("H_loaded");
    Matrix<FZ> empty1,empty2;MatGSOGram<FZ,FF> gso(H,empty1,empty2,GSO_INT_GRAM);checkpoint("gso_constructed");
    if(!gso.update_gso()) throw std::runtime_error("GSO failure");checkpoint("gso_updated");
    TraceEvaluator evaluator(d,gso,H,a,c,best,inner,initial);checkpoint("evaluator_constructed");
    evaluator.init_delta_def(precision,rho,true);checkpoint("error_initialized");
    FF radius=evaluator.expanded_radius();checkpoint("radius_expanded");
    Enumeration<FZ,FF> enumeration(gso,evaluator);checkpoint("enumerator_constructed");
    if(argc>2 && std::string(argv[2])=="--enumerate") {
      enumeration.enumerate(0,d,radius,0);checkpoint("enumeration_complete");
    }
    return 0;
  } catch(const std::exception &e) {std::cerr<<"EXCEPTION "<<e.what()<<std::endl;return 2;}
}
