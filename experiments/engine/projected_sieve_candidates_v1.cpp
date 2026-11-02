// Projected SVP candidate search. Uses the GPL-2.0-or-later G6K kernel at the pinned
// c71e084 commit. Every output vector is rechecked exactly by the Python caller.
#define main unused_complete_enumerator_main
#include "fast_lift.cpp"
#undef main
#include "kernel/siever.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>

int main(int argc,char **argv) {
  try {
    if(argc!=9) throw std::runtime_error("Usage: projected_sieve_candidates INPUT CODES DB_SIZE THREADS ALPHA CONTEXT SEED ALGORITHM");
    const std::string algorithm=argv[8];
    if(algorithm!="bgj1" && algorithm!="hk3" && algorithm!="bdgl") throw std::runtime_error("Unknown sieve algorithm");
    std::ifstream input(argv[1]);std::string magic;int d;Z a,c,best;
    input>>magic>>d>>a>>c>>best;
    if(magic!="FAST-LIFT-1" || d<2 || d>MAX_SIEVING_DIM) throw std::runtime_error("Bad dimension/header");
    std::vector<Z> inner(d),initial(d);for(auto &z:inner) input>>z;for(auto &z:initial) input>>z;
    FF::set_prec(512);Matrix<FZ> H;H.resize(d,d);
    for(int i=0;i<d;++i) for(int j=0;j<d;++j) {Z z;input>>z;mpz_set(H[i][j].get_data(),z.get_mpz_t());}
    if(!input) throw std::runtime_error("Truncated input");
    Matrix<FZ> empty1,empty2;MatGSOGram<FZ,FF> gso(H,empty1,empty2,GSO_INT_GRAM);
    if(!gso.update_gso()) throw std::runtime_error("GSO failed");
    std::vector<std::vector<Z>> exact_gram(d,std::vector<Z>(d));
    for(int i=0;i<d;++i) for(int j=0;j<d;++j)
      mpz_set(exact_gram[i][j].get_mpz_t(),H[i][j].get_data());
    FF value;long exponent,normalizer=0;
    for(int i=0;i<d;++i) {value=gso.get_r_exp(i,i,exponent);normalizer=std::max(normalizer,exponent+value.exponent());}
    std::vector<double> coefficients(d*d,0);
    for(int i=0;i<d;++i) {
      value=gso.get_r_exp(i,i,exponent);value.mul_2si(value,exponent-normalizer);
      coefficients[i*d+i]=value.get_d();
      if(!(coefficients[i*d+i]>0) || !std::isfinite(coefficients[i*d+i])) throw std::runtime_error("Bad normalized GSO");
      for(int j=0;j<i;++j) {gso.get_mu(value,i,j);coefficients[i*d+j]=value.get_d();}
    }
    size_t target=std::stoull(argv[3]);int context=std::min(d,std::stoi(argv[6]));int left=d-context;
    if(context<2 || target<100 || target>2000000) throw std::runtime_error("Resource limit");
    SieverParams parameters;parameters.threads=std::stoul(argv[4]);parameters.reserved_n=d;
    parameters.reserved_db_size=target+100;parameters.simhash_codes_basedir=argv[2];
    parameters.otf_lift=true;
    Siever kernel(d,coefficients.data(),parameters,std::stoul(argv[7]));
    int current_left=std::max(left,d-40);
    kernel.initialize_local(0,current_left,d);
    auto began=std::chrono::steady_clock::now();
    std::cerr<<"HEURISTIC_SIEVE_START "<<d<<' '<<context<<' '<<target<<' '<<algorithm<<std::endl;
    std::vector<Z> best_y=initial;unsigned improvements=0,checked=0;
    auto evaluate=[&](const std::vector<long> &coordinates) {
      bool nonzero=false;for(auto z:coordinates) nonzero=nonzero || z!=0;if(!nonzero) return;
      std::vector<Z> y;for(auto z:coordinates) y.emplace_back(z);
      Z q=0,dot=0;
      for(int i=0;i<d;++i) {
        if(y[i]==0) continue;
        q+=y[i]*y[i]*exact_gram[i][i];dot+=y[i]*inner[i];
        for(int j=0;j<i;++j) if(y[j]!=0) q+=2*y[i]*y[j]*exact_gram[i][j];
      }
      const Z k=floor_div(2*dot+a,2*a);const Z residue=dot-k*a;Z numerator=c*q+residue*residue;
      if(!mpz_divisible_p(numerator.get_mpz_t(),a.get_mpz_t())) throw std::runtime_error("Noninteger norm");
      const Z norm=numerator/a;++checked;
      if(norm>0 && norm<best) {best=norm;best_y=std::move(y);++improvements;}
    };
    auto collect=[&]() {
    std::vector<size_t> order(kernel.db.size());std::iota(order.begin(),order.end(),0);
    // This is heuristic candidate selection, not an exclusion proof. Exact
    // checking of the shortest 64 database entries plus all best lifts avoids
    // spending most of the runtime rechecking thousands of long vectors.
    const size_t count=std::min<size_t>(64,order.size());
    std::partial_sort(order.begin(),order.begin()+count,order.end(),[&](size_t i,size_t j){return kernel.db[i].len<kernel.db[j].len;});
    for(size_t index=0;index<count;++index) {
      std::vector<long> y(d,0);const auto &entry=kernel.db[order[index]];
      for(int j=current_left;j<d;++j) y[j]=entry.x[j-current_left];
      bool usable=true;
      for(int i=current_left-1;i>=0;--i) {
        long double center=0;for(int j=i+1;j<d;++j) center-=static_cast<long double>(coefficients[j*d+i])*y[j];
        if(std::abs(center)>1e18L) {usable=false;break;}y[i]=std::llround(center);
      }
      if(usable) evaluate(y);
    }
    std::vector<long> lifts((current_left+1)*d,0);std::vector<double> lengths(current_left+1,0);
    kernel.best_lifts(lifts.data(),lengths.data());
    for(int i=0;i<=current_left;++i) if(lengths[i]>0 && std::isfinite(lengths[i])) {
      std::vector<long> y(lifts.begin()+i*d,lifts.begin()+(i+1)*d);evaluate(y);
      std::cout<<"HEURISTIC_BASIS_ROW "<<i<<' '<<std::setprecision(17)<<lengths[i];
      for(long z:y) std::cout<<' '<<z;std::cout<<'\n';
    }
    auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
    std::cout<<"HEURISTIC_CANDIDATE "<<best<<' '<<kernel.db_size()<<' '<<seconds<<' '<<checked<<' '<<improvements<<'\n';
    for(const auto &z:best_y) std::cout<<z<<' ';std::cout<<std::endl;
    };
    // Begin with a small context and retain its reduced database while growing.
    // A full-dimensional random database can take far longer to initialize.
    // Completed rounds emit independently checkable witnesses immediately.
    while(true) {
      const int n=d-current_left;
      const size_t round_target=std::min(target,static_cast<size_t>(std::ceil(3.2*std::pow(4./3.,n/2.))));
      kernel.grow_db(round_target);
      collect();
      const double alpha=std::min(.8,std::stod(argv[5])*std::sqrt(static_cast<double>(context)/n));
      if(algorithm=="bgj1") kernel.bgj1_sieve(alpha);
      else if(algorithm=="hk3") kernel.hk3_sieve(alpha);
      else {
        const size_t blocks=std::min<size_t>(2,std::max<size_t>(1,n/28));
        const size_t multi_hash=4;
        double buckets=.3*std::pow(2.,(blocks-1.)/(blocks+1.))*std::pow(static_cast<double>(multi_hash),2.*blocks/(blocks+1.))*std::pow(static_cast<double>(round_target),blocks/(1.+blocks));
        buckets=std::max(std::pow(2.,blocks-1.),std::min(buckets,multi_hash*round_target/128.));
        kernel.bdgl_sieve(static_cast<size_t>(buckets),blocks,multi_hash);
      }
      collect();
      std::cerr<<"HEURISTIC_ROUND_COMPLETE "<<n<<' '<<round_target<<std::endl;
      if(current_left==left) break;
      kernel.extend_left(1);--current_left;
    }
    return 0;
  } catch(const std::exception &error) {std::cerr<<error.what()<<std::endl;return 2;}
}
