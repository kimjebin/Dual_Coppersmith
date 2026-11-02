// Exact LLL-contract verification for small Gram-only lower-bound problems.
// All arithmetic, including Bareiss divisions, is integral GMP arithmetic.
#include <gmpxx.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using Z=mpz_class;
static Z ceildiv(const Z&a,const Z&b) {Z q;mpz_cdiv_q(q.get_mpz_t(),a.get_mpz_t(),b.get_mpz_t());return q;}
static Z ceilsqrt(const Z&x) {Z q;mpz_sqrt(q.get_mpz_t(),x.get_mpz_t());if(q*q<x)++q;return q;}
int main(int argc,char**argv) {
 try {
  if(argc!=2)throw std::runtime_error("Usage: verify_bound_gram INPUT");
  std::ifstream in(argv[1]);std::string tag;int n;Z a,c,initial,x;
  in>>tag>>n>>a>>c>>initial;
  if(tag!="FAST-LIFT-1"||n<1||a!=1||c!=1||initial<=0)throw std::runtime_error("Bad header");
  for(int i=0;i<n;++i){in>>x;if(x!=0)throw std::runtime_error("Nonzero root inner product");}
  std::vector<Z> y(n);for(auto&v:y)in>>v;
  std::vector<std::vector<Z>> H(n,std::vector<Z>(n)),ell(n,std::vector<Z>(n));
  for(auto&row:H)for(auto&v:row)in>>v;
  if(!in)throw std::runtime_error("Truncated input");
  Z norm=0;
  for(int i=0;i<n;++i)for(int j=0;j<n;++j){if(H[i][j]!=H[j][i])throw std::runtime_error("Not symmetric");norm+=y[i]*H[i][j]*y[j];}
  if(norm!=initial)throw std::runtime_error("Initial witness mismatch");
  std::vector<Z> delta(n+1);delta[0]=1;
  for(int k=0;k<n;++k){
   const Z pivot=H[k][k];if(pivot<=0)throw std::runtime_error("Not positive definite");delta[k+1]=pivot;
   for(int i=k+1;i<n;++i)ell[i][k]=H[i][k];
   for(int i=k+1;i<n;++i)for(int j=i;j<n;++j){
    Z numerator=pivot*H[i][j]-H[i][k]*H[k][j];
    if(!mpz_divisible_p(numerator.get_mpz_t(),delta[k].get_mpz_t()))throw std::runtime_error("Bad Bareiss division");
    mpz_divexact(H[i][j].get_mpz_t(),numerator.get_mpz_t(),delta[k].get_mpz_t());H[j][i]=H[i][j];
   }
  }
  for(int i=1;i<n;++i){
   for(int j=0;j<i;++j)if(100*abs(ell[i][j])>51*delta[j+1])throw std::runtime_error("Size reduction fails");
   if(100*(delta[i+1]*delta[i-1]+ell[i][i-1]*ell[i][i-1])<99*delta[i]*delta[i])throw std::runtime_error("Lovasz fails");
  }
  const Z den=Z(1)<<64,radius=initial-1;std::vector<Z> bounds(n);Z largest=0;
  for(int i=n-1;i>=0;--i){
   Z weighted=0;for(int j=i+1;j<n;++j)weighted+=ceildiv(abs(ell[j][i])*den,delta[i+1])*bounds[j];
   bounds[i]=ceilsqrt(ceildiv(radius*delta[i],delta[i+1]))+ceildiv(weighted,den);
   if(bounds[i]>largest)largest=bounds[i];
  }
  int minimum=0;for(int i=1;i<n;++i)if(delta[i+1]*delta[minimum]<delta[minimum+1]*delta[i])minimum=i;
  std::cout<<"GRAM_LLL_OK "<<n<<' '<<largest<<'\n';for(auto&v:bounds)std::cout<<v<<' ';std::cout<<'\n';
  std::cout<<delta[minimum+1]<<' '<<delta[minimum]<<'\n';
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
}
