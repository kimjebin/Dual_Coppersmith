import os, sys, json, time
os.environ["PATH"]=(""+os.environ.get("PATH",""))
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("SAGE_NUM_THREADS", "1")
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import warnings, logging; warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
from itertools import product as iproduct
from math import comb
from sage.all import (ZZ, PolynomialRing, random_prime, randint, set_random_seed, matrix, diagonal_matrix, prod, gcd, inverse_mod, lcm)
from shared import dual_moment as dm, small_roots
def gen(pbits,qbits,n,Xbits,seed=0):
    set_random_seed(seed)
    p=random_prime(ZZ(2)**pbits-1,lbound=ZZ(2)**(pbits-1)); q=random_prime(ZZ(2)**qbits-1,lbound=ZZ(2)**(qbits-1)); N=p*q
    X=ZZ(2)**Xbits; P=PolynomialRing(ZZ,["x%d"%i for i in range(n)]); xs=P.gens()
    r=[randint(-int(X)+1,int(X)-1) for _ in range(n)]; a=[randint(1,int(N)-1) for _ in range(n)]
    while gcd(a[0],N)!=1: a[0]=randint(1,int(N)-1)
    S=sum(a[i]*r[i] for i in range(n)); a0=-S+int(p)*randint(1,int(q)); f=a0+sum(a[i]*xs[i] for i in range(n))
    return int(p),int(q),int(N),f,[int(v) for v in r],int(X),P,xs,a
def hm_shifts(f,N,t,P,xs,a):
    n=len(xs); fm=P({e:(ZZ(c)*inverse_mod(ZZ(a[0]),N))%N for e,c in zip(f.exponents(),f.coefficients())})
    fpw=[P.one()]
    for _ in range(t): fpw.append(fpw[-1]*fm)
    sh=[]
    for k in range(t+1):
        for iv in iproduct(*[range(t-k+1) for _ in range(n-1)]):
            if sum(iv)<=t-k: sh.append((ZZ(N)**(t-k))*fpw[k]*prod(xs[1+j]**iv[j] for j in range(n-1)))
    return sh
def hm_dual(f,N,t,P,xs,a,X,use_flatter=False):
    sh=hm_shifts(f,N,t,P,xs,a); n=len(xs)
    monos=dm._monomials_up_to(P,t); keys=[dm._exp_vec(mo) for mo in monos]; idx={kk:i for i,kk in enumerate(keys)}; nn=len(monos)
    rows=[]
    for g in sh:
        row=[0]*nn
        for e,c in g.dict().items():
            kk=tuple(e) if hasattr(e,"__iter__") else (e,)
            if kk in idx: row[idx[kk]]=c
        rows.append(row)
    A=matrix(ZZ,rows); L=dm.congruence_kernel(A,ZZ(N)**t,method="snf")
    def bpk(k): return prod((ZZ(X))**k[i] for i in range(n))
    rho=lcm([bpk(k) for k in keys]); scales=[rho//bpk(k) for k in keys]
    LBr=dm._reduce(L*diagonal_matrix(ZZ,scales),use_flatter=use_flatter)
    one=(0,)*n
    for w in LBr.rows():
        w0=ZZ(w[idx[one]])
        if w0==0: continue
        cand=[]; ok=True
        for i in range(n):
            e=tuple(1 if a2==i else 0 for a2 in range(n)); num=ZZ(X)*ZZ(w[idx[e]])
            if num%w0!=0: ok=False; break
            cand.append(int(num//w0))
        g = ZZ(f(*cand)).gcd(N) if ok else ZZ(1)
        if ok and all(abs(c)<X for c in cand) and 1 < g < N:
            return cand
    return None
def hm_primal(f,N,t,P,xs,a,X,p):
    sh=hm_shifts(f,N,t,P,xs,a); n=len(xs)
    B,monos=small_roots.create_lattice(P,sh,[X]*n)
    B=small_roots.reduce_lattice(B)
    M=ZZ(p)**t
    polys=small_roots.reconstruct_polynomials(B,None,M,monos,[X]*n)
    for roots in small_roots.find_roots(P,polys,method="groebner"):
        cand=[int(roots[xi]) for xi in xs]
        if all(-X<c<X for c in cand) and ZZ(f(*cand)).gcd(N)>1: return cand
    return None
A=[v for v in sys.argv[1:] if v!="--"]
pb,qb,n,Xb,t=int(A[0]),int(A[1]),int(A[2]),int(A[3]),int(A[4]); method=A[5]
p,q,N,f,r,X,P,xs,a=gen(pb,qb,n,Xb,0)
t0=time.time()
d = hm_primal(f,N,t,P,xs,a,X,p) if method=="primal" else hm_dual(f,N,t,P,xs,a,X,False)
print("RESULT "+json.dumps({"n":n,"dim":comb(n+t,t),"method":method,"time_s":float(time.time()-t0),"ok":bool(d==r)}),flush=True)
