"""Original ACD/HM lattices, first-jet geometry, and exact membership helpers.

Only the original lattice B is supported by this public artifact.
All bases use row vectors. First-jet derivatives describe the prediction model.
"""
from __future__ import annotations
import argparse,csv,hashlib,itertools,json,math,os,re,sys,time
from pathlib import Path
from paths import ROOT
REPO=ROOT
sys.path.insert(0,str(REPO))
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ.setdefault(key,"1")
from sage.all import ZZ,QQ,PolynomialRing,diagonal_matrix,matrix,vector,prod,gcd,lcm,identity_matrix
from sage.env import SAGE_VERSION
from second_minimum.instances import acd_basis,hm_basis
from shared.dual_moment import snf_congruence_kernel,_monomials_up_to,_exp_vec
VARIANTS={"B":(0,)}
GRID=([("acd",2,1,ell) for ell in range(10,14)]
      +[("acd",2,2,ell) for ell in range(10,13)]
      +[("acd",3,m,ell) for m in (2,3) for ell in range(5,8)]
      +[("hm",3,3,ell) for ell in range(2,10)])


def dump(path, doc):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(doc, ensure_ascii=False, sort_keys=True) + '\n',
                         encoding='utf-8')
    temporary.replace(path)


def rows(a):
    return [[str(c) for c in row] for row in a.rows()]


def norm2(v):
    return ZZ(v.dot_product(v))


def coefficients(g):
    # Sage ETuple and Python tuple do not share a dictionary hash.
    return {tuple(k): c for k, c in g.dict().items()}


def ln2norm(v):
    return math.log2(int(norm2(v))) / 2


def name(cfg):
    app, t, m, ell = cfg
    return f'{app}-t{t}-m{m}-ell{ell}-seed0'


def build(cfg):
    app, t, m, ell = cfg
    base, meta = acd_basis(t, m, ell) if app == 'acd' else hm_basis(ell)
    keys = [tuple(k) for k in meta['monomial_keys']]
    N, R = ZZ(meta['modulus']), ZZ(meta['bound'])
    ring = PolynomialRing(ZZ, [f'x{i}' for i in range(ell)])
    xs = ring.gens()
    if app == 'acd':
        ys = [xs[i] - ZZ(meta['public_samples'][i]) for i in range(ell)]
        p = gcd([N] + [ZZ(a) - r for a, r in
                       zip(meta['public_samples'], meta['planted_root'])])
    else:
        ys = [ZZ(meta['normalized_constant']) + sum(
            ZZ(a) * x for a, x in zip(meta['normalized_coefficients'], xs))]
        ys.extend(xs[1:])
        p = gcd(N, ZZ(meta['public_constant']) + sum(
            ZZ(a) * r for a, r in zip(meta['public_coefficients'], meta['planted_root'])))
        assert gcd(N, 6) == 1  # Ordinary derivative factorials are units modulo N.
    assert p.is_prime() and 1 < p < N
    # For ACD ell >> t, unused variables give derivative coefficient 1 even
    # if N has factors 2 or 3. No factorial is discarded without this reason.
    assert app != 'acd' or ell >= t + 2
    polys = [prod(ys[i] ** key[i] for i in range(ell)) for key in keys]
    forward = matrix(ZZ, [[coefficients(g).get(key, 0) for key in keys] for g in polys])
    scales = [R ** (t - sum(key)) for key in keys]
    original_exponents = [min(m, sum(k)) if app == 'acd' else k[0] for k in keys]
    raw = matrix(ZZ, [[base[i, j] // scales[j] for j in range(len(keys))]
                      for i in range(len(keys))])
    assert raw * forward.transpose() == diagonal_matrix(ZZ, [N ** a for a in original_exponents])
    root = list(map(ZZ, meta['planted_root']))
    root_moments = vector(ZZ, [prod(root[i] ** k[i] for i in range(ell)) for k in keys])
    yroot = vector(ZZ, [g(*root) for g in polys])
    jet = [list(root_moments)]
    for j in range(ell):
        jet.append([ZZ(0) if not k[j] else
                    k[j] * prod(root[i] ** (k[i] - (i == j)) for i in range(ell))
                    for k in keys])
    jet = matrix(ZZ, jet)
    return dict(cfg=cfg, base=base, meta=meta, keys=keys, N=N, R=R, p=p,
                ring=ring, xs=xs, polys=polys, forward=forward, scales=scales,
                original_exponents=original_exponents, root=root,
                root_moments=root_moments, yroot=yroot, jet=jet)


def on_root(v, data):
    # Full Veronese membership, not just degree-one coordinate ratios.
    target = [z * s for z, s in zip(data['root_moments'], data['scales'])]
    return bool(v[0]) and all(v[j] * target[0] == v[0] * target[j] for j in range(len(v)))


def predictions(data, label, exponents, root_vector):
    t, m, ell = data['cfg'][1:]
    n = len(exponents)
    logdet = sum(exponents) * math.log2(int(data['N'])) + sum(t-sum(k) for k in data['keys']) * math.log2(int(data['R']))
    root_length = ln2norm(root_vector)
    d = n - 1
    log_ball = d / 2 * math.log2(math.pi) - math.lgamma(d / 2 + 1) / math.log(2)
    gh = (logdet - root_length - log_ball) / d
    # Integer coefficients in the jet span follow from coordinates 0,e_i.
    # Diagonal affine constraints define their coefficient lattice exactly.
    Q = data['forward'] * data['jet'].transpose()
    Q = diagonal_matrix(ZZ, [data['N'] ** (m-a) for a in exponents]) * Q
    C = snf_congruence_kernel(Q, data['N'] ** m)
    index = abs(C.det())
    balanced_jet = data['jet'] * diagonal_matrix(ZZ, data['scales'])
    gramdet = (balanced_jet * balanced_jet.transpose()).det()
    covol = math.log2(int(index)) + math.log2(int(gramdet)) / 2
    jet = 0.5 * math.log2(ell / (2 * math.pi * math.e)) + (covol-root_length) / ell
    return dict(log2_root=root_length, log2_gh=gh, log2_jet=jet,
                log2_pred=min(gh, jet), pred_minus_root=min(gh, jet)-root_length,
                jet_coefficient_index=str(index), determinant_N_exponent=sum(exponents),
                determinant_R_exponent=sum(t-sum(k) for k in data['keys']))


def variant(data,label="B"):
    if label!="B":raise ValueError("This artifact supports the original lattice B only")
    exponents=list(data["original_exponents"])
    scalar=lcm([data["N"]**a//gcd(data["N"]**a,value)
                for a,value in zip(exponents,data["yroot"])])
    root_vector=vector(ZZ,[scalar*z*s for z,s in zip(data["root_moments"],data["scales"])])
    return data["base"],exponents,scalar,root_vector


def original_constraints(data):
    app,t,m,ell=data["cfg"]
    shifts=[data["N"]**(max(m-sum(k),0) if app=="acd" else m-k[0])*g
            for k,g in zip(data["keys"],data["polys"])]
    modulus=data["N"]**m;index={k:i for i,k in enumerate(data["keys"])}
    unique=set()
    for shift in shifts:
        entries=tuple(sorted((index[tuple(k)],int(c%modulus))
                             for k,c in shift.dict().items() if c%modulus))
        if entries:unique.add(entries)
    return {0:sorted(unique)}


def member(v,data,constraints,orders=(0,)):
    if tuple(orders)!=(0,):raise ValueError("Only original-lattice membership is supported")
    if any(z%s for z,s in zip(v,data["scales"])):return False
    raw=[z//s for z,s in zip(v,data["scales"])]
    modulus=data["N"]**data["cfg"][2]
    return all(sum(c*raw[i] for i,c in row)%modulus==0 for row in constraints[0])
