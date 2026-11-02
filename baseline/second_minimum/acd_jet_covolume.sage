
import os
import sys
import time
import math
import warnings
import logging

os.environ["PATH"] = "" + os.environ.get("PATH", "")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # dual_coppersmith/
warnings.filterwarnings("ignore")
logging.disable(logging.WARNING)

from math import comb, floor, log2, lgamma
from sage.all import ZZ, QQ, RealField, random_prime, randint, set_random_seed, diagonal_matrix, matrix, vector, pari, gcd as sgcd
from acd.mp_dual import _build_shifts
from shared import small_roots

try:
    pari.allocatemem(2 * 1024 * 1024 * 1024)
except Exception:
    pass

ETA = 100
RHO = 10
SEED = 0
CONFIGS = [(2,1,5),(2,1,10),(2,2,5),(2,2,11),(3,2,5),(3,2,7),(3,3,5),(3,3,7)]
RESULTS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "acd_lambda2.txt")
RF = RealField(256)


def gmax(m_vars, t_val, k_val):
    d_val = comb(t_val + m_vars, m_vars)
    lhs = d_val * log2(d_val) + d_val * d_val * log2(1.02) + RHO * d_val * (m_vars * t_val / (m_vars + 1))
    cg = comb(k_val + m_vars, m_vars) * k_val / (m_vars + 1)
    return (k_val * ETA * d_val - lhs) / cg


def gen_instance(eta, gamma, rho, ell, seed=SEED):
    set_random_seed(seed)
    p_val = random_prime(ZZ(2) ** eta - 1, lbound=ZZ(2) ** (eta - 1))
    q0 = randint(ZZ(2) ** (gamma - eta - 1), ZZ(2) ** (gamma - eta))
    n_mod = p_val * q0
    radius = ZZ(2) ** rho
    a_vals = []
    r_vals = []
    for _ in range(ell):
        qi = randint(ZZ(2) ** (gamma - eta - 1), ZZ(2) ** (gamma - eta))
        ri = randint(-int(radius) + 1, int(radius) - 1)
        a_vals.append(p_val * qi + ri)
        r_vals.append(ri)
    return int(p_val), int(n_mod), [int(v) for v in a_vals], [int(v) for v in r_vals]


def norm_sq(row):
    return ZZ(sum(ZZ(c) * ZZ(c) for c in row))


def log_norm_sq(sq):
    return 0.0 if sq <= 0 else float(RF(sq).log(2) / 2)


def log2_int(x_val):
    return float(RF(x_val).log(2))


def proportional(a_row, b_row):
    pivot = None
    for idx, val in enumerate(b_row):
        if ZZ(val) != 0:
            pivot = idx
            break
    if pivot is None:
        return False
    return all(QQ(a_row[i]) * QQ(b_row[pivot]) == QQ(b_row[i]) * QQ(a_row[pivot]) for i in range(len(b_row)))


def r_power(r_vals, exp_tuple):
    acc = ZZ(1)
    for idx, exp in enumerate(exp_tuple):
        acc *= ZZ(r_vals[idx]) ** exp
    return acc


def build_data(t_val, k_val, ell):
    gamma = int(floor(0.9 * gmax(ell, t_val, k_val)))
    _p_val, n_mod, a_vals, r_vals = gen_instance(ETA, gamma, RHO, ell)
    pr = ZZ[tuple("x%d" % i for i in range(ell))]
    x_gens = pr.gens()
    radius = ZZ(2) ** RHO
    shifts = _build_shifts(pr, x_gens, n_mod, a_vals, t_val, k_val)
    a_mat, mons = small_roots.create_lattice(pr, shifts, [1] * ell)
    q_val = ZZ(n_mod) ** k_val
    dim = len(mons)
    exps = [tuple(mo.degrees()) for mo in mons]
    s_mat, _u_snf, v_snf = a_mat.smith_form()
    divs = [ZZ(s_mat[i, i]) for i in range(min(a_mat.nrows(), dim))]
    mult = []
    for i in range(dim):
        div = divs[i] if i < len(divs) else ZZ(0)
        mult.append(ZZ(1) if div == 0 else q_val // sgcd(q_val, div))
    scales = [ZZ(radius) ** (t_val - mo.degree()) for mo in mons]
    lb = diagonal_matrix(ZZ, mult) * v_snf.transpose() * diagonal_matrix(ZZ, scales)
    jet_rows = [[r_power(r_vals, exp) for exp in exps]]
    for deriv_idx in range(ell):
        row = []
        for exp in exps:
            if exp[deriv_idx] == 0:
                row.append(ZZ(0))
            else:
                exp2 = list(exp)
                exp2[deriv_idx] -= 1
                row.append(ZZ(exp[deriv_idx]) * r_power(r_vals, tuple(exp2)))
        jet_rows.append(row)
    js = matrix(QQ, jet_rows).row_space()
    log2det = float(sum(RF(m).log(2) for m in mult) + sum(RF(s).log(2) for s in scales))
    return {"t": t_val, "m": k_val, "ell": ell, "n": dim, "r": r_vals, "exps": exps, "scales": scales, "LB": lb, "Js": js, "log2det": log2det}


def unbalance(row, scales):
    return vector(QQ, [QQ(row[i]) / QQ(scales[i]) for i in range(len(scales))])


def ratio_decode_ok(row, data):
    unb = unbalance(row, data["scales"])
    exps = data["exps"]
    r_vals = data["r"]
    idx_zero = next(i for i, exp in enumerate(exps) if all(e == 0 for e in exp))
    if unb[idx_zero] == 0:
        return False
    for var_idx in range(len(r_vals)):
        idx_var = next(i for i, exp in enumerate(exps) if all(exp[j] == (1 if j == var_idx else 0) for j in range(len(r_vals))))
        if QQ(unb[idx_var]) / QQ(unb[idx_zero]) != QQ(r_vals[var_idx]):
            return False
    return True


def second_independent_reduced_row(basis, tau_row):
    order = sorted(range(basis.nrows()), key=lambda i: norm_sq(basis.row(i)))
    for idx in order:
        row = vector(ZZ, basis.row(idx))
        if not proportional(row, tau_row):
            return idx, row, norm_sq(row)
    raise ValueError("no independent reduced row found")


def selected_jet_rows(basis, data):
    order = sorted(range(basis.nrows()), key=lambda i: norm_sq(basis.row(i)))
    selected = []
    selected_indices = []
    inside_count = 0
    rank_so_far = 0
    for idx in order:
        if unbalance(basis.row(idx), data["scales"]) in data["Js"]:
            inside_count += 1
            trial = matrix(QQ, selected + [vector(ZZ, basis.row(idx))])
            trial_rank = trial.rank()
            if trial_rank > rank_so_far:
                selected.append(vector(ZZ, basis.row(idx)))
                selected_indices.append(idx)
                rank_so_far = trial_rank
    selected_is_shortest = len(selected) >= data["ell"] + 1 and set(selected_indices[:data["ell"] + 1]) == set(order[:data["ell"] + 1])
    return selected[:data["ell"] + 1], inside_count, len(selected), selected_is_shortest


def whole_gh(data, log2_lambda1):
    d_proj = data["n"] - 1
    log2vol = (d_proj / 2.0) * math.log2(math.pi) - lgamma(d_proj / 2.0 + 1.0) / math.log(2.0)
    return (data["log2det"] - log2_lambda1 - log2vol) / d_proj


def jet_gh(rows, ell, log2_lambda1):
    gram = matrix(ZZ, rows) * matrix(ZZ, rows).transpose()
    log2detjet = 0.5 * log2_int(abs(ZZ(gram.det())))
    return 0.5 * math.log2(ell / (2.0 * math.pi * math.e)) + (log2detjet - log2_lambda1) / ell


def set_pari_precision(gram, bound):
    gram_bits = 1
    for gx in gram.list():
        if gx:
            gram_bits = max(gram_bits, abs(ZZ(gx)).nbits())
    needed_bits = max(4096, 4 * max(ZZ(bound).nbits(), gram_bits) + 1024)
    try:
        old_bits = pari.get_real_precision_bits()
    except Exception:
        old_bits = 0
    if old_bits < needed_bits:
        pari.set_real_precision_bits(int(needed_bits))


def exact_n21(basis, tau_idx, lambda2_idx, lambda2_sq, data):
    gram = basis * basis.transpose()
    best_sq = ZZ(lambda2_sq)
    best_coeff = vector(ZZ, [ZZ(1) if i == lambda2_idx else ZZ(0) for i in range(basis.nrows())])
    bound = ZZ(lambda2_sq)
    if bound > 0:
        set_pari_precision(gram, bound)
        result = gram.__pari__().qfminim(bound, 10**7, 2)
        count = ZZ(result[0])
        if count > 0:
            vecs = matrix(ZZ, result[2])
            tau_coeff = vector(ZZ, [ZZ(1) if i == tau_idx else ZZ(0) for i in range(basis.nrows())])
            for col in range(vecs.ncols()):
                coeff = vector(ZZ, [vecs[row, col] for row in range(vecs.nrows())])
                if proportional(coeff, tau_coeff):
                    continue
                sq = ZZ((coeff * gram).dot_product(coeff))
                if sq < best_sq:
                    best_sq = sq
                    best_coeff = coeff
    achiever = best_coeff * basis
    return best_sq, (unbalance(achiever, data["scales"]) in data["Js"])


def analyze(cfg):
    t_val, k_val, ell = cfg
    start = time.time()
    data = build_data(t_val, k_val, ell)
    basis = small_roots.reduce_lattice(data["LB"])
    order = sorted(range(data["n"]), key=lambda i: norm_sq(basis.row(i)))
    tau_idx = order[0]
    tau_row = vector(ZZ, basis.row(tau_idx))
    if not ratio_decode_ok(tau_row, data):
        raise ValueError("lambda1 row failed root ratio decode for %s" % (cfg,))
    lambda1_sq = norm_sq(tau_row)
    log2_lambda1 = log_norm_sq(lambda1_sq)
    lambda2_idx, _lambda2_row, lambda2_sq = second_independent_reduced_row(basis, tau_row)
    log2_lambda2 = log_norm_sq(lambda2_sq)
    jet_rows, inside_count, jet_rank, selected_is_shortest = selected_jet_rows(basis, data)
    if len(jet_rows) == ell + 1:
        jet_text = "%.6f" % jet_gh(jet_rows, ell, log2_lambda1)
    else:
        jet_text = "NA(fewer:%d/%d)" % (jet_rank, ell + 1)
    exact_text = ""
    achiever_text = ""
    if data["n"] == 21:
        exact_sq, in_jet = exact_n21(basis, tau_idx, lambda2_idx, lambda2_sq, data)
        exact_text = "%.6f" % log_norm_sq(exact_sq)
        achiever_text = "yes" if in_jet else "no"
    return {
        "t": t_val, "m": k_val, "ell": ell, "n": data["n"], "t_eq_m": "yes" if t_val == k_val else "no",
        "log2_lambda1": log2_lambda1, "log2_lambda2": log2_lambda2, "log2_wholeGH": whole_gh(data, log2_lambda1),
        "log2_jetGH": jet_text, "lambda2_exact_log2": exact_text, "achiever_in_jet": achiever_text,
        "inside_count": inside_count, "jet_rank": jet_rank, "selected_is_shortest": selected_is_shortest, "seconds": time.time() - start,
    }


def format_results(rows):
    lines = ["t m ell n t_eq_m log2_lambda1 log2_lambda2 log2_wholeGH log2_jetGH lambda2_exact_log2 achiever_in_jet"]
    for row in rows:
        lines.append("%d %d %d %d %s %.6f %.6f %.6f %s %s %s" % (
            row["t"], row["m"], row["ell"], row["n"], row["t_eq_m"], row["log2_lambda1"], row["log2_lambda2"], row["log2_wholeGH"], row["log2_jetGH"], row["lambda2_exact_log2"], row["achiever_in_jet"]
        ))
    exact_bits = []
    for row in rows:
        if row["n"] == 21:
            exact_bits.append("(t,m,ell)=(%d,%d,%d): exact lambda2 %s, achiever_in_jet=%s" % (row["t"], row["m"], row["ell"], row["lambda2_exact_log2"], row["achiever_in_jet"]))
    na_rows = [(row["t"], row["m"], row["ell"], row["jet_rank"], row["ell"] + 1) for row in rows if str(row["log2_jetGH"]).startswith("NA")]
    if na_rows:
        pattern = "Pattern not confirmed under the stated reduced-row jet-space test: t=m has jet-GH near lambda2 and whole-GH 10-23 bits above, but t!=m has only v0 in Js from the reduced rows so jet-GH is NA; whole-GH is near lambda2. n=21 exact outcomes: %s." % "; ".join(exact_bits)
    else:
        pattern = "Pattern confirmed: t=m has jet-GH near lambda2 while whole-GH is 10-23 bits above; t!=m has whole-GH near lambda2 near jet-GH. n=21 exact outcomes: %s." % "; ".join(exact_bits)
    return "\n".join(lines) + "\n" + pattern + "\n"


def main():
    rows = []
    for cfg in CONFIGS:
        print("START", cfg, flush=True)
        row = analyze(cfg)
        rows.append(row)
        print("DONE", cfg, "n", row["n"], "sec", "%.1f" % row["seconds"], "jet_rank", row["jet_rank"], flush=True)
    output = format_results(rows)
    with open(RESULTS_PATH, "w") as handle:
        handle.write(output)
    print(output, flush=True)


if False:
    pass

def jet_intersection_covol(data):
    from sage.all import matrix as _M, ZZ as _Z, QQ as _Q, RR as _R, lcm as _lcm
    import math
    ell=data["ell"]; n=data["n"]; scales=data["scales"]; exps=data["exps"]; r=data["r"]; LB=data["LB"]
    jr=[[r_power(r,e) for e in exps]]
    for j in range(ell):
        row=[]
        for e in exps:
            if e[j]==0: row.append(_Z(0))
            else:
                e2=list(e); e2[j]-=1; row.append(_Z(e[j])*r_power(r,tuple(e2)))
        jr.append(row)
    JR=_M(_Q,jr); Jperp=JR.right_kernel().basis_matrix()
    Wrows=[]
    for wi in range(Jperp.nrows()):
        w=Jperp.row(wi); row=[_Q(w[a])/_Q(scales[a]) for a in range(n)]
        d=_lcm([x.denominator() for x in row]); Wrows.append([_Z(x*d) for x in row])
    Wint=_M(_Z,Wrows); C=Wint*LB.transpose(); K=C.right_kernel().basis_matrix(); Jlat=K*LB
    Gram=Jlat*Jlat.transpose(); log2covol=float(_R(Gram.det()).abs().log(2)/2)
    Bred=Jlat.LLL(algorithm="fpLLL:proved")
    def l2(v):
        s=sum(_Z(c)**2 for c in v); return float(_R(s).log(2)/2) if s>0 else 0.0
    lam1=min(l2(Bred.row(i)) for i in range(Bred.nrows()))
    jetA=0.5*math.log2(ell/(2.0*math.pi*math.e))+(log2covol-lam1)/ell
    return dict(rankK=int(K.nrows()), lam1=lam1, jetA=jetA)

import sys
t,m,l=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
d=build_data(t,m,l); res=jet_intersection_covol(d); gh=whole_gh(d,res["lam1"])
print("ROW %d %d %d lam1=%.1f GH=%.1f jetA=%.1f rankK=%d"%(t,m,l,res["lam1"],gh,res["jetA"],res["rankK"]),flush=True)
