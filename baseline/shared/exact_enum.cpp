// Exact, resumable sphere enumeration. No floating point exclusion decisions.
// Build: c++ -O3 -std=c++17 exact_enum.cpp -lgmpxx -lgmp -o exact_enum
#include <gmpxx.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using Z = mpz_class;

static Z floor_div(const Z &a, const Z &b) {
  if (b <= 0) throw std::runtime_error("Nonpositive divisor");
  Z q;
  mpz_fdiv_q(q.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return q;
}

static Z exact_div(const Z &a, const Z &b) {
  if (b <= 0 || !mpz_divisible_p(a.get_mpz_t(), b.get_mpz_t()))
    throw std::runtime_error("Nonexact integer division");
  Z q;
  mpz_divexact(q.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return q;
}

struct Data {
  int n;
  Z root_sq, content;
  std::vector<Z> inner, delta;
  std::vector<std::vector<Z>> lower;
};

static void read_header(std::istream &in, Data &d, const std::string &expected) {
  std::string magic;
  in >> magic >> d.n >> d.root_sq >> d.content;
  if (!in || magic != expected || d.n < 1 || d.n > 4096 ||
      d.root_sq <= 0 || d.content <= 0)
    throw std::runtime_error("Invalid preparation header");
  d.inner.resize(d.n);
  for (Z &v : d.inner) in >> v;
  if (!in) throw std::runtime_error("Truncated preparation");
}

static void prepare(const char *input, const char *output, int threads) {
#ifdef _OPENMP
  omp_set_num_threads(threads);
#else
  (void) threads;
#endif
  std::ifstream in(input);
  Data d;
  read_header(in, d, "EXACT-GMP-INPUT-1");
  std::vector<std::vector<Z>> a(d.n, std::vector<Z>(d.n));
  for (auto &row : a) for (Z &v : row) in >> v;
  if (!in) throw std::runtime_error("Truncated Gram matrix");
  for (int i = 0; i < d.n; ++i)
    for (int j = 0; j < i; ++j)
      if (a[i][j] != a[j][i]) throw std::runtime_error("Nonsymmetric Gram matrix");
  d.lower.assign(d.n, std::vector<Z>(d.n));
  d.delta.push_back(1);
  for (int k = 0; k < d.n; ++k) {
    const Z pivot = a[k][k];
    if (pivot <= 0) throw std::runtime_error("Gram matrix is not positive definite");
    d.delta.push_back(pivot);
    for (int i = k + 1; i < d.n; ++i) d.lower[i][k] = a[i][k];
    std::atomic<bool> division_failed(false);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) if(d.n - k > 32)
#endif
    for (int i = k + 1; i < d.n; ++i)
      for (int j = i; j < d.n; ++j) {
        // Each upper-triangular cell and its mirror have one writer. The
        // pivot row/column are immutable throughout this parallel region.
        Z numerator = pivot * a[i][j] - a[i][k] * a[k][j];
        if (!mpz_divisible_p(numerator.get_mpz_t(), d.delta[k].get_mpz_t())) {
          division_failed.store(true);
          continue;
        }
        Z value;
        mpz_divexact(value.get_mpz_t(), numerator.get_mpz_t(), d.delta[k].get_mpz_t());
        a[i][j] = value;
        a[j][i] = std::move(value);
      }
    if (division_failed.load()) throw std::runtime_error("Nonexact Bareiss division");
    if (k % 16 == 0) std::cerr << "Bareiss " << k + 1 << '/' << d.n << '\n';
  }
  std::ofstream out(output);
  out << "EXACT-GMP-PREP-1 " << d.n << ' ' << d.root_sq << ' ' << d.content << '\n';
  for (const Z &v : d.inner) out << v << ' ';
  out << '\n';
  for (const Z &v : d.delta) out << v << ' ';
  out << '\n';
  for (int i = 1; i < d.n; ++i) {
    for (int j = 0; j < i; ++j) out << d.lower[i][j] << ' ';
    out << '\n';
  }
  out.flush();
  if (!out) throw std::runtime_error("Could not save preparation");
  std::cout << "PREPARED " << d.n << '\n';
}

static Data load(const char *path) {
  std::ifstream in(path);
  Data d;
  read_header(in, d, "EXACT-GMP-PREP-1");
  d.delta.resize(d.n + 1);
  for (Z &v : d.delta) { in >> v; if (v <= 0) throw std::runtime_error("Invalid pivot"); }
  if (d.delta[0] != 1) throw std::runtime_error("Invalid initial pivot");
  d.lower.assign(d.n, std::vector<Z>(d.n));
  for (int i = 1; i < d.n; ++i)
    for (int j = 0; j < i; ++j) in >> d.lower[i][j];
  if (!in) throw std::runtime_error("Truncated Bareiss decomposition");
  return d;
}

struct Frame {
  int level;
  Z partial, center, lo, hi, left, right;
  std::vector<Z> tail;
  Z tightened_radius;
  bool tightened = false;
};

static Z radius_for(const Data &d, const Z &best) {
  return floor_div((best - 1) * d.root_sq, d.content);
}

static bool tighten(const Data &d, Frame &f, const Z &radius) {
  if (f.tightened && f.tightened_radius == radius)
    return f.lo <= f.hi && (f.left >= f.lo || f.right <= f.hi);
  Z budget = d.delta[f.level] * (radius * d.delta[f.level + 1] - f.partial);
  if (budget < 0) return false;
  Z t;
  mpz_sqrt(t.get_mpz_t(), budget.get_mpz_t());
  Z low = -floor_div(t + f.center, d.delta[f.level + 1]);
  Z high = floor_div(t - f.center, d.delta[f.level + 1]);
  f.lo = std::max(f.lo, low);
  f.hi = std::min(f.hi, high);
  f.left = std::min(f.left, f.hi);
  f.right = std::max(f.right, f.lo);
  f.tightened_radius = radius;
  f.tightened = true;
  return f.lo <= f.hi && (f.left >= f.lo || f.right <= f.hi);
}

static bool make_frame(const Data &d, int level, std::vector<Z> tail,
                       const Z &partial, const Z &radius, bool positive, Frame &f) {
  f.level = level;
  f.partial = partial;
  f.tail = std::move(tail);
  f.center = 0;
  for (int j = level + 1; j < d.n; ++j)
    f.center += d.lower[j][level] * f.tail[j - level - 1];
  Z budget = d.delta[level] * (radius * d.delta[level + 1] - partial);
  if (budget < 0) return false;
  Z t;
  mpz_sqrt(t.get_mpz_t(), budget.get_mpz_t());
  f.lo = -floor_div(t + f.center, d.delta[level + 1]);
  f.hi = floor_div(t - f.center, d.delta[level + 1]);
  if (positive) f.lo = std::max(f.lo, Z(1));
  if (f.lo > f.hi) return false;
  Z midpoint = floor_div(-f.center, d.delta[level + 1]);
  f.left = std::min(midpoint, f.hi);
  f.right = std::max(Z(midpoint + 1), f.lo);
  f.tightened_radius = radius;
  f.tightened = true;
  return true;
}

static void write_frame(std::ostream &out, const Frame &f) {
  out << f.level << ' ' << f.partial << ' ' << f.center << ' '
      << f.lo << ' ' << f.hi << ' ' << f.left << ' ' << f.right;
  for (const Z &v : f.tail) out << ' ' << v;
  out << '\n';
}

static Frame read_frame(std::istream &in, const Data &d) {
  Frame f;
  in >> f.level >> f.partial >> f.center >> f.lo >> f.hi >> f.left >> f.right;
  if (!in || f.level < 0 || f.level >= d.n || f.partial < 0)
    throw std::runtime_error("Invalid frontier frame");
  f.tail.resize(d.n - f.level - 1);
  Z center = 0;
  for (int j = f.level + 1; j < d.n; ++j) {
    in >> f.tail[j - f.level - 1];
    center += d.lower[j][f.level] * f.tail[j - f.level - 1];
  }
  if (!in || center != f.center) throw std::runtime_error("Frontier center mismatch");
  return f;
}

static void serve(const char *path) {
  const Data d = load(path);
  std::cout << "READY " << d.n << std::endl;
  std::string command;
  while (std::cin >> command) {
    if (command == "QUIT") return;
    if (command == "INIT") {
      Z best;
      std::cin >> best;
      if (best <= 0) throw std::runtime_error("Invalid incumbent");
      std::vector<Frame> frames;
      for (int i = d.n - 1; i >= 0; --i) {
        Frame f;
        if (make_frame(d, i, std::vector<Z>(d.n - i - 1), Z(0), radius_for(d, best), true, f))
          frames.push_back(std::move(f));
      }
      std::cout << "FRAMES " << frames.size() << '\n';
      for (const Frame &f : frames) write_frame(std::cout, f);
      std::cout << "END" << std::endl;
      continue;
    }
    if (command != "JOB") throw std::runtime_error("Unknown worker command");
    Z best;
    unsigned long long limit, nodes = 0, leaves = 0;
    double seconds;
    size_t count;
    std::cin >> best >> limit >> seconds >> count;
    if (!std::cin || best <= 0 || limit == 0 || seconds <= 0 || count > 1000000)
      throw std::runtime_error("Invalid job header");
    std::vector<Frame> stack;
    for (size_t i = 0; i < count; ++i) stack.push_back(read_frame(std::cin, d));
    std::vector<Z> best_y;
    Z radius = radius_for(d, best);
    const auto began = std::chrono::steady_clock::now();
    while (!stack.empty() && nodes < limit) {
      if (nodes % 256 == 0 && std::chrono::duration<double>(
          std::chrono::steady_clock::now() - began).count() >= seconds) break;
      Frame &f = stack.back();
      if (!tighten(d, f, radius)) { stack.pop_back(); continue; }
      Z x;
      if (f.left < f.lo) x = f.right++;
      else if (f.right > f.hi) x = f.left--;
      else {
        Z lt = d.delta[f.level + 1] * f.left + f.center;
        Z rt = d.delta[f.level + 1] * f.right + f.center;
        if (abs(lt) <= abs(rt)) x = f.left--;
        else x = f.right++;
      }
      ++nodes;
      Z t = d.delta[f.level + 1] * x + f.center;
      Z partial = exact_div(d.delta[f.level] * f.partial + t * t, d.delta[f.level + 1]);
      if (partial > radius * d.delta[f.level]) continue;
      std::vector<Z> y;
      y.reserve(f.tail.size() + 1);
      y.push_back(x);
      y.insert(y.end(), f.tail.begin(), f.tail.end());
      if (f.level > 0) {
        Frame child;
        if (make_frame(d, f.level - 1, std::move(y), partial, radius, false, child))
          stack.push_back(std::move(child));
      } else {
        ++leaves;
        Z dot = 0;
        for (int i = 0; i < d.n; ++i) dot += y[i] * d.inner[i];
        Z k = floor_div(2 * dot + d.root_sq, 2 * d.root_sq);
        Z residue = dot - k * d.root_sq;
        Z value = exact_div(d.content * partial + residue * residue, d.root_sq);
        if (value < best) {
          best = value;
          best_y = std::move(y);
          radius = radius_for(d, best);
        }
      }
    }
    std::cout << "RESULT " << best << ' ' << nodes << ' ' << leaves << ' ' << stack.size()
              << ' ' << best_y.size() << '\n';
    if (!best_y.empty()) {
      for (const Z &v : best_y) std::cout << v << ' ';
      std::cout << '\n';
    }
    for (const Frame &f : stack) write_frame(std::cout, f);
    std::cout << "END" << std::endl;
  }
}

int main(int argc, char **argv) {
  try {
    if ((argc == 4 || argc == 5) && std::string(argv[1]) == "prepare") {
      int threads = argc == 5 ? std::stoi(argv[4]) : 1;
      if (threads < 1) throw std::runtime_error("Invalid preparation thread count");
      prepare(argv[2], argv[3], threads);
    }
    else if (argc == 3 && std::string(argv[1]) == "serve") serve(argv[2]);
    else throw std::runtime_error("Usage: exact_enum prepare INPUT OUTPUT [THREADS] | serve PREPARED");
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 2;
  }
}
