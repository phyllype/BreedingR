#include "estruturas.h"
#include "mme.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// O BLAS/LAPACK do proprio R: zero dependencia nova, e o usuario que trocar o BLAS do R
// acelera isto de graca. USE_FC_LEN_T + FCONE e o protocolo moderno de chamada Fortran.
#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
# define FCONE
#endif

namespace br {

Csc de_triplos(std::size_t nlin, std::size_t ncol,
               const std::vector<std::uint32_t>& li,
               const std::vector<std::uint32_t>& cj,
               const std::vector<double>& v) {
  if (li.size() != cj.size() || li.size() != v.size())
    throw Erro("triplets with different lengths");
  const std::size_t nz = li.size();
  for (std::size_t k = 0; k < nz; k++)
    if (li[k] >= nlin || cj[k] >= ncol)
      throw Erro("triplet outside the matrix");

  // passagem 1: estabiliza por LINHA
  std::vector<std::size_t> cnt(nlin + 1, 0);
  for (std::size_t k = 0; k < nz; k++) cnt[li[k] + 1]++;
  for (std::size_t r = 0; r < nlin; r++) cnt[r + 1] += cnt[r];
  std::vector<std::uint32_t> ordem(nz);
  {
    std::vector<std::size_t> pos = cnt;
    for (std::size_t k = 0; k < nz; k++) ordem[pos[li[k]]++] = static_cast<std::uint32_t>(k);
  }

  // passagem 2: estabiliza por COLUNA
  std::vector<std::size_t> cc(ncol + 1, 0);
  for (std::size_t k = 0; k < nz; k++) cc[cj[k] + 1]++;
  for (std::size_t c = 0; c < ncol; c++) cc[c + 1] += cc[c];
  const std::vector<std::size_t> lim = cc;
  std::vector<std::uint32_t> porcol(nz);
  {
    std::vector<std::size_t> pos = cc;
    for (std::size_t t = 0; t < nz; t++) {
      const std::uint32_t k = ordem[t];
      porcol[pos[cj[k]]++] = k;
    }
  }

  Csc out(nlin, ncol);
  out.linha.reserve(nz);
  out.valor.reserve(nz);
  for (std::size_t c = 0; c < ncol; c++) {
    std::size_t k = lim[c];
    while (k < lim[c + 1]) {
      const std::uint32_t r = li[porcol[k]];
      double acc = 0.0;
      while (k < lim[c + 1] && li[porcol[k]] == r) acc += v[porcol[k++]];
      out.linha.push_back(r);
      out.valor.push_back(acc);
    }
    out.colptr[c + 1] = out.linha.size();
  }
  return out;
}

// ------------------------------------------------------------------------ densa
//
// Os kernels densos chamam o LAPACK que o PROPRIO R carrega (R_ext/Lapack.h). Isso nao
// quebra a regra de dependencia: o BLAS/LAPACK vem com toda instalacao do R, e quem
// trocar o BLAS do R (OpenBLAS, MKL) acelera o pacote sem recompilar nada. O truque de
// layout: a Densa e linha-major, entao o triangulo INFERIOR daqui e o triangulo
// superior ('U') na convencao coluna-major do Fortran — mesma memoria, sem transpor.

bool chol_densa(Densa& s) {
  const std::size_t t = s.nlin;
  if (s.ncol != t) return false;
  if (t == 0) return true;
  // NaN dentro da matriz pode atravessar o dpotrf sem levantar info em algumas
  // implementacoes; o contrato daqui sempre foi devolver false, entao a varredura fica
  for (const double v : s.dados)
    if (!std::isfinite(v)) return false;
  int n = static_cast<int>(t), info = 0;
  F77_CALL(dpotrf)("U", &n, s.dados.data(), &n, &info FCONE);
  return info == 0;
}

// (L L')^-1 a partir de L, os dois no triangulo inferior COMPACTADO por coluna (coluna j =
// linhas j..n-1 em sequencia): L^-1 por blocos de 64 colunas, substituicao para frente com
// as 64 de uma vez, e Z = L^-T L^-1 por produtos internos de colunas contiguas de L^-1, um
// ladrilho de saida por thread. E o mesmo nucleo da inversa densa grande e da forma fechada
// do bloco denso final da inversa seletiva, onde o fator da cauda JA esta guardado assim.
void inversa_empacotada(const double* l, std::size_t n, double* z, int nth) {
  std::vector<std::size_t> off(n + 1, 0);
  for (std::size_t j = 0; j < n; j++) off[j + 1] = off[j] + (n - j);
  const std::size_t LB = 64;
  const std::size_t nbl = (n + LB - 1) / LB;
  std::vector<double> li(off[n], 0.0);
#ifdef _OPENMP
#pragma omp parallel num_threads(nth)
#endif
  {
    std::vector<double> x;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long bb = 0; bb < static_cast<long>(nbl); bb++) {
      const std::size_t j0 = static_cast<std::size_t>(bb) * LB, j1 = std::min(n, j0 + LB);
      const std::size_t w = j1 - j0;
      // X (linhas j0..n-1, w colunas, linha-major): L X = as colunas j0..j1-1 da identidade
      x.assign((n - j0) * w, 0.0);
      for (std::size_t c = 0; c < w; c++) x[c * w + c] = 1.0;
      for (std::size_t k = j0; k < n; k++) {
        const double* lk = l + off[k];          // coluna k de L, linhas k..n-1
        double* xk = &x[(k - j0) * w];
        const double dk = lk[0];
        for (std::size_t c = 0; c < w; c++) xk[c] /= dk;
        for (std::size_t r = k + 1; r < n; r++) {
          const double lrk = lk[r - k];
          if (lrk == 0.0) continue;
          double* xr = &x[(r - j0) * w];
          for (std::size_t c = 0; c < w; c++) xr[c] -= lrk * xk[c];
        }
      }
      for (std::size_t c = 0; c < w; c++) {
        const std::size_t j = j0 + c;
        double* dst = &li[off[j]];
        for (std::size_t r = j; r < n; r++) dst[r - j] = x[(r - j0) * w + c];
      }
    }
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
  for (long t = 0; t < static_cast<long>(nbl * (nbl + 1) / 2); t++) {
    std::size_t it = static_cast<std::size_t>((std::sqrt(8.0 * static_cast<double>(t) + 1.0) - 1.0) / 2.0);
    while (it * (it + 1) / 2 > static_cast<std::size_t>(t)) it--;
    while ((it + 1) * (it + 2) / 2 <= static_cast<std::size_t>(t)) it++;
    const std::size_t jt = static_cast<std::size_t>(t) - it * (it + 1) / 2;
    const std::size_t i0 = it * LB, i1 = std::min(n, i0 + LB);
    const std::size_t j0 = jt * LB, j1 = std::min(n, j0 + LB);
    for (std::size_t i = i0; i < i1; i++) {
      const double* ci = &li[off[i]];                    // L^-1(i.., i)
      const std::size_t jf = std::min(j1, i + 1);
      for (std::size_t j = j0; j < jf; j++) {
        const double* cj = &li[off[j]] + (i - j);        // L^-1(i.., j)
        // quatro acumuladores independentes: a soma unica e uma cadeia presa na latencia da
        // adicao; a ordem continua fixa, a mesma com qualquer numero de threads
        const std::size_t len = n - i;
        double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
        std::size_t k = 0;
        for (; k + 4 <= len; k += 4) {
          a0 += ci[k] * cj[k];
          a1 += ci[k + 1] * cj[k + 1];
          a2 += ci[k + 2] * cj[k + 2];
          a3 += ci[k + 3] * cj[k + 3];
        }
        for (; k < len; k++) a0 += ci[k] * cj[k];
        z[off[j] + (i - j)] = (a0 + a1) + (a2 + a3);
      }
    }
  }
}

// Inversa densa SPD em ladrilhos, em paralelo e determinista, para as matrizes grandes (a
// G* e a A22 do passo unico). O LAPACK de referencia que o R traz no Windows roda numa
// thread a ~1.5 GFlops, e as duas inversas n^3 dominavam o preparo com 3.000 genotipados.
// Cholesky pelo mesmo nucleo da cauda do fator esparso, no triangulo compactado por coluna,
// e depois inversa_empacotada(). Cada numero tem um dono e soma em ordem fixa: bit a bit o
// mesmo com qualquer numero de threads.
static Densa inv_pd_ladrilhos(const Densa& s) {
  const std::size_t n = s.nlin;
  std::vector<std::size_t> off(n + 1, 0);
  for (std::size_t j = 0; j < n; j++) off[j + 1] = off[j] + (n - j);
  std::vector<double> l(off[n]);
  for (std::size_t j = 0; j < n; j++)
    for (std::size_t i = j; i < n; i++) l[off[j] + (i - j)] = s.at(i, j);
  const int nth = threads();
  if (!cholesky_empacotada(l, n, nth)) throw Erro("the matrix is not positive-definite");
  std::vector<double> zp(off[n]);
  inversa_empacotada(l.data(), n, zp.data(), nth);
  Densa out(n, n);
  for (std::size_t j = 0; j < n; j++)
    for (std::size_t i = j; i < n; i++) {
      const double v = zp[off[j] + (i - j)];
      out.at(i, j) = v;
      out.at(j, i) = v;
    }
  return out;
}

Densa inv_pd(const Densa& s) {
  if (s.nlin >= 256 && s.ncol == s.nlin && !denso_lapack()) return inv_pd_ladrilhos(s);
  Densa l = s;
  if (!chol_densa(l)) throw Erro("the matrix is not positive-definite");
  int n = static_cast<int>(l.nlin), info = 0;
  if (n == 0) return l;
  // dpotri sobre o fator: devolve S^-1 no mesmo triangulo ('U' coluna-major = inferior
  // daqui); espelha para a matriz cheia, que e o contrato de sempre
  F77_CALL(dpotri)("U", &n, l.dados.data(), &n, &info FCONE);
  if (info != 0) throw Erro("the matrix is not positive-definite");
  for (std::size_t i = 0; i < l.nlin; i++)
    for (std::size_t j = i + 1; j < l.ncol; j++) l.at(i, j) = l.at(j, i);
  return l;
}

Densa inv_geral(const Densa& a) {
  const std::size_t t = a.nlin;
  if (a.ncol != t) throw Erro("matrix is not square");
  Densa m = a;
  if (t == 0) return m;
  for (const double v : m.dados)
    if (!std::isfinite(v)) throw Erro("singular matrix");
  // dgetrf + dgetri sobre a transposta implicita (linha-major lida como coluna-major da
  // pela transposta, e inv(A') = inv(A)' — a leitura linha-major do resultado ja e a
  // inversa certa)
  int n = static_cast<int>(t), info = 0;
  std::vector<int> piv(t);
  F77_CALL(dgetrf)(&n, &n, m.dados.data(), &n, piv.data(), &info);
  if (info != 0) throw Erro("singular matrix");
  int lwork = -1;
  double wq = 0.0;
  F77_CALL(dgetri)(&n, m.dados.data(), &n, piv.data(), &wq, &lwork, &info);
  lwork = static_cast<int>(wq);
  std::vector<double> work(static_cast<std::size_t>(std::max(lwork, 1)));
  F77_CALL(dgetri)(&n, m.dados.data(), &n, piv.data(), work.data(), &lwork, &info);
  if (info != 0) throw Erro("singular matrix");
  return m;
}

// Decomposicao espectral de uma simetrica pequena por rotacoes de Jacobi. Serve para
// quando a matriz pode ser SINGULAR e ainda assim precisa ser usada: Cholesky recusa,
// LU devolve lixo amplificado, e a pseudo-inversa truncada e a resposta certa. E a mesma
// escolha ja feita para o Gamma dos metafundadores.
void jacobi_sim(const Densa& a, std::vector<double>& ev, Densa& u) {
  const std::size_t n = a.nlin;
  Densa w = a;
  u = Densa(n, n);
  for (std::size_t i = 0; i < n; i++) u.at(i, i) = 1.0;
  for (int varr = 0; varr < 100; varr++) {
    double fora = 0.0;
    for (std::size_t p = 0; p < n; p++)
      for (std::size_t q = p + 1; q < n; q++) fora += w.at(p, q) * w.at(p, q);
    if (fora < 1e-30) break;
    for (std::size_t p = 0; p < n; p++)
      for (std::size_t q = p + 1; q < n; q++) {
        if (std::fabs(w.at(p, q)) < 1e-300) continue;
        const double th = 0.5 * (w.at(q, q) - w.at(p, p)) / w.at(p, q);
        const double tt = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1.0));
        const double c = 1.0 / std::sqrt(tt * tt + 1.0), s = tt * c;
        for (std::size_t k = 0; k < n; k++) {
          const double ak = w.at(p, k), bk = w.at(q, k);
          w.at(p, k) = c * ak - s * bk;
          w.at(q, k) = s * ak + c * bk;
        }
        for (std::size_t k = 0; k < n; k++) {
          const double ka = w.at(k, p), kb = w.at(k, q);
          w.at(k, p) = c * ka - s * kb;
          w.at(k, q) = s * ka + c * kb;
          const double ua = u.at(k, p), ub = u.at(k, q);
          u.at(k, p) = c * ua - s * ub;
          u.at(k, q) = s * ua + c * ub;
        }
      }
  }
  ev.assign(n, 0.0);
  for (std::size_t i = 0; i < n; i++) ev[i] = w.at(i, i);
}

double logdet_pd(const Densa& a) {
  Densa l = a;
  if (!chol_densa(l)) return std::nan("");
  double s = 0.0;
  for (std::size_t i = 0; i < l.nlin; i++) s += std::log(l.at(i, i));
  return 2.0 * s;
}

}  // namespace br
