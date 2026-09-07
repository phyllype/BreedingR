#include "estruturas.h"

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

Densa inv_do_fator(const Densa& l) {
  const std::size_t n = l.nlin;
  if (l.ncol != n) throw Erro("factor is not square");
  if (n == 0) return Densa();
  for (std::size_t i = 0; i < n; i++)
    if (l.at(i, i) == 0.0) throw Erro("singular factor: zero diagonal");

  // L^-1, coluna a coluna em rascunho contiguo
  Densa inv(n, n);
  std::vector<double> col(n, 0.0);
  for (std::size_t j = 0; j < n; j++) {
    for (std::size_t k = 0; k < j; k++) col[k] = 0.0;
    col[j] = 1.0 / l.at(j, j);
    for (std::size_t i = j + 1; i < n; i++) {
      const double* ri = l.linha(i);
      double acc = 0.0;
      for (std::size_t k = j; k < i; k++) acc += ri[k] * col[k];
      col[i] = -acc / ri[i];
    }
    for (std::size_t i = j; i < n; i++) inv.at(i, j) = col[i];
  }

  // S^-1 = (L^-1)' (L^-1), LADRILHADO.
  //
  // A versao linha a linha le o triangulo inferior inteiro da saida a cada linha de L^-1,
  // o que da n^3/2 doubles de trafego. A aritmetica e a mesma n^3/6 nos dois casos, logo o
  // custo nunca foi de conta: e de memoria. Tres ladrilhos ficam em L2 enquanto o laco
  // interno corre.
  const std::size_t B = 96;
  Densa out(n, n);
  for (std::size_t kb = 0; kb < n; kb += B) {
    const std::size_t kf = std::min(kb + B, n);
    for (std::size_t ib = 0; ib <= kb; ib += B) {
      const std::size_t iff = std::min(ib + B, n);
      for (std::size_t jb = 0; jb <= ib; jb += B) {
        const std::size_t jf = std::min(jb + B, n);
        for (std::size_t k = kb; k < kf; k++) {
          const double* rk = inv.linha(k);
          const std::size_t ihi = std::min(iff, k + 1);
          for (std::size_t i = ib; i < ihi; i++) {
            const double a = rk[i];
            if (a == 0.0) continue;
            const std::size_t jhi = std::min(jf, i + 1);
            double* dst = out.linha(i);
            for (std::size_t j = jb; j < jhi; j++) dst[j] += a * rk[j];
          }
        }
      }
    }
  }
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < i; j++) out.at(j, i) = out.at(i, j);
  return out;
}

Densa inv_pd(const Densa& s) {
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
