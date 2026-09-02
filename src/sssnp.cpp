// ssSNPBLUP: o passo unico SEM G — os marcadores como equacoes (Liu et al. 2014).
//
// O modelo equivalente, escrito por extenso porque cada bloco abaixo vem de um termo dele:
//
//   u_g = Z g + p_g,   g ~ N(0, I (1-w) s2u / kd),   p_g ~ N(0, w s2u A22),
//   u_n | u_g ~ o condicional do pedigree,           kd = 2 sum p(1-p).
//
// A precisao conjunta de (u_n, u_g, g) sai da fatoracao p(u_n|u_g) p(u_g|g) p(g), e o
// termo A^gn (A^nn)^-1 A^ng = A^gg - A22^-1 (Schur) transforma tudo em blocos de A^-1:
//
//   Q s2u = [ A^nn   A^ng                          0
//             A^gn   A^gg + ((1-w)/w) A22^-1      -(1/w) A22^-1 Z
//             0     -(1/w) Z' A22^-1               (1/w) Z' A22^-1 Z + (kd/(1-w)) I ]
//
// A parte A^nn/A^ng/A^gg ja E a penalidade que monta_mme escreve com o A^-1 esparso; o
// resto sao os acrescimos deste arquivo. G nunca e formada nem invertida: o sistema e
// resolvido por gradientes conjugados precondicionados (Vandenplas et al. 2018, 2019), e
// cada aplicacao de A22^-1 usa a identidade de Masuda et al. (2017),
//
//   A22^-1 v = A^22 v - A^21 (A^11)^-1 A^12 v,
//
// com UMA fatoracao esparsa do bloco nao-genotipado A^11 feita no comeco e uma resolucao
// triangular por aplicacao. w -> 1 desliga os marcadores e colapsa no BLUP de pedigree —
// e esse colapso e um dos gates.
//
// Limites declarados: theta e DADO (isto e um resolvedor, como o BLUP com componentes
// fixas da pratica; a REML continua nos caminhos exatos), o termo genomico e um grupo
// escalar, e G* implicita = (1-w) Z Z'/kd + w A22, SEM o ajuste afim do caminho
// genotypes= — em populacoes fora do equilibrio os dois caminhos diferem por construcao.

#include "mme.h"

#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#ifndef FCONE
# define FCONE
#endif

namespace br {

namespace {

// Aplicacao simetrica de uma CSC triangular (superior OU inferior): cada entrada guardada
// (i, j, v) contribui v x_j na linha i e, fora da diagonal, v x_i na linha j.
void symv_tri(const Csc& a, const double* x, double* y) {
  for (std::size_t j = 0; j < a.ncol; j++)
    for (std::size_t k = a.colptr[j]; k < a.colptr[j + 1]; k++) {
      const std::size_t i = a.linha[k];
      const double v = a.valor[k];
      y[i] += v * x[j];
      if (i != j) y[j] += v * x[i];
    }
}

}  // namespace

SnpBlup snp_blup(const Desenho& d, Densa& mg, const std::vector<std::string>& geno_ids,
                 const std::vector<double>& theta, double w, std::size_t maxiter,
                 double tol, bool verboso) {
  SnpBlup r;
  if (!(w > 0.0 && w < 1.0))
    throw Erro("rpg (the residual polygenic proportion) must be in (0, 1): 0 leaves the "
               "polygenic residual without a distribution, 1 turns the markers off");

  // o grupo genomico: exatamente um grupo com parentesco, escalar
  std::size_t gpar = d.modelo.grupos.size();
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++)
    if (d.modelo.grupos[g].estrutura == Estrutura::Parentesco) {
      if (gpar != d.modelo.grupos.size())
        throw Erro("snp_blup expects ONE relationship group (declared limit of this version)");
      gpar = g;
    }
  if (gpar == d.modelo.grupos.size())
    throw Erro("there is no relationship term to receive the markers");
  if (d.modelo.grupos[gpar].dim != 1)
    throw Erro("the genomic term must be a scalar group in snp_blup (declared limit of "
               "this version)");
  const double var_u = theta[d.modelo.grupos[gpar].offset];
  if (!(var_u > 0.0)) throw Erro("the additive variance in theta is not positive");

  const DesenhoTermo* dt = nullptr;
  for (const DesenhoTermo& a : d.aleatorios)
    if (a.termo == d.modelo.grupos[gpar].termos[0]) { dt = &a; break; }
  if (!dt) throw Erro("relationship group without an incidence");
  const std::size_t nl = dt->n_niveis;

  // genotipados -> posicao de nivel (== posicao no A^-1 do grupo)
  std::vector<std::size_t> pos;
  pos.reserve(geno_ids.size());
  {
    std::vector<std::pair<std::string, std::size_t>> tmp;
    tmp.reserve(nl);
    for (std::size_t i = 0; i < nl; i++) tmp.push_back({dt->niveis[i], i});
    std::sort(tmp.begin(), tmp.end());
    for (const std::string& g : geno_ids) {
      auto it = std::lower_bound(tmp.begin(), tmp.end(), std::make_pair(g, std::size_t(0)),
          [](const auto& a, const auto& b){ return a.first < b.first; });
      if (it == tmp.end() || it->first != g)
        throw Erro("genotyped animal '" + g + "' is not in the pedigree");
      pos.push_back(it->second);
    }
  }
  const std::size_t ng = pos.size();

  // Z centrada (imputacao pela media, monomorficos fora), coluna-major para o dgemv
  const std::size_t mm = mg.ncol;
  if (mg.nlin != ng) throw Erro("genotype rows and ids with different lengths");
  r.usa_marcador.assign(mm, 1);
  std::vector<double> pfreq(mm, 0.0);
  double kd = 0.0;
  std::vector<std::size_t> usados;
  for (std::size_t j = 0; j < mm; j++) {
    double soma = 0.0;
    std::size_t k = 0;
    for (std::size_t i = 0; i < ng; i++) {
      const double x = mg.at(i, j);
      if (std::isfinite(x)) { soma += x; k++; }
    }
    if (k == 0) { r.usa_marcador[j] = 0; r.n_monomorficos++; continue; }
    const double media = soma / static_cast<double>(k);
    pfreq[j] = media / 2.0;
    if (pfreq[j] <= 0.0 || pfreq[j] >= 1.0) { r.usa_marcador[j] = 0; r.n_monomorficos++; continue; }
    kd += 2.0 * pfreq[j] * (1.0 - pfreq[j]);
    usados.push_back(j);
  }
  const std::size_t mu = usados.size();
  if (!(kd > 0.0)) throw Erro("2 sum p(1-p) is not positive: the markers do not vary");
  std::vector<double> zc(ng * mu);
  for (std::size_t jj = 0; jj < mu; jj++) {
    const std::size_t j = usados[jj];
    const double dp = 2.0 * pfreq[j];
    for (std::size_t i = 0; i < ng; i++) {
      double x = mg.at(i, j);
      if (!std::isfinite(x)) { x = dp; r.n_imputados++; }
      zc[i + ng * jj] = x - dp;
    }
  }

  // blocos de A^-1 pela mascara genotipado/nao: A11 esparsa fatorada uma vez, A12 e A22
  // como triplos para as aplicacoes
  const Csc& ainv = d.kinv[gpar];
  if (ainv.ncol != nl) throw Erro("A^-1 and the term with different sizes");
  std::vector<char> eh_geno(nl, 0);
  for (std::size_t i : pos) eh_geno[i] = 1;
  std::vector<std::size_t> loc(nl);
  std::size_t n1 = 0;
  for (std::size_t i = 0; i < nl; i++) if (!eh_geno[i]) loc[i] = n1++;
  for (std::size_t k = 0; k < ng; k++) loc[pos[k]] = k;

  std::vector<std::uint32_t> i11, j11;
  std::vector<double> v11;
  std::vector<std::uint32_t> i12, j12;      // linha = nao-genotipado local, coluna = genotipado local
  std::vector<double> v12;
  std::vector<std::uint32_t> i22, j22;
  std::vector<double> v22;
  std::vector<double> diag22(ng, 0.0);      // diagonal do BLOCO A^22, o substituto do precondicionador
  for (std::size_t c = 0; c < nl; c++)
    for (std::size_t k = ainv.colptr[c]; k < ainv.colptr[c + 1]; k++) {
      const std::size_t rr = ainv.linha[k];
      const double x = ainv.valor[k];
      const bool rg = eh_geno[rr], cg = eh_geno[c];
      if (!rg && !cg) {
        std::size_t i = loc[rr], j = loc[c];
        if (i < j) std::swap(i, j);
        i11.push_back(static_cast<std::uint32_t>(i));
        j11.push_back(static_cast<std::uint32_t>(j));
        v11.push_back(x);
      } else if (rg && cg) {
        i22.push_back(static_cast<std::uint32_t>(loc[rr]));
        j22.push_back(static_cast<std::uint32_t>(loc[c]));
        v22.push_back(x);
        if (rr == c) diag22[loc[rr]] += x;
      } else {
        i12.push_back(static_cast<std::uint32_t>(rg ? loc[c] : loc[rr]));
        j12.push_back(static_cast<std::uint32_t>(rg ? loc[rr] : loc[c]));
        v12.push_back(x);
      }
    }

  Csc a22blk = de_triplos(ng, ng, i22, j22, v22);
  Csc l11;
  std::vector<std::size_t> perm11;
  if (n1 > 0) {
    Csc a11 = de_triplos(n1, n1, i11, j11, v11);
    perm11 = grau_minimo(a11);
    Csc pa = permuta_sim(a11, perm11);
    Simbolica sb = simbolica(pa);
    if (!cholesky(pa, sb, l11))
      throw Erro("the non-genotyped block of A^-1 is not positive-definite");
  }

  // A22^-1 v = A^22 v - A^21 (A^11)^-1 A^12 v, sem nunca formar A22^-1
  std::vector<double> t1(n1), t2(n1);
  auto a22inv_vezes = [&](const std::vector<double>& x, std::vector<double>& y) {
    y.assign(ng, 0.0);
    symv_tri(a22blk, x.data(), y.data());
    if (n1 == 0) return;
    std::fill(t1.begin(), t1.end(), 0.0);
    for (std::size_t k = 0; k < v12.size(); k++) t1[i12[k]] += v12[k] * x[j12[k]];
    for (std::size_t i = 0; i < n1; i++) t2[i] = t1[perm11[i]];
    std::vector<double> px = resolve(l11, t2);
    for (std::size_t i = 0; i < n1; i++) t1[perm11[i]] = px[i];
    for (std::size_t k = 0; k < v12.size(); k++) y[j12[k]] -= v12[k] * t1[i12[k]];
  };

  // as equacoes de base, nas unidades de s2e como todo o resto do motor
  Montado M = monta_mme(d, theta);
  if (!M.ok) throw Erro("theta INADMISSIBLE: some covariance is not positive-definite");
  const double fac = M.s2e / var_u;
  std::vector<std::size_t> col_geno(ng);
  for (std::size_t k = 0; k < ng; k++) col_geno[k] = M.offset_grupo[gpar] + pos[k];

  const std::size_t N = M.total + mu;
  std::vector<double> rhs(N, 0.0);
  for (std::size_t k = 0; k < M.total; k++) rhs[k] = M.rhs[k];

  // precondicionador diagonal: diag(C) mais os acrescimos, com diag(A^22) como
  // substituto de diag(A22^-1) — e um majorante (o Schur so subtrai), entao o
  // escalonamento fica do lado conservador e a correcao vem das iteracoes
  std::vector<double> prec(N, 0.0);
  for (std::size_t j = 0; j < M.c.ncol; j++)
    for (std::size_t k = M.c.colptr[j]; k < M.c.colptr[j + 1]; k++)
      if (M.c.linha[k] == j) prec[j] += M.c.valor[k];
  for (std::size_t k = 0; k < ng; k++)
    prec[col_geno[k]] += fac * (1.0 - w) / w * diag22[k];
  for (std::size_t jj = 0; jj < mu; jj++) {
    double s = 0.0;
    for (std::size_t i = 0; i < ng; i++) {
      const double z = zc[i + ng * jj];
      s += diag22[i] * z * z;
    }
    prec[M.total + jj] = fac * (s / w + kd / (1.0 - w));
  }
  for (double& p : prec) if (!(p > 0.0)) p = 1.0;

  const int ngi = static_cast<int>(ng), mui = static_cast<int>(mu), inc1 = 1;
  const double um = 1.0, zero = 0.0;
  std::vector<double> vg(ng), zv(ng), q3(ng), q4(ng), gsum(ng), gm(mu);
  auto aplica = [&](const std::vector<double>& v, std::vector<double>& y) {
    y.assign(N, 0.0);
    symv_tri(M.c, v.data(), y.data());
    for (std::size_t k = 0; k < ng; k++) vg[k] = v[col_geno[k]];
    std::fill(zv.begin(), zv.end(), 0.0);
    if (mu > 0)
      F77_CALL(dgemv)("N", &ngi, &mui, &um, zc.data(), &ngi, v.data() + M.total, &inc1,
                      &zero, zv.data(), &inc1 FCONE);
    a22inv_vezes(vg, q3);
    a22inv_vezes(zv, q4);
    for (std::size_t k = 0; k < ng; k++)
      y[col_geno[k]] += fac * ((1.0 - w) * q3[k] - q4[k]) / w;
    if (mu > 0) {
      for (std::size_t k = 0; k < ng; k++) gsum[k] = q3[k] - q4[k];
      F77_CALL(dgemv)("T", &ngi, &mui, &um, zc.data(), &ngi, gsum.data(), &inc1,
                      &zero, gm.data(), &inc1 FCONE);
      for (std::size_t jj = 0; jj < mu; jj++)
        y[M.total + jj] += fac * (-gm[jj] / w + kd / (1.0 - w) * v[M.total + jj]);
    }
  };

  // PCG classico com o precondicionador diagonal
  std::vector<double> x(N, 0.0), res(rhs), z(N), pdir(N), q(N);
  double nrhs = 0.0;
  for (double v : rhs) nrhs += v * v;
  nrhs = std::sqrt(nrhs);
  if (nrhs == 0.0) { r.ok = true; r.convergiu = true; r.solucao.assign(N, 0.0); return r; }
  for (std::size_t k = 0; k < N; k++) z[k] = res[k] / prec[k];
  double rho = 0.0;
  for (std::size_t k = 0; k < N; k++) rho += res[k] * z[k];
  pdir = z;
  std::size_t it = 0;
  double rel = 1.0;
  for (; it < maxiter; it++) {
    R_CheckUserInterrupt();
    double nr = 0.0;
    for (double v : res) nr += v * v;
    rel = std::sqrt(nr) / nrhs;
    if (verboso && it > 0 && it % 100 == 0)
      Rprintf("PCG iter %d  relative residual %.3e\n", (int) it, rel);
    if (rel < tol) break;
    aplica(pdir, q);
    double pq = 0.0;
    for (std::size_t k = 0; k < N; k++) pq += pdir[k] * q[k];
    if (!(pq > 0.0)) throw Erro("the ssSNPBLUP system lost positive-definiteness in the PCG");
    const double alfa = rho / pq;
    for (std::size_t k = 0; k < N; k++) { x[k] += alfa * pdir[k]; res[k] -= alfa * q[k]; }
    for (std::size_t k = 0; k < N; k++) z[k] = res[k] / prec[k];
    double rho2 = 0.0;
    for (std::size_t k = 0; k < N; k++) rho2 += res[k] * z[k];
    const double beta = rho2 / rho;
    rho = rho2;
    for (std::size_t k = 0; k < N; k++) pdir[k] = z[k] + beta * pdir[k];
  }

  r.ok = true;
  r.convergiu = rel < tol;
  r.iters = it;
  r.residuo = rel;
  r.solucao.assign(x.begin(), x.begin() + M.total);
  r.efeitos.assign(mm, std::nan(""));
  for (std::size_t jj = 0; jj < mu; jj++) r.efeitos[usados[jj]] = x[M.total + jj];
  r.n_fixo = M.n_fixo;
  r.offset_grupo = M.offset_grupo;
  return r;
}

}  // namespace br
