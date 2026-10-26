// PEGS: predicao genomica MULTICARACTER por Gauss-Seidel aleatorizado (Xavier e Habier,
// 2022, Genet. Sel. Evol. 54:45; implementacao de referencia bWGR::mrr). Porte do motor
// validado em Julia deste projeto, que por sua vez foi conferido contra a listagem Rcpp dos
// autores.
//
// Modelo (SNP-BLUP multivariado): para o caracter t, y_t = 1 mu_t + X b_t + e_t, os k efeitos
// de cada marcador conjuntamente b_j ~ N(0, Sb) e os residuos e_t ~ N(0, I s2e_t) SEM
// correlacao entre caracteres. Em vez de fatorar o sistema de ordem p k, os k efeitos de UM
// marcador sao resolvidos juntos, com atualizacao do residuo, varrendo os marcadores numa
// ordem SORTEADA a cada passada:
//
//   [Sb^-1 + diag_t(x_j'x_j / s2e_t)] b_j = diag_t(1/s2e_t) (x_j'e + (x_j'x_j) b_j^antigo)
//   e <- e - x_j (b_j^novo - b_j^antigo)          (so nas celulas observadas)
//
// As variancias sao reestimadas DENTRO do laco pela pseudo-expectativa (PE),
// s2b_tt = (b_t' X'y_t) / tr(X'M_t X), covariancias simetrizadas, e Sb volta ao cone das
// positivas-definidas por "bending" (Hayes e Hill, 1981) quando sai dele. Estruturas de
// covariancia opcionais (Xavier et al., 2025): simetria composta heterogenea (hcs) e fator
// analitico estendido (xfa), aplicadas depois da PE e antes do bending.
//
// O LIMITE que viaja junto (medido em dado real neste projeto): residuos NAO correlacionados
// entre caracteres. E um metodo para o mesmo caracter em AMBIENTES diferentes. Em caracteres
// do MESMO registro a correlacao residual nao tem onde entrar e vaza para a genetica (numero
// nascido x peso da leitegada: rg 0.83 aqui contra 0.33 +- 0.08 do REML bivariado).

#include "mme.h"
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

namespace br {

namespace {

// Cholesky de um sistema k x k (k = caracteres) e resolucao, rhs sobrescrito pela solucao.
// Pivo nao positivo vira 1e-12, como na referencia.
void resolve_pequeno(const std::vector<double>& A, std::vector<double>& rhs, std::size_t k,
                     std::vector<double>& L) {
  for (std::size_t j = 0; j < k; j++) {
    double s = A[j * k + j];
    for (std::size_t c = 0; c < j; c++) s -= L[j * k + c] * L[j * k + c];
    s = s > 0.0 ? s : 1e-12;
    const double ljj = std::sqrt(s);
    L[j * k + j] = ljj;
    for (std::size_t i = j + 1; i < k; i++) {
      double t = A[i * k + j];
      for (std::size_t c = 0; c < j; c++) t -= L[i * k + c] * L[j * k + c];
      L[i * k + j] = t / ljj;
    }
  }
  for (std::size_t i = 0; i < k; i++) {
    double t = rhs[i];
    for (std::size_t c = 0; c < i; c++) t -= L[i * k + c] * rhs[c];
    rhs[i] = t / L[i * k + i];
  }
  for (std::size_t i = k; i-- > 0;) {
    double t = rhs[i];
    for (std::size_t c = i + 1; c < k; c++) t -= L[c * k + i] * rhs[c];
    rhs[i] = t / L[i * k + i];
  }
}

bool positiva_definida(const Densa& a) {
  Densa c = a;
  return chol_densa(c);
}

}  // namespace

// Estrutura de covariancia sobre a estimativa de Sb, no lugar (0 livre, 1 hcs, 2 xfa).
void estrutura_pegs(Densa& vb, int tipo, std::size_t nfat) {
  const std::size_t k = vb.nlin;
  if (tipo == 0 || k < 2) return;
  if (tipo == 1) {
    // uma correlacao comum a todos os pares, cada caracter com a propria variancia
    std::vector<double> d(k);
    for (std::size_t t = 0; t < k; t++) d[t] = std::sqrt(std::max(vb.at(t, t), 1e-12));
    double s = 0.0;
    std::size_t cnt = 0;
    for (std::size_t a = 0; a < k; a++)
      for (std::size_t c = a + 1; c < k; c++) { s += vb.at(a, c) / (d[a] * d[c]); cnt++; }
    const double rho = cnt ? std::min(0.999, std::max(-0.999, s / static_cast<double>(cnt))) : 0.0;
    for (std::size_t a = 0; a < k; a++)
      for (std::size_t c = 0; c < k; c++)
        if (a != c) vb.at(a, c) = rho * d[a] * d[c];
    return;
  }
  if (tipo == 2) {
    // q fatores pela decomposicao espectral truncada, e a diagonal original de volta: a
    // variancia de cada caracter (e o h2) fica, e o corte cai todo nas covariancias
    const std::size_t q = std::min(k, std::max<std::size_t>(1, nfat));
    std::vector<double> dg(k);
    for (std::size_t t = 0; t < k; t++) dg[t] = vb.at(t, t);
    std::vector<double> val;
    Densa u;
    jacobi_sim(vb, val, u);                          // vetores nas colunas de u
    std::vector<std::size_t> ord(k);
    for (std::size_t f = 0; f < k; f++) ord[f] = f;
    std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t c) { return val[a] > val[c]; });
    Densa g(k, k);
    for (std::size_t a = 0; a < k; a++)
      for (std::size_t c = 0; c < k; c++) {
        double s = 0.0;
        for (std::size_t f = 0; f < q; f++)
          s += u.at(a, ord[f]) * std::max(val[ord[f]], 0.0) * u.at(c, ord[f]);
        g.at(a, c) = s;
      }
    for (std::size_t t = 0; t < k; t++) g.at(t, t) = dg[t];
    for (std::size_t a = 0; a < k; a++)
      for (std::size_t c = 0; c < k; c++) vb.at(a, c) = 0.5 * (g.at(a, c) + g.at(c, a));
    return;
  }
  throw Erro("cov_structure must be unstructured, hcs or xfa");
}

ResultadoPegs pegs(const Densa& y_in, const Genotipos& gt, std::size_t maxit, double tol,
                   double deflate_min, bool atualiza_vc, const Densa* vb0,
                   const std::vector<double>* ve0, int estrutura, std::size_t nfat) {
  const std::size_t n = y_in.nlin, k = y_in.ncol, p = gt.m;
  if (gt.n != n) throw Erro("genotype rows and phenotype rows with different counts");
  if (k < 1) throw Erro("pegs needs at least one trait");
  ResultadoPegs R;
  // X lida do objeto do R, sem copia: o marcador j e a coluna j, o ausente na media dela. Com
  // a matriz double e o marcador sem ausente, o ponteiro vai direto para a memoria do R.
  std::vector<double> media(p, 0.0);
  std::vector<char> inteira(p, 0);
  {
    std::vector<double> col(n);
    for (std::size_t j = 0; j < p; j++) {
      gt.coluna(j, col.data());
      double soma = 0.0;
      std::size_t c = 0;
      for (std::size_t i = 0; i < n; i++)
        if (std::isfinite(col[i])) { soma += col[i]; c++; }
      if (c == 0) throw Erro("a marker has no observed genotype");
      media[j] = soma / static_cast<double>(c);
      inteira[j] = c == n;
      R.n_imputados += n - c;
    }
  }
  std::vector<double> xbuf(n);
  auto coluna_x = [&](std::size_t j) -> const double* {
    if (gt.d && inteira[j]) return gt.d + j * n;
    gt.coluna(j, xbuf.data());
    for (double& v : xbuf) if (!std::isfinite(v)) v = media[j];
    return xbuf.data();
  };
  // Z (1 observado, 0 ausente), y centrado e mascarado, coluna-major por caracter
  std::vector<double> z(n * k, 0.0), y(n * k, 0.0), nv(k, 0.0), mu(k, 0.0);
  for (std::size_t t = 0; t < k; t++) {
    for (std::size_t i = 0; i < n; i++) {
      const double v = y_in.at(i, t);
      if (std::isfinite(v)) { z[t * n + i] = 1.0; mu[t] += v; nv[t] += 1.0; }
    }
    if (nv[t] < 2.0) throw Erro("each trait needs at least 2 observed records");
    mu[t] /= nv[t];
    for (std::size_t i = 0; i < n; i++)
      if (z[t * n + i] > 0.0) y[t * n + i] = y_in.at(i, t) - mu[t];
  }
  // somas de quadrados dos marcadores dentro dos observados de cada caracter, e tr(X'M X)
  std::vector<double> xx(p * k), trxsx(k, 0.0), msx(k, 0.0);
  for (std::size_t j = 0; j < p; j++) {
    const double* xj = coluna_x(j);
    for (std::size_t t = 0; t < k; t++) {
      double s2 = 0.0, s1 = 0.0;
      const double* zt = &z[t * n];
      for (std::size_t i = 0; i < n; i++) { s2 += xj[i] * xj[i] * zt[i]; s1 += xj[i] * zt[i]; }
      xx[j * k + t] = s2;
      msx[t] += s2 / nv[t] - (s1 / nv[t]) * (s1 / nv[t]);
    }
  }
  for (std::size_t t = 0; t < k; t++) {
    trxsx[t] = nv[t] * msx[t];
    if (!(trxsx[t] > 0.0)) throw Erro("a trait has no marker variation among its records");
  }
  std::vector<double> vy(k, 0.0), ve(k), ive(k);
  for (std::size_t t = 0; t < k; t++) {
    for (std::size_t i = 0; i < n; i++) vy[t] += y[t * n + i] * y[t * n + i];
    vy[t] /= (nv[t] - 1.0);
    ve[t] = ve0 ? (*ve0)[t] : 0.5 * vy[t];
    ive[t] = 1.0 / ve[t];
  }
  Densa vb(k, k);
  if (vb0) vb = *vb0;
  else for (std::size_t t = 0; t < k; t++) vb.at(t, t) = ve[t] / msx[t];
  Densa ig = inv_pd(vb);
  // X'y (p x k)
  std::vector<double> til(p * k, 0.0);
  for (std::size_t j = 0; j < p; j++) {
    const double* xj = coluna_x(j);
    for (std::size_t t = 0; t < k; t++) {
      double s = 0.0;
      for (std::size_t i = 0; i < n; i++) s += xj[i] * y[t * n + i];
      til[j * k + t] = s;
    }
  }
  std::vector<double> b(p * k, 0.0), e = y;
  std::vector<std::size_t> ordem(p);
  for (std::size_t j = 0; j < p; j++) ordem[j] = j;
  std::vector<double> lhs(k * k), lbuf(k * k, 0.0), rhs(k), b0(k);
  double deflate = 1.0;
  const double logtol = std::log10(tol);
  std::size_t iters = 0;
  bool convergiu = false;
  GetRNGstate();
  try {
  while (iters < maxit) {
    // Fisher-Yates com o gerador do R: set.seed() governa, e o sorteio e o mesmo em toda
    // plataforma (o std::shuffle nao e: libstdc++ e libc++ sorteiam diferente)
    for (std::size_t j = p; j > 1; j--) {
      const std::size_t r = static_cast<std::size_t>(unif_rand() * static_cast<double>(j));
      std::swap(ordem[j - 1], ordem[std::min(r, j - 1)]);
    }
    double sumsq = 0.0;
    for (std::size_t jj = 0; jj < p; jj++) {
      const std::size_t J = ordem[jj];
      const double* xj = coluna_x(J);
      for (std::size_t t = 0; t < k; t++) b0[t] = b[J * k + t];
      for (std::size_t r = 0; r < k; r++)
        for (std::size_t c = 0; c < k; c++) lhs[r * k + c] = ig.at(r, c);
      for (std::size_t t = 0; t < k; t++) lhs[t * k + t] += xx[J * k + t] * ive[t];
      for (std::size_t t = 0; t < k; t++) {
        double s = 0.0;
        const double* et = &e[t * n];
        for (std::size_t i = 0; i < n; i++) s += xj[i] * et[i];
        rhs[t] = (s + xx[J * k + t] * b0[t]) * ive[t];
      }
      resolve_pequeno(lhs, rhs, k, lbuf);
      for (std::size_t t = 0; t < k; t++) {
        const double d = rhs[t] - b0[t];
        sumsq += d * d;
        b[J * k + t] = rhs[t];
        if (d != 0.0) {
          double* et = &e[t * n];
          const double* zt = &z[t * n];
          for (std::size_t i = 0; i < n; i++) et[i] -= xj[i] * d * zt[i];
        }
      }
    }
    iters++;
    if (atualiza_vc) {
      for (std::size_t t = 0; t < k; t++) {
        double s = 0.0;
        for (std::size_t i = 0; i < n; i++) s += e[t * n + i] * y[t * n + i];
        ve[t] = s / (nv[t] - 1.0);
        ive[t] = 1.0 / ve[t];
      }
      Densa th(k, k);                      // b' X'y
      for (std::size_t a = 0; a < k; a++)
        for (std::size_t c = 0; c < k; c++) {
          double s = 0.0;
          for (std::size_t j = 0; j < p; j++) s += b[j * k + a] * til[j * k + c];
          th.at(a, c) = s;
        }
      for (std::size_t t = 0; t < k; t++) vb.at(t, t) = th.at(t, t) / trxsx[t];
      for (std::size_t a = 0; a < k; a++)
        for (std::size_t c = 0; c < k; c++)
          if (a != c) vb.at(a, c) = (th.at(a, c) + th.at(c, a)) / (trxsx[a] + trxsx[c]);
      estrutura_pegs(vb, estrutura, nfat);
      Densa a = vb;
      for (std::size_t r = 0; r < k; r++)
        for (std::size_t c = 0; c < k; c++) if (r != c) a.at(r, c) *= deflate;
      while (!positiva_definida(a) && deflate > deflate_min) {
        deflate -= 0.01;
        for (std::size_t r = 0; r < k; r++)
          for (std::size_t c = 0; c < k; c++) if (r != c) a.at(r, c) = vb.at(r, c) * deflate;
      }
      ig = inv_pd(a);
    }
    if (std::log10(sumsq + 1e-300) < logtol) { convergiu = true; break; }
    if (iters % 16 == 0) R_CheckUserInterrupt();
  }
  } catch (...) {
    PutRNGstate();                     // o gerador do R guarda o quanto andou, mesmo no erro
    throw;
  }
  PutRNGstate();

  R.mu = mu;
  R.b = Densa(p, k);
  for (std::size_t j = 0; j < p; j++)
    for (std::size_t t = 0; t < k; t++) R.b.at(j, t) = b[j * k + t];
  R.gebv = Densa(n, k);
  for (std::size_t t = 0; t < k; t++)
    for (std::size_t j = 0; j < p; j++) {
      const double bj = b[j * k + t];
      if (bj == 0.0) continue;
      const double* xj = coluna_x(j);
      for (std::size_t i = 0; i < n; i++) R.gebv.at(i, t) += xj[i] * bj;
    }
  R.h2.resize(k);
  for (std::size_t t = 0; t < k; t++) R.h2[t] = 1.0 - ve[t] / vy[t];
  R.vb = vb;
  R.ve = ve;
  R.deflate = deflate;
  R.iters = iters;
  R.convergiu = convergiu;
  return R;
}

}  // namespace br
