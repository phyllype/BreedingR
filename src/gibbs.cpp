// Gibbs (Geman e Geman, 1984; Wang, Rutledge e Gianola, 1993, 1994) para o modelo
// animal (uma caracteristica): o lado bayesiano do pacote.
//
// ## O desenho
//
// Bloco UNICO de localizacao: (b, u) | theta, y e Gaussiana com media na solucao das MME
// e covariancia C^-1 — amostrada exata via a MESMA Cholesky esparsa do REML (fator L da
// C permutada; media + sqrt(s2e) P L'^-1 z, z ~ N(0, I)). Um refactor numerico por
// iteracao, simbolica em cache. Single-site misturaria pior e jogaria fora a maquinaria
// ja validada; o bloco unico faz da media amostral o proprio BLUP, que e o gate.
//
// Condicionais conjugadas com prior PLANO (documentado; nu0 = 0, S0 = 0):
//   C_g | u  ~ InvWishart(nl_g, U' K^-1 U)          (Bartlett para a Wishart)
//   s2e | e  ~ e'e / chisq(n)
//
// ## RNG
//
// O do PROPRIO R (GetRNGstate / Rf_rnorm / Rf_rchisq): respeita set.seed() do usuario e
// nao adiciona dependencia. Por isso este arquivo, ao contrario do resto da numerica,
// inclui R.h — ele E a fronteira.

#include "mme.h"
#include <R_ext/Print.h>
#include <R_ext/Utils.h>

#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>

namespace br {

static std::size_t vech_idx(std::size_t i, std::size_t j, std::size_t d) {
  if (i < j) std::swap(i, j);
  return (j * (2 * d - j + 1)) / 2 + (i - j);
}

// x tal que L' x = z, com L triangular inferior em CSC (acesso por coluna e natural:
// a coluna j de L e a linha j de L').
static std::vector<double> resolve_lt(const Csc& L, const std::vector<double>& z) {
  const std::size_t n = L.ncol;
  std::vector<double> x(z);
  for (std::size_t jj = n; jj-- > 0;) {
    double s = x[jj];
    double diag = 0.0;
    for (std::size_t p = L.colptr[jj]; p < L.colptr[jj + 1]; p++) {
      const std::size_t i = L.linha[p];
      if (i == jj) diag = L.valor[p];
      else s -= L.valor[p] * x[i];
    }
    x[jj] = s / diag;
  }
  return x;
}

// InvWishart(nu, S) via Bartlett na Wishart(nu, S^-1): devolve C = (A A')^-1 com
// A = L T, S^-1 = L L'. Exige nu >= d.
static Densa rinvwishart(double nu, const Densa& S) {
  const std::size_t d = S.nlin;
  Densa Sinv = inv_pd(S);
  Densa Lf = Sinv;
  if (!chol_densa(Lf)) throw Erro("the matrix is not positive-definite");
  for (std::size_t i = 0; i < d; i++)
    for (std::size_t j = i + 1; j < d; j++) Lf.at(i, j) = 0.0;
  Densa T(d, d);
  for (std::size_t i = 0; i < d; i++) {
    T.at(i, i) = std::sqrt(Rf_rchisq(nu - static_cast<double>(i)));
    for (std::size_t j = 0; j < i; j++) T.at(i, j) = Rf_rnorm(0.0, 1.0);
  }
  Densa A(d, d);
  for (std::size_t i = 0; i < d; i++)
    for (std::size_t j = 0; j <= i; j++) {
      double s = 0.0;
      for (std::size_t k = j; k <= i; k++) s += Lf.at(i, k) * T.at(k, j);
      A.at(i, j) = s;
    }
  Densa W(d, d);
  for (std::size_t i = 0; i < d; i++)
    for (std::size_t j = 0; j <= i; j++) {
      double s = 0.0;
      for (std::size_t k = 0; k <= std::min(i, j); k++) s += A.at(i, k) * A.at(j, k);
      W.at(i, j) = s;
      W.at(j, i) = s;
    }
  return inv_pd(W);
}

GibbsSaida gibbs(const Desenho& d, std::size_t n_iter, std::size_t burnin,
                 std::size_t thin, bool loc_fixa, const std::vector<double>* theta_fixo,
                 bool verboso) {
  GibbsSaida S;
  const std::size_t ntheta = d.modelo.ntheta;

  // partida identica a do REML: var(y) repartida
  std::vector<double> theta = partida(d);
  if (theta_fixo) theta = *theta_fixo;

  GetRNGstate();
  CacheSimbolica cs;

  const std::size_t n_grupos = d.modelo.grupos.size();
  std::vector<std::size_t> nl_g(n_grupos, 0);
  for (std::size_t g = 0; g < n_grupos; g++)
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == d.modelo.grupos[g].termos[0]) { nl_g[g] = a.n_niveis; break; }

  std::size_t n_amostras = 0;
  std::size_t total_cols = 0;
  std::vector<double> media_u, m2_u;   // Welford (1962) para media e variancia das localizacoes

  for (std::size_t it = 0; it < n_iter; it++) {
    R_CheckUserInterrupt();
    if (verboso && n_iter >= 10 && (it + 1) % (n_iter / 10) == 0)
      Rprintf("gibbs %d/%d\n", (int) (it + 1), (int) n_iter);
    // ---- 1. localizacoes: (b, u) ~ N(solucao, s2e C_s^-1)
    Montado M = monta_mme(d, theta);
    if (!M.ok) { S.mensagem = "theta INADMISSIBLE inside the chain"; break; }
    if (total_cols == 0) {
      total_cols = M.total;
      media_u.assign(M.total, 0.0);
      m2_u.assign(M.total, 0.0);
    }
    if (!cs.pronto) {
      cs.perm = grau_minimo(M.c);
      Csc pc0 = permuta_sim(M.c, cs.perm);
      cs.sb = simbolica(pc0);
      cs.pronto = true;
    }
    Csc pc = permuta_sim(M.c, cs.perm);
    Csc L;
    if (!cholesky(pc, cs.sb, L)) { S.mensagem = "C is not positive-definite inside the chain"; break; }

    std::vector<double> pb(M.total);
    for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[cs.perm[k]];
    std::vector<double> sol_p = resolve(L, pb);

    std::vector<double> z(M.total);
    for (std::size_t k = 0; k < M.total; k++) z[k] = Rf_rnorm(0.0, 1.0);
    std::vector<double> ruido_p = resolve_lt(L, z);

    const double s2e = theta[d.modelo.offset_residual];
    std::vector<double> loc(M.total);
    for (std::size_t k = 0; k < M.total; k++)
      loc[cs.perm[k]] = sol_p[k] + std::sqrt(s2e) * ruido_p[k];

    // ---- 2. residuo: e = y - W loc; s2e ~ e'e / chisq(n)   (prior plano)
    double ee = 0.0;
    std::size_t nreg = 0;
    {
      std::vector<double> wl(d.nlin, 0.0);
      for (std::size_t r = 0; r < d.nlin; r++) {
        if (!d.usa[r]) continue;
        double s = 0.0;
        for (std::size_t j = 0; j < d.x.ncol; j++) s += d.x.at(r, j) * loc[j];
        wl[r] = s;
      }
      std::size_t col0 = M.n_fixo;
      for (std::size_t g = 0; g < n_grupos; g++) {
        for (std::size_t tm : d.modelo.grupos[g].termos)
          for (const DesenhoTermo& a : d.aleatorios)
            if (a.termo == tm) {
              for (std::size_t c = 0; c < a.z.ncol; c++) {
                const double bc = loc[col0 + c];
                if (bc == 0.0) continue;
                for (std::size_t p = a.z.colptr[c]; p < a.z.colptr[c + 1]; p++)
                  if (d.usa[a.z.linha[p]]) wl[a.z.linha[p]] += a.z.valor[p] * bc;
              }
              col0 += a.z.ncol;
            }
      }
      for (std::size_t r = 0; r < d.nlin; r++) {
        if (!d.usa[r]) continue;
        const double e = d.y[r] - wl[r];
        ee += e * e;
        nreg++;
      }
    }
    if (!loc_fixa)
      theta[d.modelo.offset_residual] = ee / Rf_rchisq(static_cast<double>(nreg));

    // ---- 3. grupos: C_g ~ InvWishart(nl, U' K^-1 U)   (prior plano)
    if (!loc_fixa) {
      for (std::size_t g = 0; g < n_grupos; g++) {
        const Grupo& gr = d.modelo.grupos[g];
        const std::size_t dim = gr.dim;
        const std::size_t nl = nl_g[g];
        const std::size_t off = M.offset_grupo[g];
        Densa Su(dim, dim);
        std::vector<std::vector<double>> ku(dim);
        for (std::size_t b = 0; b < dim; b++) {
          ku[b].assign(nl, 0.0);
          const double* u = loc.data() + off + b * nl;
          if (d.kinv[g].ncol > 0) {
            const Csc& k = d.kinv[g];
            for (std::size_t c = 0; c < nl; c++)
              for (std::size_t p = k.colptr[c]; p < k.colptr[c + 1]; p++) {
                const std::size_t r = k.linha[p];
                ku[b][r] += k.valor[p] * u[c];
                if (r != c) ku[b][c] += k.valor[p] * u[r];
              }
          } else {
            for (std::size_t l = 0; l < nl; l++) ku[b][l] = u[l];
          }
        }
        for (std::size_t a = 0; a < dim; a++)
          for (std::size_t b = 0; b <= a; b++) {
            double s = 0.0;
            const double* ua = loc.data() + off + a * nl;
            for (std::size_t l = 0; l < nl; l++) s += ua[l] * ku[b][l];
            Su.at(a, b) = s;
            Su.at(b, a) = s;
          }
        // estabiliza contra Su singular no comeco da cadeia
        for (std::size_t a = 0; a < dim; a++) Su.at(a, a) += 1e-10;
        Densa cg = rinvwishart(static_cast<double>(nl), Su);
        for (std::size_t j = 0; j < dim; j++)
          for (std::size_t i = j; i < dim; i++)
            theta[gr.offset + vech_idx(i, j, dim)] = cg.at(i, j);
      }
    }

    // ---- 4. guarda a amostra
    if (it >= burnin && ((it - burnin) % thin == 0)) {
      n_amostras++;
      for (std::size_t k = 0; k < ntheta; k++) S.amostras.push_back(theta[k]);
      const double na = static_cast<double>(n_amostras);
      for (std::size_t k = 0; k < total_cols; k++) {
        const double delta = loc[k] - media_u[k];
        media_u[k] += delta / na;
        m2_u[k] += delta * (loc[k] - media_u[k]);
      }
    }
  }

  PutRNGstate();
  S.n_amostras = n_amostras;
  S.ntheta = ntheta;
  S.media_loc = media_u;
  S.var_loc.assign(total_cols, std::nan(""));
  if (n_amostras > 1)
    for (std::size_t k = 0; k < total_cols; k++)
      S.var_loc[k] = m2_u[k] / static_cast<double>(n_amostras - 1);
  S.ok = S.mensagem.empty();
  return S;
}

}  // namespace br
