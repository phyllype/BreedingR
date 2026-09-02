// AI-REML (Gilmour, Thompson e Cullis, 1995): score analitico, informacao media, EM
// (Dempster, Laird e Rubin, 1977) de aquecimento e passo amortecido.
//
// As regras que nao se quebram, todas vindas de defeito medido no projeto de origem:
//
//  1. convergencia RELATIVA nas componentes, sqrt(soma delta^2 / soma theta^2) < tol. Nunca
//     um limiar absoluto na norma do score, que e O(n_registros): a mesma tolerancia que
//     converge com 500 registros declara nao-convergido com 200 mil, com as estimativas
//     paradas no otimo. O certificado final soma o DECREMENTO DE NEWTON dos componentes
//     livres, g' AI^-1 g — normalizado pela curvatura, ele mede o gap em -2logL e nao a
//     escala do score, entao nao reintroduz o defeito de O(n_registros).
//  2. theta inadmissivel nao e resultado. C_g nao definida interrompe o passo; nunca vira
//     resposta.
//  3. a partida vem de var(y), nunca de valores de um cartao ja ajustado — partir das
//     estimativas finais de outro programa mede o quanto o estimador se afasta da resposta,
//     nao se ele a encontra.
//  4. um ajuste que nao convergiu DIZ isso.
//
// O score e a AI vem da INVERSA SELETIVA, nunca de uma P densa (n_registros^2, impagavel).
// As pecas, por grupo g com Var(u_g) = C_g (x) K:
//
//   Q[a,b] = u_a' K^-1 u_b                    (quadraticas da solucao)
//   T[a,b] = s2e . tr(K^-1 [Cs^-1]_{a,b})     (tracos da inversa seletiva)
//   score_g = (nl C_g^-1 - C_g^-1 (Q + T) C_g^-1) / ... na parametrizacao usual
//
// e o residual usa e'e do proprio ajuste. O EM de aquecimento e de graca: a atualizacao
// C_g <- (Q + T) / nl precisa exatamente do que o score ja computou.

#include "mme.h"

#include <cstdio>
#include <limits>

// so para o TRACE opcional e a interrupcao: nenhuma conta usa o R
#include <R_ext/Print.h>
#include <R_ext/Utils.h>

namespace br {

// ---------------------------------------------------------------- as pecas de uma avaliacao

// K^-1 u_b para um grupo: percorre a esparsa do triangulo inferior espelhando.
static std::vector<double> kinv_vezes(const Csc& k, const double* u, std::size_t nl) {
  std::vector<double> out(nl, 0.0);
  if (k.ncol == 0) {          // grupo diagonal: K = I
    for (std::size_t i = 0; i < nl; i++) out[i] = u[i];
    return out;
  }
  for (std::size_t c = 0; c < nl; c++)
    for (std::size_t p = k.colptr[c]; p < k.colptr[c + 1]; p++) {
      const std::size_t r = k.linha[p];
      out[r] += k.valor[p] * u[c];
      if (r != c) out[c] += k.valor[p] * u[r];
    }
  return out;
}

Avaliacao avalia(const Desenho& d, const std::vector<double>& theta,
                 CacheSimbolica* cache) {
  Avaliacao A;
  Montado M = monta_mme(d, theta);
  if (!M.ok) return A;

  CacheSimbolica local;
  CacheSimbolica* cs = cache ? cache : &local;
  if (!cs->pronto) {
    cs->perm = grau_minimo(M.c);
    Csc pc0 = permuta_sim(M.c, cs->perm);
    cs->sb = simbolica(pc0);
    cs->pronto = true;
  }
  const std::vector<std::size_t>& perm = cs->perm;
  const Simbolica& sb = cs->sb;
  Csc pc = permuta_sim(M.c, perm);
  Csc L;
  if (!cholesky(pc, sb, L)) return A;

  std::vector<double> pb(M.total);
  for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[perm[k]];
  std::vector<double> px = resolve(L, pb);
  A.solucao.assign(M.total, 0.0);
  for (std::size_t k = 0; k < M.total; k++) A.solucao[perm[k]] = px[k];

  // -2logL
  double yy = 0.0;
  std::size_t nreg = 0;
  for (std::size_t i = 0; i < d.nlin; i++)
    if (d.usa[i]) { yy += d.y[i] * d.y[i]; nreg++; }
  double bry = 0.0;
  for (std::size_t k = 0; k < M.total; k++) bry += A.solucao[k] * M.rhs[k];
  const double n = static_cast<double>(nreg), pfixo = static_cast<double>(M.n_fixo);
  A.neg2logl = (n - pfixo) * std::log(M.s2e) + M.logdet_g + logdet(L) + (yy - bry) / M.s2e;

  // inversa seletiva, lida na numeracao original atraves da permutacao
  SelInv z = inversa_seletiva(L, 0);
  std::vector<std::size_t> inv_perm(M.total);
  for (std::size_t k = 0; k < M.total; k++) inv_perm[perm[k]] = k;
  auto z_orig = [&](std::size_t i, std::size_t j, double& out) {
    return z.get(inv_perm[i], inv_perm[j], out);
  };

  const std::size_t ntheta = d.modelo.ntheta;
  A.score.assign(ntheta, 0.0);
  A.em_theta.assign(ntheta, 0.0);
  A.ai = Densa(ntheta, ntheta);

  // residuo e^ = y - W b, que o score residual e a AI precisam
  std::vector<double> ehat(d.nlin, 0.0);
  {
    // reconstruir W b por linha: X b_f + soma Z u
    for (std::size_t i = 0; i < d.nlin; i++) {
      if (!d.usa[i]) continue;
      double wb = 0.0;
      for (std::size_t j = 0; j < d.x.ncol; j++) wb += d.x.at(i, j) * A.solucao[j];
      ehat[i] = d.y[i] - wb;
    }
    for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
      std::size_t col0 = M.offset_grupo[g];
      for (std::size_t t : d.modelo.grupos[g].termos)
        for (const DesenhoTermo& a : d.aleatorios)
          if (a.termo == t) {
            for (std::size_t c = 0; c < a.z.ncol; c++) {
              const double bc = A.solucao[col0 + c];
              if (bc == 0.0) continue;
              for (std::size_t p = a.z.colptr[c]; p < a.z.colptr[c + 1]; p++)
                if (d.usa[a.z.linha[p]]) ehat[a.z.linha[p]] -= a.z.valor[p] * bc;
            }
            col0 += a.z.ncol;
          }
    }
  }
  double ete = 0.0;
  for (std::size_t i = 0; i < d.nlin; i++)
    if (d.usa[i]) ete += ehat[i] * ehat[i];

  // por grupo: Q, T, score e EM; e os vetores v_s = soma_d C^-1[s,d] u_d para a AI
  const double s2e = M.s2e;
  double trpen_total = 0.0;   // tr(C_s^-1 Pen), para o EM do residual
  std::size_t q_ordem = 0;    // colunas aleatorias totais

  std::vector<std::vector<std::vector<double>>> vgrupo(d.modelo.grupos.size());

  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    const Grupo& gr = d.modelo.grupos[g];
    const std::size_t dim = gr.dim;
    const std::size_t off = M.offset_grupo[g];
    std::size_t nl = 0;
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == gr.termos[0]) { nl = a.n_niveis; break; }
    q_ordem += dim * nl;

    Densa cg = cov_grupo(d.modelo, theta, g);
    Densa cinv = inv_geral(cg);
    Densa cgs = cg;
    for (double& v : cgs.dados) v /= s2e;
    Densa cinv_s = inv_geral(cgs);

    // Q[a,b] = u_a' K^-1 u_b
    Densa Q(dim, dim);
    std::vector<std::vector<double>> ku(dim);
    for (std::size_t b = 0; b < dim; b++)
      ku[b] = kinv_vezes(d.kinv[g], A.solucao.data() + off + b * nl, nl);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        const double* ua = A.solucao.data() + off + a * nl;
        for (std::size_t l = 0; l < nl; l++) s += ua[l] * ku[b][l];
        Q.at(a, b) = s;
      }

    // T[a,b] = s2e tr(K^-1 [Cs^-1]_{a,b}): le a inversa seletiva SO onde K^-1 tem nao-zero.
    // Todo nao-zero de K^-1 esta na penalidade, a penalidade esta em C, e todo nao-zero de C
    // esta no padrao do fator: um pedido fora do padrao significa que esse fechamento quebrou,
    // e e CONTADO em vez de ignorado — uma soma que pula termos em silencio nao e um traco.
    Densa T(dim, dim);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        if (d.kinv[g].ncol > 0) {
          const Csc& k = d.kinv[g];
          for (std::size_t c = 0; c < nl; c++)
            for (std::size_t p = k.colptr[c]; p < k.colptr[c + 1]; p++) {
              const std::size_t r = k.linha[p];
              double zv;
              if (z_orig(off + a * nl + r, off + b * nl + c, zv)) s += k.valor[p] * zv;
              else A.fora_do_padrao++;
              if (r != c) {
                if (z_orig(off + a * nl + c, off + b * nl + r, zv)) s += k.valor[p] * zv;
                else A.fora_do_padrao++;
              }
            }
        } else {
          for (std::size_t l = 0; l < nl; l++) {
            double zv;
            if (z_orig(off + a * nl + l, off + b * nl + l, zv)) s += zv;
            else A.fora_do_padrao++;
          }
        }
        T.at(a, b) = s2e * s;
      }

    // tr(C_s^-1 Pen) deste grupo, para o EM do residual:
    // tr(kron(Cs^-1,K^-1) Z) somado — igual a soma_ab Cs^-1[a,b] tr(K^-1 Z_ba) = soma T/s2e * Cs^-1
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++)
        trpen_total += cinv_s.at(a, b) * (T.at(b, a) / s2e);

    // score do grupo: dL/dC = (nl C^-1 - C^-1 (Q + T) C^-1), montado por elemento de theta
    Densa QT(dim, dim);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) QT.at(a, b) = Q.at(a, b) + T.at(a, b);
    // Mg = nl C^-1 - C^-1 QT C^-1
    Densa CiQT(dim, dim);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        for (std::size_t k = 0; k < dim; k++) s += cinv.at(a, k) * QT.at(k, b);
        CiQT.at(a, b) = s;
      }
    Densa Mg(dim, dim);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        for (std::size_t k = 0; k < dim; k++) s += CiQT.at(a, k) * cinv.at(k, b);
        Mg.at(a, b) = static_cast<double>(nl) * cinv.at(a, b) - s;
      }
    // derivada de -2logL em theta_k: para C simetrica, o elemento (i,j) com i!=j aparece
    // duas vezes em C, entao o score soma Mg[i,j] + Mg[j,i]
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const std::size_t k = gr.theta_idx(i, j);
        A.score[k] = (i == j) ? Mg.at(i, i) : Mg.at(i, j) + Mg.at(j, i);
      }

    // EM: C_g <- (Q + T) / nl
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++)
        A.em_theta[gr.theta_idx(i, j)] =
            (QT.at(i, j) + QT.at(j, i)) / (2.0 * static_cast<double>(nl));

    // v_s = soma_d C^-1[s,d] u_d, os vetores da AI
    vgrupo[g].assign(dim, std::vector<double>(nl, 0.0));
    for (std::size_t s = 0; s < dim; s++)
      for (std::size_t dd = 0; dd < dim; dd++) {
        const double c = cinv.at(s, dd);
        if (c == 0.0) continue;
        const double* ud = A.solucao.data() + off + dd * nl;
        for (std::size_t l = 0; l < nl; l++) vgrupo[g][s][l] += c * ud[l];
      }
  }

  // Score do residual, e o termo que faltou na primeira versao deste arquivo.
  //
  //   d(-2logL)/ds2e = tr(P) - y'PPy = (n - p - q + tr(Cs^-1 Pen))/s2e - e'e/s2e^2
  //
  // A primeira versao escreveu so (n - p)/s2e - e'e/s2e^2, sem o -q + trpen: o traco de P
  // nao e (n - p)/s2e porque os efeitos ALEATORIOS tambem absorvem graus de liberdade, na
  // medida em que o parentesco os deixa. O gate de diferencas finitas pegou na hora — e
  // com o score errado o passo amortecido ainda ACEITA (a verossimilhanca esta certa),
  // so que converge para o lugar errado. Um estimador com score errado e verossimilhanca
  // certa e o pior dos mundos: converge, reporta convergencia, e estima outra coisa.
  const double nreg_d = static_cast<double>(nreg);
  A.score[d.modelo.offset_residual] =
      (nreg_d - pfixo - static_cast<double>(q_ordem) + trpen_total) / s2e
      - ete / (s2e * s2e);

  // EM do residual: s2e <- (e'e + s2e (q - trpen)) / n
  A.em_theta[d.modelo.offset_residual] =
      (ete + s2e * (static_cast<double>(q_ordem) - trpen_total)) / nreg_d;

  // ---------------- AI: f_k = dV/dtheta_k . Py, e AI[i,j] = (f_i'f_j - (W'f_i)' Cs^-1 (W'f_j)) / s2e
  //
  // Py = e^/s2e. Para o grupo g e theta_k = C[i,j], f_k = Z_g dC (x) K Z_g' Py, que com a
  // solucao em mao e Z_g aplicado a combinacoes dos u. Montamos f_k por coluna de Z.
  const std::size_t ntot = M.total;
  std::vector<std::vector<double>> f(ntheta, std::vector<double>(d.nlin, 0.0));

  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    const Grupo& gr = d.modelo.grupos[g];
    const std::size_t dim = gr.dim;
    std::size_t nl = 0;
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == gr.termos[0]) { nl = a.n_niveis; break; }

    // dC/dtheta_k para (i,j): E_ij + E_ji (i!=j) ou E_ii. f_k = Z (dC (x) K) kron... com
    // Var = C (x) K e u = (C (x) K) Z' Py, tem-se (dC (x) K) Z'Py = (dC C^-1 (x) I) u.
    // Entao o coeficiente do slot s e: soma_d [dC C^-1]_{s,d} u_d = dC aplicado aos v.
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const std::size_t k = gr.theta_idx(i, j);
        // w_s = [dC v]_s: dC = E_ij + E_ji para i!=j, E_ii para i==j
        std::vector<std::vector<double>> w(dim, std::vector<double>(nl, 0.0));
        if (i == j) {
          w[i] = vgrupo[g][i];
        } else {
          w[i] = vgrupo[g][j];
          w[j] = vgrupo[g][i];
        }
        // f_k = soma_s Z_s w_s
        std::size_t col0 = M.offset_grupo[g];
        std::size_t slot = 0;
        for (std::size_t t : d.modelo.grupos[g].termos)
          for (const DesenhoTermo& a : d.aleatorios)
            if (a.termo == t) {
              for (std::size_t cc = 0; cc < a.n_coef; cc++) {
                const std::vector<double>& ws = w[slot];
                for (std::size_t nv = 0; nv < nl; nv++) {
                  const double x = ws[nv];
                  if (x == 0.0) continue;
                  const std::size_t col = cc * a.n_niveis + nv;
                  for (std::size_t p = a.z.colptr[col]; p < a.z.colptr[col + 1]; p++)
                    if (d.usa[a.z.linha[p]]) f[k][a.z.linha[p]] += a.z.valor[p] * x;
                }
                slot++;
              }
              col0 += a.z.ncol;
            }
      }
  }
  // residual: f = Py = e^/s2e
  for (std::size_t i = 0; i < d.nlin; i++)
    if (d.usa[i]) f[d.modelo.offset_residual][i] = ehat[i] / s2e;

  // AI[i,j] = (f_i' f_j - (W' f_i)' Cs^-1 (W' f_j)) / s2e, com W'f resolvido pelas MME
  std::vector<std::vector<double>> wtf(ntheta, std::vector<double>(ntot, 0.0));
  for (std::size_t k = 0; k < ntheta; k++) {
    // W'f: X'f e Z'f
    for (std::size_t j = 0; j < d.x.ncol; j++) {
      double s = 0.0;
      for (std::size_t i = 0; i < d.nlin; i++)
        if (d.usa[i]) s += d.x.at(i, j) * f[k][i];
      wtf[k][j] = s;
    }
    for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
      std::size_t col0 = M.offset_grupo[g];
      for (std::size_t t : d.modelo.grupos[g].termos)
        for (const DesenhoTermo& a : d.aleatorios)
          if (a.termo == t) {
            for (std::size_t c = 0; c < a.z.ncol; c++) {
              double s = 0.0;
              for (std::size_t p = a.z.colptr[c]; p < a.z.colptr[c + 1]; p++)
                if (d.usa[a.z.linha[p]]) s += a.z.valor[p] * f[k][a.z.linha[p]];
              wtf[k][col0 + c] = s;
            }
            col0 += a.z.ncol;
          }
    }
  }
  // Cs^-1 (W'f) via a fatoracao ja feita
  std::vector<std::vector<double>> cwf(ntheta);
  for (std::size_t k = 0; k < ntheta; k++) {
    std::vector<double> pb2(ntot);
    for (std::size_t q2 = 0; q2 < ntot; q2++) pb2[q2] = wtf[k][perm[q2]];
    std::vector<double> px2 = resolve(L, pb2);
    cwf[k].assign(ntot, 0.0);
    for (std::size_t q2 = 0; q2 < ntot; q2++) cwf[k][perm[q2]] = px2[q2];
  }
  for (std::size_t a = 0; a < ntheta; a++)
    for (std::size_t b = a; b < ntheta; b++) {
      double ff = 0.0;
      for (std::size_t i = 0; i < d.nlin; i++)
        if (d.usa[i]) ff += f[a][i] * f[b][i];
      double wcw = 0.0;
      for (std::size_t q2 = 0; q2 < ntot; q2++) wcw += wtf[a][q2] * cwf[b][q2];
      const double v = (ff - wcw) / s2e;
      A.ai.at(a, b) = v;
      A.ai.at(b, a) = v;
    }

  A.ok = true;
  return A;
}

// ------------------- log-Cholesky coordinates (Pinheiro and Bates, 1996) for the STEP
//
// The damped AI step used to walk in theta itself, and a measured failure showed why that
// cannot work: with two terms in one group and a strong direct-social correlation, the REML
// optimum sits ON the boundary det(C_g) = 0. In theta coordinates that boundary is a curved
// wall; the admissible trench around the ridge gets thinner than ANY additive step, so every
// candidate lands outside the cone, the damping ratchets up by a factor of 10 per few
// iterations, and the fit crawls on single EM steps until it declares a false convergence
// 5-13 -2logL units above the optimum (seeds 11 and 51 of the gate in test-aireml-boundary.R).
//
// The step therefore walks in z: per group the lower Cholesky factor of C_g with the log of
// its diagonal, plus log(s2e). In z EVERY point maps to an admissible theta, the singular
// boundary is pushed to -infinity, and the ridge uncurls into a direction the damped step can
// follow. The score, the AI matrix and the likelihood are untouched: the chain rule through
// J = dtheta/dz converts them, and the answer is still reported in theta.

// Cholesky of a small dense C_g (group dimension is 1-3 in practice); false when not PD.
// Shared through mme.h: the multi-trait and AR(1) certificates use it to recognize a
// covariance group at the singularity boundary.
bool chol_pequena(const Densa& c, Densa& L) {
  const std::size_t n = c.nlin;
  L = Densa(n, n);
  for (std::size_t j = 0; j < n; j++) {
    double s = c.at(j, j);
    for (std::size_t k = 0; k < j; k++) s -= L.at(j, k) * L.at(j, k);
    if (!(s > 0.0) || !std::isfinite(s)) return false;
    L.at(j, j) = std::sqrt(s);
    for (std::size_t i = j + 1; i < n; i++) {
      double t = c.at(i, j);
      for (std::size_t k = 0; k < j; k++) t -= L.at(i, k) * L.at(j, k);
      L.at(i, j) = t / L.at(j, j);
    }
  }
  return true;
}

// z shares theta's indexing: z[theta_idx(i,j)] = log L_ii on the diagonal, L_ij below it,
// and z[offset_residual] = log s2e. False when some C_g is not PD (nothing to map).
static bool z_de_theta(const Modelo& mo, const std::vector<double>& theta,
                       std::vector<double>& z) {
  z.assign(mo.ntheta, 0.0);
  for (const Grupo& gr : mo.grupos) {
    Densa cg(gr.dim, gr.dim);
    for (std::size_t j = 0; j < gr.dim; j++)
      for (std::size_t i = j; i < gr.dim; i++) {
        const double v = theta[gr.theta_idx(i, j)];
        cg.at(i, j) = v;
        cg.at(j, i) = v;
      }
    Densa L;
    if (!chol_pequena(cg, L)) return false;
    for (std::size_t j = 0; j < gr.dim; j++)
      for (std::size_t i = j; i < gr.dim; i++)
        z[gr.theta_idx(i, j)] = (i == j) ? std::log(L.at(i, i)) : L.at(i, j);
  }
  if (!(theta[mo.offset_residual] > 0.0)) return false;
  z[mo.offset_residual] = std::log(theta[mo.offset_residual]);
  return true;
}

static std::vector<double> theta_de_z(const Modelo& mo, const std::vector<double>& z) {
  std::vector<double> theta(mo.ntheta, 0.0);
  for (const Grupo& gr : mo.grupos) {
    Densa L(gr.dim, gr.dim);
    for (std::size_t j = 0; j < gr.dim; j++)
      for (std::size_t i = j; i < gr.dim; i++) {
        const double v = z[gr.theta_idx(i, j)];
        L.at(i, j) = (i == j) ? std::exp(v) : v;
      }
    for (std::size_t j = 0; j < gr.dim; j++)
      for (std::size_t i = j; i < gr.dim; i++) {
        double s = 0.0;
        for (std::size_t k = 0; k <= j; k++) s += L.at(i, k) * L.at(j, k);
        theta[gr.theta_idx(i, j)] = s;
      }
  }
  theta[mo.offset_residual] = std::exp(z[mo.offset_residual]);
  return theta;
}

// J = dtheta/dz, block diagonal by group. For z_m = L entry (a,b) with chain factor
// s_m (= L_aa on the diagonal, 1 below), and theta_k = C_ij (i >= j):
//   dC_ij/dz_m = s_m (delta_ia L_jb + delta_ja L_ib)
static Densa jacobiano_z(const Modelo& mo, const std::vector<double>& z) {
  Densa J(mo.ntheta, mo.ntheta);
  for (const Grupo& gr : mo.grupos) {
    Densa L(gr.dim, gr.dim);
    for (std::size_t j = 0; j < gr.dim; j++)
      for (std::size_t i = j; i < gr.dim; i++) {
        const double v = z[gr.theta_idx(i, j)];
        L.at(i, j) = (i == j) ? std::exp(v) : v;
      }
    for (std::size_t b = 0; b < gr.dim; b++)
      for (std::size_t a = b; a < gr.dim; a++) {
        const std::size_t m = gr.theta_idx(a, b);
        const double sm = (a == b) ? L.at(a, a) : 1.0;
        for (std::size_t j = 0; j < gr.dim; j++)
          for (std::size_t i = j; i < gr.dim; i++) {
            double v = 0.0;
            if (i == a && b <= j) v += L.at(j, b);
            if (j == a && b <= i) v += L.at(i, b);
            J.at(gr.theta_idx(i, j), m) = sm * v;
          }
      }
  }
  J.at(mo.offset_residual, mo.offset_residual) = std::exp(z[mo.offset_residual]);
  return J;
}

// ------------------------------------------------------------------------------- o ajuste

std::vector<double> partida(const Desenho& d) {
  double soma = 0.0, soma2 = 0.0;
  std::size_t n = 0;
  for (std::size_t i = 0; i < d.nlin; i++)
    if (d.usa[i]) { soma += d.y[i]; soma2 += d.y[i] * d.y[i]; n++; }
  const double media = soma / static_cast<double>(n);
  double var = (soma2 - static_cast<double>(n) * media * media) / std::max(1.0, static_cast<double>(n - 1));
  if (!(var > 0.0) || !std::isfinite(var)) var = 1.0;

  std::vector<double> theta(d.modelo.ntheta, 0.0);
  const double frac = 0.5 / static_cast<double>(d.modelo.grupos.size());
  for (const Grupo& g : d.modelo.grupos)
    for (std::size_t i = 0; i < g.dim; i++) theta[g.theta_idx(i, i)] = frac * var;
  theta[d.modelo.offset_residual] = 0.5 * var;
  return theta;
}

Ajuste ajusta(const Desenho& d, const std::vector<double>* theta0,
              std::size_t n_em, std::size_t maxiter, double tol, bool verboso) {
  Ajuste R;
  std::vector<double> theta = theta0 ? *theta0 : partida(d);
  if (theta.size() != d.modelo.ntheta) {
    R.mensagem = "theta inicial com tamanho errado";
    return R;
  }

  CacheSimbolica cs;
  Avaliacao cur = avalia(d, theta, &cs);
  if (verboso)
    Rprintf("AI-REML: %d record(s), %d column(s), %d component(s)\n",
            (int) d.n_usadas(), (int) d.total_colunas(), (int) d.modelo.ntheta);
  if (!cur.ok) {
    R.mensagem = "o theta inicial e INADMISSIVEL: alguma covariancia nao e positiva-definida";
    return R;
  }

  // EM de aquecimento: monotono, fica no cone, e de graca. A checagem extra de
  // z_de_theta garante que o ponto de onde o passo amortecido parte tem fator de
  // Cholesky — um EM que encoste numericamente na singularidade nao vira partida.
  for (std::size_t e = 0; e < n_em; e++) {
    std::vector<double> zt;
    if (!z_de_theta(d.modelo, cur.em_theta, zt)) break;
    Avaliacao prox = avalia(d, cur.em_theta, &cs);
    if (!prox.ok) break;
    theta = cur.em_theta;
    cur = std::move(prox);
  }

  if (verboso && n_em > 0)
    Rprintf("EM warm-up (%d round(s)): -2logL %.6f\n", (int) n_em, cur.neg2logl);
  // CONJUNTO ATIVO, agora em z. Duas fronteiras REAIS do espaco de parametros passam
  // pelo mesmo mecanismo: uma variancia com piso em zero e um C_g de grupo encostando na
  // singularidade (correlacao em +/-1, o regime medido no modelo direto-social). Nos dois
  // casos a diagonal do fator de Cholesky tem um piso — absoluto para a variancia, RELATIVO
  // a maior diagonal do grupo para a singularidade — e a entrada presa no piso com o score
  // empurrando para fora e congelada, com o resto otimizado CONDICIONALMENTE a isso, que e
  // a solucao correta do problema com restricao.
  double lambda = 1e-3;
  std::vector<double> z;
  if (!z_de_theta(d.modelo, theta, z)) {
    R.mensagem = "o theta inicial e INADMISSIVEL: alguma covariancia nao e positiva-definida";
    return R;
  }

  // floors on the diagonal of z, recomputed at the current point: relative to the group's
  // largest Cholesky diagonal (so a one-term group never binds on it) and absolute against
  // the residual scale (the old variance floor, 1e-10 max(1, s2e), read through l = sqrt v)
  auto pisos_z = [&](const std::vector<double>& zc, std::vector<double>& piso,
                     std::vector<char>& relativo) {
    piso.assign(zc.size(), -std::numeric_limits<double>::infinity());
    relativo.assign(zc.size(), 0);
    const double s2e = std::exp(zc[d.modelo.offset_residual]);
    const double piso_abs = 0.5 * std::log(1e-10 * std::max(1.0, s2e));
    for (const Grupo& gr : d.modelo.grupos) {
      double zmax = -std::numeric_limits<double>::infinity();
      for (std::size_t i = 0; i < gr.dim; i++)
        zmax = std::max(zmax, zc[gr.theta_idx(i, i)]);
      const double piso_rel = zmax + std::log(1e-3);
      for (std::size_t i = 0; i < gr.dim; i++) {
        const std::size_t k = gr.theta_idx(i, i);
        piso[k] = std::max(piso_rel, piso_abs);
        relativo[k] = piso_rel > piso_abs ? 1 : 0;
      }
    }
  };

  // The pieces of the step at a point, shared by the walker and by the final certificate:
  // score and AI mapped to z by the chain rule, the floors, and two active sets over the
  // SAME floors and the SAME outward-score rule, differing only in reach. `congelado` is
  // the step's (strictly at the clamp), `na_parede` is the certificate's and the boundary
  // report's (within a factor of 2 of the floor): a component pinned at a floor with the
  // score pushing outward never zeroes its gradient — it points out of the cone by
  // construction — so the certificate must exclude it, whether the walker left it clamped
  // or resting a hair above the clamp. Excluding anything more would certify a point with
  // a live direction ignored.
  struct PecasZ {
    std::vector<double> sz, piso;
    Densa az;
    std::vector<char> congelado, na_parede, piso_rel;
  };
  auto pecas_z = [&](const std::vector<double>& zc, const Avaliacao& av) {
    PecasZ P;
    pisos_z(zc, P.piso, P.piso_rel);
    Densa J = jacobiano_z(d.modelo, zc);
    const std::size_t nt = d.modelo.ntheta;
    P.sz.assign(nt, 0.0);
    for (std::size_t m = 0; m < nt; m++)
      for (std::size_t k = 0; k < nt; k++) P.sz[m] += J.at(k, m) * av.score[k];
    P.az = Densa(nt, nt);
    for (std::size_t a = 0; a < nt; a++)
      for (std::size_t b = 0; b < nt; b++) {
        double s = 0.0;
        for (std::size_t k = 0; k < nt; k++) {
          if (J.at(k, a) == 0.0) continue;
          double t = 0.0;
          for (std::size_t l = 0; l < nt; l++) t += av.ai.at(k, l) * J.at(l, b);
          s += J.at(k, a) * t;
        }
        P.az.at(a, b) = s;
      }
    P.congelado.assign(nt, 0);
    P.na_parede.assign(nt, 0);
    for (const Grupo& gr : d.modelo.grupos)
      for (std::size_t i = 0; i < gr.dim; i++) {
        const std::size_t k = gr.theta_idx(i, i);
        P.congelado[k] = (zc[k] <= P.piso[k] + 1e-9 && P.sz[k] > 0.0) ? 1 : 0;
        // the certificate's active set, and it is WIDER than the step's on purpose: the
        // walker rests a hair above the clamp without ever being frozen there (the same
        // measured fact the boundary report handles with the same factor-2 tolerance).
        // A hair above the floor the z-gradient of the pinned diagonal is ~L^2 and its
        // AI_z diagonal ~L^4, so the quadratic form returns that direction's FULL
        // theta-scale share — measured 1.1499 on a true boundary optimum (seed 11 of
        // test-aireml-boundary.R) where the free directions were flat. On the floor
        // with the score pointing outward is pinned, whether clamped or resting.
        P.na_parede[k] = (zc[k] <= P.piso[k] + std::log(2.0) && P.sz[k] > 0.0) ? 1 : 0;
      }
    return P;
  };

  // THE CONVERGENCE CERTIFICATE, measured missing on real data (docs/CHECKLIST.md,
  // 2026-09-02): a direct+indirect fit ended converged TRUE with relDelta 5.1e-09 and
  // group scores up to -506050, sitting 6.9 -2logL units ABOVE the optimum. The damped
  // AI step and the EM step can stall together against the singularity wall with the
  // gradient far from zero, and neither the small step nor the EM confirmation sees it.
  // The evidence of an optimum is the Newton decrement restricted to the FREE components,
  //
  //   dec = g_free' [AI_free]^-1 g_free
  //
  // computed in the z coordinates of the step (the decrement is invariant to the
  // parametrization; the active set is defined in z). Near the optimum, -2logL exceeds
  // its minimum by ~dec/2 (second-order Taylor with the AI for the Hessian), so dec is a
  // GAP ON THE -2logL SCALE: the tolerance 2e-4 certifies the fit sits within ~1e-4 of
  // its optimum — four orders of magnitude below the ~4 units a likelihood-ratio test
  // calls a difference — while the measured defect (gap 6.9, dec ~14) fails it by five.
  // An AI too singular to solve on the free block certifies nothing and counts as a
  // refusal.
  const double tol_dec = 2e-4;
  auto decremento_livre = [&](const PecasZ& P) -> double {
    Densa m = P.az;
    std::vector<double> sc = P.sz;
    bool algum_livre = false;
    for (std::size_t i = 0; i < m.nlin; i++) {
      if (P.na_parede[i]) {
        sc[i] = 0.0;
        for (std::size_t j = 0; j < m.ncol; j++) { m.at(i, j) = 0.0; m.at(j, i) = 0.0; }
        m.at(i, i) = 1.0;
      } else {
        algum_livre = true;
      }
    }
    if (!algum_livre) return 0.0;
    Densa minv;
    try { minv = inv_geral(m); } catch (const Erro&) {
      return std::numeric_limits<double>::infinity();
    }
    double dec = 0.0;
    for (std::size_t i = 0; i < sc.size(); i++) {
      if (sc[i] == 0.0) continue;
      for (std::size_t j = 0; j < sc.size(); j++) dec += sc[i] * minv.at(i, j) * sc[j];
    }
    return std::fabs(dec);   // the AI is PSD; a negative here is rounding, not curvature
  };

  // EM rescue: multiplicative, stays in the cone, and does not shrink at the boundary, so
  // it walks exactly where the damped step stalls. Each EM candidate is CLAMPED through z
  // (an EM that lands numerically ON the singularity would evaluate to garbage: the MME
  // -2logL loses its meaning when det(C_g) underflows). Up to 50 steps, while each pays
  // more than 1e-8.
  auto resgate_em = [&]() -> bool {
    bool ganhou = false;
    for (int e = 0; e < 50; e++) {
      std::vector<double> ze;
      if (!z_de_theta(d.modelo, cur.em_theta, ze)) break;
      std::vector<double> piso; std::vector<char> rel;
      pisos_z(ze, piso, rel);
      for (const Grupo& gr : d.modelo.grupos)
        for (std::size_t i = 0; i < gr.dim; i++) {
          const std::size_t k = gr.theta_idx(i, i);
          ze[k] = std::max(ze[k], piso[k]);
        }
      std::vector<double> cand = theta_de_z(d.modelo, ze);
      Avaliacao prox = avalia(d, cand, &cs);
      if (!prox.ok || prox.neg2logl >= cur.neg2logl - 1e-8) break;
      double num = 0.0, den = 0.0;
      for (std::size_t i = 0; i < cand.size(); i++) {
        const double dlt = cand[i] - theta[i];
        num += dlt * dlt;
        den += cand[i] * cand[i];
      }
      R.reldelta = std::sqrt(num / std::max(den, 1e-300));
      theta = cand;
      z = ze;
      cur = std::move(prox);
      ganhou = true;
    }
    return ganhou;
  };

  for (std::size_t it = 1; it <= maxiter; it++) {
    R.iters = it;
    R_CheckUserInterrupt();
    // chain rule (score_z = J' score, AI_z = J' AI J), floors and active set, all at
    // the current point — the same pieces the certificate reads at the exits
    PecasZ P = pecas_z(z, cur);
    const std::vector<double>& piso = P.piso;
    const std::vector<double>& sz = P.sz;
    const Densa& az = P.az;
    const std::vector<char>& congelado = P.congelado;
    const std::size_t nt = theta.size();
    bool aceitou = false;
    std::vector<double> zc;
    for (int tent = 0; tent < 30 && !aceitou; tent++) {
      Densa m = az;
      for (std::size_t i = 0; i < m.nlin; i++) {
        const double di = m.at(i, i);
        m.at(i, i) = (di == 0.0) ? lambda : di * (1.0 + lambda);
      }
      // congelado: linha e coluna viram identidade e o score sai — e o sistema REDUZIDO
      // dos que ainda podem andar, sem que o bloqueado contamine a direcao dos outros
      std::vector<double> sc = sz;
      for (std::size_t i = 0; i < m.nlin; i++)
        if (congelado[i]) {
          sc[i] = 0.0;
          for (std::size_t j = 0; j < m.ncol; j++) { m.at(i, j) = 0.0; m.at(j, i) = 0.0; }
          m.at(i, i) = 1.0;
        }
      Densa minv;
      try { minv = inv_geral(m); } catch (const Erro&) { lambda *= 10.0; continue; }
      std::vector<double> passo(nt, 0.0);
      for (std::size_t i = 0; i < nt; i++) {
        if (congelado[i]) continue;
        for (std::size_t j = 0; j < nt; j++) passo[i] += minv.at(i, j) * sc[j];
      }
      // step halving on the SAME direction: the damped direction can be right while the
      // full length overshoots the curved valley, and rejecting it only to re-damp is what
      // used to ratchet lambda
      for (const double alpha : {1.0, 0.5, 0.25}) {
        zc = z;
        for (std::size_t i = 0; i < nt; i++)
          if (!congelado[i]) zc[i] -= alpha * passo[i];
        for (const Grupo& gr : d.modelo.grupos)
          for (std::size_t i = 0; i < gr.dim; i++) {
            const std::size_t k = gr.theta_idx(i, i);
            zc[k] = std::max(zc[k], piso[k]);
          }
        std::vector<double> cand = theta_de_z(d.modelo, zc);
        Avaliacao prox = avalia(d, cand, &cs);
        if (prox.ok && prox.neg2logl <= cur.neg2logl + 1e-9) {
          double num = 0.0, den = 0.0;
          for (std::size_t i = 0; i < cand.size(); i++) {
            const double dlt = cand[i] - theta[i];
            num += dlt * dlt;
            den += cand[i] * cand[i];
          }
          R.reldelta = std::sqrt(num / std::max(den, 1e-300));
          theta = cand;
          z = zc;
          cur = std::move(prox);
          lambda = std::max(lambda / 10.0, 1e-10);
          aceitou = true;
          break;
        }
      }
      if (!aceitou) lambda *= 10.0;
    }
    if (!aceitou) {
      // before giving up, hand the point to the EM rescue: measured on the direct-social
      // gate, the damped step can be stuck while EM still walks — and after a rescue the
      // damping RESTARTS, or the next AI step would inherit a lambda pumped to 1e20+ by
      // the failed attempts and stay parked forever ("amortecimento preso")
      if (resgate_em()) {
        lambda = 1e-3;
        if (verboso)
          Rprintf("      EM rescue: -2logL %.6f (no damped step improved)\n", cur.neg2logl);
        continue;
      }
      // In z every candidate is admissible, so this exit means the likelihood itself
      // refused 90 second-order variants AND the EM step. That is where the walker
      // rests — but resting is not an optimum: on the measured direct+indirect regime
      // the AI and the EM stall TOGETHER against the wall with the gradient far from
      // zero (the certificate defect above). Only the Newton decrement tells a numerical
      // optimum from a stall; the code before this stage declared converged TRUE here
      // unconditionally.
      R.decremento = decremento_livre(P);
      if (R.decremento < tol_dec) {
        R.convergiu = true;
        R.mensagem = "stopped where neither the damped AI step nor the EM step improves "
            "the likelihood (a numerical optimum tighter than the step tolerance)";
      } else {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.3g", R.decremento);
        R.mensagem = std::string("did NOT converge: neither the damped AI step nor the ") +
            "EM step improves the likelihood, but the Newton decrement of the free "
            "components is " + buf + " (tolerance 2e-4 on the -2logL scale), so the "
            "point is a stall, not a certified optimum. Try a different start=, and "
            "profile any component the message reports at a boundary";
      }
      break;
    }
    if (verboso)
      Rprintf("iter %3d  -2logL %.6f  relDelta %.3e\n",
              (int) it, cur.neg2logl, R.reldelta);
    if (R.reldelta < tol) {
      // O PASSO PEQUENO NAO PROVA OTIMO. Com uma componente encostada numa fronteira a AI
      // fica quase singular naquela direcao, o amortecimento cresce, o passo encolhe — e
      // um criterio que olha so o tamanho do passo declararia convergencia com o score
      // longe de zero. Antes de declarar, o EM tem de confirmar que nao anda mais; se
      // ainda melhora, nao havia convergencia nenhuma, e o amortecimento reinicia.
      if (resgate_em()) {
        lambda = 1e-3;
        if (verboso)
          Rprintf("      EM rescue: -2logL %.6f (the AI step had stalled)\n", cur.neg2logl);
        continue;
      }
      // ... e um EM parado tambem nao prova nada quando AI e EM travam JUNTOS (o defeito
      // do certificado, medido): converged exige ainda o decremento de Newton dos
      // componentes livres. A recusa NAO encerra o ajuste — o passo pode voltar a andar —
      // e se maxiter chegar, a mensagem final carrega o decremento reprovado.
      R.decremento = decremento_livre(pecas_z(z, cur));
      if (R.decremento < tol_dec) {
        R.convergiu = true;
        break;
      }
    }
  }
  // the decrement at the final point, whatever the exit: the certificate the user can
  // read next to fit$score, and the number the maxiter message reports
  if (std::isnan(R.decremento)) R.decremento = decremento_livre(pecas_z(z, cur));
  // norma do score no ponto final, para quem quiser conferir que ele de fato zerou:
  // e a evidencia direta de otimo, que o tamanho do passo nao da.
  R.score = cur.score;
  if (verboso)
    Rprintf("%s at iter %d, -2logL %.6f\n",
            R.convergiu ? "converged" : "STOPPED", (int) R.iters, cur.neg2logl);

  // erros padrao: 2 [AI^-1]_kk. A matriz inteira tambem sai: 2 AI^-1 e a covariancia
  // amostral das componentes, e sem ela nao ha metodo delta para h2, r ou T.
  R.se.assign(d.modelo.ntheta, std::nan(""));
  R.vcov.assign(d.modelo.ntheta * d.modelo.ntheta, std::nan(""));
  try {
    Densa inv = inv_geral(cur.ai);
    for (std::size_t k = 0; k < d.modelo.ntheta; k++) {
      const double v = 2.0 * inv.at(k, k);
      if (v > 0.0) R.se[k] = std::sqrt(v);
      for (std::size_t j = 0; j < d.modelo.ntheta; j++)
        R.vcov[k * d.modelo.ntheta + j] = 2.0 * inv.at(k, j);
    }
  } catch (const Erro&) {}

  if (!R.convergiu && R.mensagem.empty()) {
    // The maxiter exit, and it comes BEFORE the boundary notes so those append to the
    // explanation instead of silencing it. Measured by the user on the 2x2 group
    // warm-started from the reduced model: 100 iterations ended with relDelta still
    // 1.6e-04 — not a defect, a model that walks slowly. The message says exactly
    // that, and how to ask for more.
    char buf[128];
    std::snprintf(buf, sizeof(buf), " (relDelta %.3g, Newton decrement of the free "
                  "components %.3g against the 2e-4 tolerance)", R.reldelta, R.decremento);
    R.mensagem = "stopped at " + std::to_string(maxiter) + " iteration(s) without a "
        "certified optimum" + buf + ": this model asks for more iterations. Raise "
        "maxiter=; a start= from a reduced model shortens the walk";
  }
  {
    // The boundary report reads the FINAL point, not the freeze flags of the last
    // iteration: near a boundary optimum the walker can rest a hair above the floor
    // without ever being frozen there, and the user still needs to know the estimate
    // is on the edge (fit$score is NOT ~0 in the boundary direction, by design).
    std::size_t nzero = 0, nsing = 0;
    std::vector<double> zf;
    if (z_de_theta(d.modelo, theta, zf)) {
      std::vector<double> piso; std::vector<char> rel;
      pisos_z(zf, piso, rel);
      for (const Grupo& gr : d.modelo.grupos)
        for (std::size_t i = 0; i < gr.dim; i++) {
          const std::size_t k = gr.theta_idx(i, i);
          // within a factor of 2 of the floor counts as ON it: the walker rests where
          // the trench flattens, a hair above the clamp, and the reader still needs
          // to know the estimate is an edge case
          if (zf[k] <= piso[k] + std::log(2.0)) (rel[k] ? nsing : nzero)++;
        }
    }
    if (nzero > 0)
      R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") + std::to_string(nzero) +
          " component(s) at the zero boundary, fixed there while the others were "
          "optimized conditional on that";
    if (nsing > 0)
      R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") + std::to_string(nsing) +
          " covariance group direction(s) at the singularity boundary (a correlation at "
          "+/-1), held there while the others were optimized conditional on that; the "
          "boundary is real, but a delta-method interval there is not — profile the "
          "component instead";
  }
  if (cur.fora_do_padrao > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(cur.fora_do_padrao) +
        " leitura(s) fora do padrao do fator: o resultado NAO e de confianca";

  R.theta = theta;
  R.neg2logl = cur.neg2logl;
  R.solucao = cur.solucao;
  R.fora_do_padrao = cur.fora_do_padrao;

  // PEV: a diagonal da inversa seletiva das MME no otimo, vezes s2e.
  //
  // As MME estao em unidades de s2e (a penalidade e kron(Cs^-1, K^-1) com Cs = C/s2e),
  // entao [C_s^-1]_ii ja e PEV_i / s2e — a multiplicacao devolve a escala absoluta. Sem
  // esse fator a acuracia sai sistematicamente errada e ainda parece plausivel.
  {
    Montado M = monta_mme(d, theta);
    if (M.ok) {
      if (!cs.pronto) {
        cs.perm = grau_minimo(M.c);
        Csc pc0 = permuta_sim(M.c, cs.perm);
        cs.sb = simbolica(pc0);
        cs.pronto = true;
      }
      const std::vector<std::size_t>& perm = cs.perm;
      Csc pc = permuta_sim(M.c, perm);
      Csc L;
      if (cholesky(pc, cs.sb, L)) {
        SelInv z = inversa_seletiva(L, 0);
        std::vector<double> diag = z.diagonal();
        R.pev.assign(M.total, std::nan(""));
        for (std::size_t k = 0; k < M.total; k++)
          R.pev[perm[k]] = diag[k] * M.s2e;
      }
    }
  }
  return R;
}

}  // namespace br
