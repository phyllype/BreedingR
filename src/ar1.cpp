// Residuo AR(1)/CAR(1), agora para t caracteristicas: correlacao rho^|dt| dentro de
// sujeito NO TEMPO, e R0 cheia ENTRE caracteristicas, separaveis:
//
//   R = diag_s { Gamma_s (x) R0 },   Gamma_s[i,j] = rho^|t_i - t_j|
//
// t = 1 colapsa exatamente no caminho antigo (R0 1x1 = s2e), e os gates antigos viram a
// regressao disto. A separabilidade e o que paga a conta:
//
//   R^-1 = Gamma^-1 (x) R0^-1     (tridiagonal (x) pequena)
//   log|R| = soma_s [ t log|Gamma_s| + m_s log|R0| ]
//
// e TODA derivada residual tem a forma dR = D (x) S com D no espaco dos registros e S no
// das caracteristicas:
//
//   vech(R0)_ab :  D = Gamma_s      S = E_ab (simetrica)
//   rho         :  D = dGamma/drho  S = R0
//
// o que da um score unico para os dois tipos:
//
//   score = soma_s tr(Gamma^-1 D) tr(R0^-1 S)  -  tr(C^-1 W' M2 W)  -  quad
//   M2    = (Gamma^-1 D Gamma^-1) (x) (R0^-1 S R0^-1)
//   quad  = soma_ij D[i,j] . (R^-1 e)_i' S (R^-1 e)_j
//
// O clique COMPLETO do sujeito (zeros explicitos) continua no padrao de C, porque o traco
// de rho le a inversa seletiva em pares de registros distantes; fora_do_padrao denuncia.
// Empate de tempo no mesmo sujeito continua erro declarado.

#include "mme.h"
#include <unordered_map>
#include <algorithm>

namespace br {

static std::size_t vech_idx(std::size_t i, std::size_t j, std::size_t d) {
  if (i < j) std::swap(i, j);
  return (j * (2 * d - j + 1)) / 2 + (i - j);
}

// Gamma_s^-1 tridiagonal e log|Gamma_s|, pela forma fechada de Markov.
struct GammaInv {
  std::vector<double> diag, sub;   // sub[k] liga k e k+1
  double logdet = 0.0;
};

static GammaInv gamma_inv(const std::vector<double>& t, double rho) {
  const std::size_t m = t.size();
  GammaInv g;
  g.diag.assign(m, 1.0);
  g.sub.assign(m > 0 ? m - 1 : 0, 0.0);
  for (std::size_t k = 0; k + 1 < m; k++) {
    const double r = std::pow(std::fabs(rho), t[k + 1] - t[k]) *
                     ((rho < 0.0 && std::fmod(t[k + 1] - t[k], 2.0) != 0.0) ? -1.0 : 1.0);
    const double um = 1.0 - r * r;
    g.diag[k]     += r * r / um;
    g.diag[k + 1] += r * r / um;
    g.sub[k]       = -r / um;
    g.logdet      += std::log(um);
  }
  return g;
}

// Gamma_s densa (para dGamma/drho e para a forma V de referencia).
static Densa gamma_densa(const std::vector<double>& t, double rho) {
  const std::size_t m = t.size();
  Densa g(m, m);
  for (std::size_t i = 0; i < m; i++)
    for (std::size_t j = 0; j < m; j++) {
      const double dt = std::fabs(t[i] - t[j]);
      double v = std::pow(std::fabs(rho), dt);
      if (rho < 0.0 && std::fmod(dt, 2.0) != 0.0) v = -v;
      g.at(i, j) = (i == j) ? 1.0 : v;
    }
  return g;
}

// dGamma/drho, densa por sujeito. d(rho^dt)/drho = dt * rho^(dt-1).
static Densa dgamma_densa(const std::vector<double>& t, double rho) {
  const std::size_t m = t.size();
  Densa g(m, m);
  for (std::size_t i = 0; i < m; i++)
    for (std::size_t j = 0; j < m; j++) {
      if (i == j) continue;
      const double dt = std::fabs(t[i] - t[j]);
      double v;
      if (dt == 1.0) v = 1.0;
      else if (rho == 0.0) v = 0.0;                      // dt > 1 em rho = 0
      else {
        v = dt * std::pow(std::fabs(rho), dt - 1.0);
        if (rho < 0.0 && std::fmod(dt - 1.0, 2.0) != 0.0) v = -v;
      }
      g.at(i, j) = v;
    }
  return g;
}

// R0 lida do fim de theta (t = 1: o proprio s2e).
static Densa r0_de(const DesenhoAR& d, const std::vector<double>& theta) {
  Densa r(d.t, d.t);
  for (std::size_t j = 0; j < d.t; j++)
    for (std::size_t i = j; i < d.t; i++) {
      const double v = theta[d.offset_s2e + vech_idx(i, j, d.t)];
      r.at(i, j) = v;
      r.at(j, i) = v;
    }
  return r;
}

struct MontadoAR {
  Csc c;
  std::vector<double> rhs;
  double logdet_g = 0.0;
  double logdet_r = 0.0;            // soma_s t log|Gamma_s| + n log|R0|
  double yry = 0.0;
  Densa r0, r0inv;
  double rho = 0.0;
  std::size_t n_fixo = 0, total = 0;
  std::vector<std::size_t> offset_grupo;
  bool ok = false;
};

MontadoAR monta_mme_ar1(const DesenhoAR& d, const std::vector<double>& theta) {
  MontadoAR M;
  const std::size_t t = d.t;
  M.rho = theta[d.offset_rho];
  if (!(std::fabs(M.rho) < 1.0)) return M;
  M.r0 = r0_de(d, theta);
  {
    Densa l = M.r0;
    if (!chol_densa(l)) return M;   // R0 inadmissivel
  }
  M.r0inv = inv_pd(M.r0);
  const double ld_r0 = logdet_pd(M.r0);

  M.n_fixo = d.x.ncol * t;
  M.offset_grupo.resize(d.modelo.grupos.size());
  std::size_t acc = M.n_fixo;
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    M.offset_grupo[g] = acc;
    acc += d.largura(g);
  }
  M.total = acc;

  std::vector<std::uint32_t> ti, tj;
  std::vector<double> tv;
  M.rhs.assign(M.total, 0.0);
  M.logdet_r = static_cast<double>(d.n_usadas()) * ld_r0;

  const auto& lw = d.lw;

  for (std::size_t si = 0; si < d.sujeitos.size(); si++) {
    const std::vector<std::size_t>& regs = d.sujeitos[si];
    std::vector<double> ts;
    for (std::size_t r : regs) ts.push_back(d.tempo[r]);
    GammaInv gi = gamma_inv(ts, M.rho);
    M.logdet_r += static_cast<double>(t) * gi.logdet;

    const std::size_t m = regs.size();
    auto peso = [&](std::size_t i, std::size_t j) -> double {
      if (i == j) return gi.diag[i];
      if (i + 1 == j) return gi.sub[i];
      if (j + 1 == i) return gi.sub[j];
      return 0.0;
    };

    // clique do sujeito acumulado localmente (colunas ja expandidas por caracteristica)
    const std::vector<std::uint32_t>& cols_s = d.cols_suj[si];
    const std::size_t nc = cols_s.size();
    std::vector<std::uint32_t> loc(nc ? cols_s.back() + 1 : 0, 0);
    for (std::size_t k = 0; k < nc; k++) loc[cols_s[k]] = static_cast<std::uint32_t>(k);
    Densa acc_s(nc, nc);   // triangulo inferior local, denso

    for (std::size_t i = 0; i < m; i++) {
      const auto& wa = lw[regs[i]];
      for (std::size_t j = i; j < m; j++) {
        const double w = peso(i, j);
        const auto& wb = lw[regs[j]];
        for (const auto& ea : wa)
          for (const auto& eb : wb) {
            std::uint32_t p2 = loc[ea.col], q2 = loc[eb.col];
            double v = ea.val * w * M.r0inv.at(ea.trait, eb.trait) * eb.val;
            if (i == j) {
              if (p2 < q2) continue;               // so o triangulo inferior local
            } else {
              // par i != j: fora da diagonal os dois sentidos surgem e o fold soma; NA
              // diagonal so um sentido e percorrido e a contribuicao verdadeira e dupla
              if (p2 == q2) v *= 2.0;
              if (p2 < q2) std::swap(p2, q2);
            }
            acc_s.at(p2, q2) += v;
          }
        // rhs e yry
        for (const auto& ea : wa) {
          double ry = 0.0;
          for (std::size_t tb = 0; tb < t; tb++)
            ry += M.r0inv.at(ea.trait, tb) * d.y.at(regs[j], tb);
          M.rhs[ea.col] += ea.val * w * ry;
        }
        if (i != j)
          for (const auto& eb : wb) {
            double ry = 0.0;
            for (std::size_t tb = 0; tb < t; tb++)
              ry += M.r0inv.at(eb.trait, tb) * d.y.at(regs[i], tb);
            M.rhs[eb.col] += eb.val * w * ry;
          }
        double yy = 0.0;
        for (std::size_t ta = 0; ta < t; ta++)
          for (std::size_t tb = 0; tb < t; tb++)
            yy += d.y.at(regs[i], ta) * M.r0inv.at(ta, tb) * d.y.at(regs[j], tb);
        M.yry += (i == j ? 1.0 : 2.0) * w * yy;
      }
    }
    // um push por posicao do clique — ZEROS INCLUIDOS, que e o que garante o padrao
    for (std::size_t a2 = 0; a2 < nc; a2++)
      for (std::size_t b2 = 0; b2 <= a2; b2++) {
        std::uint32_t p2 = cols_s[a2], q2 = cols_s[b2];
        if (p2 > q2) std::swap(p2, q2);
        ti.push_back(p2);
        tj.push_back(q2);
        tv.push_back(acc_s.at(a2, b2));
      }
  }

  // penalidade por grupo, unidades absolutas: kron(C_g^-1, K^-1); dim ja inclui t
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    const Grupo& gr = d.modelo.grupos[g];
    const std::size_t dim = gr.dim;
    Densa cg(dim, dim);
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const double v = theta[gr.offset + vech_idx(i, j, dim)];
        cg.at(i, j) = v;
        cg.at(j, i) = v;
      }
    {
      Densa l = cg;
      if (!chol_densa(l)) return M;
    }
    Densa cinv = inv_pd(cg);
    std::size_t nl = 0;
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == gr.termos[0]) { nl = a.n_niveis; break; }
    M.logdet_g += static_cast<double>(nl) * logdet_pd(cg)
                - static_cast<double>(dim) * d.kinv_logdet[g];

    const std::size_t off = M.offset_grupo[g];
    const bool com_k = d.kinv[g].ncol > 0;
    auto poe = [&](std::size_t gi2, std::size_t gj2, double val) {
      if (gi2 > gj2) return;
      ti.push_back(static_cast<std::uint32_t>(gi2));
      tj.push_back(static_cast<std::uint32_t>(gj2));
      tv.push_back(val);
    };
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        const double f = cinv.at(a, b);
        // f == 0 NAO pula: padrao invariante em theta e o contrato do cache da simbolica
        if (com_k) {
          const Csc& k = d.kinv[g];
          for (std::size_t col = 0; col < k.ncol; col++)
            for (std::size_t p = k.colptr[col]; p < k.colptr[col + 1]; p++) {
              const std::size_t rk = k.linha[p];
              poe(off + a * nl + rk, off + b * nl + col, f * k.valor[p]);
              if (rk != col) poe(off + a * nl + col, off + b * nl + rk, f * k.valor[p]);
            }
        } else {
          for (std::size_t l = 0; l < nl; l++) poe(off + a * nl + l, off + b * nl + l, f);
        }
      }
  }

  M.c = de_triplos(M.total, M.total, ti, tj, tv);
  M.ok = true;
  return M;
}

AvaliacaoAR avalia_ar1(const DesenhoAR& d, const std::vector<double>& theta,
                       CacheSimbolica* cache) {
  AvaliacaoAR A;
  const std::size_t t = d.t;
  MontadoAR M = monta_mme_ar1(d, theta);
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

  double bry = 0.0;
  for (std::size_t k = 0; k < M.total; k++) bry += A.solucao[k] * M.rhs[k];
  A.neg2logl = M.logdet_r + M.logdet_g + logdet(L) + (M.yry - bry);

  // e^ = y - W b, por registro e caracteristica, com o esqueleto do desenho
  const auto& lw = d.lw;
  Densa ehat(d.nlin, t);
  for (std::size_t r = 0; r < d.nlin; r++) {
    if (!d.usa[r]) continue;
    for (std::size_t tau = 0; tau < t; tau++) ehat.at(r, tau) = d.y.at(r, tau);
    for (const auto& e : lw[r]) ehat.at(r, e.trait) -= e.val * A.solucao[e.col];
  }

  SelInv z = inversa_seletiva(L, 0);
  std::vector<std::size_t> inv_perm(M.total);
  for (std::size_t k = 0; k < M.total; k++) inv_perm[perm[k]] = k;
  auto z_orig = [&](std::size_t i, std::size_t j, double& out) {
    return z.get(inv_perm[i], inv_perm[j], out);
  };

  // PEV: a diagonal de C^-1 na numeracao original, unidades absolutas
  {
    std::vector<double> dz = z.diagonal();
    A.pev.assign(M.total, std::nan(""));
    for (std::size_t k = 0; k < M.total; k++) A.pev[perm[k]] = dz[k];
  }

  const std::size_t ntheta = d.modelo.ntheta;
  A.score.assign(ntheta, 0.0);
  A.ai = Densa(ntheta, ntheta);

  // ---- grupos: nl C^-1 - C^-1 (Q + T) C^-1, exatamente como na multi
  std::vector<std::vector<std::vector<double>>> vgrupo(d.modelo.grupos.size());
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    const Grupo& gr = d.modelo.grupos[g];
    const std::size_t dim = gr.dim;
    const std::size_t off = M.offset_grupo[g];
    std::size_t nl = 0;
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == gr.termos[0]) { nl = a.n_niveis; break; }

    Densa cg(dim, dim);
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const double v = theta[gr.offset + vech_idx(i, j, dim)];
        cg.at(i, j) = v;
        cg.at(j, i) = v;
      }
    Densa cinv = inv_pd(cg);

    Densa Q(dim, dim), T(dim, dim);
    std::vector<std::vector<double>> ku(dim);
    for (std::size_t b = 0; b < dim; b++) {
      ku[b].assign(nl, 0.0);
      const double* u = A.solucao.data() + off + b * nl;
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
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        const double* ua = A.solucao.data() + off + a * nl;
        for (std::size_t l = 0; l < nl; l++) s += ua[l] * ku[b][l];
        Q.at(a, b) = s;
      }
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
        T.at(a, b) = s;
      }

    Densa QT(dim, dim), CiQT(dim, dim), Mg(dim, dim);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) QT.at(a, b) = Q.at(a, b) + T.at(a, b);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        for (std::size_t k = 0; k < dim; k++) s += cinv.at(a, k) * QT.at(k, b);
        CiQT.at(a, b) = s;
      }
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        for (std::size_t k = 0; k < dim; k++) s += CiQT.at(a, k) * cinv.at(k, b);
        Mg.at(a, b) = static_cast<double>(nl) * cinv.at(a, b) - s;
      }
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const std::size_t k = gr.offset + vech_idx(i, j, dim);
        A.score[k] = (i == j) ? Mg.at(i, i) : Mg.at(i, j) + Mg.at(j, i);
      }

    vgrupo[g].assign(dim, std::vector<double>(nl, 0.0));
    for (std::size_t s2 = 0; s2 < dim; s2++)
      for (std::size_t dd = 0; dd < dim; dd++) {
        const double cq = cinv.at(s2, dd);
        if (cq == 0.0) continue;
        const double* ud = A.solucao.data() + off + dd * nl;
        for (std::size_t l = 0; l < nl; l++) vgrupo[g][s2][l] += cq * ud[l];
      }
  }

  // R^-1 e^ por registro e caracteristica: (Gamma^-1 (x) R0^-1) e
  Densa re(d.nlin, t);
  for (std::size_t si = 0; si < d.sujeitos.size(); si++) {
    const std::vector<std::size_t>& regs = d.sujeitos[si];
    std::vector<double> ts;
    for (std::size_t r : regs) ts.push_back(d.tempo[r]);
    GammaInv gi = gamma_inv(ts, M.rho);
    const std::size_t m = regs.size();
    for (std::size_t i = 0; i < m; i++) {
      for (std::size_t a = 0; a < t; a++) {
        double s = gi.diag[i] * ehat.at(regs[i], a);
        if (i > 0) s += gi.sub[i - 1] * ehat.at(regs[i - 1], a);
        if (i + 1 < m) s += gi.sub[i] * ehat.at(regs[i + 1], a);
        re.at(regs[i], a) = s;
      }
    }
  }
  // multiplica por R0^-1 nas caracteristicas
  for (std::size_t r = 0; r < d.nlin; r++) {
    if (!d.usa[r]) continue;
    double tmp[64];
    for (std::size_t a = 0; a < t; a++) {
      double s = 0.0;
      for (std::size_t b = 0; b < t; b++) s += M.r0inv.at(a, b) * re.at(r, b);
      tmp[a] = s;
    }
    for (std::size_t a = 0; a < t; a++) re.at(r, a) = tmp[a];
  }

  // ---- parametros residuais, unificados: dR = D (x) S
  //   vech(R0)_ab: D = Gamma (denso por sujeito), S = E_ab simetrica
  //   rho:         D = dGamma,                    S = R0
  const std::size_t n_res = t * (t + 1) / 2 + 1;
  {
    // S de cada parametro e R0^-1 S R0^-1
    std::vector<Densa> S(n_res), RSR(n_res);
    {
      std::size_t k = 0;
      for (std::size_t jb = 0; jb < t; jb++)
        for (std::size_t ib = jb; ib < t; ib++) {
          Densa e(t, t);
          e.at(ib, jb) = 1.0;
          e.at(jb, ib) = 1.0;
          S[k] = e;
          k++;
        }
      S[k] = M.r0;   // rho
    }
    for (std::size_t k = 0; k < n_res; k++) {
      Densa a1(t, t), a2(t, t);
      for (std::size_t i = 0; i < t; i++)
        for (std::size_t j = 0; j < t; j++) {
          double s = 0.0;
          for (std::size_t q = 0; q < t; q++) s += M.r0inv.at(i, q) * S[k].at(q, j);
          a1.at(i, j) = s;
        }
      for (std::size_t i = 0; i < t; i++)
        for (std::size_t j = 0; j < t; j++) {
          double s = 0.0;
          for (std::size_t q = 0; q < t; q++) s += a1.at(i, q) * M.r0inv.at(q, j);
          a2.at(i, j) = s;
        }
      RSR[k] = a2;
    }
    std::vector<double> tr_rinv(n_res, 0.0), meio(n_res, 0.0), quad(n_res, 0.0);

    for (std::size_t si = 0; si < d.sujeitos.size(); si++) {
      const std::vector<std::size_t>& regs = d.sujeitos[si];
      std::vector<double> ts;
      for (std::size_t r : regs) ts.push_back(d.tempo[r]);
      const std::size_t m = regs.size();
      Densa G = gamma_densa(ts, M.rho);
      Densa D = dgamma_densa(ts, M.rho);
      Densa Gi = inv_pd(G);
      GammaInv gi = gamma_inv(ts, M.rho);

      // tr(Gamma^-1 D_k): para vech(R0), D = Gamma -> traco = m; para rho, tr(Gi D)
      double trGiD_rho = 0.0;
      for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < m; j++) trGiD_rho += Gi.at(i, j) * D.at(j, i);
      for (std::size_t k = 0; k + 1 < n_res; k++) {
        double trRS = 0.0;
        for (std::size_t a = 0; a < t; a++)
          for (std::size_t b = 0; b < t; b++) trRS += M.r0inv.at(a, b) * S[k].at(b, a);
        tr_rinv[k] += static_cast<double>(m) * trRS;
      }
      tr_rinv[n_res - 1] += trGiD_rho * static_cast<double>(t);

      // M2 no espaco dos registros: para vech(R0), Gi Gamma Gi = Gi (tridiagonal!);
      // para rho, Gi D Gi denso
      Densa M2r(m, m);
      for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < m; j++) {
          double s = 0.0;
          for (std::size_t a = 0; a < m; a++)
            for (std::size_t b = 0; b < m; b++)
              s += Gi.at(i, a) * D.at(a, b) * Gi.at(b, j);
          M2r.at(i, j) = s;
        }
      auto giat = [&](std::size_t i, std::size_t j) -> double {
        if (i == j) return gi.diag[i];
        if (i + 1 == j) return gi.sub[i];
        if (j + 1 == i) return gi.sub[j];
        return 0.0;
      };

      // quad_k = soma_ij D_k[i,j] . re_i' S_k re_j   (D = Gamma para R0; dGamma para rho)
      for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < m; j++) {
          const double dG = G.at(i, j);
          const double dD = D.at(i, j);
          if (dG == 0.0 && dD == 0.0) continue;
          for (std::size_t k = 0; k < n_res; k++) {
            const double dk = (k + 1 < n_res) ? dG : dD;
            if (dk == 0.0) continue;
            double s = 0.0;
            const Densa& Sk = S[k];
            for (std::size_t a = 0; a < t; a++)
              for (std::size_t b = 0; b < t; b++)
                s += re.at(regs[i], a) * Sk.at(a, b) * re.at(regs[j], b);
            quad[k] += dk * s;
          }
        }

      // meio_k = tr(C^-1 W' [M2^rec_k (x) RSR_k] W), lido na inversa seletiva.
      // Para vech(R0) o M2 de registros e o proprio Gamma^-1 (TRIDIAGONAL); para rho e
      // denso. Um laco so, com o peso certo por parametro.
      for (std::size_t i = 0; i < m; i++) {
        const auto& wi = lw[regs[i]];
        for (std::size_t j = 0; j < m; j++) {
          const double m2_r0 = giat(i, j);          // tridiagonal
          const double m2_rho = M2r.at(i, j);       // denso
          if (m2_r0 == 0.0 && m2_rho == 0.0) continue;
          const auto& wj = lw[regs[j]];
          for (const auto& ea : wi)
            for (const auto& eb : wj) {
              double zv;
              if (!z_orig(ea.col, eb.col, zv)) { A.fora_do_padrao++; continue; }
              const double base = ea.val * eb.val * zv;
              for (std::size_t k = 0; k < n_res; k++) {
                const double m2 = (k + 1 < n_res) ? m2_r0 : m2_rho;
                if (m2 == 0.0) continue;
                meio[k] += base * m2 * RSR[k].at(ea.trait, eb.trait);
              }
            }
        }
      }
    }
    for (std::size_t k = 0; k < n_res; k++)
      A.score[d.offset_s2e + k] = tr_rinv[k] - meio[k] - quad[k];
  }

  // ---- AI: f_k' R^-1 f_l - (W'R^-1 f_k)' C^-1 (W'R^-1 f_l)
  {
    std::vector<Densa> f(ntheta, Densa(d.nlin, t));
    // grupos (o mapeamento coeficiente->coluna e o da multi: por termo, tau-major)
    for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
      const Grupo& gr = d.modelo.grupos[g];
      const std::size_t dim = gr.dim;
      std::size_t nl = 0;
      for (const DesenhoTermo& a : d.aleatorios)
        if (a.termo == gr.termos[0]) { nl = a.n_niveis; break; }
      // ncoef total (sem o fator t) e a lista de termos na ordem
      std::size_t ncoef_tot = 0;
      for (std::size_t tm : gr.termos)
        for (const DesenhoTermo& a : d.aleatorios)
          if (a.termo == tm) ncoef_tot += a.n_coef;
      for (std::size_t j = 0; j < dim; j++)
        for (std::size_t i = j; i < dim; i++) {
          const std::size_t k = gr.offset + vech_idx(i, j, dim);
          std::vector<std::vector<double>> w(dim, std::vector<double>(nl, 0.0));
          if (i == j) w[i] = vgrupo[g][i];
          else { w[i] = vgrupo[g][j]; w[j] = vgrupo[g][i]; }
          // coeficiente b (tau-major sobre a lista de termos): coluna do nivel nv =
          // colbase(termo) + (tau * n_coef + ct) * nl + nv, e o registro contribui na
          // CARACTERISTICA tau
          for (std::size_t b = 0; b < dim; b++) {
            const std::vector<double>& ws = w[b];
            bool tem = false;
            for (std::size_t l = 0; l < nl; l++) if (ws[l] != 0.0) { tem = true; break; }
            if (!tem) continue;
            const std::size_t tau = b / ncoef_tot;
            std::size_t cpos = b % ncoef_tot;
            std::size_t cbase = 0;
            for (std::size_t tm : gr.termos) {
              for (const DesenhoTermo& a : d.aleatorios)
                if (a.termo == tm) {
                  if (cpos < a.n_coef) {
                    const std::size_t ct = cpos;
                    for (std::size_t nv = 0; nv < nl; nv++) {
                      const double xv = ws[nv];
                      if (xv == 0.0) continue;
                      const std::size_t colz = ct * a.n_niveis + nv;
                      for (std::size_t p = a.z.colptr[colz]; p < a.z.colptr[colz + 1]; p++)
                        if (d.usa[a.z.linha[p]])
                          f[k].at(a.z.linha[p], tau) += a.z.valor[p] * xv;
                    }
                    cpos = static_cast<std::size_t>(-1);
                  } else {
                    cpos -= a.n_coef;
                  }
                  (void) cbase;
                }
              if (cpos == static_cast<std::size_t>(-1)) break;
            }
          }
        }
    }
    // residuais: f_k = (D_k (x) S_k) (R^-1 e)
    for (std::size_t si = 0; si < d.sujeitos.size(); si++) {
      const std::vector<std::size_t>& regs = d.sujeitos[si];
      std::vector<double> ts;
      for (std::size_t r : regs) ts.push_back(d.tempo[r]);
      const std::size_t m = regs.size();
      Densa G = gamma_densa(ts, M.rho);
      Densa D = dgamma_densa(ts, M.rho);
      for (std::size_t i = 0; i < m; i++) {
        // acumula por caracteristica: soma_j D[i,j] re_j primeiro
        std::vector<double> sG(t, 0.0), sD(t, 0.0);
        for (std::size_t j = 0; j < m; j++) {
          const double dg = G.at(i, j), dd = D.at(i, j);
          for (std::size_t b = 0; b < t; b++) {
            sG[b] += dg * re.at(regs[j], b);
            sD[b] += dd * re.at(regs[j], b);
          }
        }
        std::size_t k = 0;
        for (std::size_t jb = 0; jb < t; jb++)
          for (std::size_t ib = jb; ib < t; ib++) {
            // S = E_ab simetrica
            f[d.offset_s2e + k].at(regs[i], ib) += sG[jb];
            if (ib != jb) f[d.offset_s2e + k].at(regs[i], jb) += sG[ib];
            k++;
          }
        for (std::size_t a = 0; a < t; a++) {
          double s = 0.0;
          for (std::size_t b = 0; b < t; b++) s += M.r0.at(a, b) * sD[b];
          f[d.offset_rho].at(regs[i], a) = s;
        }
      }
    }

    // R^-1 f e W'R^-1 f
    std::vector<Densa> rf(ntheta, Densa(d.nlin, t));
    for (std::size_t k = 0; k < ntheta; k++) {
      for (std::size_t si = 0; si < d.sujeitos.size(); si++) {
        const std::vector<std::size_t>& regs = d.sujeitos[si];
        std::vector<double> ts;
        for (std::size_t r : regs) ts.push_back(d.tempo[r]);
        GammaInv gi = gamma_inv(ts, M.rho);
        const std::size_t m = regs.size();
        for (std::size_t i = 0; i < m; i++)
          for (std::size_t a = 0; a < t; a++) {
            double s = gi.diag[i] * f[k].at(regs[i], a);
            if (i > 0) s += gi.sub[i - 1] * f[k].at(regs[i - 1], a);
            if (i + 1 < m) s += gi.sub[i] * f[k].at(regs[i + 1], a);
            rf[k].at(regs[i], a) = s;
          }
      }
      for (std::size_t r = 0; r < d.nlin; r++) {
        if (!d.usa[r]) continue;
        double tmp[64];
        for (std::size_t a = 0; a < t; a++) {
          double s = 0.0;
          for (std::size_t b = 0; b < t; b++) s += M.r0inv.at(a, b) * rf[k].at(r, b);
          tmp[a] = s;
        }
        for (std::size_t a = 0; a < t; a++) rf[k].at(r, a) = tmp[a];
      }
    }
    std::vector<std::vector<double>> wtf(ntheta, std::vector<double>(M.total, 0.0));
    for (std::size_t r = 0; r < d.nlin; r++) {
      if (!d.usa[r]) continue;
      for (std::size_t k = 0; k < ntheta; k++)
        for (const auto& e : lw[r]) wtf[k][e.col] += e.val * rf[k].at(r, e.trait);
    }
    std::vector<std::vector<double>> cwf(ntheta);
    for (std::size_t k = 0; k < ntheta; k++) {
      std::vector<double> pb2(M.total);
      for (std::size_t q2 = 0; q2 < M.total; q2++) pb2[q2] = wtf[k][perm[q2]];
      std::vector<double> px2 = resolve(L, pb2);
      cwf[k].assign(M.total, 0.0);
      for (std::size_t q2 = 0; q2 < M.total; q2++) cwf[k][perm[q2]] = px2[q2];
    }
    for (std::size_t a = 0; a < ntheta; a++)
      for (std::size_t b = a; b < ntheta; b++) {
        double frf = 0.0;
        for (std::size_t r = 0; r < d.nlin; r++)
          if (d.usa[r])
            for (std::size_t tau = 0; tau < t; tau++)
              frf += f[a].at(r, tau) * rf[b].at(r, tau);
        double wcw = 0.0;
        for (std::size_t q2 = 0; q2 < M.total; q2++) wcw += wtf[a][q2] * cwf[b][q2];
        const double v = frf - wcw;
        A.ai.at(a, b) = v;
        A.ai.at(b, a) = v;
      }
  }

  A.ok = true;
  return A;
}

// -2logL pela forma V densa, com R = diag_s(Gamma_s (x) R0) explicita. So para os gates.
double neg2logl_densa_V_ar1(const DesenhoAR& d, const std::vector<double>& theta) {
  const std::size_t t = d.t;
  const double rho = theta[d.offset_rho];
  if (!(std::fabs(rho) < 1.0)) return std::nan("");
  Densa r0(t, t);
  for (std::size_t j = 0; j < t; j++)
    for (std::size_t i = j; i < t; i++) {
      const double v = theta[d.offset_s2e + vech_idx(i, j, t)];
      r0.at(i, j) = v;
      r0.at(j, i) = v;
    }
  {
    Densa l = r0;
    if (!chol_densa(l)) return std::nan("");
  }

  std::vector<std::size_t> linhas;
  for (std::size_t i = 0; i < d.nlin; i++)
    if (d.usa[i]) linhas.push_back(i);
  const std::size_t Nr = linhas.size();
  const std::size_t N = Nr * t;              // (registro, caracteristica), registro-major
  std::vector<std::size_t> pos(d.nlin, static_cast<std::size_t>(-1));
  for (std::size_t k = 0; k < Nr; k++) pos[linhas[k]] = k;

  Densa V(N, N);
  for (std::size_t si = 0; si < d.sujeitos.size(); si++) {
    const std::vector<std::size_t>& regs = d.sujeitos[si];
    std::vector<double> ts;
    for (std::size_t r : regs) ts.push_back(d.tempo[r]);
    Densa G = gamma_densa(ts, rho);
    for (std::size_t i = 0; i < regs.size(); i++)
      for (std::size_t j = 0; j < regs.size(); j++)
        for (std::size_t a = 0; a < t; a++)
          for (std::size_t b = 0; b < t; b++)
            V.at(pos[regs[i]] * t + a, pos[regs[j]] * t + b) = G.at(i, j) * r0.at(a, b);
  }

  std::vector<std::size_t> off_grupo(d.modelo.grupos.size());
  {
    std::size_t acc = d.x.ncol * t;
    for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
      off_grupo[g] = acc;
      acc += d.largura(g);
    }
  }

  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    const Grupo& gr = d.modelo.grupos[g];
    const std::size_t dim = gr.dim;
    std::size_t nl = 0;
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == gr.termos[0]) { nl = a.n_niveis; break; }
    Densa cg(dim, dim);
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const double v = theta[gr.offset + vech_idx(i, j, dim)];
        cg.at(i, j) = v;
        cg.at(j, i) = v;
      }
    Densa K;
    if (d.kinv[g].ncol > 0) K = inv_pd(d.kinv[g].densa_simetrica());
    else { K = Densa(nl, nl); for (std::size_t i = 0; i < nl; i++) K.at(i, i) = 1.0; }

    // Z (N x dim*nl) com o mapeamento da montagem: coeficiente b tau-major sobre termos
    std::size_t ncoef_tot = dim / t;
    const std::size_t larg = dim * nl;
    Densa Z(N, larg);
    for (std::size_t b = 0; b < dim; b++) {
      const std::size_t tau = b / ncoef_tot;
      std::size_t cpos = b % ncoef_tot;
      for (std::size_t tm : gr.termos) {
        bool feito = false;
        for (const DesenhoTermo& a : d.aleatorios)
          if (a.termo == tm) {
            if (cpos < a.n_coef) {
              Densa zd = a.z.densa();
              for (std::size_t k = 0; k < Nr; k++)
                for (std::size_t nv = 0; nv < nl; nv++) {
                  const double v = zd.at(linhas[k], cpos * nl + nv);
                  if (v != 0.0) Z.at(k * t + tau, b * nl + nv) = v;
                }
              feito = true;
            } else {
              cpos -= a.n_coef;
            }
          }
        if (feito) break;
      }
    }
    Densa CK(larg, larg);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++)
        for (std::size_t i = 0; i < nl; i++)
          for (std::size_t j = 0; j < nl; j++)
            CK.at(a * nl + i, b * nl + j) = cg.at(a, b) * K.at(i, j);
    Densa ZC(N, larg);
    for (std::size_t r2 = 0; r2 < N; r2++)
      for (std::size_t c2 = 0; c2 < larg; c2++) {
        double s = 0.0;
        for (std::size_t k = 0; k < larg; k++) s += Z.at(r2, k) * CK.at(k, c2);
        ZC.at(r2, c2) = s;
      }
    for (std::size_t r2 = 0; r2 < N; r2++)
      for (std::size_t c2 = 0; c2 < N; c2++) {
        double s = 0.0;
        for (std::size_t k = 0; k < larg; k++) s += ZC.at(r2, k) * Z.at(c2, k);
        V.at(r2, c2) += s;
      }
  }

  const double ldV = logdet_pd(V);
  if (std::isnan(ldV)) return std::nan("");
  Densa Vi = inv_pd(V);

  const std::size_t p = d.x.ncol * t;
  Densa X(N, p);
  std::vector<double> y(N);
  for (std::size_t k = 0; k < Nr; k++)
    for (std::size_t a = 0; a < t; a++) {
      y[k * t + a] = d.y.at(linhas[k], a);
      for (std::size_t j = 0; j < d.x.ncol; j++)
        X.at(k * t + a, j * t + a) = d.x.at(linhas[k], j);
    }
  Densa XtVi(p, N);
  for (std::size_t a = 0; a < p; a++)
    for (std::size_t r2 = 0; r2 < N; r2++) {
      double s = 0.0;
      for (std::size_t k = 0; k < N; k++) s += X.at(k, a) * Vi.at(k, r2);
      XtVi.at(a, r2) = s;
    }
  Densa XtViX(p, p);
  std::vector<double> XtViy(p, 0.0);
  for (std::size_t a = 0; a < p; a++) {
    for (std::size_t b = 0; b < p; b++) {
      double s = 0.0;
      for (std::size_t r2 = 0; r2 < N; r2++) s += XtVi.at(a, r2) * X.at(r2, b);
      XtViX.at(a, b) = s;
    }
    for (std::size_t r2 = 0; r2 < N; r2++) XtViy[a] += XtVi.at(a, r2) * y[r2];
  }
  const double ldX = logdet_pd(XtViX);
  if (std::isnan(ldX)) return std::nan("");
  Densa Xi = inv_pd(XtViX);
  double yViy = 0.0;
  for (std::size_t r2 = 0; r2 < N; r2++) {
    double s = 0.0;
    for (std::size_t k = 0; k < N; k++) s += Vi.at(r2, k) * y[k];
    yViy += y[r2] * s;
  }
  double quad = 0.0;
  for (std::size_t a = 0; a < p; a++)
    for (std::size_t b = 0; b < p; b++) quad += XtViy[a] * Xi.at(a, b) * XtViy[b];
  return ldV + ldX + (yViy - quad);
}

}  // namespace br
