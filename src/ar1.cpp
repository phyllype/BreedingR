// Residuo AR(1)/CAR(1) (Wade e Quaas, 1993), agora para t caracteristicas: correlacao rho^|dt| dentro de
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
// No padrao de C entram so os pares de colunas do mesmo registro e de registros vizinhos
// (Gamma^-1 tridiagonal), e o score, a AI e os residuos so leem esses pares da inversa
// seletiva; fora_do_padrao denuncia o contrario. Empate de tempo no mesmo sujeito continua
// erro declarado. O tempo que chega aqui e o do eixo de eixo_tempo_ar1 (ar1b.cpp): na grade,
// o numero de passos, inteiro exato.

#include "mme.h"
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <R_ext/Random.h>

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

// O SINAL de r' (antes em dgamma_densa, que saiu):
// Gamma(dt) = s(dt) |rho|^dt, com s(dt) = -1 quando rho < 0 e fmod(dt, 2) != 0, que e a
// convencao de gamma_densa. Derivando essa MESMA expressao, para rho < 0 vale
// |rho| = -rho, logo dGamma/drho = -s(dt) dt |rho|^(dt-1). O sinal e portanto uma funcao de
// fmod(dt, 2), e nao de fmod(dt - 1, 2): as duas coincidem para dt INTEIRO e divergem fora
// dele. Com dt = 2.5 e rho < 0, a regra antiga devolvia o sinal trocado, e a
// verossimilhanca e o gradiente passavam a descrever Gammas diferentes: o caminhante
// descia uma superficie que nao estava medindo.
// dGamma^-1/drho, TRIDIAGONAL como a propria Gamma^-1, e dlog|Gamma|/drho. Por intervalo
// com r = Gamma(dt) e r' = dr/drho (a convencao de sinal acima), u = 1 - r^2:
// a diagonal dos dois registros ganha 2 r r'/u^2, o elo entre eles -r'(1 + r^2)/u^2, e
// log|Gamma| = soma log u da -2 r r'/u. Com isto o score e a AI do rho saem em O(m) por
// sujeito; antes montavam Gamma e dGamma densas e Gamma^-1 dGamma Gamma^-1 num laco
// quadruplo, O(m^4) por sujeito e por avaliacao.
struct DGammaInv {
  std::vector<double> ddiag, dsub;
  double dlogdet = 0.0;
};

static DGammaInv dgamma_inv(const std::vector<double>& t, double rho) {
  const std::size_t m = t.size();
  DGammaInv g;
  g.ddiag.assign(m, 0.0);
  g.dsub.assign(m > 0 ? m - 1 : 0, 0.0);
  for (std::size_t k = 0; k + 1 < m; k++) {
    const double dt = t[k + 1] - t[k];
    const double r = std::pow(std::fabs(rho), dt) *
                     ((rho < 0.0 && std::fmod(dt, 2.0) != 0.0) ? -1.0 : 1.0);
    double rp;
    if (dt == 1.0) rp = 1.0;
    else if (rho == 0.0) rp = 0.0;
    else {
      rp = dt * std::pow(std::fabs(rho), dt - 1.0);
      if (rho < 0.0 && std::fmod(dt, 2.0) == 0.0) rp = -rp;
    }
    const double u = 1.0 - r * r;
    g.ddiag[k]     += 2.0 * r * rp / (u * u);
    g.ddiag[k + 1] += 2.0 * r * rp / (u * u);
    g.dsub[k]       = -rp * (1.0 + r * r) / (u * u);
    g.dlogdet      += -2.0 * r * rp / u;
  }
  return g;
}

// x com Gamma^-1 x = v (Gamma^-1 tridiagonal simetrica positiva-definida, Thomas): e o
// produto Gamma v sem formar Gamma
static std::vector<double> resolve_tri(const GammaInv& gi, const std::vector<double>& v) {
  const std::size_t m = v.size();
  std::vector<double> c(m, 0.0), d(m, 0.0), x(m, 0.0);
  if (m == 0) return x;
  double b0 = gi.diag[0];
  c[0] = m > 1 ? gi.sub[0] / b0 : 0.0;
  d[0] = v[0] / b0;
  for (std::size_t i = 1; i < m; i++) {
    const double den = gi.diag[i] - gi.sub[i - 1] * c[i - 1];
    c[i] = (i + 1 < m) ? gi.sub[i] / den : 0.0;
    d[i] = (v[i] - gi.sub[i - 1] * d[i - 1]) / den;
  }
  x[m - 1] = d[m - 1];
  for (std::size_t i = m - 1; i-- > 0;) x[i] = d[i] - c[i] * x[i + 1];
  return x;
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
  M.total = acc + d.n_mv;   // as colunas mv vem depois dos grupos, sem penalidade

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

    // W' R^-1 W do sujeito. Gamma^-1 e TRIDIAGONAL, entao so entram os pares de colunas do
    // MESMO registro e de registros VIZINHOS no tempo, e e so isso que vai para o padrao de C
    // (sempre os mesmos pares, com zero incluido quando rho = 0: o padrao nao depende de
    // theta, que e o contrato do cache da simbolica). Antes o padrao levava o clique COMPLETO
    // do sujeito, com zeros explicitos, herdado do tempo em que o traco do rho lia a inversa
    // seletiva em pares de registros distantes; hoje o score, a AI e os residuos so leem
    // pares vizinhos (fora_do_padrao denuncia o contrario). Numa serie longa com um CG de
    // tempo o clique juntava todos os niveis num bloco denso: 1 sujeito x 20000 dias com 2858
    // CGs passava de 169 s e 2 GB (medido em f2c0d69). As contribuicoes de cada par de
    // colunas somam na ordem do laco, como no acumulador denso de antes.
    std::vector<std::uint64_t> chave;
    std::vector<double> valor;
    for (std::size_t i = 0; i < m; i++) {
      const auto& wa = lw[regs[i]];
      const std::size_t j1 = std::min(m - 1, i + 1);
      for (std::size_t j = i; j <= j1; j++) {
        const double w = (i == j) ? gi.diag[i] : gi.sub[i];
        const auto& wb = lw[regs[j]];
        for (const auto& ea : wa)
          for (const auto& eb : wb) {
            std::uint32_t p2 = ea.col, q2 = eb.col;
            double v = ea.val * w * M.r0inv.at(ea.trait, eb.trait) * eb.val;
            if (i == j) {
              if (p2 > q2) continue;               // so o triangulo superior
            } else {
              // par i != j: fora da diagonal os dois sentidos surgem e o fold soma; NA
              // diagonal so um sentido e percorrido e a contribuicao verdadeira e dupla
              if (p2 == q2) v *= 2.0;
              if (p2 > q2) std::swap(p2, q2);
            }
            chave.push_back((static_cast<std::uint64_t>(q2) << 32) | p2);
            valor.push_back(v);
          }
        // rhs (y'R^-1 y nao e mais acumulado: a forma quadratica sai dos residuos, em
        // avalia_ar1)
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
      }
    }
    // um push por par de colunas do sujeito, ZEROS INCLUIDOS, que e o que garante o padrao;
    // a ordenacao estavel soma as contribuicoes de um par na ordem em que vieram
    std::vector<std::size_t> ord(chave.size());
    std::iota(ord.begin(), ord.end(), static_cast<std::size_t>(0));
    std::stable_sort(ord.begin(), ord.end(),
                     [&](std::size_t a, std::size_t b) { return chave[a] < chave[b]; });
    for (std::size_t k = 0; k < ord.size();) {
      const std::uint64_t c = chave[ord[k]];
      double s = 0.0;
      while (k < ord.size() && chave[ord[k]] == c) s += valor[ord[k++]];
      ti.push_back(static_cast<std::uint32_t>(c & 0xffffffffu));
      tj.push_back(static_cast<std::uint32_t>(c >> 32));
      tv.push_back(s);
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
  // rho NEGATIVO exige grade inteira, e a razao esta no campo tempo_inteiro (mme.h): fora
  // dela Gamma(dt) = rho^dt com rho < 0 nao e uma funcao de correlacao valida, a
  // multiplicatividade quebra e as rotas densa e esparsa deixam de concordar. Marcar o
  // theta como inadmissivel e o mesmo tratamento que uma covariancia fora do cone recebe: o
  // passo rejeita e o amortecimento sobe, sem projetar nada em silencio.
  if (!d.tempo_inteiro && theta[d.offset_rho] < 0.0) return A;
  MontadoAR M = monta_mme_ar1(d, theta);
  if (!M.ok) return A;

  CacheSimbolica local;
  CacheSimbolica* cs = cache ? cache : &local;
  if (!cs->pronto) {
    cs->perm = grau_minimo(M.c);
    cs->sb = simbolica(permuta_sim(M.c, cs->perm));
    cs->pronto = true;
  }
  const std::vector<std::size_t>& perm = cs->perm;
  const Simbolica& sb = cs->sb;
  const Csc& pc = permuta_cache(M.c, *cs);
  Csc L;
  if (!cholesky(pc, sb, L)) return A;

  std::vector<double> pb(M.total);
  for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[perm[k]];
  std::vector<double> px = resolve(L, pb);
  A.solucao.assign(M.total, 0.0);
  for (std::size_t k = 0; k < M.total; k++) A.solucao[perm[k]] = px[k];

  // -2logL: os termos que nao dependem de y aqui; y'Py entra depois de R^-1 e, pelos
  // residuos, e'R^-1 e + u'G^-1 u, e nao por y'R^-1 y - b'r. O racional (invariancia a
  // uma constante somada a y, medida no colapso rho = 0 contra o model()) esta em avalia(),
  // src/aireml.cpp. As celulas ausentes entram com y = 0 e a coluna mv delas: a equacao da
  // coluna mv zera (R^-1 e) na celula, e a identidade e'R^-1 e + b'Pen b = y'R^-1 y - b'r
  // continua exata.
  const double parte_fixa = M.logdet_r + M.logdet_g + logdet(L);
  double upen = 0.0;   // soma_g soma_ab C_g^-1[a,b] Q[a,b] = u'G^-1 u

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
  // proposta EM: comeca no theta corrente, os grupos sao sobrescritos abaixo, e R0 e
  // rho ficam onde estao, o rho nao tem passo M fechado e o R0 multivariado ainda
  // nao tem o seu. Quem chama so aceita a proposta se ela baixar a verossimilhanca.
  A.em_theta = theta;
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
      for (std::size_t b = 0; b < dim; b++) upen += cinv.at(a, b) * Q.at(a, b);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) {
        double s = 0.0;
        if (d.kinv[g].ncol > 0) {
          s = traco_kinv(d.kinv[g], nl, off + a * nl, off + b * nl, z_orig, A.fora_do_padrao);
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
    // EM do grupo: C_g <- (Q + T) / nl, com o QT que o score ja montou acima
    for (std::size_t j2 = 0; j2 < dim; j2++)
      for (std::size_t i2 = j2; i2 < dim; i2++)
        A.em_theta[gr.offset + vech_idx(i2, j2, dim)] =
            (QT.at(i2, j2) + QT.at(j2, i2)) / (2.0 * static_cast<double>(nl));
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
  {
    double ere = 0.0;
    for (std::size_t r = 0; r < d.nlin; r++) {
      if (!d.usa[r]) continue;
      for (std::size_t a = 0; a < t; a++) ere += ehat.at(r, a) * re.at(r, a);
    }
    A.neg2logl = parte_fixa + ere + upen;
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
      GammaInv gi = gamma_inv(ts, M.rho);
      DGammaInv dq = dgamma_inv(ts, M.rho);

      // tr(R^-1 dR): para vech(R0), m tr(R0^-1 S); para rho, t dlog|Gamma|/drho
      for (std::size_t k = 0; k + 1 < n_res; k++) {
        double trRS = 0.0;
        for (std::size_t a = 0; a < t; a++)
          for (std::size_t b = 0; b < t; b++) trRS += M.r0inv.at(a, b) * S[k].at(b, a);
        tr_rinv[k] += static_cast<double>(m) * trRS;
      }
      tr_rinv[n_res - 1] += dq.dlogdet * static_cast<double>(t);

      // R^-1 dR R^-1 no espaco dos registros, TRIDIAGONAL nos dois casos: para vech(R0)
      // e Gamma^-1 (com R0^-1 S R0^-1 nas caracteristicas), para rho e
      // Gamma^-1 dGamma Gamma^-1 = -dGamma^-1/drho (com R0^-1)
      auto giat = [&](std::size_t i, std::size_t j) -> double {
        if (i == j) return gi.diag[i];
        if (i + 1 == j) return gi.sub[i];
        if (j + 1 == i) return gi.sub[j];
        return 0.0;
      };
      auto dqat = [&](std::size_t i, std::size_t j) -> double {
        if (i == j) return -dq.ddiag[i];
        if (i + 1 == j) return -dq.dsub[i];
        if (j + 1 == i) return -dq.dsub[j];
        return 0.0;
      };

      // quad_k = e' (W_k (x) RSR_k) e, com W_k tridiagonal: re'(D (x) S)re reescrito em e
      for (std::size_t i = 0; i < m; i++) {
        const std::size_t j0 = i > 0 ? i - 1 : 0, j1 = std::min(m - 1, i + 1);
        for (std::size_t j = j0; j <= j1; j++) {
          const double w_r0 = giat(i, j), w_rho = dqat(i, j);
          for (std::size_t k = 0; k < n_res; k++) {
            const double wk = (k + 1 < n_res) ? w_r0 : w_rho;
            if (wk == 0.0) continue;
            double s2 = 0.0;
            for (std::size_t a = 0; a < t; a++)
              for (std::size_t b = 0; b < t; b++)
                s2 += ehat.at(regs[i], a) * RSR[k].at(a, b) * ehat.at(regs[j], b);
            quad[k] += wk * s2;
          }
        }
      }

      // meio_k = tr(C^-1 W' [W_k (x) RSR_k] W), lido na inversa seletiva so nos pares de
      // registros VIZINHOS, que estao no padrao do fator
      for (std::size_t i = 0; i < m; i++) {
        const auto& wi = lw[regs[i]];
        const std::size_t j0 = i > 0 ? i - 1 : 0, j1 = std::min(m - 1, i + 1);
        for (std::size_t j = j0; j <= j1; j++) {
          const double m2_r0 = giat(i, j);
          const double m2_rho = dqat(i, j);
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
      // f_k = (D_k (x) S_k) R^-1 e sem Gamma densa: para vech(R0) (Gamma (x) S)R^-1 e =
      // (I (x) S R0^-1) e, local ao registro; para rho (dGamma (x) R0) R^-1 e =
      // -(Gamma dGamma^-1 (x) I) e, uma solucao tridiagonal por caracteristica
      GammaInv gi = gamma_inv(ts, M.rho);
      DGammaInv dq = dgamma_inv(ts, M.rho);
      for (std::size_t i = 0; i < m; i++) {
        std::vector<double> sG(t, 0.0);
        for (std::size_t b = 0; b < t; b++)
          for (std::size_t c2 = 0; c2 < t; c2++) sG[b] += M.r0inv.at(b, c2) * ehat.at(regs[i], c2);
        std::size_t k = 0;
        for (std::size_t jb = 0; jb < t; jb++)
          for (std::size_t ib = jb; ib < t; ib++) {
            f[d.offset_s2e + k].at(regs[i], ib) += sG[jb];
            if (ib != jb) f[d.offset_s2e + k].at(regs[i], jb) += sG[ib];
            k++;
          }
      }
      for (std::size_t a = 0; a < t; a++) {
        std::vector<double> v(m, 0.0);
        for (std::size_t i = 0; i < m; i++) {
          double s2 = dq.ddiag[i] * ehat.at(regs[i], a);
          if (i > 0) s2 += dq.dsub[i - 1] * ehat.at(regs[i - 1], a);
          if (i + 1 < m) s2 += dq.dsub[i] * ehat.at(regs[i + 1], a);
          v[i] = s2;
        }
        std::vector<double> x = resolve_tri(gi, v);
        for (std::size_t i = 0; i < m; i++) f[d.offset_rho].at(regs[i], a) = -x[i];
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

  // SO AS CELULAS OBSERVADAS: a verossimilhanca marginal delas, pela submatriz de V (a rota
  // do SAS), que nao passa pelas colunas mv da rota esparsa. Sem celula ausente e a V de
  // sempre.
  std::vector<std::size_t> cel;
  for (std::size_t k = 0; k < Nr; k++)
    for (std::size_t a = 0; a < t; a++)
      if (d.obs.empty() || d.obs[linhas[k] * t + a]) cel.push_back(k * t + a);
  if (cel.size() != N) {
    Densa Vo(cel.size(), cel.size());
    for (std::size_t i2 = 0; i2 < cel.size(); i2++)
      for (std::size_t j2 = 0; j2 < cel.size(); j2++) Vo.at(i2, j2) = V.at(cel[i2], cel[j2]);
    V = Vo;
  }
  const double ldV = logdet_pd(V);
  if (std::isnan(ldV)) return std::nan("");
  Densa Vi = inv_pd(V);

  const std::size_t p = d.x.ncol * t;
  Densa Xc(N, p);
  std::vector<double> yc(N);
  for (std::size_t k = 0; k < Nr; k++)
    for (std::size_t a = 0; a < t; a++) {
      yc[k * t + a] = d.y.at(linhas[k], a);
      for (std::size_t ex = d.x.linptr[linhas[k]]; ex < d.x.linptr[linhas[k] + 1]; ex++)
        Xc.at(k * t + a, d.x.coluna[ex] * t + a) = d.x.valor[ex];
    }
  const std::size_t No = cel.size();
  Densa X(No, p);
  std::vector<double> y(No);
  for (std::size_t i2 = 0; i2 < No; i2++) {
    y[i2] = yc[cel[i2]];
    for (std::size_t a = 0; a < p; a++) X.at(i2, a) = Xc.at(cel[i2], a);
  }
  Densa XtVi(p, No);
  for (std::size_t a = 0; a < p; a++)
    for (std::size_t r2 = 0; r2 < No; r2++) {
      double s = 0.0;
      for (std::size_t k = 0; k < No; k++) s += X.at(k, a) * Vi.at(k, r2);
      XtVi.at(a, r2) = s;
    }
  Densa XtViX(p, p);
  std::vector<double> XtViy(p, 0.0);
  for (std::size_t a = 0; a < p; a++) {
    for (std::size_t b = 0; b < p; b++) {
      double s = 0.0;
      for (std::size_t r2 = 0; r2 < No; r2++) s += XtVi.at(a, r2) * X.at(r2, b);
      XtViX.at(a, b) = s;
    }
    for (std::size_t r2 = 0; r2 < No; r2++) XtViy[a] += XtVi.at(a, r2) * y[r2];
  }
  const double ldX = logdet_pd(XtViX);
  if (std::isnan(ldX)) return std::nan("");
  Densa Xi = inv_pd(XtViX);
  double yViy = 0.0;
  for (std::size_t r2 = 0; r2 < No; r2++) {
    double s = 0.0;
    for (std::size_t k = 0; k < No; k++) s += Vi.at(r2, k) * y[k];
    yViy += y[r2] * s;
  }
  double quad = 0.0;
  for (std::size_t a = 0; a < p; a++)
    for (std::size_t b = 0; b < p; b++) quad += XtViy[a] * Xi.at(a, b) * XtViy[b];
  return ldV + ldX + (yViy - quad);
}

// ---------------------------------------------------------------- residuos no otimo
//
// As MME daqui estao em unidades absolutas (C = W'R^-1 W + G^-1), entao
//   Var(e^_{r,tau}) = R0[tau,tau] - w_{r,tau}' C^-1 w_{r,tau}
// com w_{r,tau} as entradas da linha r na caracteristica tau (lw). Os pares de colunas de uma
// linha estao no padrao de C (monta_mme_ar1 poe todo par do mesmo registro). Uma montagem e
// um fator a mais no fim do ajuste, no theta final.
ResiduosMT residuos_ar1(const DesenhoAR& d, const std::vector<double>& theta) {
  ResiduosMT S;
  const std::size_t t = d.t;
  S.t = t;
  S.nlin = d.nlin;
  if (theta.size() != d.modelo.ntheta) return S;
  MontadoAR M = monta_mme_ar1(d, theta);
  if (!M.ok) return S;
  const std::vector<std::size_t> perm = grau_minimo(M.c);
  const Csc pc = permuta_sim(M.c, perm);
  const Simbolica sb = simbolica(pc);
  Csc L;
  if (!cholesky(pc, sb, L)) return S;
  std::vector<double> pb(M.total);
  for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[perm[k]];
  const std::vector<double> px = resolve(L, pb);
  std::vector<double> s(M.total);
  for (std::size_t k = 0; k < M.total; k++) s[perm[k]] = px[k];
  SelInv z = inversa_seletiva(L, 0);
  std::vector<std::size_t> inv_perm(M.total);
  for (std::size_t k = 0; k < M.total; k++) inv_perm[perm[k]] = k;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  S.residuo.assign(d.nlin * t, nan);
  S.ajustado.assign(d.nlin * t, nan);
  S.ajustado_fixo.assign(d.nlin * t, nan);
  S.variancia.assign(d.nlin * t, nan);
  const std::size_t n_fixo = d.x.ncol * t;
  for (std::size_t r = 0; r < d.nlin; r++) {
    if (!d.usa[r]) continue;
    for (std::size_t tau = 0; tau < t; tau++) {
      if (!d.obs.empty() && !d.obs[r * t + tau]) continue;
      double fa = 0.0, fx = 0.0, q = 0.0;
      const auto& lr = d.lw[r];
      for (std::size_t p = 0; p < lr.size(); p++) {
        if (lr[p].trait != tau) continue;
        const double v = lr[p].val * s[lr[p].col];
        fa += v;
        if (lr[p].col < n_fixo) fx += v;
        for (std::size_t p2 = p; p2 < lr.size(); p2++) {
          if (lr[p2].trait != tau) continue;
          double c;
          if (!z.get(inv_perm[lr[p].col], inv_perm[lr[p2].col], c)) { S.fora++; continue; }
          q += (p == p2 ? 1.0 : 2.0) * lr[p].val * lr[p2].val * c;
        }
      }
      S.ajustado[r * t + tau] = fa;
      S.ajustado_fixo[r * t + tau] = fx;
      S.residuo[r * t + tau] = d.y.at(r, tau) - fa;
      S.variancia[r * t + tau] = M.r0.at(tau, tau) - q;
    }
  }
  S.ok = true;
  return S;
}

// Bootstrap parametrico do AR(1)/CAR(1) univariado com theta fixo, preparado uma vez: u*_g
// como no model(), e* pela recorrencia do processo dentro do sujeito, e_1 = sd z,
// e_k = r e_{k-1} + sd sqrt(1 - r^2) z com r = Gamma(t_k - t_{k-1}) (a convencao de sinal de
// gamma_inv), y* = X b^ + Z u* + e*, e o residuo pelo mesmo fator de C em theta^. Unidades
// absolutas. O preparo guarda o desenho, o fator, os elos r e a Gamma^-1 de cada sujeito.
namespace {

class PreparoAR : public PreparoNulos {
 public:
  DesenhoAR d;
  std::vector<std::size_t> perm, usadas;
  Csc L;
  SorteioAleatorio SA;
  std::vector<double> xb;
  std::vector<std::vector<double>> elo;
  std::vector<GammaInv> gis;
  double s2e = 0.0;
  std::size_t q_tot = 0;

  // rhs = W' R^-1 y, R^-1 = Gamma^-1 / s2e por sujeito
  void rhs_de(const std::vector<double>& ys, std::vector<double>& rhs) const {
    std::fill(rhs.begin(), rhs.end(), 0.0);
    for (std::size_t si = 0; si < d.sujeitos.size(); si++) {
      const std::vector<std::size_t>& regs = d.sujeitos[si];
      const GammaInv& gi = gis[si];
      const std::size_t m = regs.size();
      for (std::size_t k = 0; k < m; k++) {
        double v = gi.diag[k] * ys[regs[k]];
        if (k > 0) v += gi.sub[k - 1] * ys[regs[k - 1]];
        if (k + 1 < m) v += gi.sub[k] * ys[regs[k + 1]];
        v /= s2e;
        for (const auto& en : d.lw[regs[k]]) rhs[en.col] += en.val * v;
      }
    }
  }

  NulosSaida lote(std::size_t nsim, bool com_u, double familia) override {
    NulosSaida S;
    const std::size_t nu = usadas.size(), n_fixo = d.x.ncol;
    const double sd_e = std::sqrt(s2e);
    S.n_usadas = nu;
    S.nsim = nsim;
    S.residuo.assign(nu * nsim, 0.0);
    S.fixo.assign(nu * nsim, 0.0);
    if (com_u) S.u.assign((q_tot - n_fixo) * nsim, 0.0);
    std::vector<double> u(q_tot, 0.0), ys(d.nlin, 0.0), rhs(q_tot), pb(q_tot), s(q_tot);
    GetRNGstate();
    try {
      for (std::size_t rep = 0; rep < nsim; rep++) {
        checa_interrupcao();
        sorteia_aleatorios(SA, u);
        for (std::size_t si = 0; si < d.sujeitos.size(); si++) {
          const std::vector<std::size_t>& regs = d.sujeitos[si];
          double ant = 0.0;
          for (std::size_t k = 0; k < regs.size(); k++) {
            const double zk = erro_padrao(familia);
            const double r = elo[si][k];
            const double e = (k == 0) ? sd_e * zk : r * ant + sd_e * std::sqrt(1.0 - r * r) * zk;
            ant = e;
            const std::size_t i = regs[k];
            double v = xb[i];
            for (const auto& en : d.lw[i])
              if (en.col >= n_fixo) v += en.val * u[en.col];
            ys[i] = v + e;
          }
        }
        rhs_de(ys, rhs);
        for (std::size_t k = 0; k < q_tot; k++) pb[k] = rhs[perm[k]];
        const std::vector<double> px = resolve(L, pb);
        for (std::size_t k = 0; k < q_tot; k++) s[perm[k]] = px[k];
        double* er = S.residuo.data() + rep * nu;
        double* fr = S.fixo.data() + rep * nu;
        for (std::size_t k = 0; k < nu; k++) {
          const std::size_t i = usadas[k];
          double fa = 0.0, fx = 0.0;
          for (const auto& en : d.lw[i]) {
            const double v = en.val * s[en.col];
            fa += v;
            if (en.col < n_fixo) fx += v;
          }
          er[k] = ys[i] - fa;
          fr[k] = fx;
        }
        if (com_u)
          for (std::size_t k = n_fixo; k < q_tot; k++)
            S.u[rep * (q_tot - n_fixo) + (k - n_fixo)] = u[k];
      }
    } catch (...) {
      PutRNGstate();
      throw;
    }
    PutRNGstate();
    return S;
  }
};

}  // namespace

std::unique_ptr<PreparoNulos> prepara_nulos_ar1(DesenhoAR&& d0, const std::vector<double>& theta,
                                                const std::vector<double>& b) {
  if (d0.t != 1)
    throw Erro("the parametric bootstrap of a multi-trait AR(1) fit is not available yet: "
               "diagnose each trait with its own model_ar1()");
  std::unique_ptr<PreparoAR> P(new PreparoAR());
  P->d = std::move(d0);
  const DesenhoAR& d = P->d;
  if (b.size() != d.x.ncol)
    throw Erro("b with " + std::to_string(b.size()) + " entries for " +
               std::to_string(d.x.ncol) + " fixed columns");
  double rho = 0.0;
  std::vector<std::size_t> offset_grupo;
  {
    MontadoAR M = monta_mme_ar1(d, theta);
    if (!M.ok) throw Erro("theta INADMISSIBLE (s2e <= 0 or |rho| >= 1)");
    P->perm = grau_minimo(M.c);
    {
      const Csc pc = permuta_sim(M.c, P->perm);
      const Simbolica sb = simbolica(pc);
      if (!cholesky(pc, sb, P->L))
        throw Erro("the mixed model equations are not positive-definite");
    }
    offset_grupo = M.offset_grupo;
    P->s2e = M.r0.at(0, 0);
    rho = M.rho;
    P->q_tot = M.total;
    for (std::size_t i = 0; i < d.nlin; i++)
      if (d.usa[i]) P->usadas.push_back(i);
    // o residuo dos dados no fator de agora (o rhs da montagem e o W' R^-1 y dos dados)
    std::vector<double> pb(M.total), s(M.total);
    for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[P->perm[k]];
    const std::vector<double> px = resolve(P->L, pb);
    for (std::size_t k = 0; k < M.total; k++) s[P->perm[k]] = px[k];
    P->residuo_dados.assign(P->usadas.size(), 0.0);
    for (std::size_t k = 0; k < P->usadas.size(); k++) {
      const std::size_t i = P->usadas[k];
      double fa = 0.0;
      for (const auto& en : d.lw[i]) fa += en.val * s[en.col];
      P->residuo_dados[k] = d.y.at(i, 0) - fa;
    }
  }
  // o sorteio (as K^-1 fatoradas) depois que C saiu de cena
  P->SA = prepara_sorteio(d.modelo, d.kinv, d.aleatorios, theta, offset_grupo);
  P->xb.assign(d.nlin, 0.0);
  for (std::size_t i : P->usadas)
    for (const auto& e : d.lw[i])
      if (e.col < d.x.ncol) P->xb[i] += e.val * b[e.col];
  // r de cada elo e o Gamma^-1 de cada sujeito, fixos em todas as replicas
  P->elo.resize(d.sujeitos.size());
  P->gis.resize(d.sujeitos.size());
  for (std::size_t si = 0; si < d.sujeitos.size(); si++) {
    const std::vector<std::size_t>& regs = d.sujeitos[si];
    std::vector<double> ts;
    for (std::size_t r : regs) ts.push_back(d.tempo[r]);
    P->gis[si] = gamma_inv(ts, rho);
    P->elo[si].assign(regs.size(), 0.0);
    for (std::size_t k = 1; k < regs.size(); k++) {
      const double dt = ts[k] - ts[k - 1];
      P->elo[si][k] = std::pow(std::fabs(rho), dt) *
                      ((rho < 0.0 && std::fmod(dt, 2.0) != 0.0) ? -1.0 : 1.0);
    }
  }
  P->diag_k.resize(d.modelo.grupos.size());
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++)
    if (P->SA.lk[g].ncol > 0) P->diag_k[g] = diagonal_de_k(P->SA.lk[g], P->SA.perm_k[g]);
  P->niveis = niveis_grupos(d);
  P->usa = d.usa;
  return std::unique_ptr<PreparoNulos>(P.release());
}

}  // namespace br
