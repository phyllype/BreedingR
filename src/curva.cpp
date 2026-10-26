// Curva de crescimento nao linear ajustada direto nas pesagens: o motor de gibbs_growth().
//
// ## O modelo
//
// Pesagem j do animal i, na idade t_ij:
//
//   y_ij = f(t_ij; phi_i, nu) + e_ij                   escala crua
//   log y_ij = log f(t_ij; phi_i, nu) + e_ij           escala log (erro multiplicativo)
//   e_ij ~ N(0, s2e_c), c a classe da pesagem
//
//   f(t; phi, nu) = A (1 + nu e)^(-1/nu),  e = exp(-k (t - ti)),  phi = (log A, ti, log k)
//
// nu = 1 logistica, nu -> 0 Gompertz, nu = -1/3 von Bertalanffy, nu livre Richards (dado aqui;
// estimado so na fase B). ti e a idade da inflexao para todo nu: d2f/dt2 = -A k^2 e
// (1 + nu e)^(-1/nu - 2) (1 - e), zero em e = 1. A conta e sempre em log:
//
//   log f = log A - g(nu, e),  g(nu, e) = log1p(nu e) / nu,  g(0, e) = e
//   g pela serie e (1 - nu e / 2 + (nu e)^2 / 3) com |nu e| < 1e-6
//   nu > 0 e nu e > 1e8: log(1 + nu e) = log nu + log e + log1p(1 / (nu e))
//   1 + nu e <= 0 (nu < 0, abaixo da borda do suporte): log f = -inf, verossimilhanca 0
//
//   d log f / d log A = 1,  d log f / d ti = -k e / (1 + nu e),
//   d log f / d log k = k (t - ti) e / (1 + nu e);  na crua, d f = f d log f
//
// O nivel 2 e o modelo multicaracter do model_mt() com os tres parametros como
// caracteristicas latentes, uma linha por animal:
//
//   phi_i = X_i b + a_i + soma_r Z_ri c_r + p_i,  a ~ N(0, G0 (x) K),  c_r ~ N(0, C_r (x) I),
//   p_i ~ N(0, Sigma_p)
//
// com um phi por animal, p_i e o termo permanente e o residuo do nivel 2 ao mesmo tempo: e o
// R0 do desenho multicaracter (o motor o chama de res; o R o renomeia pe).
//
// ## A varredura (Metropolis dentro de Gibbs)
//
// 1. As curvas, animal a animal e em paralelo. Com mu_i = W_i loc, Sigma_p, s2e e nu fixos, o
//    alvo e log pi_i(phi) = soma_j log N(y_ij; f_ij(phi), s2e) - (phi - mu_i)' Sigma_p^-1
//    (phi - mu_i) / 2. (a) O modo pelo Gauss-Newton amortecido a partir de mu_i, NUNCA do phi
//    corrente: H = J'DJ + Sigma_p^-1, passo H^-1 g cortado pela metade ate log pi subir (30
//    vezes no maximo), parada com o decremento g'H^-1 g < 1e-10 ou 25 passos. O modo m_i e
//    H(m_i) sao funcoes so dos outros blocos, e por isso a proposta independente abaixo e uma
//    proposta valida de Metropolis-Hastings. Com nu < 0 e mu_i pondo alguma pesagem fora do
//    suporte, a partida move ti para min_j t_ij + log(-1/(2 nu)) / k, onde 1 + nu e = 1/2 na
//    pesagem mais cedo (regra deterministica, tambem so dos outros blocos). (b) Proposta
//    independente t multivariada com 5 graus de liberdade, centro m_i e escala H^-1, aceita
//    com a razao pi(phi*) q(phi) / (pi(phi) q(phi*)). (c) Passeio aleatorio phi + s_i L z com
//    L L' = H^-1, aceito com pi(phi*) / pi(phi); s_i parte de 2.38 / sqrt(3) e, SO no
//    aquecimento, segue log s_i += (it + 1)^-0.6 (aceito - 0.30), congelado depois (Andrieu e
//    Thoms 2008: a cadeia guardada e homogenea). Os numeros aleatorios de cada animal (3
//    normais, 1 qui-quadrado, 1 uniforme, 3 normais, 1 uniforme) saem do gerador do R em
//    serie, na ordem dos animais, para um buffer antes do laco paralelo: a cadeia e a mesma
//    com qualquer numero de threads. As somas por classe que o passo 4 le saem por pesagem e
//    sao somadas em serie, na ordem fixa.
// 3. As localizacoes do nivel 2 num bloco: d.y recebe phi, as MME absolutas do multicaracter
//    em theta = (grupos, Sigma_p) (monta_mme_mt), a Cholesky esparsa com a simbolica em cache
//    e loc = solucao + P L^-T z, como no gibbs() (sem o fator s2e: o sistema e absoluto).
//    prior_fixed soma 1/dp_tau^2 na diagonal de cada equacao fixa do parametro tau (nada no
//    parametro de priori plana).
// 4. As covariancias, conjugadas: Sigma_p | phi, loc ~ InvWishart(n - 4, E'E) na plana, cada
//    grupo C_g | U ~ InvWishart(nl - 4, U'K^-1 U) e cada s2e_c ~ SSE_c / chisq(n_c - 2); a
//    propria soma gl0 diag(escala) e gl0 graus de liberdade, a uniforme no desvio usa
//    chisq(n_c - 1). rinvwishart() e gl_priori() sao as do gibbs().
// 5. Nas iteracoes guardadas: os componentes, as medias e variancias (Welford) de phi, das
//    localizacoes e dos derivados, os acumuladores do WAIC e os residuos pivotais.
//
// As opcoes internas dos portoes chegam em OpcoesCurva: `variante` 2 tira a razao q da
// proposta independente (a lei estacionaria vira pi q) e 3 tira a priori do nivel 2 da
// aceitacao (o alvo vira a verossimilhanca so); `movimentos` roda um movimento so;
// `segura` 1 prende phi e os componentes e pula os passos 1 e 4 (portao E3), 2 prende so phi
// (portao E8); CurvaConf::linear troca f pelo polinomio phi1 + phi2 t + phi3 t^2 (portao E4,
// contra o gibbs() na regressao aleatoria equivalente).

#include "mme.h"
#include <algorithm>
#include <R_ext/Print.h>
#include <R_ext/Utils.h>

#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>

namespace br {

static const double MENOS_INF = -std::numeric_limits<double>::infinity();

static std::size_t vech3(std::size_t i, std::size_t j) {
  if (i < j) std::swap(i, j);
  return (j * (2 * 3 - j + 1)) / 2 + (i - j);
}

// ---- algebra 3 x 3 (por linhas)

// Cholesky inferior; false se nao positiva-definida
static bool chol3(const double* a, double* l) {
  const double d0 = a[0];
  if (!(d0 > 0.0) || !std::isfinite(d0)) return false;
  const double l00 = std::sqrt(d0);
  const double l10 = a[3] / l00, l20 = a[6] / l00;
  const double d1 = a[4] - l10 * l10;
  if (!(d1 > 0.0)) return false;
  const double l11 = std::sqrt(d1);
  const double l21 = (a[7] - l20 * l10) / l11;
  const double d2 = a[8] - l20 * l20 - l21 * l21;
  if (!(d2 > 0.0)) return false;
  const double l22 = std::sqrt(d2);
  l[0] = l00; l[1] = 0.0; l[2] = 0.0;
  l[3] = l10; l[4] = l11; l[5] = 0.0;
  l[6] = l20; l[7] = l21; l[8] = l22;
  return true;
}

// x = (L L')^-1 b
static void resolve_chol3(const double* l, const double* b, double* x) {
  const double w0 = b[0] / l[0];
  const double w1 = (b[1] - l[3] * w0) / l[4];
  const double w2 = (b[2] - l[6] * w0 - l[7] * w1) / l[8];
  x[2] = w2 / l[8];
  x[1] = (w1 - l[7] * x[2]) / l[4];
  x[0] = (w0 - l[3] * x[1] - l[6] * x[2]) / l[0];
}

// x = L^-T z (covariancia (L L')^-1 quando z ~ N(0, I))
static void resolve_lt3(const double* l, const double* z, double* x) {
  x[2] = z[2] / l[8];
  x[1] = (z[1] - l[7] * x[2]) / l[4];
  x[0] = (z[0] - l[3] * x[1] - l[6] * x[2]) / l[0];
}

// ||L' v||^2 = v' (L L') v
static double forma_lt3(const double* l, const double* v) {
  const double a0 = l[0] * v[0] + l[3] * v[1] + l[6] * v[2];
  const double a1 = l[4] * v[1] + l[7] * v[2];
  const double a2 = l[8] * v[2];
  return a0 * a0 + a1 * a1 + a2 * a2;
}

// inversa de uma 3 x 3 positiva-definida pela Cholesky
static bool inv3(const double* a, double* out) {
  double l[9];
  if (!chol3(a, l)) return false;
  for (int c = 0; c < 3; c++) {
    double e[3] = {0.0, 0.0, 0.0}, x[3];
    e[c] = 1.0;
    resolve_chol3(l, e, x);
    for (int r = 0; r < 3; r++) out[r * 3 + c] = x[r];
  }
  for (int r = 0; r < 3; r++)
    for (int c = 0; c < r; c++) {
      const double m = 0.5 * (out[r * 3 + c] + out[c * 3 + r]);
      out[r * 3 + c] = m;
      out[c * 3 + r] = m;
    }
  return true;
}

static double quad3(const double* p, const double* v) {
  double s = 0.0;
  for (int r = 0; r < 3; r++) {
    double q = 0.0;
    for (int c = 0; c < 3; c++) q += p[r * 3 + c] * v[c];
    s += v[r] * q;
  }
  return s;
}

// a covariancia 3 x 3 guardada em vech (por colunas) a partir de off
static void cov3_de(const std::vector<double>& theta, std::size_t off, double* c) {
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 3; j++) c[i * 3 + j] = theta[off + vech3(i, j)];
}

// ---- a familia de curvas

// log f(t; phi, nu) e, com dl, as derivadas em phi. false fora do suporte ou com
// transbordamento (lf = -inf)
static bool log_curva(double t, const double* phi, double nu, double& lf, double* dl) {
  const double k = std::exp(phi[2]);
  const double dt = t - phi[1];
  const double le = -k * dt;   // log e, guardado como esta
  double g, r;                 // g(nu, e) e r = e / (1 + nu e)
  if (nu == 0.0) {
    const double e = std::exp(le);
    g = e;
    r = e;
  } else {
    const double e = std::exp(le);
    const double ne = nu * e;
    if (std::fabs(ne) < 1e-6) {
      g = e * (1.0 - ne / 2.0 + ne * ne / 3.0);
      r = e / (1.0 + ne);
    } else if (nu > 0.0 && ne > 1e8) {
      const double lne = std::log(nu) + le;
      g = (lne + std::log1p(std::exp(-lne))) / nu;
      r = 1.0 / (nu + std::exp(-le));
    } else {
      if (!(1.0 + ne > 0.0)) {
        lf = MENOS_INF;
        return false;
      }
      g = std::log1p(ne) / nu;
      r = e / (1.0 + ne);
    }
  }
  if (!std::isfinite(g) || !std::isfinite(r)) {
    lf = MENOS_INF;
    return false;
  }
  lf = phi[0] - g;
  if (dl) {
    dl[0] = 1.0;
    dl[1] = -k * r;
    dl[2] = k * dt * r;
  }
  return std::isfinite(lf);
}

bool curva_analise(const CurvaConf& c, double t, const double* phi, double& f, double* jac) {
  if (c.linear) {
    f = phi[0] + phi[1] * t + phi[2] * t * t;
    if (jac) {
      jac[0] = 1.0;
      jac[1] = t;
      jac[2] = t * t;
    }
    return std::isfinite(f);
  }
  double lf;
  double dl[3];
  if (!log_curva(t, phi, c.nu, lf, jac ? dl : nullptr)) {
    f = MENOS_INF;
    return false;
  }
  if (c.log_escala) {
    f = lf;
    if (jac)
      for (int q = 0; q < 3; q++) jac[q] = dl[q];
  } else {
    f = std::exp(lf);
    if (jac)
      for (int q = 0; q < 3; q++) jac[q] = f * dl[q];
  }
  return true;
}

double curva_peso(const CurvaConf& c, double t, const double* phi) {
  if (c.linear) return phi[0] + phi[1] * t + phi[2] * t * t;
  double lf;
  if (!log_curva(t, phi, c.nu, lf, nullptr)) return 0.0;
  return std::exp(lf);
}

// c_A(nu) = log1p(nu) / nu, 1 em nu = 0 (serie perto de zero). O peso na inflexao e
// W_inf = A (1 + nu)^(-1/nu), log W_inf = log A - c_A(nu); a taxa maxima, na inflexao, e
// R_max = A k (1 + nu)^(-1/nu - 1), log R_max = log A + log k - (1 + nu) c_A(nu)
static double c_inflexao(double nu) {
  if (std::fabs(nu) < 1e-8) return 1.0 - nu / 2.0 + nu * nu / 3.0;
  return std::log1p(nu) / nu;
}
static double peso_inflexao(const CurvaConf& c, const double* phi) {
  if (c.linear) return std::nan("");
  return std::exp(phi[0] - c_inflexao(c.nu));
}
static double taxa_maxima(const CurvaConf& c, const double* phi) {
  if (c.linear) return std::nan("");
  return std::exp(phi[0] + phi[2] - (1.0 + c.nu) * c_inflexao(c.nu));
}

// ---- as pesagens por animal

PesagensCurva agrupa_pesagens(const std::vector<int>& animal, const std::vector<double>& t,
                              const std::vector<double>& y, const std::vector<double>& logpeso,
                              const std::vector<int>& classe, std::size_t n_anim,
                              std::size_t n_classes) {
  const std::size_t n = animal.size();
  if (t.size() != n || y.size() != n || classe.size() != n || logpeso.size() != n)
    throw Erro("weighings: age, weight and class with different lengths");
  if (n_classes == 0) throw Erro("weighings: no residual class");
  PesagensCurva p;
  p.n_anim = n_anim;
  p.n_classes = n_classes;
  p.ini.assign(n_anim + 1, 0);
  for (std::size_t k = 0; k < n; k++) {
    if (animal[k] < 0 || static_cast<std::size_t>(animal[k]) >= n_anim)
      throw Erro("weighings: an animal index outside the table");
    if (classe[k] < 0 || static_cast<std::size_t>(classe[k]) >= n_classes)
      throw Erro("weighings: a residual class outside the classes");
    if (!std::isfinite(t[k]) || !std::isfinite(y[k]))
      throw Erro("weighings: a non-finite age or weight reached the engine");
    p.ini[static_cast<std::size_t>(animal[k]) + 1]++;
  }
  for (std::size_t i = 0; i < n_anim; i++) p.ini[i + 1] += p.ini[i];
  std::vector<std::size_t> pos(p.ini.begin(), p.ini.end() - 1);
  p.ordem.assign(n, 0);
  p.t.assign(n, 0.0);
  p.y.assign(n, 0.0);
  p.logpeso.assign(n, 0.0);
  p.classe.assign(n, 0);
  for (std::size_t k = 0; k < n; k++) {
    const std::size_t q = pos[static_cast<std::size_t>(animal[k])]++;
    p.ordem[q] = k;
    p.t[q] = t[k];
    p.y[q] = y[k];
    p.logpeso[q] = logpeso[k];
    p.classe[q] = classe[k];
  }
  for (std::size_t i = 0; i < n_anim; i++)
    if (p.ini[i + 1] == p.ini[i])
      throw Erro("weighings: an animal of the level-2 table has no weighing");
  return p;
}

// ---- o alvo de um animal

struct AlvoAnimal {
  const CurvaConf* c = nullptr;
  const double* t = nullptr;
  const double* y = nullptr;
  const int* cl = nullptr;
  std::size_t n = 0;
  const double* mu = nullptr;     // 3
  const double* prec = nullptr;   // 3 x 3, P = Sigma_p^-1
  const double* uprec = nullptr;  // 3 x 3 triangular superior, P = U'U
  const double* s2e = nullptr;    // por classe
};

static AlvoAnimal alvo_de(const CurvaConf& c, const PesagensCurva& p, std::size_t i,
                          const double* mu, const double* prec, const double* uprec,
                          const double* s2e) {
  AlvoAnimal a;
  a.c = &c;
  const std::size_t ini = p.ini[i];
  a.t = p.t.data() + ini;
  a.y = p.y.data() + ini;
  a.cl = p.classe.data() + ini;
  a.n = p.ini[i + 1] - ini;
  a.mu = mu;
  a.prec = prec;
  a.uprec = uprec;
  a.s2e = s2e;
  return a;
}

// o termo da priori, -(phi - mu)' P (phi - mu) / 2, a mesma aritmetica em todo lugar
static double termo_priori(const AlvoAnimal& a, const double* phi) {
  const double d[3] = {phi[0] - a.mu[0], phi[1] - a.mu[1], phi[2] - a.mu[2]};
  return -0.5 * quad3(a.prec, d);
}

static double log_alvo(const AlvoAnimal& a, const double* phi, bool com_priori) {
  double s = 0.0;
  for (std::size_t j = 0; j < a.n; j++) {
    double f;
    if (!curva_analise(*a.c, a.t[j], phi, f, nullptr)) return MENOS_INF;
    const double r = a.y[j] - f;
    s -= 0.5 * r * r / a.s2e[a.cl[j]];
  }
  if (!std::isfinite(s)) return MENOS_INF;
  if (com_priori) s += termo_priori(a, phi);
  return s;
}

// Junta a linha (a, beta) ao fator R (3 x 3 triangular superior, por linhas) e a c por
// rotacoes de Givens: a QR sequencial de [sqrt(D) J; U] sem formar J'DJ. Pela equacao normal,
// com uma pesagem de residuo pequeno e a curva longe dela (e = exp(-k (t - ti)) grande), J'DJ
// tem condicao da ordem de 1e19 e a Cholesky de J'DJ + P perde a positividade (medido no caso
// bimodal do portao E2); a QR e estavel em relacao a norma de A, e os valores singulares
// pequenos (os da priori) saem com erro relativo da ordem de 1e-7 ali
static void givens_linha(double* r, double* c, double* a, double beta) {
  for (int k = 0; k < 3; k++) {
    if (a[k] == 0.0) continue;
    const double rk = r[k * 3 + k];
    const double h = std::hypot(rk, a[k]);
    const double cs = rk / h, sn = a[k] / h;
    r[k * 3 + k] = h;
    for (int j = k + 1; j < 3; j++) {
      const double t1 = cs * r[k * 3 + j] + sn * a[j];
      a[j] = -sn * r[k * 3 + j] + cs * a[j];
      r[k * 3 + j] = t1;
    }
    const double t2 = cs * c[k] + sn * beta;
    beta = -sn * c[k] + cs * beta;
    c[k] = t2;
    a[k] = 0.0;
  }
}

// Num ponto: log pi, o fator L (inferior) de H = J'DJ + P = L L' e c com o gradiente g = L c.
// Assim o passo de Gauss-Newton e H^-1 g = L^-T c e o decremento g'H^-1 g e c'c. P = U'U, com
// U = chol(P)' na partida da QR e c = -U (phi - mu)
static bool gn_avalia(const AlvoAnimal& a, const double* phi, double& lp, double* l, double* c) {
  double r[9];
  for (int q = 0; q < 9; q++) r[q] = a.uprec[q];
  const double d[3] = {phi[0] - a.mu[0], phi[1] - a.mu[1], phi[2] - a.mu[2]};
  for (int k = 0; k < 3; k++)
    c[k] = -(a.uprec[k * 3] * d[0] + a.uprec[k * 3 + 1] * d[1] + a.uprec[k * 3 + 2] * d[2]);
  double s = 0.0;
  for (std::size_t j = 0; j < a.n; j++) {
    double f, jac[3];
    if (!curva_analise(*a.c, a.t[j], phi, f, jac)) {
      lp = MENOS_INF;
      return false;
    }
    const double res = a.y[j] - f;
    s -= 0.5 * res * res / a.s2e[a.cl[j]];
    const double sw = 1.0 / std::sqrt(a.s2e[a.cl[j]]);
    double linha[3] = {sw * jac[0], sw * jac[1], sw * jac[2]};
    givens_linha(r, c, linha, sw * res);
  }
  if (!std::isfinite(s)) {
    lp = MENOS_INF;
    return false;
  }
  lp = s + termo_priori(a, phi);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) l[i * 3 + j] = j <= i ? r[j * 3 + i] : 0.0;
  for (int q = 0; q < 9; q++)
    if (!std::isfinite(l[q])) return false;
  for (int k = 0; k < 3; k++)
    if (!(l[k * 3 + k] > 0.0) || !std::isfinite(c[k])) return false;
  return std::isfinite(lp);
}

// A partida do Gauss-Newton: mu, com ti movido (so com nu < 0) quando mu poe alguma pesagem
// fora do suporte
static void partida_modo(const AlvoAnimal& a, double* phi) {
  for (int q = 0; q < 3; q++) phi[q] = a.mu[q];
  const CurvaConf& c = *a.c;
  if (c.linear || !(c.nu < 0.0) || a.n == 0) return;
  double tmin = a.t[0];
  bool fora = false;
  for (std::size_t j = 0; j < a.n; j++) {
    tmin = std::min(tmin, a.t[j]);
    double lf;
    if (!log_curva(a.t[j], phi, c.nu, lf, nullptr)) fora = true;
  }
  if (fora) phi[1] = tmin + std::log(-1.0 / (2.0 * c.nu)) / std::exp(phi[2]);
}

// O modo pelo Gauss-Newton amortecido e o fator L de H nele (H = L L'). false quando a
// partida ja nao avalia: a proposta fica entao a priori do nivel 2, centro na partida e
// precisao P (L = chol(P)), ainda uma funcao so dos outros blocos
static bool modo_animal(const AlvoAnimal& a, double* m, double* l, int& passos) {
  double phi[3];
  partida_modo(a, phi);
  double lp, c[3], lc[9];
  passos = 0;
  if (!gn_avalia(a, phi, lp, lc, c)) {
    for (int q = 0; q < 3; q++) m[q] = phi[q];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) l[i * 3 + j] = j <= i ? a.uprec[j * 3 + i] : 0.0;
    return false;
  }
  for (;; passos++) {
    const double dec = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    if (!(dec >= 1e-10) || passos >= 25) break;
    double delta[3];
    resolve_lt3(lc, c, delta);
    double passo = 1.0, cand[3], lpc = MENOS_INF;
    bool subiu = false;
    for (int metade = 0; metade <= 30; metade++) {
      for (int q = 0; q < 3; q++) cand[q] = phi[q] + passo * delta[q];
      lpc = log_alvo(a, cand, true);
      if (lpc > lp) {
        subiu = true;
        break;
      }
      passo *= 0.5;
    }
    if (!subiu) break;
    double cn[3], ln[9], lpn;
    if (!gn_avalia(a, cand, lpn, ln, cn)) break;
    for (int q = 0; q < 3; q++) phi[q] = cand[q];
    lp = lpn;
    for (int q = 0; q < 3; q++) c[q] = cn[q];
    for (int q = 0; q < 9; q++) lc[q] = ln[q];
  }
  for (int q = 0; q < 3; q++) m[q] = phi[q];
  for (int q = 0; q < 9; q++) l[q] = lc[q];
  return true;
}

// o estado de um animal na cadeia
struct EstadoAnimal {
  double phi[3] = {0.0, 0.0, 0.0};
  double log_s = 0.0;
  double tentativas = 0.0, ac_ind = 0.0, ac_rw = 0.0;
  double falhas = 0.0;   // varreduras em que o modo nao avaliou (a proposta foi a priori)
};

static const double LOG_S0 = std::log(2.38 / std::sqrt(3.0));

// Os dois movimentos do passo 1 de um animal, com o modo m e o fator L de H nele (L L' = H)
// ja calculados e os 9 numeros aleatorios dele em `u`: z(3), w ~ chisq(5), uniforme, z'(3),
// uniforme. `it` e a varredura (a adaptacao usa it < burnin)
static void movimentos_animal(const AlvoAnimal& a, EstadoAnimal& e, const double* u,
                              int variante, int movimentos, std::size_t it, std::size_t burnin,
                              const double* m, const double* l) {
  const bool com_priori = variante != 3;
  double lp = log_alvo(a, e.phi, com_priori);
  const bool conta = it >= burnin;
  if (conta) e.tentativas += 1.0;
  // (b) independente: t com 5 graus de liberdade, centro m e escala H^-1
  if (movimentos != 1) {
    const double esc = std::sqrt(5.0 / u[3]);
    double zz[3] = {u[0] * esc, u[1] * esc, u[2] * esc}, x[3], cand[3];
    resolve_lt3(l, zz, x);
    for (int q = 0; q < 3; q++) cand[q] = m[q] + x[q];
    const double lpc = log_alvo(a, cand, com_priori);
    double razao = lpc - lp;
    if (variante != 2) {
      const double qc = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) * 5.0 / u[3];
      const double dv[3] = {e.phi[0] - m[0], e.phi[1] - m[1], e.phi[2] - m[2]};
      const double qa = forma_lt3(l, dv);
      // log q(atual) - log q(proposta), q ~ (1 + Q/5)^(-(5 + 3)/2)
      razao += -4.0 * std::log1p(qa / 5.0) + 4.0 * std::log1p(qc / 5.0);
    }
    if (std::isfinite(lpc) && std::log(u[4]) < razao) {
      for (int q = 0; q < 3; q++) e.phi[q] = cand[q];
      lp = lpc;
      if (conta) e.ac_ind += 1.0;
    }
  }
  // (c) passeio aleatorio com a forma da Laplace, escala adaptada so no aquecimento
  if (movimentos != 2) {
    const double s = std::exp(e.log_s);
    double zz[3] = {u[5] * s, u[6] * s, u[7] * s}, x[3], cand[3];
    resolve_lt3(l, zz, x);
    for (int q = 0; q < 3; q++) cand[q] = e.phi[q] + x[q];
    const double lpc = log_alvo(a, cand, com_priori);
    const bool aceito = std::isfinite(lpc) && std::log(u[8]) < lpc - lp;
    if (aceito) {
      for (int q = 0; q < 3; q++) e.phi[q] = cand[q];
      if (conta) e.ac_rw += 1.0;
    }
    if (it < burnin)
      e.log_s += std::pow(static_cast<double>(it) + 1.0, -0.6) * ((aceito ? 1.0 : 0.0) - 0.30);
  }
}

// O passo 1 de um animal: o modo e H a partir de mu (so dos outros blocos) e os dois
// movimentos
static void passo_animal(const AlvoAnimal& a, EstadoAnimal& e, const double* u, int variante,
                         int movimentos, std::size_t it, std::size_t burnin) {
  double m[3], l[9];
  int passos = 0;
  if (!modo_animal(a, m, l, passos)) e.falhas += 1.0;
  movimentos_animal(a, e, u, variante, movimentos, it, burnin, m, l);
}

// os 9 numeros de um animal, do gerador do R (so fora de regiao paralela)
static void sorteia_animal(double* u) {
  u[0] = Rf_rnorm(0.0, 1.0);
  u[1] = Rf_rnorm(0.0, 1.0);
  u[2] = Rf_rnorm(0.0, 1.0);
  u[3] = Rf_rchisq(5.0);
  u[4] = Rf_runif(0.0, 1.0);
  u[5] = Rf_rnorm(0.0, 1.0);
  u[6] = Rf_rnorm(0.0, 1.0);
  u[7] = Rf_rnorm(0.0, 1.0);
  u[8] = Rf_runif(0.0, 1.0);
}

// P = Sigma^-1 e U = chol(P)' (superior, P = U'U) de uma covariancia 3 x 3
static bool precisao3(const double* sigma, double* prec, double* uprec) {
  if (!inv3(sigma, prec)) return false;
  double l[9];
  if (!chol3(prec, l)) return false;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) uprec[i * 3 + j] = j >= i ? l[j * 3 + i] : 0.0;
  return true;
}

static void densa3(const Densa& s, double* out) {
  if (s.nlin != 3 || s.ncol != 3) throw Erro("the level-2 covariance must be 3 x 3");
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 3; j++) out[i * 3 + j] = s.at(i, j);
}

ModosCurva modos_curva(const CurvaConf& c, const PesagensCurva& p, const std::vector<double>& mu,
                       const Densa& sigma, const std::vector<double>& s2e) {
  const std::size_t n = p.n_anim;
  if (mu.size() != 3 && mu.size() != 3 * n) throw Erro("mu: one row, or one row per animal");
  if (s2e.size() != p.n_classes) throw Erro("s2e: one value per residual class");
  for (double v : s2e)
    if (!(v > 0.0) || !std::isfinite(v)) throw Erro("s2e must be positive and finite");
  double sg[9], prec[9], uprec[9];
  densa3(sigma, sg);
  if (!precisao3(sg, prec, uprec)) throw Erro("the level-2 covariance is not positive-definite");
  ModosCurva M;
  M.m.assign(3 * n, 0.0);
  M.h.assign(9 * n, 0.0);
  M.passos.assign(n, 0);
  M.ok.assign(n, 0);
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads())
  for (std::size_t i = 0; i < n; i++) {
    const double* mi = mu.size() == 3 ? mu.data() : mu.data() + 3 * i;
    const AlvoAnimal a = alvo_de(c, p, i, mi, prec, uprec, s2e.data());
    int passos = 0;
    double l[9];
    const bool ok = modo_animal(a, M.m.data() + 3 * i, l, passos);
    for (int u = 0; u < 3; u++)
      for (int v = 0; v < 3; v++)
        M.h[9 * i + u * 3 + v] = l[u * 3] * l[v * 3] + l[u * 3 + 1] * l[v * 3 + 1] +
                                 l[u * 3 + 2] * l[v * 3 + 2];
    M.passos[i] = passos;
    M.ok[i] = ok && std::isfinite(log_alvo(a, M.m.data() + 3 * i, true)) ? 1 : 0;
  }
  return M;
}

// O gerador do R aberto num bloco: GetRNGstate() na construcao e PutRNGstate() em devolve()
// ou, quando uma excecao (uma interrupcao, um erro do motor) sai antes, no destrutor. Assim o
// .Random.seed guarda o quanto a cadeia andou mesmo no erro, como em pegs.cpp; no caminho
// normal as chamadas ficam nos mesmos pontos de antes e nenhum sorteio muda
class GeradorR {
 public:
  GeradorR() { GetRNGstate(); }
  ~GeradorR() { devolve(); }
  GeradorR(const GeradorR&) = delete;
  GeradorR& operator=(const GeradorR&) = delete;
  void devolve() {
    if (aberto_) {
      aberto_ = false;
      PutRNGstate();
    }
  }
 private:
  bool aberto_ = true;
};

PassoCurvaSaida passo_curva_isolado(const CurvaConf& c, const PesagensCurva& p,
                                    const double* mu, const Densa& sigma,
                                    const std::vector<double>& s2e, const double* phi0,
                                    std::size_t burnin, std::size_t n, int variante,
                                    int movimentos) {
  if (p.n_anim != 1) throw Erro("the single-animal step takes one animal");
  if (s2e.size() != p.n_classes) throw Erro("s2e: one value per residual class");
  if (variante < 1 || variante > 3) throw Erro("unknown variant");
  if (movimentos < 0 || movimentos > 2) throw Erro("unknown move set");
  double sg[9], prec[9], uprec[9];
  densa3(sigma, sg);
  if (!precisao3(sg, prec, uprec)) throw Erro("the level-2 covariance is not positive-definite");
  PassoCurvaSaida S;
  const AlvoAnimal a = alvo_de(c, p, 0, mu, prec, uprec, s2e.data());
  // mu, Sigma, s2e e nu nao mudam aqui, entao o modo e H (funcoes deles so) sao os mesmos em
  // toda varredura: calculados uma vez, a cadeia e a que passo_animal() daria recalculando
  double lm[9];
  S.ok = modo_animal(a, S.m, lm, S.passos);
  for (int u = 0; u < 3; u++)
    for (int v = 0; v < 3; v++)
      S.h[u * 3 + v] = lm[u * 3] * lm[v * 3] + lm[u * 3 + 1] * lm[v * 3 + 1] +
                       lm[u * 3 + 2] * lm[v * 3 + 2];
  EstadoAnimal e;
  e.log_s = LOG_S0;
  for (int q = 0; q < 3; q++) e.phi[q] = phi0 ? phi0[q] : S.m[q];
  if (!std::isfinite(log_alvo(a, e.phi, variante != 3)))
    throw Erro("the starting curve is outside the support of the target");
  S.cadeia.reserve(3 * n);
  double u[9];
  GeradorR gerador;
  for (std::size_t it = 0; it < burnin + n; it++) {
    if (it % 4096 == 0) checa_interrupcao();
    sorteia_animal(u);
    movimentos_animal(a, e, u, variante, movimentos, it, burnin, S.m, lm);
    if (it >= burnin)
      for (int q = 0; q < 3; q++) S.cadeia.push_back(e.phi[q]);
  }
  gerador.devolve();
  S.ac_ind = e.tentativas > 0.0 ? e.ac_ind / e.tentativas : std::nan("");
  S.ac_rw = e.tentativas > 0.0 ? e.ac_rw / e.tentativas : std::nan("");
  S.escala = std::exp(e.log_s);
  return S;
}

// ---- a cadeia

// uma entrada de uma linha de W do nivel 2: coluna global, caracteristica e valor
struct EntW {
  std::uint32_t col;
  std::uint32_t tr;
  double v;
};

// A folga contra uma escala singular no comeco da cadeia, como no gibbs(), mas RELATIVA A CADA
// DIAGONAL: 1e-10 vezes a propria entrada. A do gibbs() e 1e-10 vezes a media da diagonal, que
// aqui mistura as unidades de log A, ti e log k (o ti em tempo ao quadrado), e com ela a
// cadeia deixava de ser equivariante a unidade do tempo (portao E5: 2.5e-3 de erro relativo
// num componente depois de 300 varreduras). Uma diagonal nula fica com a media das outras
static void folga_diagonal(Densa& s) {
  double tr = 0.0;
  for (std::size_t a = 0; a < s.nlin; a++) tr += s.at(a, a);
  for (std::size_t a = 0; a < s.nlin; a++) {
    const double v = s.at(a, a);
    s.at(a, a) += 1e-10 * (v > 0.0 ? v : (tr > 0.0 ? tr / static_cast<double>(s.nlin) : 1.0));
  }
}

// Welford (1962) de um vetor
static void welford(std::vector<double>& m, std::vector<double>& m2, const double* x,
                    std::size_t n, double na) {
  for (std::size_t k = 0; k < n; k++) {
    const double delta = x[k] - m[k];
    m[k] += delta / na;
    m2[k] += delta * (x[k] - m[k]);
  }
}

SaidaCurva gibbs_curva(DesenhoMT& d, const PesagensCurva& p, const OpcoesCurva& o) {
  SaidaCurva S;
  const std::size_t T = 3;
  const std::size_t n = p.n_anim;
  const CurvaConf& cc = o.curva;
  if (d.t != T) throw Erro("the level-2 design must have the three curve parameters");
  if (d.nlin != n) throw Erro("the level-2 table must have one row per animal with weighings");
  for (std::size_t r = 0; r < n; r++)
    if (!d.usa[r])
      throw Erro("an animal with weighings does not enter the level-2 model (row " +
                 std::to_string(r + 1) + " of the per-animal table)");
  const Modelo& mo = d.modelo;
  const std::size_t ng = mo.grupos.size();
  const std::size_t ntheta = mo.ntheta;
  const std::size_t nc = p.n_classes;
  if (o.priori_grupos.size() != ng) throw Erro("one prior per level-2 group");
  if (o.theta0.size() != ntheta)
    throw Erro("start: " + std::to_string(o.theta0.size()) + " level-2 components, the "
               "design has " + std::to_string(ntheta));
  if (o.s2e0.size() != nc) throw Erro("start: one residual variance per class");
  if (o.phi0.size() != 3 * n) throw Erro("start: one curve per animal");
  if (o.thin == 0) throw Erro("thin must be positive");
  if (o.variante < 1 || o.variante > 4) throw Erro("unknown variant");
  if (o.movimentos < 0 || o.movimentos > 2) throw Erro("unknown move set");
  // os grupos: um termo de um coeficiente cada, entao dim = 3
  std::vector<std::size_t> nl_g(ng, 0);
  std::vector<const DesenhoTermo*> termo_g(ng, nullptr);
  for (std::size_t g = 0; g < ng; g++) {
    const Grupo& gr = mo.grupos[g];
    if (gr.dim != T || gr.termos.size() != 1)
      throw Erro("level-2 group '" + gr.nome + "': one term with one coefficient per group");
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == gr.termos[0]) { nl_g[g] = a.n_niveis; termo_g[g] = &a; break; }
    if (!termo_g[g] || termo_g[g]->n_coef != 1)
      throw Erro("level-2 group '" + gr.nome + "' without its incidence");
    if (gr.estrutura == Estrutura::Declarada && d.kinv[g].ncol != nl_g[g])
      throw Erro("a declared K with a number of rows other than the levels of its term");
  }
  if (o.grupo_genetico >= static_cast<int>(ng)) throw Erro("unknown genetic group");
  std::vector<double> theta = o.theta0;
  std::vector<double> s2e = o.s2e0;
  std::vector<double> phi = o.phi0;
  for (std::size_t g = 0; g < ng; g++) {
    double cg[9], l[9];
    cov3_de(theta, mo.grupos[g].offset, cg);
    if (!chol3(cg, l))
      throw Erro("start: the covariance of group '" + mo.grupos[g].nome +
                 "' is not positive-definite");
  }
  {
    double sp[9], l[9];
    cov3_de(theta, mo.offset_residual, sp);
    if (!chol3(sp, l)) throw Erro("start: the permanent covariance is not positive-definite");
  }
  for (double v : s2e)
    if (!(v > 0.0) || !std::isfinite(v)) throw Erro("start: the residual variance must be positive");
  // o phi de partida tem de estar no suporte de toda pesagem do animal
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = p.ini[i]; j < p.ini[i + 1]; j++) {
      double f;
      if (!curva_analise(cc, p.t[j], phi.data() + 3 * i, f, nullptr))
        throw Erro("start: the curve of animal row " + std::to_string(i + 1) +
                   " leaves a weighing outside its support");
    }
  const std::size_t nw = p.t.size();
  std::vector<std::size_t> n_classe(nc, 0);
  for (int c : p.classe) n_classe[static_cast<std::size_t>(c)]++;

  // as prioris no formato do gibbs() para gl_priori
  auto pg_de = [](const PrioriCurva& pc) {
    PrioriGibbs q;
    q.tipo = pc.tipo;
    q.gl0 = pc.gl0;
    return q;
  };
  for (std::size_t g = 0; g < ng; g++)
    if (!(gl_priori(pg_de(o.priori_grupos[g]), nl_g[g], T) > static_cast<double>(T) - 1.0))
      throw Erro("too few levels in group '" + mo.grupos[g].nome + "' for its prior");
  if (!(gl_priori(pg_de(o.priori_pe), n, T) > static_cast<double>(T) - 1.0))
    throw Erro("too few animals with weighings for the prior of the permanent term");
  for (std::size_t c = 0; c < nc; c++)
    if (!(gl_priori(pg_de(o.priori_res), n_classe[c], 1) > 0.0))
      throw Erro("too few weighings in a residual class for its prior");

  // ---- o estado inicial das localizacoes: a solucao das MME com y = phi0 (partida S5)
  auto poe_y = [&]() {
    for (std::size_t r = 0; r < n; r++)
      for (std::size_t tau = 0; tau < T; tau++) d.y.at(r, tau) = phi[r * T + tau];
  };
  CacheSimbolica cs;
  std::vector<std::size_t> pos_diag_fixos;
  std::vector<double> prec_fixos;      // 1/dp^2 da equacao fixa k (0 no parametro plano)
  if (o.dp_fixos.size() != T) throw Erro("prior_fixed: one standard deviation per parameter");
  const bool com_priori_fixa = o.dp_fixos[0] > 0.0 || o.dp_fixos[1] > 0.0 || o.dp_fixos[2] > 0.0;
  std::vector<std::vector<EntW>> lw;   // as linhas de W, montadas na primeira montagem
  std::vector<double> loc, mu(3 * n, 0.0), xb(3 * n, 0.0);
  std::size_t n_fixo = 0, total = 0;
  std::vector<std::size_t> off_g;
  // monta, fatora e devolve a solucao; ruido = true soma P L^-T z
  auto nivel2 = [&](bool ruido, std::vector<double>* sol_out) -> bool {
    poe_y();
    MontadoMT M = monta_mme_mt(d, theta);
    if (!M.ok) return false;
    if (total == 0) {
      total = M.total;
      n_fixo = M.n_fixo;
      off_g = M.offset_grupo;
    }
    if (com_priori_fixa) {
      if (pos_diag_fixos.size() != M.n_fixo) {
        pos_diag_fixos.clear();
        prec_fixos.assign(M.n_fixo, -1.0);
        // o parametro de cada equacao fixa sai de eq_fixa[j*T + tau]
        for (std::size_t j = 0; j < d.x.ncol; j++)
          for (std::size_t tau = 0; tau < T; tau++) {
            const int e = d.eq_fixa[j * T + tau];
            if (e < 0) continue;
            const double dp = o.dp_fixos[tau];
            prec_fixos[static_cast<std::size_t>(e)] = dp > 0.0 ? 1.0 / (dp * dp) : 0.0;
          }
        for (std::size_t k = 0; k < M.n_fixo; k++) {
          if (prec_fixos[k] < 0.0)
            throw Erro("a fixed-effect equation of the level-2 equations with no parameter");
          const auto ini = M.c.linha.begin() + static_cast<std::ptrdiff_t>(M.c.colptr[k]);
          const auto fim = M.c.linha.begin() + static_cast<std::ptrdiff_t>(M.c.colptr[k + 1]);
          const auto q = std::lower_bound(ini, fim, static_cast<std::uint32_t>(k));
          if (q == fim || *q != k)
            throw Erro("a fixed-effect column with no diagonal in the level-2 equations");
          pos_diag_fixos.push_back(static_cast<std::size_t>(q - M.c.linha.begin()));
        }
      }
      for (std::size_t k = 0; k < pos_diag_fixos.size(); k++)
        M.c.valor[pos_diag_fixos[k]] += prec_fixos[k];
    }
    if (!cs.pronto) {
      cs.perm = grau_minimo(M.c);
      cs.sb = simbolica(permuta_sim(M.c, cs.perm));
      cs.pronto = true;
    }
    const Csc& pc = permuta_cache(M.c, cs);
    Csc L;
    if (!cholesky(pc, cs.sb, L)) return false;
    std::vector<double> pb(M.total);
    for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[cs.perm[k]];
    std::vector<double> sol_p = resolve(L, pb);
    if (sol_out) {
      sol_out->assign(M.total, 0.0);
      for (std::size_t k = 0; k < M.total; k++) (*sol_out)[cs.perm[k]] = sol_p[k];
    }
    loc.assign(M.total, 0.0);
    if (ruido) {
      std::vector<double> z(M.total);
      for (std::size_t k = 0; k < M.total; k++) z[k] = Rf_rnorm(0.0, 1.0);
      const std::vector<double> rp = resolve_lt(L, z);
      // variante 4: o braco errado do portao E3, o ruido vezes raiz de 2
      const double fator = o.variante == 4 ? std::sqrt(2.0) : 1.0;
      for (std::size_t k = 0; k < M.total; k++) loc[cs.perm[k]] = sol_p[k] + fator * rp[k];
    } else {
      for (std::size_t k = 0; k < M.total; k++) loc[cs.perm[k]] = sol_p[k];
    }
    return true;
  };
  auto monta_lw = [&]() {
    lw.assign(n, {});
    for (std::size_t r = 0; r < n; r++)
      for (std::size_t ex = d.x.linptr[r]; ex < d.x.linptr[r + 1]; ex++) {
        const std::size_t j = d.x.coluna[ex];
        const double v = d.x.valor[ex];
        if (v == 0.0) continue;
        for (std::size_t tau = 0; tau < T; tau++) {
          const int e = d.eq_fixa[j * T + tau];
          if (e >= 0)
            lw[r].push_back({static_cast<std::uint32_t>(e), static_cast<std::uint32_t>(tau), v});
        }
      }
    for (std::size_t g = 0; g < ng; g++) {
      const DesenhoTermo& a = *termo_g[g];
      for (std::size_t c = 0; c < a.z.ncol; c++)
        for (std::size_t k = a.z.colptr[c]; k < a.z.colptr[c + 1]; k++) {
          const std::size_t r = a.z.linha[k];
          const double v = a.z.valor[k];
          if (v == 0.0 || r >= n) continue;
          const std::size_t ct = c / a.n_niveis, nv = c % a.n_niveis;
          for (std::size_t tau = 0; tau < T; tau++)
            lw[r].push_back({static_cast<std::uint32_t>(off_g[g] + (tau * a.n_coef + ct) *
                                                                   a.n_niveis + nv),
                             static_cast<std::uint32_t>(tau), v});
        }
    }
  };
  // mu = W loc e a parte fixa X b de cada animal
  auto calcula_mu = [&]() {
    for (std::size_t r = 0; r < n; r++) {
      double m[3] = {0.0, 0.0, 0.0}, f[3] = {0.0, 0.0, 0.0};
      for (const EntW& e : lw[r]) {
        m[e.tr] += e.v * loc[e.col];
        if (e.col < n_fixo) f[e.tr] += e.v * loc[e.col];
      }
      for (std::size_t tau = 0; tau < T; tau++) {
        mu[r * T + tau] = m[tau];
        xb[r * T + tau] = f[tau];
      }
    }
  };

  GeradorR gerador;
  std::vector<double> sol_partida;
  if (!nivel2(false, &sol_partida)) {
    gerador.devolve();
    throw Erro("the level-2 equations at the start are not positive-definite");
  }
  monta_lw();
  calcula_mu();

  // ---- os acumuladores
  std::vector<EstadoAnimal> est(n);
  for (std::size_t i = 0; i < n; i++) {
    est[i].log_s = LOG_S0;
    for (std::size_t q = 0; q < 3; q++) est[i].phi[q] = phi[i * T + q];
  }
  S.ncomp = ntheta + nc;
  const std::size_t n_guardadas = o.n_iter > o.burnin ? (o.n_iter - o.burnin - 1) / o.thin + 1 : 0;
  S.amostras.reserve(n_guardadas * S.ncomp);
  std::vector<double> m_phi(3 * n, 0.0), m2_phi(3 * n, 0.0);
  std::vector<double> m_loc(total, 0.0), m2_loc(total, 0.0);
  const std::size_t nb = d.n_fixa;
  std::vector<double> m_b(nb, 0.0), m2_b(nb, 0.0), b_dado(nb, 0.0);
  S.lse.assign(nw, MENOS_INF);
  S.media_ld.assign(nw, 0.0);
  S.m2_ld.assign(nw, 0.0);
  S.media_ajuste.assign(nw, 0.0);
  const std::size_t n_pdm = std::min(o.n_pdm, n_guardadas);
  std::vector<std::size_t> quais_pdm;
  for (std::size_t k = 0; k < n_pdm; k++)
    quais_pdm.push_back(((2 * k + 1) * n_guardadas) / (2 * n_pdm));
  std::size_t prox_pdm = 0;
  S.pdm.reserve(n_pdm * nw);
  // os derivados: do grupo genetico, nas idades e na grade
  const bool com_der = o.grupo_genetico >= 0;
  const std::size_t gg = com_der ? static_cast<std::size_t>(o.grupo_genetico) : 0;
  const std::size_t nlg = com_der ? nl_g[gg] : 0;
  const std::size_t na_id = o.idades.size(), ngr = o.grade.size();
  const std::size_t nd = na_id + 2;
  if (o.classe_idade.size() != na_id) throw Erro("one residual class per age");
  if (com_der) {
    S.media_der.assign(nlg * nd, 0.0);
    S.m2_der.assign(nlg * nd, 0.0);
    S.media_gc.assign(nlg * ngr, 0.0);
    S.m2_gc.assign(nlg * ngr, 0.0);
  }
  S.media_ac.assign(n * ngr, 0.0);
  S.pop.reserve(n_guardadas * ngr);
  S.h2_idade.reserve(n_guardadas * na_id);
  S.loglik.reserve(n_guardadas);

  std::vector<double> buf(9 * n, 0.0), ef(nw, 0.0), ee(nw, 0.0);
  std::vector<double> sse(nc, 0.0);
  std::size_t n_amostras = 0;
  const double log_2pi = std::log(2.0 * M_PI);

  for (std::size_t it = 0; it < o.n_iter; it++) {
    checa_interrupcao();
    if (o.verboso && o.n_iter >= 10 && (it + 1) % (o.n_iter / 10) == 0)
      Rprintf("gibbs_growth %d/%d\n", (int) (it + 1), (int) o.n_iter);
    // ---- 1. as curvas
    if (o.segura == 0) {
      double sp[9], prec[9], uprec[9];
      cov3_de(theta, mo.offset_residual, sp);
      if (!precisao3(sp, prec, uprec)) {
        S.mensagem = "the permanent covariance drawn inside the chain is not positive-definite";
        break;
      }
      for (std::size_t i = 0; i < n; i++) sorteia_animal(buf.data() + 9 * i);
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads())
      for (std::size_t i = 0; i < n; i++) {
        const AlvoAnimal a = alvo_de(cc, p, i, mu.data() + 3 * i, prec, uprec, s2e.data());
        passo_animal(a, est[i], buf.data() + 9 * i, o.variante, o.movimentos, it, o.burnin);
      }
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t q = 0; q < 3; q++) phi[i * T + q] = est[i].phi[q];
    }
    // a curva e o residuo de cada pesagem no phi corrente (escala de analise)
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads())
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = p.ini[i]; j < p.ini[i + 1]; j++) {
        double f;
        curva_analise(cc, p.t[j], phi.data() + 3 * i, f, nullptr);
        ef[j] = f;
        ee[j] = p.y[j] - f;
      }
    // ---- 3. as localizacoes do nivel 2
    {
      std::vector<double> sol1;
      const bool primeira = o.segura == 1 && it == 0;
      if (!nivel2(true, primeira ? &sol1 : nullptr)) {
        S.mensagem = "the level-2 equations are not positive-definite inside the chain (a "
                     "variance drawn near zero)";
        break;
      }
      if (primeira) {
        S.sol_primeira = sol1;
        S.b_primeira.assign(sol1.begin(), sol1.begin() + static_cast<std::ptrdiff_t>(d.n_fixa));
        for (std::size_t tau = 0; tau < T; tau++)
          desfaz_escala_b(d.volta_x[tau], [&](std::size_t j) { return d.eq_fixa[j * T + tau]; },
                          S.b_primeira.data(), true);
      }
      calcula_mu();
    }
    // ---- 4. as covariancias e as variancias
    bool parou = false;
    if (o.segura != 1) {
      // cada grupo: C_g | U ~ InvWishart(gl, U'K^-1 U + priori)
      for (std::size_t g = 0; g < ng && !parou; g++) {
        const Grupo& gr = mo.grupos[g];
        const std::size_t nl = nl_g[g];
        const std::size_t off = off_g[g];
        Densa su(T, T);
        std::vector<std::vector<double>> ku(T);
        for (std::size_t b = 0; b < T; b++) {
          ku[b].assign(nl, 0.0);
          const double* u = loc.data() + off + b * nl;
          if (d.kinv[g].ncol > 0) {
            const Csc& k = d.kinv[g];
            for (std::size_t c = 0; c < nl; c++)
              for (std::size_t q = k.colptr[c]; q < k.colptr[c + 1]; q++) {
                const std::size_t r = k.linha[q];
                ku[b][r] += k.valor[q] * u[c];
                if (r != c) ku[b][c] += k.valor[q] * u[r];
              }
          } else {
            for (std::size_t l = 0; l < nl; l++) ku[b][l] = u[l];
          }
        }
        for (std::size_t a = 0; a < T; a++)
          for (std::size_t b = 0; b <= a; b++) {
            double s = 0.0;
            const double* ua = loc.data() + off + a * nl;
            for (std::size_t l = 0; l < nl; l++) s += ua[l] * ku[b][l];
            su.at(a, b) = s;
            su.at(b, a) = s;
          }
        const PrioriCurva& pr = o.priori_grupos[g];
        if (pr.tipo == 3)
          for (std::size_t a = 0; a < T; a++) su.at(a, a) += pr.gl0 * pr.escala[a];
        folga_diagonal(su);
        Densa cg;
        try {
          cg = rinvwishart(gl_priori(pg_de(pr), nl, T), su);
        } catch (const Erro&) {
          S.mensagem = "the scale matrix of group '" + gr.nome + "' is not positive-definite "
                       "inside the chain (its variance drawn near zero)";
          parou = true;
          break;
        }
        for (std::size_t j = 0; j < T; j++)
          for (std::size_t i = j; i < T; i++) theta[gr.offset + vech3(i, j)] = cg.at(i, j);
      }
      // Sigma_p | phi, loc ~ InvWishart(gl, E'E + priori), E = phi - W loc
      if (!parou) {
        Densa se(T, T);
        for (std::size_t r = 0; r < n; r++)
          for (std::size_t a = 0; a < T; a++)
            for (std::size_t b = 0; b <= a; b++)
              se.at(a, b) += (phi[r * T + a] - mu[r * T + a]) * (phi[r * T + b] - mu[r * T + b]);
        for (std::size_t a = 0; a < T; a++)
          for (std::size_t b = 0; b < a; b++) se.at(b, a) = se.at(a, b);
        const PrioriCurva& pr = o.priori_pe;
        if (pr.tipo == 3)
          for (std::size_t a = 0; a < T; a++) se.at(a, a) += pr.gl0 * pr.escala[a];
        folga_diagonal(se);
        Densa cp;
        try {
          cp = rinvwishart(gl_priori(pg_de(pr), n, T), se);
        } catch (const Erro&) {
          S.mensagem = "the scale matrix of the permanent term is not positive-definite "
                       "inside the chain";
          parou = true;
        }
        if (!parou)
          for (std::size_t j = 0; j < T; j++)
            for (std::size_t i = j; i < T; i++)
              theta[mo.offset_residual + vech3(i, j)] = cp.at(i, j);
      }
      // s2e_c | y, phi ~ (SSE_c + priori) / chisq(gl), somado em serie na ordem fixa
      if (!parou) {
        std::fill(sse.begin(), sse.end(), 0.0);
        for (std::size_t j = 0; j < nw; j++)
          sse[static_cast<std::size_t>(p.classe[j])] += ee[j] * ee[j];
        const PrioriCurva& pr = o.priori_res;
        for (std::size_t c = 0; c < nc; c++)
          s2e[c] = (sse[c] + (pr.tipo == 3 ? pr.gl0 * pr.escala[0] : 0.0)) /
                   Rf_rchisq(gl_priori(pg_de(pr), n_classe[c], 1));
      }
    }
    if (parou) break;

    // ---- 5. guarda
    if (it >= o.burnin && ((it - o.burnin) % o.thin == 0)) {
      n_amostras++;
      const double na = static_cast<double>(n_amostras);
      for (std::size_t k = 0; k < ntheta; k++) S.amostras.push_back(theta[k]);
      for (std::size_t c = 0; c < nc; c++) S.amostras.push_back(s2e[c]);
      welford(m_phi, m2_phi, phi.data(), 3 * n, na);
      welford(m_loc, m2_loc, loc.data(), total, na);
      for (std::size_t k = 0; k < nb; k++) b_dado[k] = loc[k];
      for (std::size_t tau = 0; tau < T; tau++)
        desfaz_escala_b(d.volta_x[tau], [&](std::size_t j) { return d.eq_fixa[j * T + tau]; },
                        b_dado.data(), true);
      welford(m_b, m2_b, b_dado.data(), nb, na);
      // o WAIC, a curva media e o log L do dado
      double ll = 0.0;
      for (std::size_t j = 0; j < nw; j++) {
        const double v = s2e[static_cast<std::size_t>(p.classe[j])];
        const double ld0 = -0.5 * (log_2pi + std::log(v)) - 0.5 * ee[j] * ee[j] / v;
        ll += ld0;
        const double ld = ld0 - p.logpeso[j];
        const double mx = std::max(S.lse[j], ld);
        S.lse[j] = mx + std::log(std::exp(S.lse[j] - mx) + std::exp(ld - mx));
        const double delta = ld - S.media_ld[j];
        S.media_ld[j] += delta / na;
        S.m2_ld[j] += delta * (ld - S.media_ld[j]);
        S.media_ajuste[j] += (ef[j] - S.media_ajuste[j]) / na;
      }
      S.loglik.push_back(ll);
      if (prox_pdm < quais_pdm.size() && quais_pdm[prox_pdm] == n_amostras - 1) {
        for (std::size_t j = 0; j < nw; j++)
          S.pdm.push_back(ee[j] / std::sqrt(s2e[static_cast<std::size_t>(p.classe[j])]));
        prox_pdm++;
      }
      // phi_ref: a media de X_i b nos animais com pesagem
      double ref[3] = {0.0, 0.0, 0.0};
      for (std::size_t r = 0; r < n; r++)
        for (std::size_t tau = 0; tau < T; tau++) ref[tau] += xb[r * T + tau];
      for (std::size_t tau = 0; tau < T; tau++) ref[tau] /= static_cast<double>(n);
      for (std::size_t g2 = 0; g2 < ngr; g2++) S.pop.push_back(curva_peso(cc, o.grade[g2], ref));
      // a curva de cada animal na grade
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads())
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t g2 = 0; g2 < ngr; g2++) {
          const double v = curva_peso(cc, o.grade[g2], phi.data() + 3 * i);
          S.media_ac[i * ngr + g2] += (v - S.media_ac[i * ngr + g2]) / na;
        }
      if (com_der) {
        // o h2 de primeira ordem de W(t) e os valores geneticos dos derivados
        double g0[9], tot[9], sp[9];
        cov3_de(theta, mo.grupos[gg].offset, g0);
        cov3_de(theta, mo.offset_residual, sp);
        for (int q = 0; q < 9; q++) tot[q] = sp[q];
        for (std::size_t g = 0; g < ng; g++) {
          double cg[9];
          cov3_de(theta, mo.grupos[g].offset, cg);
          for (int q = 0; q < 9; q++) tot[q] += cg[q];
        }
        for (std::size_t k = 0; k < na_id; k++) {
          double f, jac[3];
          if (!curva_analise(cc, o.idades[k], ref, f, jac)) {
            S.h2_idade.push_back(std::nan(""));
            continue;
          }
          const double num = quad3(g0, jac);
          const double den = quad3(tot, jac) + s2e[static_cast<std::size_t>(o.classe_idade[k])];
          S.h2_idade.push_back(num / den);
        }
        std::vector<double> base(nd), base_gc(ngr);
        for (std::size_t k = 0; k < na_id; k++) base[k] = curva_peso(cc, o.idades[k], ref);
        base[na_id] = peso_inflexao(cc, ref);
        base[na_id + 1] = taxa_maxima(cc, ref);
        const std::size_t off = off_g[gg];
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads())
        for (std::size_t l = 0; l < nlg; l++) {
          const double ph[3] = {ref[0] + loc[off + l], ref[1] + loc[off + nlg + l],
                                ref[2] + loc[off + 2 * nlg + l]};
          for (std::size_t k = 0; k < nd; k++) {
            const double x = k < na_id ? curva_peso(cc, o.idades[k], ph) - base[k]
                             : k == na_id ? peso_inflexao(cc, ph) - base[k]
                                          : taxa_maxima(cc, ph) - base[k];
            const std::size_t q = l * nd + k;
            const double delta = x - S.media_der[q];
            S.media_der[q] += delta / na;
            S.m2_der[q] += delta * (x - S.media_der[q]);
          }
          for (std::size_t g2 = 0; g2 < ngr; g2++) {
            const double x = curva_peso(cc, o.grade[g2], ph);
            const std::size_t q = l * ngr + g2;
            const double delta = x - S.media_gc[q];
            S.media_gc[q] += delta / na;
            S.m2_gc[q] += delta * (x - S.media_gc[q]);
          }
        }
      }
    }
  }
  gerador.devolve();
  S.n_amostras = n_amostras;
  S.n_pdm = prox_pdm;
  const double den = n_amostras > 1 ? static_cast<double>(n_amostras - 1) : std::nan("");
  S.media_phi = m_phi;
  S.var_phi.assign(3 * n, std::nan(""));
  S.media_loc = m_loc;
  S.var_loc.assign(total, std::nan(""));
  S.media_b = m_b;
  S.var_b.assign(nb, std::nan(""));
  if (n_amostras > 1) {
    for (std::size_t k = 0; k < 3 * n; k++) S.var_phi[k] = m2_phi[k] / den;
    for (std::size_t k = 0; k < total; k++) S.var_loc[k] = m2_loc[k] / den;
    for (std::size_t k = 0; k < nb; k++) S.var_b[k] = m2_b[k] / den;
  }
  S.ac_ind.assign(n, std::nan(""));
  S.ac_rw.assign(n, std::nan(""));
  S.escala_rw.assign(n, 0.0);
  S.tentativas.assign(n, 0.0);
  for (std::size_t i = 0; i < n; i++) {
    S.tentativas[i] = est[i].tentativas;
    if (est[i].tentativas > 0.0) {
      S.ac_ind[i] = est[i].ac_ind / est[i].tentativas;
      S.ac_rw[i] = est[i].ac_rw / est[i].tentativas;
    }
    S.escala_rw[i] = std::exp(est[i].log_s);
    S.modos_fora += static_cast<std::size_t>(est[i].falhas);
  }
  return S;
}

}  // namespace br
