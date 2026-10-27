// Gibbs (Geman e Geman, 1984; Wang, Rutledge e Gianola, 1993, 1994) para o modelo
// animal (uma caracteristica): o lado bayesiano do pacote.
//
// ## O desenho
//
// Bloco UNICO de localizacao: (b, u) | theta, y e Gaussiana com media na solucao das MME
// e covariancia C^-1, amostrada exata via a MESMA Cholesky esparsa do REML (fator L da
// C permutada; media + sqrt(s2e) P L'^-1 z, z ~ N(0, I)). Um refactor numerico por
// iteracao, simbolica em cache. Single-site misturaria pior e jogaria fora a maquinaria
// ja validada; o bloco unico faz da media amostral o proprio BLUP, que e o gate.
//
// Condicionais conjugadas. A priori de referencia de Jeffreys, p(s2) ~ 1/s2 e
// p(C) ~ |C|^-(d+1)/2 (antes este cabecalho dizia "prior PLANO", que e outra coisa: a plana
// na variancia da chisq(nl - 2)):
//   C_g | u  ~ InvWishart(nl_g, U' K^-1 U)          (Bartlett para a Wishart)
//   s2e | e  ~ e'e / chisq(n)
// As outras (PrioriGibbs) so mudam os graus de liberdade e somam nu0 S0 a escala. O padrao do
// gibbs() gaussiano (R/gibbs.R) e a plana nos grupos de covariancia (a meia-t no grupo com
// niveis de menos para ela, ExtrasGibbs::recuo_poucos_niveis), a meia-t nos termos escalares
// (ExtrasGibbs::priori_escalar) e a plana no residuo (priori_res, que o R passa); no probit, a
// meia-t com nu = 4 em todo termo. O passo ASIS em todo termo estimado no gaussiano e so nos grupos
// de covariancia no probit (ExtrasGibbs::asis_so_grupos), e, no probit, o PX-DA.
//
// A meia-t (tipo 4; Huang e Wand 2013) e hierarquica: C | a ~ InvWishart(nu + d - 1,
// 2 nu diag(1/a)), a_k ~ IG(1/2, 1/A_k^2), cada desvio-padrao meia-t(nu, A_k) e, com nu = 2,
// cada correlacao uniforme em (-1, 1). As duas condicionais sao conjugadas:
//   C_g | u, a ~ InvWishart(nl + nu + d - 1, U' K^-1 U + 2 nu diag(1/a))
//   a_k | C    ~ IG((nu + d)/2, nu (C^-1)_kk + 1/A_k^2),
//               sorteado como 2 (nu (C^-1)_kk + 1/A_k^2) / chisq(nu + d)
// A forma (nu + d)/2 sai da constante da InvWishart, |S|^(gl/2) com gl = nu + d - 1, vezes a
// IG(1/2) da priori; ela NAO depende de nl (uma forma (nu + nl)/2 faria a priori depender do
// tamanho do dado), e o portao V1 (test-gibbs-meia-t.R) confere as marginais da priori
// contra a meia-t e a densidade (1 - r^2)^(nu/2 - 1) da correlacao. A escala A_k e
// automatica (r2_meia_t, abaixo) e a cadeia inteira, partida inclusive, e equivariante na
// unidade de y, na escala de uma covariavel da base e na escala de K.
//
// Com asis (ExtrasGibbs), cada grupo ganha depois do sorteio centrado um passo de
// entrelacamento com a parametrizacao nao centrada U~ = U chol(C)^-T (Yu e Meng 2011), um
// Metropolis-Hastings cuja aceitacao e so a razao da priori vezes o jacobiano (passo_asis,
// abaixo, tem a derivacao).
//
// ## RNG
//
// O do PROPRIO R (GetRNGstate / Rf_rnorm / Rf_rchisq): respeita set.seed() do usuario e
// nao adiciona dependencia. Por isso este arquivo, ao contrario do resto da numerica,
// inclui R.h, ele E a fronteira.

#include "mme.h"
#include <algorithm>
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
std::vector<double> resolve_lt(const Csc& L, const std::vector<double>& z) {
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
// A = L T, S^-1 = L L'. Exige nu >= d. Com w_out, devolve tambem W = A A' = C^-1, a
// Wishart sorteada, que a meia-t le sem inverter C de novo.
Densa rinvwishart(double nu, const Densa& S, Densa* w_out) {
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
  if (w_out) *w_out = W;
  return inv_pd(W);
}

// a_k | C da meia-t: 2 (nu (C^-1)_kk + 1/A_k^2) / chisq(gl), com gl = nu + d. cinv e C^-1.
static void sorteia_a_meia_t(const Densa& cinv, double nu, const std::vector<double>& A,
                             double gl, std::vector<double>& a) {
  for (std::size_t k = 0; k < a.size(); k++)
    a[k] = 2.0 * (nu * cinv.at(k, k) + 1.0 / (A[k] * A[k])) / Rf_rchisq(gl);
}

// Um passo da meia-t num grupo: C | a, depois a | C. Su chega com U' K^-1 U (zero na priori
// pura do portao V1) e nl e o numero de niveis. variante = 1 sorteia a com a forma
// (nu + nl)/2 no lugar de (nu + d)/2: e a forma errada, que so o braco do portao usa.
static Densa passo_meia_t(Densa Su, double nl, double nu, const std::vector<double>& A,
                          std::vector<double>& a, int variante) {
  const std::size_t dim = Su.nlin;
  const double d = static_cast<double>(dim);
  for (std::size_t k = 0; k < dim; k++) Su.at(k, k) += 2.0 * nu / a[k];
  Densa W;
  Densa C = rinvwishart(nl + nu + d - 1.0, Su, &W);
  sorteia_a_meia_t(W, nu, A, variante == 1 ? nu + nl : nu + d, a);
  return C;
}

// O portao V1: as duas condicionais da meia-t sem dado (U' K^-1 U = 0, nl = 0), partindo de
// C = diag(A^2), `burnin` passos descartados e um de cada `thin` guardado. Devolve n amostras
// do vech de C (por colunas, a ordem de theta). O passo e o MESMO passo_meia_t da cadeia.
std::vector<double> gibbs_priori_meia_t(std::size_t dim, double nu, const std::vector<double>& A,
                                        std::size_t n, std::size_t thin, std::size_t burnin,
                                        int variante) {
  if (dim == 0 || A.size() != dim) throw Erro("A needs one entry per dimension");
  if (!(nu > 0.0)) throw Erro("nu must be positive");
  if (thin == 0) throw Erro("thin must be positive");
  for (double v : A)
    if (!(v > 0.0) || !std::isfinite(v)) throw Erro("A must be positive and finite");
  const std::size_t np = dim * (dim + 1) / 2;
  std::vector<double> out;
  out.reserve(n * np);
  GetRNGstate();
  std::vector<double> a(dim, 1.0);
  Densa c0(dim, dim);
  for (std::size_t k = 0; k < dim; k++) c0.at(k, k) = A[k] * A[k];
  sorteia_a_meia_t(inv_pd(c0), nu, A, variante == 1 ? nu : nu + static_cast<double>(dim), a);
  const Densa zero(dim, dim);
  for (std::size_t it = 0; it < burnin + n * thin; it++) {
    if (it % 4096 == 0) R_CheckUserInterrupt();
    const Densa C = passo_meia_t(zero, 0.0, nu, A, a, variante);
    if (it >= burnin && (it - burnin + 1) % thin == 0)
      for (std::size_t j = 0; j < dim; j++)
        for (std::size_t i = j; i < dim; i++) out.push_back(C.at(i, j));
  }
  PutRNGstate();
  return out;
}

// quem fala nos erros da meia-t: o usuario que nao escreveu prior = "half_t" a recebe no
// padrao (termos escalares no gaussiano, todo termo no probit), e a mensagem diz isso
static const std::string quem_meia_t =
    "the half-t prior (prior = \"half_t\", or the default of a scalar random term, of a "
    "covariance group with too few levels for \"flat\" and of family = \"probit\")";

// A escala automatica da meia-t, por componente k de cada grupo (os coeficientes dos termos
// do grupo, na ordem de theta): r_k^2 = (media, nas linhas usadas, de ||z_ik||^2, a norma ao
// quadrado da linha do bloco de incidencia do componente) x (media da diagonal de K nos
// niveis do componente com incidencia numa linha usada). r_k e 1 no animal() sem
// endogamia, cerca de sqrt(media(n - 1)) num indirect() sem diluicao e o valor quadratico
// medio da covariavel numa inclinacao de regressao aleatoria; com o fator de K, um K
// reescalado por c leva r_k por sqrt(c). A diagonal de K vem da inversa seletiva da K^-1 do
// grupo (a de H no passo unico, a da APY, a da precisao declarada): e a K que a cadeia usa.
// Nos niveis com registro porque numa precisao aumentada (dominance_inverse()) os niveis de
// subclasse nao tem registro e a priori deles (F/4) nao e a escala dos efeitos que o dado
// ve. Devolve r_k^2 por grupo (vazio no grupo preso, que nao anda, e no grupo que nao usa a
// meia-t, quais[g] = 0: no padrao gaussiano so os termos escalares a usam).
static std::vector<std::vector<double>> r2_meia_t(const Desenho& d,
                                                  const std::vector<std::size_t>& nl_g,
                                                  const std::vector<char>& quais) {
  const std::size_t n_grupos = d.modelo.grupos.size();
  std::vector<std::vector<double>> r2(n_grupos);
  const std::vector<std::string> nomes = d.modelo.nomes_theta();
  const double n_reg = static_cast<double>(d.n_usadas());
  for (std::size_t g = 0; g < n_grupos; g++) {
    const Grupo& gr = d.modelo.grupos[g];
    if (!quais[g] || d.modelo.preso(gr.offset)) continue;
    const std::size_t nl = nl_g[g];
    std::vector<double> dk;
    if (g < d.kinv.size() && d.kinv[g].ncol > 0) {
      if (d.kinv[g].ncol != nl)
        throw Erro("a declared K with a number of rows other than the levels of its term");
      const std::vector<std::size_t> perm = grau_minimo(d.kinv[g]);
      const Csc pk = permuta_sim(d.kinv[g], perm);
      const Simbolica sb = simbolica(pk);
      Csc L;
      if (!cholesky(pk, sb, L))
        throw Erro(quem_meia_t + ": the precision K^-1 of group '" + gr.nome +
                   "' is not positive-definite, so the mean diagonal of K that sets the "
                   "scale of its prior cannot be read");
      dk = diagonal_de_k(L, perm);
    }
    std::size_t b = 0;
    for (std::size_t tm : gr.termos)
      for (const DesenhoTermo& a : d.aleatorios) {
        if (a.termo != tm) continue;
        for (std::size_t c = 0; c < a.n_coef; c++, b++) {
          double sz = 0.0, sk = 0.0;
          std::size_t nk = 0;
          for (std::size_t l = 0; l < a.n_niveis; l++) {
            const std::size_t col = c * a.n_niveis + l;
            bool tem = false;
            for (std::size_t p = a.z.colptr[col]; p < a.z.colptr[col + 1]; p++)
              if (d.usa[a.z.linha[p]] && a.z.valor[p] != 0.0) {
                sz += a.z.valor[p] * a.z.valor[p];
                tem = true;
              }
            if (tem) {
              sk += dk.empty() ? 1.0 : dk[l];
              nk++;
            }
          }
          const std::string nome = b < gr.dim ? nomes[gr.theta_idx(b, b)] : gr.nome;
          if (nk == 0 || !(sz > 0.0))
            throw Erro(quem_meia_t + ": '" + nome + "' has no incidence on the records, "
                       "so the scale of its prior (sd(y) / r, with r^2 the mean squared row "
                       "of its incidence) is undefined");
          const double v = (sz / n_reg) * (sk / static_cast<double>(nk));
          if (!(v > 0.0) || !std::isfinite(v))
            throw Erro(quem_meia_t + ": the mean diagonal of K over the levels of '" +
                       nome + "' is not positive, so the scale of its prior is undefined");
          r2[g].push_back(v);
        }
        break;
      }
    if (r2[g].size() != gr.dim)
      throw Erro(quem_meia_t + ": the incidence of group '" + gr.nome +
                 "' does not match its dimension");
  }
  return r2;
}

// graus de liberdade e acrescimo de escala da condicional de um bloco dim x dim com nl
// niveis (ou do residuo, dim = 1 e nl = n registros), pela priori escolhida
double gl_priori(const PrioriGibbs& pr, std::size_t nl, std::size_t dim) {
  const double n = static_cast<double>(nl), d = static_cast<double>(dim);
  switch (pr.tipo) {
    case 1: return n - d - 1.0;
    case 2: return n - 1.0;
    case 3: return n + pr.gl0;
    default: return n;
  }
}

// ---- ASIS: entrelacamento das parametrizacoes centrada e nao centrada de cada grupo (Yu e
// Meng 2011; o molde e o de Kastner e Fruhwirth-Schnatter 2014 para a volatilidade
// estocastica). A cadeia de sempre e a CENTRADA: (b, u) | C e depois C | U. Quando o dado
// informa pouco sobre U, C | U fica preso a U e a cadeia anda devagar; na parametrizacao
// NAO centrada, U~ = U L^-T com L = chol(C), a priori de U~ (cada coluna N(0, K)) nao depende
// de C, e C | U~ fica preso ao dado. O passo leva (U, C) a (U~, C), sorteia C de novo dado
// U~ e o resto, e volta a U = U~ L*'.
//
// A conta (conferida aqui antes de escrever o codigo). Com vec(U) ~ N(0, C (x) K) e C = L L',
// vec(U~) = (L^-1 (x) I) vec(U) ~ N(0, I (x) K), livre de C; a mudanca de variavel de U para
// U~ tem jacobiano |L|^nl = |C|^(nl/2), que cancela o |C|^(-nl/2) da densidade de U, entao
// p(U~, C) = N(U~; 0, I (x) K) p(C) e, dado U~ (e a, na meia-t), p(C | U~, resto) e
// proporcional a verossimilhanca vezes p(C | a). A parte do grupo no preditor e
//   Z_g u_g = soma_j Z_j U_j = soma_j soma_{m <= j} L_jm Z_j U~_m,
// LINEAR em vech(L): a regressao do residuo parcial r = y - X b - (os outros grupos) nas
// d(d+1)/2 colunas X_jm = Z_j U~_m, com ruido s2e (1 no probit, na liability). A densidade
// alvo em vech(L), no conjunto L_ii > 0, e
//   pi(L) ~ N(r; X vech(L), s2e I) p(C = L L' | a) J(L),
// com J(L) = 2^d prod_i L_ii^(d - i + 1) (i de 1 a d) o jacobiano de vech(L) -> vech(L L')
// para L triangular inferior (no d = 2: as derivadas em (L11, L21, L22) dao uma matriz
// triangular de diagonal 2 L11, L11, 2 L22, determinante 4 L11^2 L22). A proposta e a
// posteriori gaussiana de vech(L) sob a priori plana, q(L) = N(vech(L); bhat, s2e (X'X)^-1),
// bhat = (X'X)^-1 X'r; como ||r - X l||^2 = ||r - X bhat||^2 + (l - bhat)' X'X (l - bhat),
// q e a verossimilhanca vezes uma constante. Proposta independente, entao a aceitacao de
// Metropolis-Hastings e
//   min(1, pi(L*) q(L) / (pi(L) q(L*))) = min(1, p(L* L*' | a) J(L*) / (p(L L' | a) J(L))),
// a verossimilhanca cancela; uma proposta com algum L*_ii <= 0 tem pi = 0 e e recusada (a
// cadeia fica). A priori de cada tipo entra pela mesma forma: a condicional centrada e
// InvWishart(gl, U'K^-1 U + S0), e a verossimilhanca de U da |C|^(-nl/2) exp(-tr(U'K^-1 U
// C^-1)/2), entao p(C) ~ |C|^(-kappa/2) exp(-tr(S0 C^-1)/2) com kappa = gl - nl + d + 1:
// jeffreys d + 1, plana 0, uniforme no desvio d (= 1), c(df, scale) df + d + 1 com S0 = df
// scale I, meia-t nu + 2d com S0 = 2 nu diag(1/a). Com log|C| = 2 soma log L_ii,
//   log p(C) + log J(L) = soma_i (d - i + 1 - kappa) log L_ii - tr(S0 C^-1)/2 + constante.
// O alvo precisa integrar perto de L_ii = 0, onde a verossimilhanca em vech(L) fica positiva
// (U = U~ L' com L_ii = 0 ainda tem efeito pelas outras entradas, e no escalar o efeito zero
// tem verossimilhanca positiva): o expoente d - i + 1 - kappa tem de passar de -1, e a linha
// que manda e i = d, kappa < 2. A plana (0) e a uniforme no desvio (1, escalar) passam; a
// propria c(df, scale) e a meia-t passam pelo termo exp(-s0_k (C^-1)_kk / 2), que com
// (C^-1)_dd = 1/L_dd^2 vai a zero mais depressa que qualquer potencia. Jeffreys (kappa = d + 1)
// da prod_i L_ii^(-i): L^-1 no escalar e L11^-1 L22^-2 no 2 x 2, sem integral. O
// Metropolis-Hastings independente entao aceita a ida para o polo e quase nao volta (a razao
// de aceitacao sai de uma L pequena com probabilidade da ordem de L / L*); por isso o gibbs()
// (R/gibbs.R) recusa asis = TRUE num grupo estimado sob Jeffreys, e validation/gibbs_asis.R
// (modo jeffreys) mede o braco recusado.
// A folga relativa de 1e-10 da InvWishart centrada (fora da meia-t) nao entra aqui: e uma
// estabilizacao numerica do sorteio, nao parte da priori. A conta e equivariante: com c y,
// L vai a c L, U~ nao muda, bhat e o sorteio vao por c e a razao de aceitacao nao muda; com
// uma covariavel da base por c, L vai a D L com D diagonal, e o mesmo vale. Custo O(nnz(Z) d)
// por grupo, sem fatorar nada de novo.
struct CompAsis {
  const DesenhoTermo* a = nullptr;
  std::size_t c = 0;   // o coeficiente do termo (as colunas c * n_niveis + l de a->z)
};

// chol_pequena (aireml.cpp, declarada em mme.h): L triangular inferior com zeros acima

// log p(L L') + log J(L), sem as constantes: soma_i (d - i - kappa) log L_ii, com i de 0, menos
// soma_k s0_k (C^-1)_kk / 2, (C^-1)_kk = ||L^-1 e_k||^2 pela substituicao para frente
// (sem_jacobiano: o braco errado do portao, so a priori)
static double log_alvo_asis(const Densa& L, double kappa, const std::vector<double>& s0,
                            bool sem_jacobiano) {
  const std::size_t n = L.nlin;
  double v = 0.0;
  for (std::size_t i = 0; i < n; i++)
    v += ((sem_jacobiano ? 0.0 : static_cast<double>(n - i)) - kappa) * std::log(L.at(i, i));
  std::vector<double> col(n);
  for (std::size_t k = 0; k < n; k++) {
    if (s0[k] == 0.0) continue;
    std::fill(col.begin(), col.end(), 0.0);
    col[k] = 1.0 / L.at(k, k);
    double q = col[k] * col[k];
    for (std::size_t i = k + 1; i < n; i++) {
      double t = 0.0;
      for (std::size_t m = k; m < i; m++) t -= L.at(i, m) * col[m];
      col[i] = t / L.at(i, i);
      q += col[i] * col[i];
    }
    v -= 0.5 * s0[k] * q;
  }
  return v;
}

// Um passo ASIS no grupo de componentes `comps` (nl niveis, localizacoes em loc[off + j nl +
// l], vech de C em theta[th_off + ...]). e_iter e o residuo realizado da iteracao (y - W loc;
// no probit, a liability menos W loc) e sai atualizado com o loc novo. Devolve 1 aceito, 0
// recusado, 2 proposta fora (algum L*_ii <= 0), -1 nao rodou (C ou X'X sem Cholesky). x, ut e
// parcial sao espacos de trabalho reaproveitados entre iteracoes. variante: 1 o passo, 2 e 3
// os bracos errados do portao (sem o jacobiano; aceita toda proposta de dentro)
static int passo_asis(const Desenho& d, const std::vector<CompAsis>& comps, std::size_t nl,
                      std::size_t off, std::size_t th_off, double s2e, double kappa,
                      const std::vector<double>& s0, std::vector<double>& theta,
                      std::vector<double>& loc, std::vector<double>& e_iter,
                      const std::vector<std::size_t>& usadas,
                      std::vector<std::vector<double>>& x, std::vector<double>& ut,
                      std::vector<double>& parcial, int variante) {
  const std::size_t dim = comps.size();
  const std::size_t nq = dim * (dim + 1) / 2;
  Densa C(dim, dim);
  for (std::size_t j = 0; j < dim; j++)
    for (std::size_t i = j; i < dim; i++) {
      C.at(i, j) = theta[th_off + vech_idx(i, j, dim)];
      C.at(j, i) = C.at(i, j);
    }
  Densa L;
  if (!chol_pequena(C, L)) return -1;
  // U~ = U L^-T, nivel a nivel: L u~_l = u_l (as linhas de U e de U~ como colunas)
  ut.assign(dim * nl, 0.0);
  for (std::size_t l = 0; l < nl; l++)
    for (std::size_t i = 0; i < dim; i++) {
      double s = loc[off + i * nl + l];
      for (std::size_t m = 0; m < i; m++) s -= L.at(i, m) * ut[m * nl + l];
      ut[i * nl + l] = s / L.at(i, i);
    }
  // as colunas X_jm = Z_j U~_m, m <= j, no indice vech_idx(j, m) (o de theta)
  x.resize(nq);
  for (std::vector<double>& v : x) v.assign(d.nlin, 0.0);
  for (std::size_t j = 0; j < dim; j++) {
    const DesenhoTermo& a = *comps[j].a;
    for (std::size_t l = 0; l < nl; l++) {
      const std::size_t col = comps[j].c * a.n_niveis + l;
      for (std::size_t p = a.z.colptr[col]; p < a.z.colptr[col + 1]; p++) {
        const std::size_t r = a.z.linha[p];
        const double v = a.z.valor[p];
        if (!d.usa[r] || v == 0.0) continue;
        for (std::size_t m = 0; m <= j; m++) x[vech_idx(j, m, dim)][r] += v * ut[m * nl + l];
      }
    }
  }
  std::vector<double> beta(nq);
  for (std::size_t j = 0; j < dim; j++)
    for (std::size_t m = 0; m <= j; m++) beta[vech_idx(j, m, dim)] = L.at(j, m);
  // o residuo parcial r = e + Z_g u_g = e + X vech(L), X'X e X'r
  Densa G(nq, nq);
  std::vector<double> h(nq, 0.0);
  parcial.assign(d.nlin, 0.0);
  for (std::size_t r : usadas) {
    double v = e_iter[r];
    for (std::size_t q = 0; q < nq; q++) v += beta[q] * x[q][r];
    parcial[r] = v;
    for (std::size_t q = 0; q < nq; q++) {
      h[q] += x[q][r] * v;
      for (std::size_t q2 = 0; q2 <= q; q2++) G.at(q, q2) += x[q][r] * x[q2][r];
    }
  }
  for (std::size_t q = 0; q < nq; q++)
    for (std::size_t q2 = 0; q2 < q; q2++) G.at(q2, q) = G.at(q, q2);
  Densa Lg;
  if (!chol_pequena(G, Lg)) return -1;
  // bhat = (X'X)^-1 X'r e o sorteio bhat + sqrt(s2e) Lg^-T z, de covariancia s2e (X'X)^-1
  std::vector<double> w(h);
  for (std::size_t i = 0; i < nq; i++) {
    for (std::size_t k = 0; k < i; k++) w[i] -= Lg.at(i, k) * w[k];
    w[i] /= Lg.at(i, i);
  }
  std::vector<double> z(nq);
  for (std::size_t q = 0; q < nq; q++) z[q] = w[q] + std::sqrt(s2e) * Rf_rnorm(0.0, 1.0);
  for (std::size_t i = nq; i-- > 0;) {
    for (std::size_t k = i + 1; k < nq; k++) z[i] -= Lg.at(k, i) * z[k];
    z[i] /= Lg.at(i, i);
  }
  Densa Ls(dim, dim);
  for (std::size_t j = 0; j < dim; j++)
    for (std::size_t m = 0; m <= j; m++) Ls.at(j, m) = z[vech_idx(j, m, dim)];
  for (std::size_t i = 0; i < dim; i++)
    if (!(Ls.at(i, i) > 0.0) || !std::isfinite(Ls.at(i, i))) return 2;
  const double la = log_alvo_asis(Ls, kappa, s0, variante == 2) -
                    log_alvo_asis(L, kappa, s0, variante == 2);
  if (variante != 3 && !(std::log(Rf_runif(0.0, 1.0)) < la)) return 0;
  // aceito: C = L* L*', U = U~ L*' e o residuo de cada linha usada com o loc novo
  for (std::size_t j = 0; j < dim; j++)
    for (std::size_t i = j; i < dim; i++) {
      double s = 0.0;
      for (std::size_t k = 0; k <= j; k++) s += Ls.at(i, k) * Ls.at(j, k);
      theta[th_off + vech_idx(i, j, dim)] = s;
    }
  for (std::size_t l = 0; l < nl; l++)
    for (std::size_t j = 0; j < dim; j++) {
      double s = 0.0;
      for (std::size_t m = 0; m <= j; m++) s += Ls.at(j, m) * ut[m * nl + l];
      loc[off + j * nl + l] = s;
    }
  for (std::size_t r : usadas) {
    double v = parcial[r];
    for (std::size_t q = 0; q < nq; q++) v -= z[q] * x[q][r];
    e_iter[r] = v;
  }
  return 1;
}

// Normal padrao truncada, pela inversa da CDF em escala LOG (Rf_pnorm5/Rf_qnorm5 com
// log_p): z > a quando acima, z <= a quando abaixo. Em escala linear a cauda longe (a = 9)
// da u * Phi(-9) ~ 1e-19 * u e o quantil sai infinito; em log nao.
static double normal_truncada(double a, bool acima) {
  const double lu = std::log(Rf_runif(0.0, 1.0));
  if (acima) return Rf_qnorm5(lu + Rf_pnorm5(a, 0.0, 1.0, 0, 1), 0.0, 1.0, 0, 1);
  return Rf_qnorm5(lu + Rf_pnorm5(a, 0.0, 1.0, 1, 1), 0.0, 1.0, 1, 1);
}

GibbsSaida gibbs(const Desenho& d_in, std::size_t n_iter, std::size_t burnin,
                 std::size_t thin, bool loc_fixa, const std::vector<double>* theta_fixo,
                 bool verboso, const PrioriGibbs& priori, bool probit, std::size_t n_pdm,
                 const PrioriGibbs* priori_res, const ExtrasGibbs* extras) {
  GibbsSaida S;
  const ExtrasGibbs sem_extras;
  const ExtrasGibbs& ex = extras ? *extras : sem_extras;
  if ((ex.dp_fixos > 0.0 || ex.px) && !probit)
    throw Erro("the normal prior of the fixed effects and the PX-DA move are for "
               "family = \"probit\"");
  if (ex.dp_fixos < 0.0 || !std::isfinite(ex.dp_fixos))
    throw Erro("the prior standard deviation of the fixed effects must be positive");
  // a priori do residuo: a mesma dos grupos quando nao vem a parte (o comportamento de
  // sempre, a mesma aritmetica)
  const PrioriGibbs& pres = priori_res ? *priori_res : priori;
  const bool meia_t_res = pres.tipo == 4 && !probit;
  // PROBIT (Albert e Chib 1993; no melhoramento, Sorensen, Andersen, Gianola e Korsgaard
  // 1995): uma liability l_i ~ N(w_i' loc, 1) truncada em (0, inf) quando y = 1 e em
  // (-inf, 0] quando y = 0 entra no lugar de y, e o resto da cadeia e a gaussiana de
  // sempre com a variancia residual FIXA em 1, que da a escala da liability. O limiar fica
  // em 0 e o intercepto implicito do desenho faz o papel dele. O desenho e copiado uma
  // vez so no probit, porque o y dele e reescrito a cada iteracao.
  Desenho d_probit;
  std::vector<char> y01;
  if (probit) {
    d_probit = d_in;
    y01.assign(d_in.nlin, 0);
    for (std::size_t r = 0; r < d_in.nlin; r++) {
      if (!d_in.usa[r]) continue;
      if (d_in.y[r] != 0.0 && d_in.y[r] != 1.0)
        throw Erro("family = \"probit\" takes a 0/1 trait");
      y01[r] = d_in.y[r] == 1.0;
      d_probit.y[r] = y01[r] ? 0.5 : -0.5;
    }
  }
  const Desenho& d = probit ? d_probit : d_in;
  const std::size_t ntheta = d.modelo.ntheta;
  // partida: var(y) repartida, como no REML; num grupo de K DECLARADA a parte dele e
  // dividida pela media geometrica dos autovalores de K (a mesma conta dos espelhos,
  // multitrait2.cpp), o que torna a cadeia EQUIVARIANTE: com cK a partida cai por c. Numa
  // precisao aumentada (kernel(Kinv = dominance_inverse())) a escala vem pronta, a de D e
  // nao a de K_aug, e a cadeia parte do mesmo ponto de kernel(K = D). As prioris agem em
  // sigma2 e nao dependem de K: a posteriori e a mesma, so a partida e a mistura (com os
  // niveis de subclasse latentes a mais) mudam.
  std::vector<double> theta = partida(d);
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    const Grupo& gr = d.modelo.grupos[g];
    if (gr.estrutura != Estrutura::Declarada || g >= d.kinv.size() || d.kinv[g].ncol == 0)
      continue;
    const double gm = escala_kernel(d.kinv_escala, d.kinv_logdet, d.kinv, g);
    if (std::isfinite(gm) && gm > 0.0)
      for (std::size_t i = 0; i < gr.dim; i++) theta[gr.theta_idx(i, i)] /= gm;
  }
  const std::size_t n_grupos = d.modelo.grupos.size();
  std::vector<std::size_t> nl_g(n_grupos, 0);
  for (std::size_t g = 0; g < n_grupos; g++)
    for (const DesenhoTermo& a : d.aleatorios)
      if (a.termo == d.modelo.grupos[g].termos[0]) { nl_g[g] = a.n_niveis; break; }
  // a priori de cada grupo: a `priori` dos grupos, ou a dos termos escalares quando ela vem a
  // parte (o padrao gaussiano: meia-t no grupo de dimensao 1, a plana nos de covariancia).
  // Sem escalar_a_parte todo grupo aponta para `priori` e a aritmetica e a de antes.
  // O recuo do padrao: a plana num grupo d x d da a condicional InvWishart(nl - d - 1, .),
  // propria so com nl - d - 1 > d - 1, nl > 2d (a conta que, sem o recuo, para a cadeia com
  // "too few levels" no passo 3), e com media finita so com nl - d - 1 > d + 1, nl > 2d + 2: a
  // media a posteriori de C e a media dessa condicional sobre a posteriori (esperanca
  // iterada), entao sem ela a media que o gibbs() reporta nao existe. Com
  // nl <= 2d + 2 o grupo fica com a priori dos escalares (a meia-t, propria em qualquer nl) e o
  // nome volta em grupos_recuo, que o R poe na mensagem. Com theta_fixed nada e sorteado e
  // nada recua
  std::vector<const PrioriGibbs*> pri_g(n_grupos, &priori);
  std::vector<char> meia_t_g(n_grupos, 0);
  bool meia_t = false;
  for (std::size_t g = 0; g < n_grupos; g++) {
    const Grupo& gr = d.modelo.grupos[g];
    if (ex.escalar_a_parte && gr.dim == 1) {
      pri_g[g] = &ex.priori_escalar;
    } else if (ex.escalar_a_parte && ex.recuo_poucos_niveis && !loc_fixa &&
               !d.modelo.preso(gr.offset) && pri_g[g]->tipo != 4 &&
               !(gl_priori(*pri_g[g], nl_g[g], gr.dim) > static_cast<double>(gr.dim) + 1.0)) {
      pri_g[g] = &ex.priori_escalar;
      S.grupos_recuo.push_back(gr.nome);
      S.grupos_recuo_nl.push_back(static_cast<int>(nl_g[g]));
      S.grupos_recuo_d.push_back(static_cast<int>(gr.dim));
    }
    meia_t_g[g] = pri_g[g]->tipo == 4;
    meia_t = meia_t || meia_t_g[g];
  }
  if (ex.escalar_a_parte && !ex.priori_escalar.nomes_mult.empty())
    throw Erro("the prior of the scalar terms takes no multiplier by name");
  // MEIA-T: a escala A_k de cada componente (r2_meia_t) vezes o multiplicador do usuario, e
  // a do residuo, sd(y) vezes o dele; no probit s = 1, a escala da liability. A partida de
  // cada variancia de grupo passa a ser a parte de var(y) dela dividida por r_k^2 (com as
  // covariancias em 0), a mesma regra da escala: assim a cadeia inteira, partida inclusive,
  // e equivariante na unidade de y, na escala de uma covariavel da base (a partida comum
  // var(y) / 2 / n_grupos nao muda quando a covariavel muda de unidade) e na escala de K.
  // Com theta_fixed nada disto entra: a cadeia nao sorteia componentes.
  const double nu_res = pres.gl0;
  std::vector<std::vector<double>> A_g(n_grupos), a_g(n_grupos);
  double A_res = 0.0, a_res = 0.0;
  S.escala_meia_t.assign(ntheta, std::nan(""));
  if (!loc_fixa && (meia_t || meia_t_res)) {
    double s = 1.0;
    if (!probit) {
      double soma = 0.0, ss = 0.0;
      std::size_t n = 0;
      for (std::size_t r = 0; r < d.nlin; r++)
        if (d.usa[r]) { soma += d.y[r]; n++; }
      const double m = soma / static_cast<double>(n);
      for (std::size_t r = 0; r < d.nlin; r++)
        if (d.usa[r]) ss += (d.y[r] - m) * (d.y[r] - m);
      s = n > 1 ? std::sqrt(ss / static_cast<double>(n - 1)) : 0.0;
      if (!(s > 0.0) || !std::isfinite(s))
        throw Erro(quem_meia_t + " takes its scale from sd(y), and the trait has no "
                   "spread over the records");
    }
    const std::vector<std::string> nomes = d.modelo.nomes_theta();
    if (meia_t) {
      const std::vector<std::vector<double>> r2 = r2_meia_t(d, nl_g, meia_t_g);
      // os nomes de half_t(scale =) sao os da `priori` dos grupos (a dos escalares nao tem)
      std::vector<char> usado(priori.nomes_mult.size(), 0);
      const double frac = 0.5 / static_cast<double>(n_grupos);
      for (std::size_t g = 0; g < n_grupos; g++) {
        const Grupo& gr = d.modelo.grupos[g];
        if (r2[g].empty()) continue;
        const PrioriGibbs& pg = *pri_g[g];
        A_g[g].assign(gr.dim, 0.0);
        a_g[g].assign(gr.dim, 0.0);
        for (std::size_t i = 0; i < gr.dim; i++)
          for (std::size_t j = 0; j <= i; j++) theta[gr.theta_idx(i, j)] = 0.0;
        for (std::size_t k = 0; k < gr.dim; k++) {
          const std::size_t ik = gr.theta_idx(k, k);
          double mult = pg.nomes_mult.empty() ? pg.escala : 1.0;
          for (std::size_t q = 0; q < pg.nomes_mult.size(); q++)
            if (pg.nomes_mult[q] == nomes[ik]) { mult = pg.mult[q]; usado[q] = 1; }
          A_g[g][k] = mult * s / std::sqrt(r2[g][k]);
          S.escala_meia_t[ik] = A_g[g][k];
          theta[ik] = frac * s * s / r2[g][k];
        }
      }
      for (std::size_t q = 0; q < usado.size(); q++)
        if (!usado[q]) {
          std::string lista;
          for (std::size_t g = 0; g < n_grupos; g++)
            for (std::size_t k = 0; !r2[g].empty() && k < d.modelo.grupos[g].dim; k++)
              lista += (lista.empty() ? "'" : ", '") +
                       nomes[d.modelo.grupos[g].theta_idx(k, k)] + "'";
          throw Erro("half_t(scale =) names '" + priori.nomes_mult[q] + "', which is not a "
                     "variance the chain samples; the names are " + lista);
        }
    }
    if (meia_t_res) {
      A_res = pres.escala * s;
      S.escala_meia_t[d.modelo.offset_residual] = A_res;
    }
  }
  // a partida pedida (start = do gibbs()) no lugar da de cima, a regra da meia-t inclusive: e
  // o que da as cadeias partidas dispersas que o R-hat entre cadeias pede. No probit o residuo
  // fica em 1 (linha de baixo); cada grupo que anda precisa partir positivo-definido
  if (!ex.partida.empty()) {
    if (theta_fixo) throw Erro("start = and theta_fixed = do not go together");
    const std::vector<std::string> nomes = d.modelo.nomes_theta();
    std::string lista;
    for (const std::string& s : nomes) lista += (lista.empty() ? "'" : ", '") + s + "'";
    if (ex.partida.size() != ntheta)
      throw Erro("start has " + std::to_string(ex.partida.size()) + " entries and the chain "
                 "has " + std::to_string(ntheta) + " components: " + lista);
    if (!ex.nomes_partida.empty() && ex.nomes_partida != nomes)
      throw Erro("the names of start are not the components of the chain, in their order: " +
                 lista);
    for (double v : ex.partida)
      if (!std::isfinite(v)) throw Erro("start must be finite");
    theta = ex.partida;
    for (std::size_t g = 0; g < n_grupos; g++) {
      if (d.modelo.preso(d.modelo.grupos[g].offset)) continue;
      Densa cg = cov_grupo(d.modelo, theta, g);
      if (!chol_densa(cg))
        throw Erro("start: the covariance of group '" + d.modelo.grupos[g].nome +
                   "' is not positive-definite");
    }
    if (!probit && !(theta[d.modelo.offset_residual] > 0.0))
      throw Erro("start: the residual variance must be positive");
  }
  if (probit) theta[d.modelo.offset_residual] = 1.0;
  if (theta_fixo) theta = *theta_fixo;
  if (probit && theta[d.modelo.offset_residual] != 1.0)
    throw Erro("with family = \"probit\" the residual variance is 1 by definition: "
               "theta_fixed must end with 1");
  // kernel(fixed = v): o componente preso comeca no valor pedido e nao e amostrado
  for (std::size_t k = 0; k < ntheta; k++)
    if (d.modelo.preso(k)) theta[k] = d.modelo.theta_fixo[k];
  // a partida vai no resultado: e o que confere que uma precisao declarada parte do mesmo
  // ponto da K densa equivalente
  S.partida = theta;
  if (priori.tipo == 2)
    for (const Grupo& gr : d.modelo.grupos)
      if (gr.dim > 1)
        throw Erro("prior = \"uniform_sd\" is for scalar variances; a covariance group "
                   "of dimension > 1 needs \"jeffreys\", \"flat\", \"half_t\" or "
                   "c(df = , scale = )");

  // O PISO DO RESIDUO (lista de correcoes, item 39). Num residuo quase nao identificado (so a
  // soma P0 + s2e separa, a correlacao dentro do animal perto de 1) a cadeia pode andar ate
  // var(residual) ~1e-13 da escala, onde as MME W'W + s2e G^-1 perdem a Cholesky e a cadeia
  // para (uma reproducao sintetica de 100 animais de 12 maes parava em 12063 de 19000
  // amostras). O piso nao cobre toda parada: no dado real a cadeia parou com o residuo em
  // 1e-7 var(y), acima dele, na crista de um termo aninhado num efeito fixo (a mensagem da
  // parada, abaixo, diz o estado). O sorteio do residuo fica TRUNCADO em
  // piso_res = 1e-10 var(y), que e o piso do residuo no REML (pisos_z, aireml.cpp): uma
  // priori truncada la, longe da massa a posteriori. O sorteio que cai abaixo do piso e
  // trocado por um da condicional truncada, pela inversa da acumulada; o par (sorteio comum,
  // troca so abaixo) da EXATAMENTE a truncada, e a cadeia que nunca toca o piso fica a mesma
  // bit a bit
  double piso_res = 0.0;
  std::size_t n_piso_res = 0;
  if (!probit) {
    double s = 0.0, s2 = 0.0, n = 0.0;
    for (std::size_t r = 0; r < d.nlin; r++)
      if (d.usa[r]) { s += d.y[r]; s2 += d.y[r] * d.y[r]; n += 1.0; }
    const double v = n > 1.0 ? (s2 - s * s / n) / (n - 1.0) : 0.0;
    piso_res = 1e-10 * (v > 0.0 && std::isfinite(v) ? v : 1.0);
  }
  // s2 = num / chisq(gl) truncado em s2 >= piso_res, isto e chisq <= num / piso_res
  auto trunca_res = [&](double s2, double num, double gl) {
    if (!(s2 < piso_res) || !(num > 0.0)) return s2;
    n_piso_res++;
    const double lim = num / piso_res;
    const double u = unif_rand() * Rf_pchisq(lim, gl, 1, 0);
    const double q = Rf_qchisq(u, gl, 1, 0);
    return q > 0.0 && std::isfinite(q) ? num / q : piso_res;
  };

  GetRNGstate();
  CacheSimbolica cs;

  // meia-t: os a partem de um sorteio do passo 2 na partida (o residuo e depois cada grupo,
  // a ordem da iteracao). So aqui a cadeia da meia-t consome o gerador antes da primeira
  // localizacao; nas outras prioris nada muda.
  if (!loc_fixa && meia_t_res)
    a_res = 2.0 * (nu_res / theta[d.modelo.offset_residual] + 1.0 / (A_res * A_res)) /
            Rf_rchisq(nu_res + 1.0);
  if (!loc_fixa && meia_t)
    for (std::size_t g = 0; g < n_grupos; g++) {
      if (A_g[g].empty()) continue;
      const Grupo& gr = d.modelo.grupos[g];
      const double nu_ht = pri_g[g]->gl0;
      Densa c0(gr.dim, gr.dim);
      for (std::size_t i = 0; i < gr.dim; i++)
        for (std::size_t j = 0; j < gr.dim; j++) c0.at(i, j) = theta[gr.theta_idx(i, j)];
      sorteia_a_meia_t(inv_pd(c0), nu_ht, A_g[g], nu_ht + static_cast<double>(gr.dim), a_g[g]);
    }

  std::size_t n_amostras = 0;
  std::size_t total_cols = 0;
  std::vector<double> media_u, m2_u;   // Welford (1962) para media e variancia das localizacoes
  std::vector<double> media_b(d.x.ncol, 0.0), m2_b(d.x.ncol, 0.0), b_dado(d.x.ncol, 0.0);

  // O residuo realizado da iteracao (e_s = y - W loc_s; no probit a liability nova menos
  // W loc_s), a media dele nas amostras guardadas, e os sorteios das discrepancias pivotais:
  // n_pdm amostras guardadas espacadas por igual, cada uma e_s / sqrt(s2e_s) com o s2e
  // sorteado DEPOIS de loc_s, o par (loc_s, s2e_s) e um sorteio da posteriori conjunta
  const std::size_t n_guardadas = n_iter > burnin ? (n_iter - burnin - 1) / thin + 1 : 0;
  const std::size_t n_pdm_ef = std::min(n_pdm, n_guardadas);
  std::vector<std::size_t> quais_pdm;
  for (std::size_t k = 0; k < n_pdm_ef; k++)
    quais_pdm.push_back(((2 * k + 1) * n_guardadas) / (2 * n_pdm_ef));
  std::size_t prox_pdm = 0;
  std::vector<std::size_t> usadas;
  for (std::size_t r = 0; r < d.nlin; r++)
    if (d.usa[r]) usadas.push_back(r);
  std::vector<double> e_iter(d.nlin, 0.0), media_e(d.nlin, 0.0);
  S.residuos_pdm.reserve(n_pdm_ef * usadas.size());
  // posicoes da diagonal do bloco fixo em M.c (o padrao vem do cache e nao muda)
  std::vector<std::size_t> pos_diag_fixos;
  // ASIS: os componentes de cada grupo que anda (termo e coeficiente, na ordem de theta, a
  // mesma das colunas do bloco do grupo) e os espacos de trabalho do passo
  if (ex.asis < 0 || ex.asis > 3) throw Erro("asis: unknown variant");
  const bool asis = ex.asis > 0 && !loc_fixa;
  std::vector<std::vector<CompAsis>> comp_asis(n_grupos);
  std::vector<std::vector<double>> asis_x;
  std::vector<double> asis_ut, asis_parcial;
  if (asis) {
    S.asis_tentativas.assign(n_grupos, 0.0);
    S.asis_aceites.assign(n_grupos, 0.0);
    S.asis_fora.assign(n_grupos, 0.0);
    for (std::size_t g = 0; g < n_grupos; g++) {
      const Grupo& gr = d.modelo.grupos[g];
      if (d.modelo.preso(gr.offset)) continue;
      // o padrao do probit deixa os termos escalares de fora (o passo nao foi medido neles)
      if (ex.asis_so_grupos && gr.dim == 1) continue;
      for (std::size_t tm : gr.termos)
        for (const DesenhoTermo& a : d.aleatorios)
          if (a.termo == tm) {
            if (a.n_niveis != nl_g[g])
              throw Erro("asis = TRUE: the terms of group '" + gr.nome + "' do not share "
                         "their levels");
            for (std::size_t c = 0; c < a.n_coef; c++) comp_asis[g].push_back({&a, c});
            break;
          }
      if (comp_asis[g].size() != gr.dim)
        throw Erro("asis = TRUE: the incidence of group '" + gr.nome +
                   "' does not match its dimension");
    }
  }

  for (std::size_t it = 0; it < n_iter; it++) {
    R_CheckUserInterrupt();
    if (verboso && n_iter >= 10 && (it + 1) % (n_iter / 10) == 0) {
      Rprintf("gibbs %d/%d\n", (int) (it + 1), (int) n_iter);
      imprime_theta(theta, d.modelo.nomes_theta(), d.modelo);
    }
    // ---- 1. localizacoes: (b, u) ~ N(solucao, s2e C_s^-1)
    Montado M = monta_mme(d, theta, &cs);
    if (!M.ok) { S.mensagem = "theta INADMISSIBLE inside the chain"; break; }
    // priori normal dos efeitos fixos, b ~ N(0, dp^2 I) na codificacao do motor: soma
    // s2e / dp^2 (o sistema esta em unidades de s2e; no probit s2e = 1) na diagonal do bloco
    // fixo. O padrao do fator nao muda (a diagonal de X'X ja esta la), entao a simbolica em
    // cache vale. Com o dado binario e os fixos planos, multiplicar todas as liabilities por
    // tau multiplica a verossimilhanca integrada em b por cerca de tau^p: e o que esta
    // priori corta
    if (ex.dp_fixos > 0.0) {
      if (pos_diag_fixos.size() != M.n_fixo) {
        pos_diag_fixos.clear();
        for (std::size_t k = 0; k < M.n_fixo; k++) {
          const auto ini = M.c.linha.begin() + static_cast<std::ptrdiff_t>(M.c.colptr[k]);
          const auto fim = M.c.linha.begin() + static_cast<std::ptrdiff_t>(M.c.colptr[k + 1]);
          const auto p = std::lower_bound(ini, fim, static_cast<std::uint32_t>(k));
          if (p == fim || *p != k)
            throw Erro("a fixed-effect column with no diagonal in the mixed-model equations");
          pos_diag_fixos.push_back(static_cast<std::size_t>(p - M.c.linha.begin()));
        }
      }
      const double prec = M.s2e / (ex.dp_fixos * ex.dp_fixos);
      for (std::size_t p : pos_diag_fixos) M.c.valor[p] += prec;
    }
    if (total_cols == 0) {
      total_cols = M.total;
      media_u.assign(M.total, 0.0);
      m2_u.assign(M.total, 0.0);
    }
    if (!cs.pronto) {
      cs.perm = grau_minimo(M.c);
      cs.sb = simbolica(permuta_sim(M.c, cs.perm));
      cs.pronto = true;
      S.bloco_denso = bloco_denso_simbolico(cs.sb);
      S.colunas_fator = cs.sb.n;
    }
    const Csc& pc = permuta_cache(M.c, cs);
    Csc L;
    // a matriz dos coeficientes das MME (nao a covariancia de um grupo): perde a Cholesky
    // quando uma variancia sorteada cai a ~1e-13 da escala, como o residuo sob Jeffreys com um
    // registro por animal (validation/indirect_gibbs_recovery.R, modo half_t)
    if (!cholesky(pc, cs.sb, L)) {
      // o ESTADO no ponto da falha (item 39 da lista de correcoes): o residuo sobre var(y) e,
      // por grupo, a menor variancia sobre var(y) e o menor autovalor sobre o maior, para que a
      // mensagem diga o que foi ao limite
      std::string estado, grande;
      char buf[160];
      if (!probit) {
        std::snprintf(buf, sizeof buf, "var(residual) %.3g var(y)",
                      theta[d.modelo.offset_residual] / (piso_res / 1e-10));
        estado = buf;
      }
      for (std::size_t g = 0; g < n_grupos; g++) {
        const Densa cg = cov_grupo(d.modelo, theta, g);
        std::vector<double> ev;
        Densa u;
        jacobi_sim(cg, ev, u);
        double mn = ev.empty() ? 0.0 : ev[0], mx = mn, dmin = cg.at(0, 0);
        for (double e : ev) { mn = std::min(mn, e); mx = std::max(mx, e); }
        for (std::size_t a = 0; a < cg.nlin; a++) dmin = std::min(dmin, cg.at(a, a));
        if (!probit && mx > 10.0 * (piso_res / 1e-10))
          grande += (grande.empty() ? "'" : ", '") + d.modelo.grupos[g].nome + "'";
        if (cg.nlin > 1)
          std::snprintf(buf, sizeof buf, "group '%s': smallest variance %.3g var(y), smallest "
                        "eigenvalue %.3g of the largest", d.modelo.grupos[g].nome.c_str(),
                        dmin / (probit ? 1.0 : piso_res / 1e-10), mx > 0.0 ? mn / mx : 0.0);
        else
          std::snprintf(buf, sizeof buf, "var(%s) %.3g var(y)", d.modelo.grupos[g].nome.c_str(),
                        dmin / (probit ? 1.0 : piso_res / 1e-10));
        estado += (estado.empty() ? "" : ", ") + std::string(buf);
      }
      S.mensagem = "the coefficient matrix of the mixed-model equations is not "
                   "positive-definite inside the chain (a variance drawn near zero; at the "
                   "failure: " + estado + ")";
      // medido no dado real (item 39): o residuo em 1e-7 var(y) e a familia em 20 var(y) com a
      // correlacao no limite, a familia aninhada na caixa fixa. A parte entre caixas desse termo
      // so a priori identifica, e a cadeia anda por la ate o fator nao sair
      if (!grande.empty())
        S.mensagem += "; the variance of group(s) " + grande + " passed 10 var(y) there: a "
                      "term whose levels each sit in one class of a fixed effect (a dam in one "
                      "box) has its part between classes identified only by the prior, and the "
                      "chain wanders there; fit that class as random() or merge the classes";
      break;
    }

    std::vector<double> pb(M.total);
    for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[cs.perm[k]];
    std::vector<double> sol_p = resolve(L, pb);

    // PX-DA de Haar (Wang e Roy 2018, Algoritmo 2; Hobert e Marchev 2008): depois da
    // liability e antes das localizacoes, todas as liabilities usadas l passam a g l, com
    // g^2 ~ Gamma(n/2, taxa q/2) e q = l'l - sol'W'l, a forma quadratica de l com (b, u)
    // integrados (o sistema do probit esta em s2e = 1). Dado C, p(l | C) e proporcional a
    // exp(-q(l)/2) nas restricoes de sinal, q(g l) = g^2 q(l) e as restricoes nao mudam com
    // g > 0 (o limiar e 0): o passo deixa p(l | C) invariante, e as localizacoes sorteadas
    // em seguida da liability nova fecham a varredura sem mudar a posteriori. A media do
    // sorteio e linear em l, entao ela vai a g sol sem resolver de novo
    if (ex.px) {
      double ll = 0.0;
      for (std::size_t r : usadas) ll += d_probit.y[r] * d_probit.y[r];
      double bry = 0.0;
      for (std::size_t k = 0; k < M.total; k++) bry += sol_p[k] * pb[k];
      const double q = ll - bry;
      if (q > 0.0 && std::isfinite(q)) {
        const double gx = std::sqrt(Rf_rchisq(static_cast<double>(usadas.size())) / q);
        for (std::size_t r : usadas) d_probit.y[r] *= gx;
        for (double& v : sol_p) v *= gx;
      }
    }

    std::vector<double> z(M.total);
    for (std::size_t k = 0; k < M.total; k++) z[k] = Rf_rnorm(0.0, 1.0);
    std::vector<double> ruido_p = resolve_lt(L, z);

    const double s2e = theta[d.modelo.offset_residual];
    std::vector<double> loc(M.total);
    for (std::size_t k = 0; k < M.total; k++)
      loc[cs.perm[k]] = sol_p[k] + std::sqrt(s2e) * ruido_p[k];

    // ---- 2. residuo: e = y - W loc; s2e ~ e'e / chisq(n) com Jeffreys (as outras prioris mudam
    // os graus de liberdade, gl_priori: n - 2 na plana, n - 1 na uniforme no desvio)
    double ee = 0.0;
    std::size_t nreg = 0;
    {
      std::vector<double> wl(d.nlin, 0.0);
      for (std::size_t r = 0; r < d.nlin; r++) {
        if (!d.usa[r]) continue;
        double s = 0.0;
        for (std::size_t e = d.x.linptr[r]; e < d.x.linptr[r + 1]; e++)
          s += d.x.valor[e] * loc[d.x.coluna[e]];
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
        e_iter[r] = e;
      }
      // probit: a liability nova, condicional as localizacoes, e o y da proxima rodada; o
      // residuo latente realizado e o dessa liability nova (o par com loc_s)
      if (probit)
        for (std::size_t r = 0; r < d.nlin; r++)
          if (d.usa[r]) {
            d_probit.y[r] = wl[r] + normal_truncada(-wl[r], y01[r] != 0);
            e_iter[r] = d_probit.y[r] - wl[r];
          }
    }
    if (!loc_fixa && !probit && meia_t_res) {
      // meia-t no residuo: s2e | e, a ~ (e'e + 2 nu / a) / chisq(n + nu), depois
      // a | s2e ~ 2 (nu / s2e + 1/A^2) / chisq(nu + 1)
      const double num = ee + 2.0 * nu_res / a_res;
      const double gl = static_cast<double>(nreg) + nu_res;
      const double s2 = trunca_res(num / Rf_rchisq(gl), num, gl);
      theta[d.modelo.offset_residual] = s2;
      a_res = 2.0 * (nu_res / s2 + 1.0 / (A_res * A_res)) / Rf_rchisq(nu_res + 1.0);
    } else if (!loc_fixa && !probit) {
      const double gl = gl_priori(pres, nreg, 1);
      if (!(gl > 0.0)) throw Erro("too few records for this prior on the residual");
      const double num = ee + (pres.tipo == 3 ? pres.gl0 * pres.escala : 0.0);
      theta[d.modelo.offset_residual] = trunca_res(num / Rf_rchisq(gl), num, gl);
    }

    // ---- 3. grupos: C_g ~ InvWishart(nl, U' K^-1 U) com Jeffreys (as outras prioris mudam os
    // graus de liberdade, gl_priori: nl - d - 1 na plana)
    bool parou_grupo = false;
    if (!loc_fixa) {
      for (std::size_t g = 0; g < n_grupos; g++) {
        const Grupo& gr = d.modelo.grupos[g];
        const PrioriGibbs& pg = *pri_g[g];
        const std::size_t dim = gr.dim;
        const std::size_t nl = nl_g[g];
        const std::size_t off = M.offset_grupo[g];
        // grupo preso por kernel(fixed =): escalar por construcao (entrada.cpp recusa o resto)
        if (d.modelo.preso(gr.offset)) continue;
        if (gr.estrutura == Estrutura::Declarada && d.kinv[g].ncol != nl)
          throw Erro("a declared K with a number of rows other than the levels of its term");
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
        // meia-t: C | a e a | C (passo_meia_t). U'K^-1 U + 2 nu diag(1/a) e positiva-definida
        // por construcao, entao a folga de 1e-10 de baixo nao entra
        if (meia_t_g[g]) {
          Densa cg;
          try {
            cg = passo_meia_t(Su, static_cast<double>(nl), pg.gl0, A_g[g], a_g[g], 0);
          } catch (const Erro&) {
            parou_grupo = true;
          }
          if (parou_grupo) {
            S.mensagem = "the scale matrix of group '" + gr.nome + "' is not positive-definite "
                         "inside the chain (its variance drawn near zero)";
            break;
          }
          for (std::size_t j = 0; j < dim; j++)
            for (std::size_t i = j; i < dim; i++)
              theta[gr.offset + vech_idx(i, j, dim)] = cg.at(i, j);
          continue;
        }
        // priori propria: soma nu0 S0 na diagonal da escala
        if (pg.tipo == 3)
          for (std::size_t a = 0; a < dim; a++) Su.at(a, a) += pg.gl0 * pg.escala;
        // estabiliza contra Su singular no comeco da cadeia. RELATIVO a escala de Su: um
        // 1e-10 absoluto era um piso escondido na priori, e nao escalava com K, entao a
        // cadeia com cK nao era a cadeia com K reescalada
        double tr = 0.0;
        for (std::size_t a = 0; a < dim; a++) tr += Su.at(a, a);
        const double eps = 1e-10 * (tr > 0.0 ? tr / static_cast<double>(dim) : 1.0);
        for (std::size_t a = 0; a < dim; a++) Su.at(a, a) += eps;
        const double gl = gl_priori(pg, nl, dim);
        if (!(gl > static_cast<double>(dim) - 1.0))
          throw Erro("too few levels in a covariance group for this prior: the inverse "
                     "Wishart needs more degrees of freedom than the dimension minus one");
        // uma variancia presa no zero (Jeffreys num termo escalar que o dado nao segura,
        // validation/gibbs_scalar_prior.R) chega a U'K^-1 U subnormal e a escala da InvWishart
        // perde a Cholesky: a cadeia para com a mensagem, como a das MME, e as outras cadeias
        // seguem. Antes o erro derrubava o gibbs() inteiro, todas as cadeias
        Densa cg;
        try {
          cg = rinvwishart(gl, Su);
        } catch (const Erro&) {
          parou_grupo = true;
        }
        if (parou_grupo) {
          S.mensagem = "the scale matrix of group '" + gr.nome + "' is not positive-definite "
                       "inside the chain (its variance drawn near zero)";
          break;
        }
        for (std::size_t j = 0; j < dim; j++)
          for (std::size_t i = j; i < dim; i++)
            theta[gr.offset + vech_idx(i, j, dim)] = cg.at(i, j);
      }
    }
    if (parou_grupo) break;

    // ---- 3b. ASIS em cada grupo que anda, depois do sorteio centrado de C (e de a, na
    // meia-t), na ordem dos grupos: o residuo realizado sai de cada passo com o loc novo e e
    // o de partida do grupo seguinte (passo_asis, acima, tem a conta)
    if (asis) {
      const double s2e_asis = theta[d.modelo.offset_residual];
      for (std::size_t g = 0; g < n_grupos; g++) {
        if (comp_asis[g].empty()) continue;
        const Grupo& gr = d.modelo.grupos[g];
        const std::size_t dim = gr.dim;
        std::vector<double> s0(dim, 0.0);
        const PrioriGibbs& pg = *pri_g[g];
        double kappa;
        if (meia_t_g[g]) {
          kappa = pg.gl0 + 2.0 * static_cast<double>(dim);
          for (std::size_t k = 0; k < dim; k++) s0[k] = 2.0 * pg.gl0 / a_g[g][k];
        } else {
          kappa = gl_priori(pg, nl_g[g], dim) - static_cast<double>(nl_g[g]) +
                  static_cast<double>(dim) + 1.0;
          if (pg.tipo == 3)
            for (std::size_t k = 0; k < dim; k++) s0[k] = pg.gl0 * pg.escala;
        }
        const int res = passo_asis(d, comp_asis[g], nl_g[g], M.offset_grupo[g], gr.offset,
                                   s2e_asis, kappa, s0, theta, loc, e_iter, usadas, asis_x,
                                   asis_ut, asis_parcial, ex.asis);
        if (it >= burnin && res >= 0) {
          S.asis_tentativas[g] += 1.0;
          if (res == 1) S.asis_aceites[g] += 1.0;
          if (res == 2) S.asis_fora[g] += 1.0;
        }
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
      // b desta amostra na unidade do dado (fixos.cpp), e o Welford dele
      for (std::size_t k = 0; k < d.x.ncol && k < total_cols; k++) b_dado[k] = loc[k];
      desfaz_escala_b(d.volta_x[0], [](std::size_t j) { return (int) j; },
                      b_dado.data(), true);
      for (std::size_t k = 0; k < d.x.ncol && k < total_cols; k++) {
        const double delta = b_dado[k] - media_b[k];
        media_b[k] += delta / na;
        m2_b[k] += delta * (b_dado[k] - media_b[k]);
      }
      for (std::size_t r : usadas) media_e[r] += (e_iter[r] - media_e[r]) / na;
      if (prox_pdm < quais_pdm.size() && quais_pdm[prox_pdm] == n_amostras - 1) {
        const double sd = std::sqrt(theta[d.modelo.offset_residual]);
        for (std::size_t r : usadas) S.residuos_pdm.push_back(e_iter[r] / sd);
        prox_pdm++;
      }
    }
  }

  PutRNGstate();
  if (n_piso_res > 0) {
    char buf[96];
    std::snprintf(buf, sizeof buf, "%.3g", piso_res);
    const std::string nota = "var(residual) was drawn below its floor " + std::string(buf) +
        " (1e-10 var(y), where the prior is truncated) " + std::to_string(n_piso_res) +
        " time(s) and was drawn again above it: the data hold almost no information to "
        "separate the residual from the other variances";
    S.mensagem = S.mensagem.empty() ? nota : S.mensagem + "; " + nota;
  }
  S.n_amostras = n_amostras;
  S.ntheta = ntheta;
  S.media_loc = media_u;
  S.var_loc.assign(total_cols, std::nan(""));
  if (n_amostras > 1)
    for (std::size_t k = 0; k < total_cols; k++)
      S.var_loc[k] = m2_u[k] / static_cast<double>(n_amostras - 1);
  if (n_amostras > 0 && total_cols >= d.x.ncol) {
    S.media_b = media_b;
    S.var_b.assign(d.x.ncol, std::nan(""));
    if (n_amostras > 1)
      for (std::size_t k = 0; k < d.x.ncol; k++)
        S.var_b[k] = m2_b[k] / static_cast<double>(n_amostras - 1);
  }
  S.media_residuo.assign(d.nlin, std::nan(""));
  if (n_amostras > 0)
    for (std::size_t r : usadas) S.media_residuo[r] = media_e[r];
  S.n_pdm = prox_pdm;
  // W E[loc] e X E[b], pela linearidade as medias a posteriori dos ajustados: as mesmas
  // contas do passo 2, com a media das localizacoes no lugar do sorteio
  S.ajustado.assign(d.nlin, std::nan(""));
  S.ajustado_fixo.assign(d.nlin, std::nan(""));
  if (n_amostras > 0 && total_cols > 0) {
    for (std::size_t r : usadas) {
      double s = 0.0;
      for (std::size_t e = d.x.linptr[r]; e < d.x.linptr[r + 1]; e++)
        s += d.x.valor[e] * media_u[d.x.coluna[e]];
      S.ajustado_fixo[r] = s;
      S.ajustado[r] = s;
    }
    std::size_t col0 = d.x.ncol;
    for (std::size_t g = 0; g < n_grupos; g++)
      for (std::size_t tm : d.modelo.grupos[g].termos)
        for (const DesenhoTermo& a : d.aleatorios)
          if (a.termo == tm) {
            for (std::size_t c = 0; c < a.z.ncol; c++) {
              const double bc = media_u[col0 + c];
              if (bc == 0.0) continue;
              for (std::size_t p = a.z.colptr[c]; p < a.z.colptr[c + 1]; p++)
                if (d.usa[a.z.linha[p]]) S.ajustado[a.z.linha[p]] += a.z.valor[p] * bc;
            }
            col0 += a.z.ncol;
          }
  }
  S.ok = S.mensagem.empty();
  return S;
}

}  // namespace br
