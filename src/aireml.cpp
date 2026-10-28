// AI-REML (Gilmour, Thompson e Cullis, 1995): score analitico, informacao media, EM
// (Dempster, Laird e Rubin, 1977) de aquecimento e passo amortecido.
//
// As regras que nao se quebram, todas vindas de defeito medido no projeto de origem:
//
//  1. convergencia RELATIVA nas componentes, sqrt(soma delta^2 / soma theta^2) < tol. Nunca
//     um limiar absoluto na norma do score, que e O(n_registros): a mesma tolerancia que
//     converge com 500 registros declara nao-convergido com 200 mil, com as estimativas
//     paradas no otimo. O certificado final soma o DECREMENTO DE NEWTON dos componentes
//     livres, g' AI^-1 g, normalizado pela curvatura, ele mede o gap em -2logL e nao a
//     escala do score, entao nao reintroduz o defeito de O(n_registros).
//  2. theta inadmissivel nao e resultado. C_g nao definida interrompe o passo; nunca vira
//     resposta.
//  3. a partida vem de var(y), nunca de valores de um cartao ja ajustado, partir das
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
#include <cstdlib>
#include <limits>

// so para o TRACE opcional e a interrupcao: nenhuma conta usa o R
#include <R_ext/Print.h>
#include <R_ext/Utils.h>

namespace br {
// MATRIZ DE INFORMACAO SINGULAR, dita pelo nome dos componentes.
//
// O erro padrao sai de 2 AI^-1, e quando a AI e singular numa direcao ele sai NaN sem que o
// ajuste diga por que. Medido num direto + indireto com 800 animais e poucas baias com
// parentes: var(animal), cov e var(indirect) com erro padrao NaN, e tres pontos bem
// diferentes -- var(indirect) 0.27, 0.16 e 0.02 -- com o MESMO -2logL. O dado nao separava
// os tres, e o numero reportado era um ponto arbitrario de uma crista plana. Quem le o
// ajuste precisa saber disso por escrito, e saber QUAIS componentes.
//
// O corte e o do certificado: autovalor abaixo de 1e-8 do maior e direcao sem curvatura. Os
// componentes com peso >= 0.2 no autovetor dessa direcao sao os que o dado nao separa.
//
// O CORTE E FEITO NAS UNIDADES PADRONIZADAS, F AI F com F = diag(fator), e nao na AI crua
// das unidades do usuario. A AI crua mistura unidades: com y em unidades c, a linha de uma
// variancia escala por 1/c^4 e a do rho do AR(1) nao escala, entao o mesmo ajuste era
// SINGULAR com y em minutos e regular com y em horas, com as MESMAS estimativas (medido no
// P1-5: 5 de 10 tracos univariados com EP NaN em unidades cruas, todos finitos em /100).
// `fator` e a unidade de cada componente sob a padronizacao por caracteristica de
// model_mt()/model_ar1(), s_a s_b com s o desvio padrao EXATO de cada caracteristica (1 para
// o rho): y e c y dao a mesma F AI F ate o arredondamento do dado, para qualquer c, e o
// mesmo veredito. O desvio exato, e nao a potencia de 2 logo acima dele com que a caminhada
// anda: aquela so e invariante a fatores potencia de 2, e com y vezes 100 a razao entre as
// duas muda por um fator em (1/2, 2) por caracteristica, que na AI vira ate 16 entre a linha
// de uma variancia e a do rho, o bastante para um autovalor perto do corte trocar de lado.
// No model() todos os componentes tem a unidade de y^2, a AI inteira escala por 1/c^4 e o
// corte relativo ao maior autovalor ja e invariante: ali `fator` fica nulo e a conta e a de
// sempre.
//
// A forma de correlacao, D^-1/2 AI D^-1/2, que tambem e invariante, foi tentada e
// descartada: ela normaliza a diagonal, e um parametro SEM informacao nenhuma deixa de ser
// visto (medido no desenho de test-grupo-iid-niveis.R, a covariancia entre dois termos de
// rotulos disjuntos com diagonal da AI 3,1e-31 contra 5 a 170 dos outros: SINGULAR e EP NaN
// com o corte cru, EP 9,6e15 e nenhum aviso na forma de correlacao).
std::string aviso_informacao(const Densa& ai_bruta, const std::vector<std::string>& nomes,
                             std::vector<double>* se, const std::vector<double>* fator) {
  const std::size_t n = ai_bruta.nlin;
  if (n == 0 || nomes.size() != n) return "";
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      if (!std::isfinite(ai_bruta.at(i, j))) return "";
  Densa ai = ai_bruta;
  if (fator && fator->size() == n)
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) ai.at(i, j) *= (*fator)[i] * (*fator)[j];
  std::vector<double> ev;
  Densa u;
  jacobi_sim(ai, ev, u);
  double maior = 0.0;
  for (double e : ev) maior = std::max(maior, std::fabs(e));
  if (!(maior > 0.0)) return "";
  const double corte = 1e-8 * maior;
  std::string grupos;
  for (std::size_t k = 0; k < ev.size(); k++) {
    if (ev[k] > corte) continue;
    std::string lista;
    for (std::size_t i = 0; i < n; i++)
      if (std::fabs(u.at(i, k)) >= 0.2) {
        lista += (lista.empty() ? "" : ", ") + nomes[i];
        // a mensagem diz NaN e o EP tem de ser NaN: invertida numericamente, uma AI
        // singular devolve NaN ou um numero enorme conforme o arredondamento (medido:
        // 6,8e5 depois de uma troca de ordenacao que nao mudou nada no ajuste)
        if (se && se->size() == n) (*se)[i] = std::nan("");
      }
    if (!lista.empty()) grupos += (grupos.empty() ? "" : "; and ") + lista;
  }
  if (grupos.empty()) return "";
  // O CONSELHO segue os componentes do veredito (lista de correcoes, item 20): o texto do
  // efeito indireto (baias de tamanhos diferentes) saia tambem num modelo materno, onde baia
  // nao existe
  const bool indireto = grupos.find("indirect") != std::string::npos;
  const bool materno = grupos.find("maternal") != std::string::npos;
  std::string mais;
  if (indireto) mais += "relatives sharing the pen, pens of DIFFERENT sizes";
  if (materno)
    mais += std::string(mais.empty() ? "" : "; ") +
        "dams with records of their own and several offspring, and granddams with records";
  if (mais.empty()) mais = "relatives sharing a level, more records per level";
  return "the information matrix is SINGULAR: the data do not separate " + grupos +
      ". Along that direction every combination of those components fits the data "
      "equally well, so their standard errors are NaN and the point estimate is one "
      "arbitrary point of a flat ridge, not an estimate: refitting from another start= "
      "lands elsewhere on the ridge with the same -2logL. Only combinations of them are "
      "estimable. This is the design, not a failed fit: it needs more information (" + mais +
      "), a simpler model, or some of those components held fixed";
}


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

// O -2logL do model() pelo sistema transformado (o racional em BlocoT, src/mme.h). As
// coordenadas do grupo transformado tem Var(u*) = s2e (I (x) K), e a forma quadratica sai como
// em avalia(), pelos residuos: e'e + u*'(I (x) K^-1)u* nos grupos transformados, e
// u'(C_g^-1 (x) K^-1)u nos outros, tudo sobre s2e
static bool neg2logl_transformada(const Desenho& d, const std::vector<double>& theta,
                                  const std::vector<Densa>& F, double& valor) {
  Montado M = monta_mme(d, theta, nullptr, &F);
  if (!M.ok) return false;
  const std::vector<std::size_t> perm = grau_minimo(M.c);
  const Csc pc = permuta_sim(M.c, perm);
  const Simbolica sb = simbolica(pc);
  Csc L;
  if (!cholesky(pc, sb, L)) return false;
  std::vector<double> pb(M.total);
  for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[perm[k]];
  const std::vector<double> px = resolve(L, pb);
  std::vector<double> sol(M.total, 0.0);
  for (std::size_t k = 0; k < M.total; k++) sol[perm[k]] = px[k];

  // upen_transformada da u'(C_g^-1 (x) K^-1)u com a C_g absoluta, que ja e o b'Pen b / s2e de
  // avalia(), e u*'(I (x) K^-1)u* / s2e nos grupos transformados
  const double s2e = M.s2e;
  double upen = 0.0;
  if (!upen_transformada(d.modelo, d.aleatorios, d.kinv, theta, M.offset_grupo, F, sol, s2e, upen))
    return false;
  destransforma_solucao(M.blocos_t, sol);

  std::vector<double> ehat(d.nlin, 0.0);
  std::size_t nreg = 0;
  for (std::size_t i = 0; i < d.nlin; i++) {
    if (!d.usa[i]) continue;
    nreg++;
    double wb = 0.0;
    for (std::size_t e = d.x.linptr[i]; e < d.x.linptr[i + 1]; e++)
      wb += d.x.valor[e] * sol[d.x.coluna[e]];
    ehat[i] = d.y[i] - wb;
  }
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    std::size_t col0 = M.offset_grupo[g];
    for (std::size_t t : d.modelo.grupos[g].termos)
      for (const DesenhoTermo& a : d.aleatorios)
        if (a.termo == t) {
          for (std::size_t cc = 0; cc < a.z.ncol; cc++) {
            const double bc = sol[col0 + cc];
            if (bc == 0.0) continue;
            for (std::size_t q = a.z.colptr[cc]; q < a.z.colptr[cc + 1]; q++)
              if (d.usa[a.z.linha[q]]) ehat[a.z.linha[q]] -= a.z.valor[q] * bc;
          }
          col0 += a.z.ncol;
        }
  }
  double ete = 0.0;
  for (std::size_t i = 0; i < d.nlin; i++)
    if (d.usa[i]) ete += ehat[i] * ehat[i];
  valor = (static_cast<double>(nreg) - static_cast<double>(M.n_fixo)) * std::log(s2e) +
          M.logdet_g + logdet(L) + ete / s2e + upen + d.logdet_escala_x;
  return std::isfinite(valor);
}

// O RESIDUO POR CLASSE (lista de correcoes, item 35). Com R = diag(s2_k(i)), as linhas
// escaladas por sqrt(w_i), w_i = s2_1 / s2_k(i), dao o modelo de Var(e~) = s2_1 I, e as MME,
// a solucao, os escores dos grupos e a AI do motor valem nessas linhas como estao: a AI no
// espaco original e (f~_a'f~_b - (W~'f~_a)' Cs^-1 (W~'f~_b)) / s2_1 com f~ = D^1/2 f, e o f~ de
// um grupo e o de sempre nas colunas escaladas. O que muda e so o bloco do residuo: o f~ da
// classe k e 1{i em k} e~_i / s2_k, o escore e (1/s2_k) soma_{i em k} [(1 - h_ii) - e~_i^2/s2_1]
// com h_ii a alavanca do registro, o EM e s2_k <- (s2_k/n_k) soma_{i em k} (e~_i^2/s2_1 + h_ii),
// e o -2logL ganha o jacobiano soma_k n_k log(s2_k / s2_1)
static Desenho escala_por_classe(const Desenho& d, const std::vector<double>& theta,
                                 std::vector<double>& w) {
  Desenho e = d;
  const std::size_t off = d.modelo.offset_residual;
  w.assign(d.nlin, 1.0);
  for (std::size_t i = 0; i < d.nlin; i++) w[i] = theta[off] / theta[off + d.classe_res[i]];
  std::vector<double> r(d.nlin);
  for (std::size_t i = 0; i < d.nlin; i++) r[i] = std::sqrt(w[i]);
  for (std::size_t i = 0; i < d.nlin; i++) {
    e.y[i] *= r[i];
    for (std::size_t q = e.x.linptr[i]; q < e.x.linptr[i + 1]; q++) e.x.valor[q] *= r[i];
  }
  for (DesenhoTermo& a : e.aleatorios)
    for (std::size_t q = 0; q < a.z.valor.size(); q++) a.z.valor[q] *= r[a.z.linha[q]];
  return e;
}

Avaliacao avalia(const Desenho& d0, const std::vector<double>& theta,
                 CacheSimbolica* cache) {
  // com residuo por classe, o resto de avalia() le as linhas escaladas; a parte dos dados das
  // MME muda com theta, entao o cache dela nao vale (a simbolica, do padrao, vale)
  const bool het = d0.modelo.n_residuos() > 1;
  for (std::size_t k = 0; het && k < d0.modelo.n_residuos(); k++)
    if (!(theta[d0.modelo.offset_residual + k] > 0.0)) return Avaliacao();
  std::vector<double> w_het;
  Desenho d_het;
  if (het) {
    d_het = escala_por_classe(d0, theta, w_het);
    if (cache) cache->mont_dono = nullptr;
  }
  const Desenho& d = het ? d_het : d0;
  Avaliacao A;
  Montado M = monta_mme(d, theta, cache);
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

  // -2logL: os termos que nao dependem de y aqui; y'Py entra depois, pelos residuos.
  //
  // A FORMA QUADRATICA SAI DOS RESIDUOS, e nao de y'y - b'r. As duas sao a mesma conta no
  // otimo das MME (b'r = b'(W'W + Pen) b), mas y'y - b'r subtrai dois numeros da ordem de
  // n (media de y)^2 e herda em PRIMEIRA ordem o erro da solucao b: somar uma constante a y
  // mudava o -2logL, que o REML nao deixa mudar. Medido no dado do portao D4
  // (test-unidade-partida-fronteira.R), em theta fixo: relativo 1,1e-8 com y + 1000 dp e
  // 1,1e-6 com y + 1e4 dp no model(), 6,6e-7 no AR(1) e 3,9e-7 no multicaracter com 1e4 dp,
  // e o colapso rho = 0 do model_ar1() contra o model() a -7,8e-11 com 100 dp e -7,9e-7 com
  // 1e4 dp; pelos residuos, tudo abaixo de 3e-14. No P1-5 (dado real) o colapso dava -9e-6.
  // e'e + b'Pen b e o minimo da soma de quadrados penalizada, estacionario em b: o erro da
  // solucao entra so em SEGUNDA ordem, e os residuos ja nao carregam a media.
  std::size_t nreg = 0;
  for (std::size_t i = 0; i < d.nlin; i++)
    if (d.usa[i]) nreg++;
  const double pfixo = static_cast<double>(M.n_fixo);
  const double parte_fixa = (static_cast<double>(nreg) - pfixo) * std::log(M.s2e) +
                            M.logdet_g + logdet(L);
  double upen = 0.0;   // soma_g soma_ab C_g^-1[a,b] Q[a,b] = b'Pen b / s2e

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
      for (std::size_t e = d.x.linptr[i]; e < d.x.linptr[i + 1]; e++)
        wb += d.x.valor[e] * A.solucao[d.x.coluna[e]];
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
    Densa cgs = cg;
    for (double& v : cgs.dados) v /= s2e;
    // Na fronteira (correlacao 1, a crista de um desenho aliasado) a Cholesky da montagem
    // pode aceitar a C_g com um pivo de 1e-17 e a LU daqui achar um pivo zero exato: o
    // arredondamento das duas difere, e qual delas vence muda de plataforma (medido no CI
    // do Linux, com o mesmo teste passando no Windows). O ponto e inadmissivel, como na
    // montagem, e nao um erro que derruba o ajuste.
    Densa cinv, cinv_s;
    try {
      cinv = inv_geral(cg);
      cinv_s = inv_geral(cgs);
    } catch (const Erro&) {
      return A;
    }

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
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++) upen += cinv.at(a, b) * Q.at(a, b);

    // T[a,b] = s2e tr(K^-1 [Cs^-1]_{a,b}): le a inversa seletiva SO onde K^-1 tem nao-zero.
    // Todo nao-zero de K^-1 esta na penalidade, a penalidade esta em C, e todo nao-zero de C
    // esta no padrao do fator: um pedido fora do padrao significa que esse fechamento quebrou,
    // e e CONTADO em vez de ignorado, uma soma que pula termos em silencio nao e um traco.
    Densa T(dim, dim);
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
        T.at(a, b) = s2e * s;
      }

    // tr(C_s^-1 Pen) deste grupo, para o EM do residual:
    // tr(kron(Cs^-1,K^-1) Z) somado, igual a soma_ab Cs^-1[a,b] tr(K^-1 Z_ba) = soma T/s2e * Cs^-1
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
  // medida em que o parentesco os deixa. O gate de diferencas finitas pegou na hora, e
  // com o score errado o passo amortecido ainda ACEITA (a verossimilhanca esta certa),
  // so que converge para o lugar errado. Um estimador com score errado e verossimilhanca
  // certa e o pior dos mundos: converge, reporta convergencia, e estima outra coisa.
  const double nreg_d = static_cast<double>(nreg);
  // y'Py = (e'e + b'Pen b) / s2e, com e e b do proprio ajuste (ver o -2logL acima)
  // + o termo da padronizacao das covariaveis: o -2logL sai na escala do dado (fixos.cpp)
  A.neg2logl = parte_fixa + ete / s2e + upen + d.logdet_escala_x;
  A.score[d.modelo.offset_residual] =
      (nreg_d - pfixo - static_cast<double>(q_ordem) + trpen_total) / s2e
      - ete / (s2e * s2e);
  // o residuo por classe: a alavanca h_ii de cada registro nas linhas escaladas, pela inversa
  // seletiva (os pares de colunas de uma linha estao no padrao do fator, porque estao em C)
  const std::size_t nres = d.modelo.n_residuos(), off_r = d.modelo.offset_residual;
  std::vector<double> alav;
  if (het) {
    const LinhasW lw = linhas_w(d, M.offset_grupo);
    alav.assign(d.nlin, 0.0);
    for (std::size_t i = 0; i < d.nlin; i++) {
      if (!d.usa[i]) continue;
      double h = 0.0;
      for (std::size_t a = lw.ptr[i]; a < lw.ptr[i + 1]; a++)
        for (std::size_t b = lw.ptr[i]; b < lw.ptr[i + 1]; b++) {
          double zv = 0.0;
          if (z_orig(lw.col[a], lw.col[b], zv)) h += lw.val[a] * lw.val[b] * zv;
          else A.fora_do_padrao++;
        }
      alav[i] = h;
    }
    std::vector<double> soma(nres, 0.0), soma_e(nres, 0.0), n_k(nres, 0.0);
    double jac = 0.0;
    for (std::size_t i = 0; i < d.nlin; i++) {
      if (!d.usa[i]) continue;
      const std::size_t k = d.classe_res[i];
      soma[k] += (1.0 - alav[i]) - ehat[i] * ehat[i] / s2e;
      soma_e[k] += ehat[i] * ehat[i] / s2e + alav[i];
      n_k[k] += 1.0;
      jac -= std::log(w_het[i]);
    }
    A.neg2logl += jac;
    for (std::size_t k = 0; k < nres; k++) {
      const double s2k = theta[off_r + k];
      A.score[off_r + k] = soma[k] / s2k;
      A.em_theta[off_r + k] = n_k[k] > 0.0 ? s2k * soma_e[k] / n_k[k] : s2k;
    }
  }

  // EM do residual: s2e <- (e'e + s2e (p + q - trpen)) / n. E[e'e | y] soma
  // s2e tr(C^-1 W'W) = s2e (p + q - trpen), e os p fixos entram ali como dado faltante de
  // priori plana (lista de correcoes, item 15). Sem o p, o ponto fixo do EM era
  // e'e / (n - q + trpen), o do residuo de ML, e nao o zero do score acima,
  // e'e / (n - p - q + trpen): com muitos niveis fixos por registro o aquecimento levava o
  // residuo para baixo do otimo REML
  if (!het)
    A.em_theta[d.modelo.offset_residual] =
        (ete + s2e * (pfixo + static_cast<double>(q_ordem) - trpen_total)) / nreg_d;

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
    if (d.usa[i]) {
      if (het) {
        const std::size_t k = d.classe_res[i];
        f[off_r + k][i] = ehat[i] / theta[off_r + k];
      } else {
        f[d.modelo.offset_residual][i] = ehat[i] / s2e;
      }
    }

  // AI[i,j] = (f_i' f_j - (W' f_i)' Cs^-1 (W' f_j)) / s2e, com W'f resolvido pelas MME
  std::vector<std::vector<double>> wtf(ntheta, std::vector<double>(ntot, 0.0));
  for (std::size_t k = 0; k < ntheta; k++) {
    // W'f: X'f e Z'f. X'f por linhas: cada entrada soma os registros na ordem crescente
    for (std::size_t i = 0; i < d.nlin; i++) {
      if (!d.usa[i]) continue;
      for (std::size_t e = d.x.linptr[i]; e < d.x.linptr[i + 1]; e++)
        wtf[k][d.x.coluna[e]] += d.x.valor[e] * f[k][i];
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

  // perto da singularidade de um grupo, o VALOR vem do sistema transformado; o score e a AI
  // continuam das MME, eles so dao a direcao, e quem aceita o passo e o valor
  {
    std::vector<Densa> cgs(d.modelo.grupos.size()), F;
    for (std::size_t g = 0; g < cgs.size(); g++) {
      cgs[g] = cov_grupo(d.modelo, theta, g);
      for (double& v : cgs[g].dados) v /= s2e;
    }
    double v = 0.0;
    if (fatores_quase_singulares(cgs, F) && neg2logl_transformada(d, theta, F, v)) {
      // o jacobiano do residuo por classe, que o valor transformado nao carrega
      if (het)
        for (std::size_t i = 0; i < d.nlin; i++)
          if (d.usa[i]) v -= std::log(w_het[i]);
      A.neg2logl = v;
    }
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
// O mapa diz QUAIS pedacos do theta andam em cada coordenada, para que os tres
// ajustadores usem esta mesma maquinaria. Um bloco e uma covariancia inteira no triangulo
// do theta (um grupo, ou o R0 de t caracteristicas); um escalar_log e uma variancia
// isolada, que e como o residuo do univariado sempre andou; e uma correlacao anda em
// atanh, o que empurra |rho| = 1 para o infinito da mesma forma que o log-Cholesky empurra
// det(C) = 0.
static std::size_t idx_tri(std::size_t off, std::size_t dim, std::size_t i, std::size_t j) {
  if (i < j) std::swap(i, j);
  return off + (j * (2 * dim - j + 1)) / 2 + (i - j);
}

bool z_de_theta(const MapaZ& mp, std::size_t ntheta, const std::vector<double>& theta,
                std::vector<double>& z) {
  z.assign(ntheta, 0.0);
  for (std::size_t b = 0; b < mp.blocos.size(); b++) {
    const std::size_t off = mp.blocos[b].first, dim = mp.blocos[b].second;
    Densa cg(dim, dim);
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const double v = theta[idx_tri(off, dim, i, j)];
        cg.at(i, j) = v;
        cg.at(j, i) = v;
      }
    Densa L;
    if (!chol_pequena(cg, L)) return false;
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++)
        z[idx_tri(off, dim, i, j)] = (i == j) ? std::log(L.at(i, i)) : L.at(i, j);
  }
  for (std::size_t k : mp.escalares_log) {
    if (!(theta[k] > 0.0)) return false;
    z[k] = std::log(theta[k]);
  }
  for (std::size_t k : mp.correlacoes) {
    if (!(std::fabs(theta[k]) < 1.0)) return false;
    z[k] = std::atanh(theta[k]);
  }
  return true;
}

// Pisos na diagonal de z, por bloco, recalculados no ponto corrente: relativo a maior
// diagonal de Cholesky do bloco (assim um bloco de uma linha so nunca encosta nele) e
// absoluto contra uma escala (piso de variancia 1e-10 escala, lido atraves de l = sqrt(v)).
// O model() passa a var(y) dos registros usados, fixa na caminhada (o racional e a medida
// estao em ajusta); model_mt() e model_ar1() caminham com cada caracteristica padronizada e
// passam max(1, variancia residual media) na unidade padronizada, que la e sempre 1. O piso relativo de 1e-3 limita a condicao do bloco a ~1e6, que
// e o que mantem a inversa do bloco utilizavel: sem ele o passo em z chega a fronteira
// exata (autovalor medido em 1.5e-14) e a solucao do MME para os efeitos fixos perde
// digitos contra a mesma GLS montada densa. Escalares em log e correlacoes em atanh nao
// tem piso: as duas transformacoes ja levam a fronteira ao infinito.
void pisos_z(const MapaZ& mp, const std::vector<double>& z, double escala_residual,
             std::vector<double>& piso, std::vector<char>& relativo) {
  piso.assign(z.size(), -std::numeric_limits<double>::infinity());
  relativo.assign(z.size(), 0);
  const double piso_abs = 0.5 * std::log(1e-10 * escala_residual);
  for (std::size_t b = 0; b < mp.blocos.size(); b++) {
    const std::size_t off = mp.blocos[b].first, dim = mp.blocos[b].second;
    double zmax = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < dim; i++)
      zmax = std::max(zmax, z[idx_tri(off, dim, i, i)]);
    const double piso_rel = zmax + std::log(1e-3);
    for (std::size_t i = 0; i < dim; i++) {
      const std::size_t k = idx_tri(off, dim, i, i);
      piso[k] = std::max(piso_rel, piso_abs);
      relativo[k] = piso_rel > piso_abs ? 1 : 0;
    }
  }
}

// TODAS as pecas do passo num ponto, numa chamada: z, o score e a AI levados para z pela
// regra da cadeia, os pisos, e os DOIS conjuntos ativos. Existe porque separar essas pecas
// entre o passo e o certificado foi um defeito medido: o passo foi portado para z e o
// certificado ficou em theta, com regra propria. O teste que ele usava, lmin < 1e-3 lmax,
// e exatamente a desigualdade que o piso do passo torna FALSA por construcao, no ponto em
// que o laco para, o grampo deixa lmin = 1e-3 lmax exato e o `<` estrito nunca dispara. A
// direcao grampeada ficava entao congelada no passo (relDelta ~ 0 sempre) e cobrada
// integralmente no certificado (decremento travado em 3.25e+03), e o ajuste girava ate o
// maxiter sem mover a verossimilhanca.
//
// Os dois conjuntos saem dos MESMOS pisos e da MESMA regra de score apontando para fora, e
// o do certificado e deliberadamente mais largo (piso + log 2 contra piso + 1e-9): o
// caminhante descansa um fio acima do grampo sem nunca ser congelado ali, e cobrar essa
// direcao certificaria um ponto com uma direcao viva ignorada.
PassoZ passo_z(const MapaZ& mp, std::size_t ntheta, const std::vector<double>& theta,
               const std::vector<double>& score, const Densa& ai, double escala_residual) {
  PassoZ P;
  P.ok = z_de_theta(mp, ntheta, theta, P.zc);
  if (!P.ok) return P;
  std::vector<char> rel;
  pisos_z(mp, P.zc, escala_residual, P.piso, rel);
  const Densa J = jacobiano_z(mp, ntheta, P.zc);
  P.sz.assign(ntheta, 0.0);
  for (std::size_t m = 0; m < ntheta; m++)
    for (std::size_t k = 0; k < ntheta; k++) P.sz[m] += J.at(k, m) * score[k];
  P.az = Densa(ntheta, ntheta);
  for (std::size_t a = 0; a < ntheta; a++)
    for (std::size_t b = 0; b < ntheta; b++) {
      double s = 0.0;
      for (std::size_t k = 0; k < ntheta; k++) {
        if (J.at(k, a) == 0.0) continue;
        double u = 0.0;
        for (std::size_t l = 0; l < ntheta; l++) u += ai.at(k, l) * J.at(l, b);
        s += J.at(k, a) * u;
      }
      P.az.at(a, b) = s;
    }
  P.congelado.assign(ntheta, 0);
  P.na_parede.assign(ntheta, 0);
  for (std::size_t i = 0; i < ntheta; i++) {
    if (!(P.piso[i] > -std::numeric_limits<double>::infinity())) continue;
    // O PASSO so congela quem esta no grampo com o score empurrando para FORA: se o score
    // aponta para dentro, a coordenada tem de poder subir e sair do piso sozinha.
    if (P.sz[i] > 0.0 && P.zc[i] <= P.piso[i] + 1e-9) P.congelado[i] = 1;
    // O CERTIFICADO nao olha o sinal. Estar descansando no piso depois de o caminhante ter
    // tido todas as iteracoes para sair dele JA e a evidencia de que a restricao esta
    // ativa; exigir score positivo ali reprova o proprio otimo restrito quando o
    // multiplicador e numericamente zero. Medido: um bloco declarado preso em posto 1 com
    // sz = -1e-04 ficava de fora da exclusao e sustentava o decremento em 0.034 num ponto
    // onde 4000 sorteios multivariados e uma busca em linha nao acharam melhora NENHUMA, e
    // onde perturbar o bloco em 1e-06 ja da theta inadmissivel.
    if (P.zc[i] <= P.piso[i] + std::log(2.0)) P.na_parede[i] = 1;
  }
  // e a parede prende o BLOCO, nao a coordenada: o piso relativo fixa a RAZAO entre as
  // diagonais de Cholesky, entao um bloco encostado nele esta confinado a uma face de
  // dimensao menor e nenhuma das coordenadas dele, a fora da diagonal inclusive, que nao
  // tem piso proprio, entra no certificado.
  //
  // TENTADO E DESFEITO no item ar1real: tirar do certificado so a diagonal na parede e
  // cobrar as outras coordenadas do bloco, as que andam sobre a face de posto 1. A ideia
  // pega o ponto da face que nao e o otimo da face (test-fixed-solutions.R, model_mt com
  // h2 = 0, semente 32: decremento 1,7e-4 com o bloco fora num ponto 9e-4 acima do melhor de
  // 30 partidas, 0,039 com a face dentro), mas o score analitico do bloco perto da
  // singularidade e lixo: na semente 33 do mesmo desenho, com L_22 / L_11 = 1e-3 (o piso) e
  // r = +1, o gradiente em (log L_11, L_21) dava (57,6; -117) contra (0,0012; 0,018) das
  // diferencas finitas da forma V densa, e o proprio -2logL esparso estava 2,8e-6 abaixo da
  // forma V. O piso relativo nao limita a condicao do bloco quando a fora da diagonal domina
  // (ali L_11 = 0,049 e L_21 = 0,98: condicao ~1e11), e o MME perde os digitos. Com a face
  // no certificado esse ajuste nunca certificava (decremento 0,02 num ponto 2e-5 acima do
  // melhor) e parava em converged = FALSE. A face fica fora; quem a anda bem e o passo, com
  // a curvatura da face (abaixo), e quem confere outra bacia e a exploracao da fronteira.
  for (std::size_t b = 0; b < mp.blocos.size(); b++) {
    const std::size_t off = mp.blocos[b].first, dim = mp.blocos[b].second;
    if (dim < 2) continue;
    bool preso = false;
    for (std::size_t i = 0; i < dim && !preso; i++)
      if (P.na_parede[idx_tri(off, dim, i, i)]) preso = true;
    if (!preso) continue;
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) P.na_parede[idx_tri(off, dim, i, j)] = 1;
  }
  // O certificado exclui mais que o passo, e por dois motivos distintos.
  //
  // O primeiro e o de sempre: uma coordenada na parede com o score apontando para fora
  // aponta para fora do cone por construcao, o gradiente dela nunca zera, e cobra-la
  // impediria qualquer ponto de fronteira de se certificar.
  //
  // O segundo e numerico e foi MEDIDO. Quando um bloco encosta na singularidade, a
  // diagonal de Cholesky presa tem dtheta/dz = 2 L^2, que colapsa junto com L: num
  // bivariado com o grupo declarado preso em posto 1, L = 4.0e-04 dava dtheta/dz = 3.2e-07
  // e a linha correspondente de AI_z ficava da ordem de 1e-13. O gradiente ali ja tinha
  // zerado (sz = -1e-04, dentro do arredondamento), entao a regra do score para fora NAO
  // a excluia, e inverter uma matriz com essa direcao quase nula devolvia decremento 14.87
  // num ponto onde uma descida por coordenada ainda achava 0.0063 unidade, mil vezes a
  // folga real. Uma direcao sem curvatura em z nao carrega informacao para o certificado:
  // a forma quadratica ali e ruido amplificado, nao gap de verossimilhanca.
  double maior = 0.0;
  for (std::size_t i = 0; i < ntheta; i++) maior = std::max(maior, std::fabs(P.az.at(i, i)));
  if (maior > 0.0)
    for (std::size_t i = 0; i < ntheta; i++)
      if (std::fabs(P.az.at(i, i)) <= 1e-10 * maior) {
        P.na_parede[i] = 1;
        // e o PASSO tambem a congela. Deixa-la no sistema amortecido nao move nada (o
        // Jacobiano dela colapsou) e ainda estraga as outras direcoes, porque a matriz que
        // se inverte para achar o passo fica com a mesma direcao quase nula que estragava
        // o certificado.
        P.congelado[i] = 1;
      }
  // A CURVATURA DA FACE. AI_z = J' AI J e a regra da cadeia de PRIMEIRA ordem; a Hessiana de
  // f(theta(z)) tem ainda soma_k score_k d2 theta_k / dz dz', que no interior some junto com o
  // score e numa fronteira nao: ali o score de theta e o multiplicador da restricao, e o termo
  // e a curvatura da face (C = v v' e uma superficie curva em theta). Sem ele a AI_z
  // subestimava a curvatura de log L_11 sobre a face de posto 1 em 2 a 4 vezes (medido por
  // diferencas finitas no test-fixed-solutions.R, model_mt com h2 = 0, semente 32: 1,23 a
  // 2,36 contra 0,56 a 0,65; com o termo, 1,15 a 1,38), o passo de Newton passava do ponto e
  // o caminhante rastejava pela face com relDelta ~1e-4 por centenas de iteracoes.
  //
  // Por bloco, f = <S, L L'> com S_ii = score_ii e S_ij = score_ij / 2, e para as coordenadas
  // a = (i, m) e b = (j, n) do triangulo de L: d2 f / dz_a dz_b = 2 S_ij D_a D_b quando m = n
  // (D = L_ii na diagonal, em log, e 1 fora dela), mais 2 (S L_.m)_i L_ii na propria diagonal
  // a = b = (i, i). So nos blocos com alguma diagonal na parede (no interior o caminho fica
  // o de sempre), e so a parte positiva-semidefinida: longe do otimo o score pode deixar o
  // termo indefinido, e somado a AI_z ele so pode encurtar o passo, nunca virar a direcao.
  for (std::size_t b = 0; b < mp.blocos.size(); b++) {
    const std::size_t off = mp.blocos[b].first, dim = mp.blocos[b].second;
    bool na = false;
    for (std::size_t i = 0; i < dim; i++)
      if (P.na_parede[idx_tri(off, dim, i, i)]) na = true;
    if (!na) continue;
    Densa S(dim, dim), L(dim, dim);
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const std::size_t k = idx_tri(off, dim, i, j);
        S.at(i, j) = (i == j) ? score[k] : 0.5 * score[k];
        S.at(j, i) = S.at(i, j);
        L.at(i, j) = (i == j) ? std::exp(P.zc[k]) : P.zc[k];
      }
    // as coordenadas do bloco na ordem do triangulo: (linha, coluna) de L
    std::vector<std::size_t> lin, col;
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) { lin.push_back(i); col.push_back(j); }
    const std::size_t nb = lin.size();
    Densa H(nb, nb);
    bool finito = true;
    for (std::size_t a = 0; a < nb; a++)
      for (std::size_t c = 0; c < nb; c++) {
        if (col[a] != col[c]) continue;
        const double da = (lin[a] == col[a]) ? L.at(lin[a], lin[a]) : 1.0;
        const double dc = (lin[c] == col[c]) ? L.at(lin[c], lin[c]) : 1.0;
        double v = 2.0 * S.at(lin[a], lin[c]) * da * dc;
        if (a == c && lin[a] == col[a]) {
          double sl = 0.0;
          for (std::size_t q = 0; q < dim; q++) sl += S.at(lin[a], q) * L.at(q, col[a]);
          v += 2.0 * sl * L.at(lin[a], lin[a]);
        }
        if (!std::isfinite(v)) finito = false;
        H.at(a, c) = v;
      }
    if (!finito) continue;
    std::vector<double> ev;
    Densa u;
    jacobi_sim(H, ev, u);
    for (std::size_t a = 0; a < nb; a++)
      for (std::size_t c = 0; c < nb; c++) {
        double v = 0.0;
        for (std::size_t q = 0; q < nb; q++)
          if (ev[q] > 0.0) v += u.at(a, q) * ev[q] * u.at(c, q);
        P.az.at(idx_tri(off, dim, lin[a], col[a]), idx_tri(off, dim, lin[c], col[c])) += v;
      }
  }
  return P;
}

// O decremento de Newton dos componentes LIVRES, g' AI^-1 g nas coordenadas z em que o
// passo anda e em que os conjuntos ativos estao definidos. Um componente na parede aponta
// para fora do cone por construcao e o gradiente dele nunca zera; cobra-lo impediria
// qualquer ponto de fronteira de se certificar. O que se certifica e o otimo CONDICIONAL
// ao que esta preso.
std::string lista_na_parede(const std::vector<char>& parede,
                            const std::vector<std::string>& nomes) {
  std::string s;
  std::size_t n = 0, mostra = 0;
  for (std::size_t k = 0; k < parede.size() && k < nomes.size(); k++) {
    if (!parede[k]) continue;
    n++;
    if (mostra < 8) {
      s += (mostra ? ", " : "") + nomes[k];
      mostra++;
    }
  }
  if (n > mostra) s += " and " + std::to_string(n - mostra) + " more";
  return s;
}

double decremento_z(const PassoZ& P) {
  if (!P.ok) return std::numeric_limits<double>::infinity();
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
  // A forma quadratica por PSEUDO-INVERSA TRUNCADA, e nao por inversa. Uma AI_z com uma
  // direcao mal condicionada nao e uma excecao rara aqui: e o estado normal perto de uma
  // fronteira de covariancia, onde dtheta/dz encolhe com a diagonal de Cholesky presa.
  // Invertendo cheio, o decremento devolvia 0.0339 num ponto que uma descida por
  // coordenada melhorava em 1.3e-05, duas mil vezes a folga real, e o certificado
  // reprovava um ajuste que estava no otimo. Truncar em 1e-8 do maior autovalor mede o gap
  // nas direcoes que existem e descarta as que nao carregam curvatura. E a mesma escolha
  // ja feita para o Gamma dos metafundadores.
  std::vector<double> ev;
  Densa u;
  jacobi_sim(m, ev, u);
  double maior_ev = 0.0;
  for (double e : ev) maior_ev = std::max(maior_ev, std::fabs(e));
  if (!(maior_ev > 0.0)) return 0.0;
  const double corte = 1e-8 * maior_ev;
  double dec = 0.0;
  for (std::size_t k = 0; k < ev.size(); k++) {
    if (ev[k] <= corte) continue;
    double proj = 0.0;
    for (std::size_t i = 0; i < sc.size(); i++) proj += u.at(i, k) * sc[i];
    dec += proj * proj / ev[k];
  }
  return std::fabs(dec);   // a AI e PSD; negativo aqui e arredondamento, nao curvatura
}

// REENTRADA NO CONE, a condicao de KKT conferida pela VEROSSIMILHANCA, para os
// ajustadores que caminham com passo_z (model_mt(), model_ar1()). Medido no P1-8 (model_mt,
// pares diarios reais): o ajuste saia com converged = TRUE numa fronteira de covariancia
// 46,8 e 16,3 unidades de -2logL ACIMA do otimo interior. Duas pecas faziam isso juntas. O
// certificado exclui da conta todo componente na parede sem olhar o sinal do score (passo_z,
// acima: a exigencia de sinal reprovava otimos restritos verdadeiros, porque perto da
// singularidade o score em z e ruido), e o passo congela a coordenada cuja curvatura em z
// colapsou (dtheta/dz = 2 L^2 vai a zero com L), tambem sem olhar o sinal. Resultado: um
// bloco parado perto da parede com o gradiente apontando PARA DENTRO do cone nunca era
// empurrado de volta, e o certificado o ignorava.
//
// A condicao de KKT de min f(C) com C no cone PSD e G >= 0 e G C = 0, com G a matriz do
// gradiente de f = -2logL no bloco (G_ii = score_ii, G_ij = score_ij / 2, pois C_ij = C_ji
// e um componente so). Um autovalor negativo de G, com autovetor w, e uma direcao C + e ww'
// que fica no cone e baixa o -2logL em primeira ordem: o ponto nao era otimo. Em z essa
// direcao some, porque o Jacobiano da diagonal presa colapsou; em theta ela continua la. O
// score perto da singularidade carrega ruido, entao o sinal dele nao decide sozinho: cada
// direcao vira um ponto, com o comprimento do modelo quadratico da AI (e = -lambda / d'AI d,
// d = vech(ww'), e metade e o dobro dele), e quem decide e o -2logL avaliado.
//
// Medido (test-unidade-partida-fronteira.R, P8-1, seed 12, a variancia da ninhada de y2
// partindo de 1e-14): o bloco parou em r = -0,9999 com score -2,25 na var(random@y2),
// apontando para dentro, e o ajuste certificava 0,041 unidade acima do otimo interior. Os
// candidatos por coordenada pioravam todos, de 0,003 a 79 unidades, porque a menor fracao
// da escala (1e-3 de 1, a escala da caracteristica padronizada) ja e 2,4 vezes a variancia
// inteira de y2 no bloco; o primeiro ponto direcional baixou 0,009, e dali o caminhante chega ao otimo.
// Os candidatos por coordenada ficam como reserva, para quando o gradiente nao aponta
// nada: a variancia condicional L_ii^2 em 1e-3, 1e-2 e 1e-1 da escala da coordenada (a
// variancia marginal C_ii dela ou, na parede do zero, a escala residual), so quando isso ao
// menos dobra L_ii, e o bloco inteiro descorrelacionado com as variancias no mesmo piso
// (levantar so uma diagonal mantem a fora da diagonal de Cholesky, e num bloco de posto 1
// isso mantem a correlacao em +/-1). E uma sonda, nao uma prova: um otimo restrito que so
// melhora longe dos candidatos passa.
CandidatosKKT candidatos_reentrada(const MapaZ& mp, std::size_t ntheta,
                                   const std::vector<double>& z,
                                   const std::vector<char>& parede, double escala,
                                   const std::vector<double>& score, const Densa& ai) {
  CandidatosKKT out;
  if (z.size() != ntheta || parede.size() != ntheta) return out;
  const std::vector<double> theta = theta_de_z(mp, ntheta, z);
  const bool com_gradiente = score.size() == ntheta && ai.nlin == ntheta && ai.ncol == ntheta;
  for (std::size_t b = 0; b < mp.blocos.size(); b++) {
    const std::size_t off = mp.blocos[b].first, dim = mp.blocos[b].second;
    bool bloco_na_parede = false;
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++)
        if (parede[idx_tri(off, dim, i, j)]) bloco_na_parede = true;
    if (!bloco_na_parede) continue;
    // as direcionais: autovetores de autovalor negativo de G
    if (com_gradiente) {
      Densa G(dim, dim);
      bool finito = true;
      for (std::size_t j = 0; j < dim; j++)
        for (std::size_t i = j; i < dim; i++) {
          const double s = score[idx_tri(off, dim, i, j)];
          if (!std::isfinite(s)) finito = false;
          G.at(i, j) = (i == j) ? s : 0.5 * s;
          G.at(j, i) = G.at(i, j);
        }
      if (finito) {
        std::vector<double> ev;
        Densa u;
        jacobi_sim(G, ev, u);
        for (std::size_t k = 0; k < ev.size(); k++) {
          if (!(ev[k] < 0.0)) continue;
          std::vector<double> dv(ntheta, 0.0);
          for (std::size_t j = 0; j < dim; j++)
            for (std::size_t i = j; i < dim; i++)
              dv[idx_tri(off, dim, i, j)] = u.at(i, k) * u.at(j, k);
          double q = 0.0;
          for (std::size_t a = off; a < off + dim * (dim + 1) / 2; a++)
            for (std::size_t c = off; c < off + dim * (dim + 1) / 2; c++)
              q += dv[a] * ai.at(a, c) * dv[c];
          if (!(q > 0.0) || !std::isfinite(q)) continue;
          const double e = -ev[k] / q;
          for (const double m : {1.0, 0.5, 2.0}) {
            std::vector<double> th = theta;
            for (std::size_t a = 0; a < ntheta; a++) th[a] += m * e * dv[a];
            std::vector<double> zc;
            if (z_de_theta(mp, ntheta, th, zc)) out.direcionais.push_back(zc);
          }
        }
      }
    }
    // as por coordenada
    Densa L(dim, dim);
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const double v = z[idx_tri(off, dim, i, j)];
        L.at(i, j) = (i == j) ? std::exp(v) : v;
      }
    std::vector<double> cii(dim, 0.0);
    for (std::size_t i = 0; i < dim; i++)
      for (std::size_t m = 0; m <= i; m++) cii[i] += L.at(i, m) * L.at(i, m);
    for (std::size_t i = 0; i < dim; i++) {
      const std::size_t k = idx_tri(off, dim, i, i);
      if (!parede[k]) continue;
      const double ref = std::max(cii[i], escala);
      const double v = L.at(i, i) * L.at(i, i);
      for (const double fr : {1e-3, 1e-2, 1e-1}) {
        const double vn = fr * ref;
        if (!(vn > 4.0 * v) || !std::isfinite(vn)) continue;
        std::vector<double> zc = z;
        zc[k] = 0.5 * std::log(vn);
        out.coordenadas.push_back(zc);
      }
    }
    if (dim > 1)
      for (const double fr : {1e-3, 1e-2, 1e-1}) {
        std::vector<double> zc = z;
        for (std::size_t j = 0; j < dim; j++)
          for (std::size_t i = j; i < dim; i++)
            zc[idx_tri(off, dim, i, j)] =
                (i == j) ? 0.5 * std::log(std::max(cii[i], fr * escala)) : 0.0;
        bool finito = true;
        for (double v : zc) if (!std::isfinite(v)) finito = false;
        if (finito) out.coordenadas.push_back(zc);
      }
  }
  return out;
}

// A SONDA DA FACE, pela verossimilhanca. O certificado tira da conta o bloco inteiro que
// encosta na parede (passo_z), porque o score analitico dele perto da singularidade e lixo;
// mas um ponto da face que nao e o otimo da face passava junto. Medido no
// test-fixed-solutions.R (model_mt com h2 = 0, semente 32): a partida dos univariados
// certificava 9e-4 unidade de -2logL acima do melhor de 30 partidas, com o bloco genetico em
// r = -1, e das 30 partidas aleatorias 10 certificavam mais de 1e-4 acima dele.
//
// Aqui a face e conferida sem o score do bloco. As coordenadas do bloco que andam sobre a
// face (todas menos as diagonais na parede e as de curvatura colapsada) entram por
// diferencas finitas centrais do -2logL, com passo h = 1e-2 em z: o gradiente, a Hessiana
// entre elas e, pelo score das OUTRAS coordenadas avaliado nos mesmos pontos, a Hessiana
// cruzada com o resto do modelo. As outras coordenadas livres entram com o score e a AI_z
// de sempre. O passo de Newton CONJUNTO desse modelo quadratico (comprimentos 1, 1/2 e 1/4)
// e os proprios pontos da sonda sao avaliados, e volta o melhor deles quando algum baixa o
// -2logL; quem chama decide pelo ganho: acima de tol_reentrada o certificado e recusado e o
// caminhante vai para ele, e um certificado aceito ainda e lapidado por ela enquanto cada
// passo ganha mais que tol_lapida (o racional esta em mme.h). O conjunto importa: com
// h2 ~ 0 a covariancia genetica e a residual trocam variancia ao longo de um vale estreito,
// e uma sonda so da face, com o resto parado, deixava 9 das 30 partidas certificarem mais de
// 1e-4 acima do melhor (o otimo condicional de cada bloco nao e o otimo conjunto).
//
// O passo h = 1e-2 e o que o ruido permite: no canto da semente 33 do mesmo desenho o
// -2logL esparso erra ~1e-6 (2,8e-6 abaixo da forma V densa) e com h = 1e-5 a diferenca
// finita ja era lixo; com 1e-2 o erro do gradiente fica perto de 1e-4. Custa 2n + n(n - 1)/2
// + 3 avaliacoes, com n as coordenadas da face, e so quando o resto do certificado ja passou.
// E uma sonda, nao uma prova: um otimo que so melhora longe do modelo quadratico passa.
bool sonda_face(const MapaZ& mp, std::size_t ntheta, const PassoZ& P, double f0,
                const std::function<bool(const std::vector<double>&, double&,
                                         std::vector<double>&)>& f,
                std::vector<double>& z_novo, double& f_novo) {
  if (!P.ok || P.zc.size() != ntheta || !std::isfinite(f0)) return false;
  const double h = 1e-2;
  // as coordenadas da face (F) e as livres do certificado (Rl)
  std::vector<std::size_t> F, Rl;
  for (std::size_t b = 0; b < mp.blocos.size(); b++) {
    const std::size_t off = mp.blocos[b].first, dim = mp.blocos[b].second;
    if (dim < 2) continue;
    bool na = false;
    for (std::size_t i = 0; i < dim; i++)
      if (P.zc[idx_tri(off, dim, i, i)] <= P.piso[idx_tri(off, dim, i, i)] + std::log(2.0))
        na = true;
    if (!na) continue;
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const std::size_t k = idx_tri(off, dim, i, j);
        if (P.congelado[k]) continue;
        if (i == j && P.zc[k] <= P.piso[k] + std::log(2.0)) continue;
        F.push_back(k);
      }
  }
  if (F.empty()) return false;
  for (std::size_t k = 0; k < ntheta; k++)
    if (!P.na_parede[k]) Rl.push_back(k);
  const std::size_t nf = F.size(), nr = Rl.size(), n = nf + nr;
  bool achou = false;
  f_novo = f0;
  std::vector<double> gz;
  auto prova = [&](const std::vector<double>& zc, double& v) -> bool {
    std::vector<double> zz = zc;
    for (std::size_t i = 0; i < ntheta; i++)
      if (zz[i] < P.piso[i]) zz[i] = P.piso[i];
    if (!f(zz, v, gz) || !std::isfinite(v)) return false;
    if (v < f_novo) {
      f_novo = v;
      z_novo = zz;
      achou = true;
    }
    return true;
  };
  std::vector<double> g(n, 0.0), mais(nf, 0.0);
  Densa H(n, n);
  for (std::size_t a = 0; a < nf; a++) {
    std::vector<double> zc = P.zc;
    double vm = 0.0, vn = 0.0;
    zc[F[a]] += h;
    if (!prova(zc, vm)) return achou;
    const std::vector<double> gm = gz;
    zc[F[a]] -= 2.0 * h;
    if (!prova(zc, vn)) return achou;
    mais[a] = vm;
    g[a] = (vm - vn) / (2.0 * h);
    H.at(a, a) = (vm - 2.0 * f0 + vn) / (h * h);
    for (std::size_t c = 0; c < nr; c++) {
      H.at(nf + c, a) = (gm[Rl[c]] - gz[Rl[c]]) / (2.0 * h);
      H.at(a, nf + c) = H.at(nf + c, a);
    }
  }
  for (std::size_t a = 0; a < nf; a++)
    for (std::size_t c = a + 1; c < nf; c++) {
      std::vector<double> zc = P.zc;
      double v = 0.0;
      zc[F[a]] += h;
      zc[F[c]] += h;
      if (!prova(zc, v)) return achou;
      H.at(a, c) = (v - mais[a] - mais[c] + f0) / (h * h);
      H.at(c, a) = H.at(a, c);
    }
  for (std::size_t c = 0; c < nr; c++) {
    g[nf + c] = P.sz[Rl[c]];
    for (std::size_t e = 0; e < nr; e++) H.at(nf + c, nf + e) = P.az.at(Rl[c], Rl[e]);
  }
  // o passo de Newton conjunto, com a Hessiana levada a positiva-definida: autovalores em
  // modulo, e os abaixo de 1e-3 do maior levantados a ele (uma direcao plana nao vira um
  // passo infinito)
  std::vector<double> ev;
  Densa u;
  jacobi_sim(H, ev, u);
  double maior = 0.0;
  for (double e : ev) maior = std::max(maior, std::fabs(e));
  if (!(maior > 0.0) || !std::isfinite(maior)) return achou;
  std::vector<double> passo(n, 0.0);
  for (std::size_t q = 0; q < n; q++) {
    const double e = std::max(std::fabs(ev[q]), 1e-3 * maior);
    double proj = 0.0;
    for (std::size_t a = 0; a < n; a++) proj += u.at(a, q) * g[a];
    for (std::size_t a = 0; a < n; a++) passo[a] -= u.at(a, q) * proj / e;
  }
  for (const double alpha : {1.0, 0.5, 0.25}) {
    std::vector<double> zc = P.zc;
    for (std::size_t a = 0; a < nf; a++) zc[F[a]] += alpha * passo[a];
    for (std::size_t c = 0; c < nr; c++) zc[Rl[c]] += alpha * passo[nf + c];
    double v = 0.0;
    prova(zc, v);
  }
  return achou;
}

// O PONTO DESCORRELACIONADO NA PAREDE, o outro candidato da lapidacao (tol_lapida, mme.h).
// Uma caracteristica cuja variancia MARGINAL no bloco esta no piso (a diagonal de Cholesky na
// parede e C_ii ate 4 vezes a variancia do piso) esta na fronteira do zero, e ali a
// correlacao dela com as outras nao tem sentido; mas o piso relativo nao e o zero, e no piso a
// correlacao ainda mexe no -2logL. Medido (o desenho 2, semente 1, da busca do D3, model_ar1
// bivariado com var(animal@y1) no piso): a partida de var(y) certificava r = -0,79 com
// var(animal@y1) = 1,0e-6, 6,7e-4 unidade de -2logL acima do ponto da partida dos
// univariados; no mesmo ponto, com a covariancia zerada, o -2logL baixa 7,3e-4. A sonda da
// face nao acha essa direcao: no piso, L_21 carrega quase toda a var(animal@y2), o vale que
// gira L_21 para zero com a variancia marginal fixa e um arco, e nenhum dos pontos das
// diferencas finitas nem o passo de Newton delas melhorava. Aqui o candidato e o proprio
// ponto do limite: as covariancias dessas caracteristicas em zero, TODAS as variancias onde
// estao, de volta a z e aos pisos. Vazio quando nenhum bloco tem caracteristica assim (a face
// r = +-1, com a diagonal condicional no piso e a variancia marginal longe dele, nao entra).
std::vector<double> descorrelada_na_parede(const MapaZ& mp, std::size_t ntheta,
                                           const PassoZ& P) {
  if (!P.ok || P.zc.size() != ntheta || P.piso.size() != ntheta) return {};
  std::vector<double> th = theta_de_z(mp, ntheta, P.zc);
  bool mudou = false;
  for (std::size_t b = 0; b < mp.blocos.size(); b++) {
    const std::size_t off = mp.blocos[b].first, dim = mp.blocos[b].second;
    if (dim < 2) continue;
    for (std::size_t i = 0; i < dim; i++) {
      const std::size_t k = idx_tri(off, dim, i, i);
      if (!(P.zc[k] <= P.piso[k] + std::log(2.0))) continue;
      if (!(th[k] <= 4.0 * std::exp(2.0 * P.piso[k]))) continue;
      for (std::size_t j = 0; j < dim; j++)
        if (j != i && th[idx_tri(off, dim, i, j)] != 0.0) {
          th[idx_tri(off, dim, i, j)] = 0.0;
          mudou = true;
        }
    }
  }
  if (!mudou) return {};
  std::vector<double> z;
  if (!z_de_theta(mp, ntheta, th, z)) return {};
  for (std::size_t i = 0; i < ntheta; i++)
    if (z[i] < P.piso[i]) z[i] = P.piso[i];
  return z;
}

// A EXPLORACAO DA FRONTEIRA, a segunda metade do conserto do P8-1. A sonda de KKT acima e
// LOCAL: ela avalia pontos perto do bloco na parede, com o score do ponto, e o score de um
// bloco quase singular nao tem sentido. Medido (a revisao do item ar1real, sim_ninhadas
// seed 3 de test-unidade-partida-fronteira.R, var(random@y2) partindo de 1e-12 com a
// covariancia zero): converged = TRUE em 13 iteracoes, decremento 9,8e-8, 4,73 unidades de
// -2logL acima do otimo interior, com o bloco da ninhada em r = -1 e var(random@y1) =
// 1,6e-9 contra 0,233 e r = +0,76 no otimo. A diagonal do gradiente do bloco ali e 4,3e8, o
// unico autovalor negativo dele preve ganho de 2e-5 (abaixo da folga), e os candidatos
// descorrelacionados da sonda pioram o -2logL NO PONTO (+0,094, +2,77 e +62); mas o segmento
// ate o otimo ja baixa 0,0079 a 0,1% do caminho e desce monotonico ate -4,73. Era outra
// bacia: o resto do modelo (var(animal@y1) 1,84 contra 1,11) tinha se ajustado ao bloco
// na parede. Numa varredura de 48 partidas assim (seeds 1 a 12, cov da ninhada +-0,25, uma
// variancia da ninhada em 1e-12), 6 certificados falsos, de 0,0056 a 4,73 unidades.
//
// Nenhum teste local ve outra bacia; uma CAMINHADA a ve. Quando um ajuste termina com algum
// bloco de covariancia na parede (certificado, ou parado sem certificar antes do teto de
// iteracoes), uma segunda caminhada parte de dentro do cone: cada bloco com coordenada na
// parede vira diagonal, com cada variancia em max(C_ii, 10% da partida de var(y)), e o resto
// fica onde a primeira deixou. Medido no caso acima: da partida assim, com 1%, 10%, 50% e
// 100% da partida de var(y), as quatro caminhadas certificam o otimo interior em 7 a 9
// iteracoes. Volta o melhor dos dois. A segunda vence quando baixa o -2logL por mais que
// tol_reentrada (1e-4: a primeira, certificada ou nao, nao era o otimo); quando certifica e
// fica abaixo da primeira por qualquer margem; quando certifica a menos de 2e-4 de uma
// primeira que nao certificou (a folga do proprio certificado); e, sem certificado nas duas,
// quando fica abaixo. So uma primeira certificada contra uma segunda sem certificado e
// pouco abaixo dela fica com a primeira. A regra anterior (a segunda so vencia por mais que
// 1e-4) devolvia o pior de dois certificados: no test-mf-genomico.R a segunda caminhada
// certificava 2,49e-5 abaixo e voltava a primeira, 2,49e-5 acima do base. Senao volta a
// primeira, com a segunda dita na mensagem. Ate duas rodadas, para a segunda caminhada que
// tambem termina na parede, e a segunda rodada so depois de uma vitoria por mais que a folga
// (o racional esta no fim de explora_fronteira). O
// orcamento e o mesmo maxiter: a segunda conta a partir da primeira, `iters` e o total, e
// sem iteracao sobrando nao ha exploracao. E uma busca, nao uma prova: um otimo melhor que
// nenhuma das duas caminhadas alcanca passa, e o custo e uma caminhada a mais em todo ajuste
// que termina na parede.
static std::vector<double> partida_reentrada(const MapaZ& mp, const std::vector<double>& theta,
                                             const std::vector<char>& parede,
                                             const std::vector<double>& auto_) {
  std::vector<double> th = theta;
  bool algum = false;
  for (std::size_t b = 0; b < mp.blocos.size(); b++) {
    const std::size_t off = mp.blocos[b].first, dim = mp.blocos[b].second;
    bool na = false;
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++)
        if (parede[idx_tri(off, dim, i, j)]) na = true;
    if (!na) continue;
    algum = true;
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const std::size_t k = idx_tri(off, dim, i, j);
        th[k] = (i == j) ? std::max(theta[k], 0.1 * auto_[k]) : 0.0;
      }
  }
  if (!algum) return {};
  for (double v : th)
    if (!std::isfinite(v)) return {};
  return th;
}

AjusteMT explora_fronteira(AjusteMT R, const MapaZ& mp, const std::vector<double>& auto_,
                           std::size_t maxiter,
                           const std::function<AjusteMT(const std::vector<double>&,
                                                        std::size_t)>& caminha,
                           bool verboso) {
  const char* como = "a second walk from inside the cone (every block at a covariance "
                     "boundary made diagonal, its variances at least 10% of the var(y) start)";
  for (int rodada = 0; rodada < 2; rodada++) {
    if (!(R.convergiu || R.parou) || R.iters >= maxiter || !std::isfinite(R.neg2logl))
      return R;
    if (R.theta.size() != auto_.size() || R.na_parede.size() != auto_.size()) return R;
    const std::vector<double> th = partida_reentrada(mp, R.theta, R.na_parede, auto_);
    if (th.empty()) return R;
    if (verboso)
      Rprintf("boundary check: %s, from iteration %d\n", como, (int) R.iters);
    AjusteMT E = caminha(th, R.iters);
    // negativo: a segunda caminhada baixou o -2logL
    const double gap = E.neg2logl - R.neg2logl;
    const bool vence = std::isfinite(E.neg2logl) &&
                       (gap < -tol_reentrada ||
                        (E.convergiu && (gap < 0.0 || (!R.convergiu && gap < 2e-4))) ||
                        (!R.convergiu && !E.convergiu && gap < 0.0));
    char buf[480];
    if (!vence) {
      R.iters = std::max(R.iters, E.iters);
      if (std::isfinite(E.neg2logl))
        std::snprintf(buf, sizeof(buf), "%s ended %.3g -2logL unit(s) %s%s", como,
                      std::fabs(gap), gap >= 0.0 ? "above this point" : "below this point, "
                      "inside the 1e-4 slack", E.convergiu ? "" : " without a certificate");
      else
        std::snprintf(buf, sizeof(buf), "%s could not be evaluated", como);
      R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") + buf;
      return R;
    }
    if (R.convergiu && gap < -tol_reentrada)
      std::snprintf(buf, sizeof(buf), "a first walk certified a covariance boundary at -2logL "
                    "%.6f that is not the optimum: %s %s this point, %.3g lower", R.neg2logl,
                    como, E.convergiu ? "certified" : "reached, without a certificate,",
                    -gap);
    else if (R.convergiu)
      std::snprintf(buf, sizeof(buf), "a first walk certified a covariance boundary at -2logL "
                    "%.6f: %s certified this point, %.3g lower (inside the 1e-4 slack, the "
                    "lower of the two certified points is returned)", R.neg2logl, como, -gap);
    else
      std::snprintf(buf, sizeof(buf), "a first walk stopped without a certificate at a "
                    "covariance boundary at -2logL %.6f: %s %s this point (%+.3g)", R.neg2logl,
                    como, E.convergiu ? "certified" : "reached, also without a certificate,",
                    gap);
    E.mensagem += std::string(E.mensagem.empty() ? "" : "; ") + buf;
    E.partida = R.partida;
    R = std::move(E);
    // outra rodada so quando a segunda caminhada mostrou que a primeira NAO era o otimo (mais
    // que a folga abaixo): a segunda tambem pode ter parado numa parede errada. Uma vitoria
    // dentro da folga encerra. Medido no caso de fronteira de validation/
    // unit_shift_start_scale.R (24 mil registros, ninhada de y2 em zero): a segunda caminhada
    // certificava 4,9e-8 abaixo da primeira, a terceira ficava 9e-8 acima, e as tres custavam
    // 23 iteracoes e 94 s contra as 8 da primeira.
    if (!(gap < -tol_reentrada)) return R;
  }
  return R;
}

// O mapa do univariado: cada grupo e um bloco, o residuo e o escalar em log. Escrito assim
// para que o passo do univariado continue IDENTICO ao que sempre foi.
static MapaZ mapa_z(const Modelo& mo) {
  MapaZ mp;
  for (const Grupo& gr : mo.grupos)
    mp.blocos.push_back(std::make_pair(gr.offset, gr.dim));
  for (std::size_t k = 0; k < mo.n_residuos(); k++) mp.escalares_log.push_back(mo.offset_residual + k);
  return mp;
}

static bool z_de_theta(const Modelo& mo, const std::vector<double>& theta,
                       std::vector<double>& z) {
  return z_de_theta(mapa_z(mo), mo.ntheta, theta, z);
}

std::vector<double> theta_de_z(const MapaZ& mp, std::size_t ntheta,
                               const std::vector<double>& z) {
  std::vector<double> theta(ntheta, 0.0);
  for (std::size_t b = 0; b < mp.blocos.size(); b++) {
    const std::size_t off = mp.blocos[b].first, dim = mp.blocos[b].second;
    Densa L(dim, dim);
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const double v = z[idx_tri(off, dim, i, j)];
        L.at(i, j) = (i == j) ? std::exp(v) : v;
      }
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        double s = 0.0;
        for (std::size_t k = 0; k <= j; k++) s += L.at(i, k) * L.at(j, k);
        theta[idx_tri(off, dim, i, j)] = s;
      }
  }
  for (std::size_t k : mp.escalares_log) theta[k] = std::exp(z[k]);
  for (std::size_t k : mp.correlacoes) theta[k] = std::tanh(z[k]);
  return theta;
}

static std::vector<double> theta_de_z(const Modelo& mo, const std::vector<double>& z) {
  return theta_de_z(mapa_z(mo), mo.ntheta, z);
}

// J = dtheta/dz, block diagonal by group. For z_m = L entry (a,b) with chain factor
// s_m (= L_aa on the diagonal, 1 below), and theta_k = C_ij (i >= j):
//   dC_ij/dz_m = s_m (delta_ia L_jb + delta_ja L_ib)
bool vcov_condicional(const Densa& ai, const Densa& J, const std::vector<std::size_t>& fz,
                      const std::vector<std::size_t>& livre, std::vector<double>& vcov,
                      std::vector<double>& se) {
  const std::size_t nth = J.nlin, nf = fz.size();
  Densa aj(nth, nf);   // AI J_F
  for (std::size_t i = 0; i < nth; i++)
    for (std::size_t b = 0; b < nf; b++) {
      double s = 0.0;
      for (std::size_t j = 0; j < nth; j++) s += ai.at(i, j) * J.at(j, fz[b]);
      aj.at(i, b) = s;
    }
  Densa A(nf, nf);     // J_F' AI J_F
  for (std::size_t a = 0; a < nf; a++)
    for (std::size_t b = 0; b < nf; b++) {
      double s = 0.0;
      for (std::size_t i = 0; i < nth; i++) s += J.at(i, fz[a]) * aj.at(i, b);
      A.at(a, b) = s;
    }
  Densa inv;
  try { inv = inv_geral(A); } catch (const Erro&) { return false; }
  Densa ji(nth, nf);   // J_F (J_F' AI J_F)^-1
  for (std::size_t i = 0; i < nth; i++)
    for (std::size_t b = 0; b < nf; b++) {
      double s = 0.0;
      for (std::size_t c = 0; c < nf; c++) s += J.at(i, fz[c]) * inv.at(c, b);
      ji.at(i, b) = s;
    }
  for (const std::size_t i : livre) {
    for (const std::size_t k : livre) {
      double s = 0.0;
      for (std::size_t b = 0; b < nf; b++) s += ji.at(i, b) * J.at(k, fz[b]);
      vcov[i * nth + k] = 2.0 * s;
    }
    const double v = vcov[i * nth + i];
    if (v > 0.0) se[i] = std::sqrt(v);
  }
  return true;
}

Densa jacobiano_z(const MapaZ& mp, std::size_t ntheta, const std::vector<double>& z) {
  Densa J(ntheta, ntheta);
  for (std::size_t bl = 0; bl < mp.blocos.size(); bl++) {
    const std::size_t off = mp.blocos[bl].first, dim = mp.blocos[bl].second;
    Densa L(dim, dim);
    for (std::size_t j = 0; j < dim; j++)
      for (std::size_t i = j; i < dim; i++) {
        const double v = z[idx_tri(off, dim, i, j)];
        L.at(i, j) = (i == j) ? std::exp(v) : v;
      }
    for (std::size_t b = 0; b < dim; b++)
      for (std::size_t a = b; a < dim; a++) {
        const std::size_t m = idx_tri(off, dim, a, b);
        const double sm = (a == b) ? L.at(a, a) : 1.0;
        for (std::size_t j = 0; j < dim; j++)
          for (std::size_t i = j; i < dim; i++) {
            double v = 0.0;
            if (i == a && b <= j) v += L.at(j, b);
            if (j == a && b <= i) v += L.at(i, b);
            J.at(idx_tri(off, dim, i, j), m) = sm * v;
          }
      }
  }
  for (std::size_t k : mp.escalares_log) J.at(k, k) = std::exp(z[k]);
  // d tanh(z)/dz = 1 - tanh(z)^2
  for (std::size_t k : mp.correlacoes) {
    const double r = std::tanh(z[k]);
    J.at(k, k) = 1.0 - r * r;
  }
  return J;
}

static Densa jacobiano_z(const Modelo& mo, const std::vector<double>& z) {
  return jacobiano_z(mapa_z(mo), mo.ntheta, z);
}

// O residuo, os ajustados e a alavanca no fator do bloco final (o mesmo theta da PEV). A
// solucao sai DESTE fator, e nao da ultima avaliacao: residuo e alavanca ficam no mesmo
// theta por construcao. Tudo na escala escalada por sqrt(w); quem chama divide.
static void residuos_no_fator(const Desenho& d, const Montado& M, const Csc& L,
                              const SelInv& z, const std::vector<std::size_t>& perm,
                              Ajuste& R) {
  const std::size_t tot = M.total;
  std::vector<double> pb(tot);
  for (std::size_t k = 0; k < tot; k++) pb[k] = M.rhs[perm[k]];
  const std::vector<double> px = resolve(L, pb);
  std::vector<double> s(tot);
  for (std::size_t k = 0; k < tot; k++) s[perm[k]] = px[k];
  std::vector<std::size_t> inv_perm(tot);
  for (std::size_t k = 0; k < tot; k++) inv_perm[perm[k]] = k;
  const LinhasW W = linhas_w(d, M.offset_grupo);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  R.residuo.assign(d.nlin, nan);
  R.ajustado.assign(d.nlin, nan);
  R.ajustado_fixo.assign(d.nlin, nan);
  R.alavanca.assign(d.nlin, nan);
  R.alavanca_fora = 0;
  for (std::size_t i = 0; i < d.nlin; i++) {
    if (!d.usa[i]) continue;
    double fa = 0.0, fx = 0.0;
    for (std::size_t p = W.ptr[i]; p < W.ptr[i + 1]; p++) {
      const double v = W.val[p] * s[W.col[p]];
      fa += v;
      if (W.col[p] < d.x.ncol) fx += v;
    }
    R.ajustado[i] = fa;
    R.ajustado_fixo[i] = fx;
    R.residuo[i] = d.y[i] - fa;
    // h_ii = soma_{j,k} w_ij w_ik C^-1[j,k], cada par fora da diagonal duas vezes
    double h = 0.0;
    for (std::size_t p = W.ptr[i]; p < W.ptr[i + 1]; p++)
      for (std::size_t q = p; q < W.ptr[i + 1]; q++) {
        double c;
        if (!z.get(inv_perm[W.col[p]], inv_perm[W.col[q]], c)) {
          R.alavanca_fora++;
          continue;
        }
        h += (p == q ? 1.0 : 2.0) * W.val[p] * W.val[q] * c;
      }
    R.alavanca[i] = h;
  }
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
  for (std::size_t k = 0; k < d.modelo.n_residuos(); k++)
    theta[d.modelo.offset_residual + k] = 0.5 * var;
  return theta;
}

Ajuste ajusta(const Desenho& d, const std::vector<double>* theta0,
              std::size_t n_em, std::size_t maxiter, double tol, bool verboso) {
  Ajuste R;
  std::vector<double> theta = theta0 ? *theta0 : partida(d);
  // o componente preso comeca NO VALOR PEDIDO. Sem isto ele congelaria onde a partida o
  // deixou, que e o oposto do que fixed= promete.
  for (std::size_t k = 0; k < d.modelo.ntheta; k++)
    if (d.modelo.preso(k)) theta[k] = d.modelo.theta_fixo[k];
  if (theta.size() != d.modelo.ntheta) {
    R.mensagem = "start= has the wrong number of components for this model";
    return R;
  }
  // A PARTIDA DO RESIDUO ABAIXO DO PISO (aberto de 2026-10-10): um residuo de partida positivo
  // abaixo de 1e-10 var(y), o piso da caminhada (escala_piso, abaixo), parte do residuo da
  // partida automatica, e a mensagem diz. Medido (SIMULACAO, 300 animais): com var(animal) 0,3
  // e var(residual) 1e-14 na partida, as MME nao fatoravam e o ajuste parava em "the starting
  // theta is INADMISSIBLE", que era falso (o 1e-14 e positivo). Subir so ate o piso era pior:
  // a caminhada parava no canto do piso com converged = TRUE, 91 unidades de -2logL acima do
  // otimo
  std::string nota_partida;
  {
    double soma = 0.0, soma2 = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < d.nlin; i++)
      if (d.usa[i]) { soma += d.y[i]; soma2 += d.y[i] * d.y[i]; n++; }
    double v = n > 1 ? (soma2 - soma * soma / static_cast<double>(n)) / static_cast<double>(n - 1) : 0.0;
    if (!(v > 0.0) || !std::isfinite(v)) v = 1.0;
    const std::vector<double> auto_ = partida(d);
    for (std::size_t k = 0; k < d.modelo.n_residuos(); k++) {
      const std::size_t kr = d.modelo.offset_residual + k;
      if (!d.modelo.preso(kr) && theta[kr] > 0.0 && theta[kr] < 1e-10 * v) {
        char buf[200];
        std::snprintf(buf, sizeof buf, "the starting residual variance %.3g was below its floor "
                      "(1e-10 var(y)), and the fit started it from %.4g instead", theta[kr],
                      auto_[kr]);
        if (nota_partida.empty()) nota_partida = buf;
        theta[kr] = auto_[kr];
      }
    }
  }

  CacheSimbolica cs;
  Avaliacao cur = avalia(d, theta, &cs);
  // O BLOCO DENSO FINAL DO FATOR, dito na hora. Ele e o k do custo k^3 de cada iteracao, e o
  // caso que motivou isto e silencioso: num passo unico sem apy_core= o bloco genotipado e
  // um clique e k e o numero de genotipados, medido 400 -> 400, 800 -> 800, 1600 -> 1600,
  // contra nucleo + 1 com a APY. Quem pagava k^3 sem querer so descobria montando o H na
  // mao. O numero sai do padrao simbolico, que ja existe, sem custo extra.
  if (cs.pronto) {
    R.bloco_denso = bloco_denso_simbolico(cs.sb);
    R.colunas_fator = cs.sb.n;
  }
  if (verboso) {
    Rprintf("AI-REML: %d record(s), %d column(s), %d component(s)\n",
            (int) d.n_usadas(), (int) d.total_colunas(), (int) d.modelo.ntheta);
    if (R.colunas_fator > 0)
      Rprintf("  dense block of the factor: %d of %d columns, about %.3g flops per "
              "factorization in it\n", (int) R.bloco_denso, (int) R.colunas_fator,
              std::pow((double) R.bloco_denso, 3.0) / 3.0);
  }
  if (!cur.ok) {
    R.mensagem = "the starting theta is INADMISSIBLE: some covariance matrix is not positive-definite";
    return R;
  }

  // EM de aquecimento: monotono, fica no cone, e de graca. A checagem extra de
  // z_de_theta garante que o ponto de onde o passo amortecido parte tem fator de
  // Cholesky, um EM que encoste numericamente na singularidade nao vira partida.
  for (std::size_t e = 0; e < n_em; e++) {
    std::vector<double> zt;
    // o EM tambem respeita o preso: ele recalcula TODOS os componentes por formula
    // fechada, entao sem esta restauracao ele desfaria o fixed= a cada rodada.
    for (std::size_t k = 0; k < d.modelo.ntheta; k++)
      if (d.modelo.preso(k)) cur.em_theta[k] = theta[k];
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
  // casos a diagonal do fator de Cholesky tem um piso, absoluto para a variancia, RELATIVO
  // a maior diagonal do grupo para a singularidade, e a entrada presa no piso com o score
  // empurrando para fora e congelada, com o resto otimizado CONDICIONALMENTE a isso, que e
  // a solucao correta do problema com restricao.
  double lambda = 1e-3;
  std::vector<double> z;
  if (!z_de_theta(d.modelo, theta, z)) {
    R.mensagem = "the starting theta is INADMISSIBLE: some covariance matrix is not positive-definite";
    return R;
  }

  // floors on the diagonal of z, recomputed at the current point: relative to the group's
  // largest Cholesky diagonal (so a one-term group never binds on it) and absolute against
  // a FIXED scale, the variance of y over the records used (read through l = sqrt v)
  //
  // A ESCALA DO PISO ABSOLUTO E FIXA NA CAMINHADA. Era 1e-10 max(1, s2e) com a s2e do ponto
  // corrente: o max(1, .) prendia o piso em 1e-10 na unidade do usuario sempre que s2e < 1
  // (com y vezes 1e-5 var(pe) ficava presa nele e o ajuste certificava var(animal) 14%
  // abaixo; com 1e-6 nem convergia, 1,685 unidades acima), e a s2e corrente faz o piso andar
  // junto com o residuo. Tirar so o max(1, .) deixou isso a mostra: uma variancia presa no
  // piso passa a andar com a s2e, o passo nao sabe disso, e um ajuste que o base certificava
  // em 19 iteracoes (test-mgs-misto.R, var(animal) na fronteira) rodou 300 com relDelta
  // 1e-5 sem convergir. Com var(y), calculada uma vez, o piso e invariante a unidade e nao
  // se move.
  double escala_piso = 1.0;
  {
    double soma = 0.0, soma2 = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < d.nlin; i++)
      if (d.usa[i]) { soma += d.y[i]; soma2 += d.y[i] * d.y[i]; n++; }
    if (n > 1) {
      const double media = soma / static_cast<double>(n);
      double ss = 0.0;
      for (std::size_t i = 0; i < d.nlin; i++)
        if (d.usa[i]) ss += (d.y[i] - media) * (d.y[i] - media);
      const double v = ss / static_cast<double>(n - 1);
      if (v > 0.0 && std::isfinite(v)) escala_piso = v;
    }
  }
  // O PISO DO RESIDUO (lista de correcoes, item 9). No univariado o residuo e um escalar em
  // log, e o mapa tem so os grupos como blocos; sem piso, o residuo que ia a zero andava em
  // log ate -infinito, nunca entrava no conjunto ativo, e o ajuste parava com o decremento em
  // 1e11 e uma mensagem que nao dizia que era o residuo (medido com dois registros iguais por
  // animal, var(residual) a 9e-14). Aqui ele tem o piso absoluto das variancias dos grupos,
  // 1e-10 da escala, e as mesmas regras: grampo, conjunto ativo, relato e EP condicional
  auto pisos_z = [&](const std::vector<double>& zc, std::vector<double>& piso,
                     std::vector<char>& relativo) {
    br::pisos_z(mapa_z(d.modelo), zc, escala_piso, piso, relativo);
    for (std::size_t k = 0; k < d.modelo.n_residuos(); k++) {
      piso[d.modelo.offset_residual + k] = std::log(1e-10 * escala_piso);
      relativo[d.modelo.offset_residual + k] = 0;
    }
  };
  // os indices de theta com piso em z: a diagonal de cada grupo e o residuo
  std::vector<std::size_t> com_piso;
  for (const Grupo& gr : d.modelo.grupos)
    for (std::size_t i = 0; i < gr.dim; i++) com_piso.push_back(gr.theta_idx(i, i));
  for (std::size_t k = 0; k < d.modelo.n_residuos(); k++)
    com_piso.push_back(d.modelo.offset_residual + k);

  // The pieces of the step at a point, shared by the walker and by the final certificate:
  // score and AI mapped to z by the chain rule, the floors, and two active sets over the
  // SAME floors and the SAME outward-score rule, differing only in reach. `congelado` is
  // the step's (strictly at the clamp), `na_parede` is the certificate's and the boundary
  // report's (within a factor of 2 of the floor): a component pinned at a floor with the
  // score pushing outward never zeroes its gradient, it points out of the cone by
  // construction, so the certificate must exclude it, whether the walker left it clamped
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
    for (const std::size_t k : com_piso) {
        // preso pelo usuario ou grampeado na fronteira: nos dois casos a direcao sai
        // do sistema, o score dela e ignorado e o resto e otimizado CONDICIONALMENTE.
        P.congelado[k] = d.modelo.preso(k) ? 1 :
            ((zc[k] <= P.piso[k] + 1e-9 && P.sz[k] > 0.0) ? 1 : 0);
        // the certificate's active set, and it is WIDER than the step's on purpose: the
        // walker rests a hair above the clamp without ever being frozen there (the same
        // measured fact the boundary report handles with the same factor-2 tolerance).
        // A hair above the floor the z-gradient of the pinned diagonal is ~L^2 and its
        // AI_z diagonal ~L^4, so the quadratic form returns that direction's FULL
        // theta-scale share, measured 1.1499 on a true boundary optimum (seed 11 of
        // test-aireml-boundary.R) where the free directions were flat. On the floor
        // with the score pointing outward is pinned, whether clamped or resting.
        // um componente PRESO tambem sai do certificado: ele nao pode andar, entao o
        // score dele nunca zera e cobra-lo impediria qualquer ajuste com fixed= de
        // convergir. O que se certifica e o otimo CONDICIONAL ao que foi preso.
        P.na_parede[k] = d.modelo.preso(k) ? 1 :
            ((zc[k] <= P.piso[k] + std::log(2.0) && P.sz[k] > 0.0) ? 1 : 0);
      }
    return P;
  };

  // THE CONVERGENCE CERTIFICATE, measured missing on real data (2026-09-02): a
  // direct+indirect fit ended converged TRUE with relDelta 5.1e-09 and
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
  // its optimum, four orders of magnitude below the ~4 units a likelihood-ratio test
  // calls a difference, while the measured defect (gap 6.9, dec ~14) fails it by five.
  // The inverse of the free block is SPECTRAL, with each eigenvalue floored at 1e-9 of the
  // largest. On a flat ridge the AI is singular to the rounding (measured on the direct +
  // indirect gate of test-grupo-multitermo.R: eigenvalues 5.3e3, 5.4e2, 54 and +-1e-12, the
  // gradient 1e-9 along the last), and the plain inverse gave a decrement of 1e-6 at one point
  // and of 1.7e7 at a point 1e-3 away with the same -2logL to 1e-8: the number followed the
  // rounding of the null eigenvalue. With the floor, a direction of numerically zero curvature
  // counts g^2 / (1e-9 lambda_max): a gradient of 1e-9 there certifies, and a large one (the
  // case above, scores of 5e5) still refuses.
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
    for (const double v : m.dados)
      if (!std::isfinite(v)) return std::numeric_limits<double>::infinity();
    std::vector<double> ev;
    Densa u;
    jacobi_sim(m, ev, u);
    double maior = 0.0;
    for (const double e : ev) maior = std::max(maior, std::fabs(e));
    if (!(maior > 0.0)) return std::numeric_limits<double>::infinity();
    // the AI is PSD; a negative eigenvalue here is rounding, not curvature
    const double piso = 1e-9 * maior;
    double dec = 0.0;
    for (std::size_t k = 0; k < ev.size(); k++) {
      double proj = 0.0;
      for (std::size_t i = 0; i < sc.size(); i++) proj += u.at(i, k) * sc[i];
      dec += proj * proj / std::max(std::fabs(ev[k]), piso);
    }
    return dec;
  };

  // EM rescue: multiplicative, stays in the cone, and does not shrink at the boundary, so
  // it walks exactly where the damped step stalls. Each EM candidate is CLAMPED through z
  // (an EM that lands numerically ON the singularity would evaluate to garbage: the MME
  // -2logL loses its meaning when det(C_g) underflows). Up to 50 steps, while each pays
  // more than 1e-8.
  // o gancho de teste do D2 (o racional esta em multitrait2.cpp): recusa todo passo AI e EM
  const bool sem_passo = std::getenv("UAITREM_TESTE_SEM_PASSO") != nullptr;
  auto resgate_em = [&]() -> bool {
    bool ganhou = false;
    if (sem_passo) return false;
    for (int e = 0; e < 50; e++) {
      std::vector<double> ze;
      // mesma restauracao do EM de aquecimento: o resgate tambem recomputa tudo
      for (std::size_t k = 0; k < d.modelo.ntheta; k++)
        if (d.modelo.preso(k)) cur.em_theta[k] = theta[k];
      if (!z_de_theta(d.modelo, cur.em_theta, ze)) break;
      std::vector<double> piso; std::vector<char> rel;
      pisos_z(ze, piso, rel);
      for (const std::size_t k : com_piso) ze[k] = std::max(ze[k], piso[k]);
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

  // SEM REENTRADA POR SONDA AQUI, e de proposito. A condicao de KKT o univariado ja confere
  // pelo sinal: o certificado so exclui uma coordenada na parede com o score em z apontando
  // para FORA (na_parede, acima), e o passo so congela nessa mesma condicao, entao uma
  // variancia perto da parede com o gradiente para DENTRO nunca e congelada: o passo em z
  // dela cresce como 1 / L^2 e o caminhante sai sozinho. A sonda pela verossimilhanca dos
  // ajustadores multicaracter (candidatos_reentrada) existe porque passo_z congela pela
  // curvatura sem olhar o sinal. A versao intermediaria deste conserto pos aqui a sonda por
  // coordenada, e ela nao disparou nenhuma vez: partidas com var(ninhada) em 1e-3, 1e-6 e
  // 1e-10, e em r = +/-(1 - 1.5e-6) e +/-(1 - 1e-8) num grupo direto + indireto, chegaram
  // ao otimo ou pararam com converged = FALSE, nunca um TRUE falso.

  // iteracoes seguidas aceitas sem ganho relevante, e se o laco parou por elas (o corte de
  // estagnacao, no fim do laco)
  int parado = 0;
  bool estagnou = false;
  for (std::size_t it = 1; it <= maxiter; it++) {
    R.iters = it;
    R_CheckUserInterrupt();
    // chain rule (score_z = J' score, AI_z = J' AI J), floors and active set, all at
    // the current point, the same pieces the certificate reads at the exits
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
      // congelado: linha e coluna viram identidade e o score sai, e o sistema REDUZIDO
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
        for (const std::size_t k : com_piso) zc[k] = std::max(zc[k], piso[k]);
        if (sem_passo) continue;
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
          const double ganho = cur.neg2logl - prox.neg2logl;
          parado = ganho > 1e-8 * std::max(1.0, std::fabs(cur.neg2logl)) ? 0 : parado + 1;
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
      // gate, the damped step can be stuck while EM still walks, and after a rescue the
      // damping RESTARTS, or the next AI step would inherit a lambda pumped to 1e20+ by
      // the failed attempts and stay parked forever ("amortecimento preso")
      if (resgate_em()) {
        lambda = 1e-3;
        parado = 0;
        if (verboso)
          Rprintf("      EM rescue: -2logL %.6f (no damped step improved)\n", cur.neg2logl);
        continue;
      }
      // In z every candidate is admissible, so this exit means the likelihood itself
      // refused 90 second-order variants AND the EM step. That is where the walker
      // rests, but resting is not an optimum: on the measured direct+indirect regime
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
        R.mensagem = std::string("did NOT converge: at iteration ") + std::to_string(it) +
            " neither the damped AI step nor the EM step improved the likelihood, but the "
            "Newton decrement of the free components is " + buf + " (tolerance 2e-4 on the "
            "-2logL scale), so the point is a stall, not a certified optimum. Try a "
            "different start=, and profile any component the message reports at a boundary";
      }
      break;
    }
    if (verboso) {
      Rprintf("iter %3d  -2logL %.6f  relDelta %.3e\n",
              (int) it, cur.neg2logl, R.reldelta);
      imprime_theta(theta, d.modelo.nomes_theta(), d.modelo);
    }
    if (R.reldelta < tol) {
      // O PASSO PEQUENO NAO PROVA OTIMO. Com uma componente encostada numa fronteira a AI
      // fica quase singular naquela direcao, o amortecimento cresce, o passo encolhe, e
      // um criterio que olha so o tamanho do passo declararia convergencia com o score
      // longe de zero. Antes de declarar, o EM tem de confirmar que nao anda mais; se
      // ainda melhora, nao havia convergencia nenhuma, e o amortecimento reinicia.
      if (resgate_em()) {
        lambda = 1e-3;
        parado = 0;
        if (verboso)
          Rprintf("      EM rescue: -2logL %.6f (the AI step had stalled)\n", cur.neg2logl);
        continue;
      }
      // ... e um EM parado tambem nao prova nada quando AI e EM travam JUNTOS (o defeito
      // do certificado, medido): converged exige ainda o decremento de Newton dos
      // componentes livres. A recusa NAO encerra o ajuste, o passo pode voltar a andar,
      // e se maxiter chegar, a mensagem final carrega o decremento reprovado.
      R.decremento = decremento_livre(pecas_z(z, cur));
      if (R.decremento < tol_dec) {
        R.convergiu = true;
        break;
      }
    }
    // ESTAGNACAO, como nos espelhos (ajusta_mt, ajusta_ar1): 25 iteracoes seguidas aceitas
    // sem ganho acima de 1e-8 do -2logL e sem EM que ande. A recusa do certificado acima nao
    // encerra o laco, porque o passo pode voltar a andar; mas um caminhante que nao anda mais
    // seguia ate maxiter repetindo o mesmo ponto. Medido num direto + indireto com o social
    // aliasado (simulacao, 192 registros): 200 iteracoes com relDelta 0 e o mesmo -2logL, e
    // num modelo grande cada uma e uma avaliacao inteira. O veredito segue o certificado,
    // como no fim do orcamento (a crista plana, abaixo).
    if (parado >= 25) {
      R.decremento = decremento_livre(pecas_z(z, cur));
      estagnou = true;
      break;
    }
  }
  // the decrement at the final point, whatever the exit: the certificate the user can
  // read next to fit$score, and the number the maxiter message reports
  if (std::isnan(R.decremento)) R.decremento = decremento_livre(pecas_z(z, cur));
  // CRISTA PLANA NO FIM DO ORCAMENTO. O certificado so e consultado depois que relDelta cai
  // abaixo da tolerancia, e isso e deliberado: parar pelo decremento a cada iteracao encerra
  // cedo demais, com o -2logL ja otimo mas theta ainda impreciso, e reprova os portoes de
  // exatidao (medido: 19 falhas). Sobra um caso em que a regra mentia. Quando a
  // verossimilhanca e plana ao longo de uma direcao pouco identificada, os componentes
  // deslizam por ela, relDelta nunca cai, e o laco esgotava maxiter dizendo converged FALSE
  // com o decremento ja abaixo da tolerancia. Medido num direto + indireto com 800 animais:
  // 300 iteracoes, relDelta 2.2e-4, decremento 1.3e-6 contra 2e-4. O model_mt() convergia no
  // mesmo dado. Aqui o veredicto passa a seguir o certificado, e a mensagem diz que a
  // direcao e plana, que e a informacao que importa para ler aqueles componentes.
  // So quando o laco andou: com maxiter = 0 nao ha caminhada nenhuma, e uma AVALIACAO no
  // start= nao e convergencia, mesmo com o decremento pequeno (o ponto pode ser um otimo
  // dado de fora, mas quem o certifica e o ajuste que chega nele).
  if (estagnou && !R.convergiu && R.decremento >= tol_dec) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.3g", R.decremento);
    R.mensagem = std::string("did NOT converge: at iteration ") + std::to_string(R.iters) +
        " the step had made no progress for 25 iterations (no gain above 1e-8 of the -2logL, "
        "and no EM step that moved), and the Newton decrement of the free components is " +
        buf + " (tolerance 2e-4 on the -2logL scale), so the point is a stall, not a "
        "certified optimum. More iterations give the same point: fit from another start=, "
        "and profile any component the message reports at a boundary";
  }
  if (!R.convergiu && maxiter > 0 && (R.iters >= maxiter || estagnou) &&
      R.decremento < tol_dec) {
    R.convergiu = true;
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.3g", R.reldelta);
    R.mensagem = std::string("converged by the Newton decrement while the components were "
        "still moving (relDelta ") + buf + "): the likelihood is flat along some direction, "
        "so the values along it are about equally supported. Read the standard errors and "
        "the correlation between the components before trusting any one of them";
  }
  // norma do score no ponto final, para quem quiser conferir que ele de fato zerou:
  // e a evidencia direta de otimo, que o tamanho do passo nao da.
  R.score = cur.score;
  if (verboso)
    Rprintf("%s at iter %d, -2logL %.6f\n",
            R.convergiu ? "converged" : "STOPPED", (int) R.iters, cur.neg2logl);

  // erros padrao: 2 [AI_ll^-1]_kk no bloco LIVRE (lista de correcoes, itens 17 e 18). Um
  // componente DADO nao foi estimado: o preso (kernel(fixed =)) e o componente escalar parado
  // na fronteira zero (a menos de 2 vezes o piso, a regra do relatorio abaixo, "fixed there
  // while the others were optimized conditional on that"). Ele fica sem EP, e os EP dos outros
  // saem da informacao CONDICIONAL a ele. Antes o componente preso ganhava EP e os outros
  // vinham da AI inteira: medido num kernel(pen, fixed = 0.36), var(kernel) com EP 0,108 e os
  // EP de var(animal) e var(residual) iguais aos do ajuste com a baia livre. A matriz inteira
  // tambem sai: 2 AI_ll^-1 e a covariancia amostral das componentes, sem ela nao ha metodo
  // delta para h2, r ou T; a linha e a coluna do componente dado sao 0, porque ali ele entra
  // como constante. O diagnostico SINGULAR (aviso_informacao, no fim) le o mesmo bloco livre.
  const std::size_t nth = d.modelo.ntheta;
  std::vector<char> dado(nth, 0);
  for (std::size_t k = 0; k < nth; k++)
    if (d.modelo.preso(k)) dado[k] = 1;
  // as coordenadas de z de um grupo de VARIAS dimensoes paradas no piso (a parede do bloco: uma
  // correlacao +-1, ou uma combinacao dos efeitos sem variancia), pela mesma regra de 2 vezes
  // o piso do relatorio abaixo
  std::vector<double> zf;
  std::vector<char> parede(nth, 0);
  bool bloco_na_parede = false;
  if (z_de_theta(d.modelo, theta, zf)) {
    std::vector<double> piso; std::vector<char> rel;
    pisos_z(zf, piso, rel);
    for (const Grupo& gr : d.modelo.grupos) {
      if (gr.dim == 1) {
        const std::size_t k = gr.theta_idx(0, 0);
        if (!rel[k] && zf[k] <= piso[k] + std::log(2.0)) dado[k] = 1;
        continue;
      }
      for (std::size_t i = 0; i < gr.dim; i++) {
        const std::size_t k = gr.theta_idx(i, i);
        if (zf[k] <= piso[k] + std::log(2.0)) { parede[k] = 1; bloco_na_parede = true; }
      }
    }
    for (std::size_t k = 0; k < d.modelo.n_residuos(); k++) {
      const std::size_t kr = d.modelo.offset_residual + k;
      if (zf[kr] <= piso[kr] + std::log(2.0)) dado[kr] = 1;
    }
  }
  std::vector<std::size_t> livre;
  for (std::size_t k = 0; k < nth; k++)
    if (!dado[k]) livre.push_back(k);
  Densa ai_livre(livre.size(), livre.size());
  for (std::size_t a = 0; a < livre.size(); a++)
    for (std::size_t b = 0; b < livre.size(); b++) ai_livre.at(a, b) = cur.ai.at(livre[a], livre[b]);
  R.se.assign(nth, std::nan(""));
  R.vcov.assign(nth * nth, 0.0);
  for (const std::size_t a : livre)
    for (const std::size_t b : livre) R.vcov[a * nth + b] = std::nan("");
  if (!bloco_na_parede) {
    try {
      Densa inv = inv_geral(ai_livre);
      for (std::size_t a = 0; a < livre.size(); a++) {
        const double v = 2.0 * inv.at(a, a);
        if (v > 0.0) R.se[livre[a]] = std::sqrt(v);
        for (std::size_t b = 0; b < livre.size(); b++)
          R.vcov[livre[a] * nth + livre[b]] = 2.0 * inv.at(a, b);
      }
    } catch (const Erro&) {}
  } else {
    // O BLOCO NA PAREDE (lista de correcoes, item 17, a metade dos grupos). O grupo foi
    // segurado na parede enquanto os outros eram otimizados condicionais a ela, entao o EP e o
    // da informacao CONDICIONAL a parede, no mapa do passo: V = 2 J_F (J_F' AI J_F)^-1 J_F', com
    // J = d theta / d z e F as coordenadas de z livres (fora as do piso e as dadas). Sem parede
    // J_F = J e V = 2 AI^-1, a regra de antes; por isso este caminho so roda com a parede, e um
    // ajuste sem ela segue bit a bit. Medido num rn() 2 x 2 de posto 1 verdadeiro (simulacao,
    // 480 animais, 4 registros cada, r = 1 no ajuste): o EP de var(pe) era 0,0518 pela AI
    // inteira, 0,0486 assim, e 0,0489 pela curvatura do perfil de -2logL (profile_theta, so
    // avaliacoes); var(residual) 0,0271, 0,02686 e 0,0268. Um componente cujo theta so depende
    // de coordenadas seguradas (uma variancia em zero dentro do grupo) fica sem EP.
    std::vector<std::size_t> fz;
    for (std::size_t m = 0; m < nth; m++)
      if (!parede[m] && !dado[m]) fz.push_back(m);
    vcov_condicional(cur.ai, jacobiano_z(d.modelo, zf), fz, livre, R.vcov, R.se);
  }

  if (!R.convergiu && maxiter == 0) {
    // uma avaliacao, nao um ajuste: dito como tal, com o decremento no ponto dado
    char buf[240];
    if (n_em > 0)
      std::snprintf(buf, sizeof(buf), "maxiter = 0: no AI iteration was run, so the "
                    "components are the start moved by %d EM round(s) and evaluated, not "
                    "estimated; the Newton decrement there is %.3g", (int) n_em, R.decremento);
    else
      std::snprintf(buf, sizeof(buf), "maxiter = 0: no AI iteration was run, so the "
                    "components are the start evaluated, not estimated; the Newton decrement "
                    "there is %.3g", R.decremento);
    R.mensagem = buf;
  } else if (!R.convergiu && R.mensagem.empty()) {
    // The maxiter exit, and it comes BEFORE the boundary notes so those append to the
    // explanation instead of silencing it. Measured by the user on the 2x2 group
    // warm-started from the reduced model: 100 iterations ended with relDelta still
    // 1.6e-04, not a defect, a model that walks slowly. The message says exactly
    // that, and how to ask for more.
    char buf[128];
    std::snprintf(buf, sizeof(buf), " (relDelta %.3g, Newton decrement of the free "
                  "components %.3g against the 2e-4 tolerance)", R.reldelta, R.decremento);
    // relDelta 0: congelou, e mais iteracoes dao o mesmo ponto (lista de correcoes, item 8)
    R.mensagem = "stopped at " + std::to_string(maxiter) + " iteration(s) without a "
        "certified optimum" + buf + (R.reldelta == 0.0
        ? ": the walk froze, no step changed the components, so more iterations give the "
          "same point. Fit from another start=, or read the components resting at a boundary"
        : ": this model asks for more iterations. Raise maxiter=; a start= from a reduced "
          "model shortens the walk");
  }
  {
    // The boundary report reads the FINAL point, not the freeze flags of the last
    // iteration: near a boundary optimum the walker can rest a hair above the floor
    // without ever being frozen there, and the user still needs to know the estimate
    // is on the edge (fit$score is NOT ~0 in the boundary direction, by design).
    // Os componentes e os grupos sao NOMEADOS (lista de correcoes, item 21), e o grupo singular
    // diz se ha uma correlacao de par em +-1: num grupo de 3 ou mais a singularidade pode vir
    // de uma combinacao dos efeitos sem correlacao extrema nenhuma, e a mensagem dizia "a
    // correlation at +/-1" mesmo assim
    std::size_t nzero = 0, nsing = 0;
    std::string zeros, singulares;
    const std::vector<std::string> nomes_f = d.modelo.nomes_theta();
    std::vector<double> zf;
    if (z_de_theta(d.modelo, theta, zf)) {
      std::vector<double> piso; std::vector<char> rel;
      pisos_z(zf, piso, rel);
      for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
        const Grupo& gr = d.modelo.grupos[g];
        std::size_t sing = 0;
        for (std::size_t i = 0; i < gr.dim; i++) {
          const std::size_t k = gr.theta_idx(i, i);
          // within a factor of 2 of the floor counts as ON it: the walker rests where
          // the trench flattens, a hair above the clamp, and the reader still needs
          // to know the estimate is an edge case
          if (!(zf[k] <= piso[k] + std::log(2.0))) continue;
          if (rel[k]) {
            sing++;
          } else {
            nzero++;
            zeros += (zeros.empty() ? "" : ", ") + nomes_f[k];
          }
        }
        if (!sing) continue;
        nsing += sing;
        const Densa c = cov_grupo(d.modelo, theta, g);
        double rmax = 0.0;
        for (std::size_t i = 0; i < gr.dim; i++)
          for (std::size_t j = i + 1; j < gr.dim; j++)
            if (c.at(i, i) > 0.0 && c.at(j, j) > 0.0)
              rmax = std::max(rmax, std::fabs(c.at(i, j)) / std::sqrt(c.at(i, i) * c.at(j, j)));
        char buf[200];
        if (rmax >= 0.999)
          std::snprintf(buf, sizeof(buf), "group %s: a correlation at +/-1", gr.nome.c_str());
        else
          std::snprintf(buf, sizeof(buf), "group %s: a combination of its effects with no "
                        "variance, and no correlation at +/-1, the largest is %.3g",
                        gr.nome.c_str(), rmax);
        singulares += (singulares.empty() ? "" : "; ") + std::string(buf);
      }
      const std::vector<std::string> nomes_r = d.modelo.nomes_theta();
      for (std::size_t k = 0; k < d.modelo.n_residuos(); k++) {
        const std::size_t kr = d.modelo.offset_residual + k;
        if (zf[kr] <= piso[kr] + std::log(2.0))
          R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") + nomes_r[kr] + " at the "
              "zero boundary (held at 1e-10 of the variance of the records while the others "
              "were optimized conditional on that): the random effects take all the variation "
              "between the records, so the model has more random effects than the data hold "
              "apart from the residual";
      }
    }
    if (nzero > 0)
      R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") + std::to_string(nzero) +
          " component(s) at the zero boundary (" + zeros + "), fixed there while the others "
          "were optimized conditional on that";
    // sem "the boundary is real": o item 41 da lista mediu um ponto do proprio modelo 3,08
    // unidades abaixo de uma parede relatada assim
    if (nsing > 0)
      R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") + std::to_string(nsing) +
          " covariance group direction(s) at the singularity boundary (" + singulares +
          "), held there while the others were optimized conditional on that; a "
          "delta-method interval there is not valid, profile the component instead";
  }
  if (cur.fora_do_padrao > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(cur.fora_do_padrao) +
        " value(s) read outside the sparsity pattern of the factor: the result is NOT reliable";

  R.theta = theta;
  R.neg2logl = cur.neg2logl;
  R.solucao = cur.solucao;
  R.fora_do_padrao = cur.fora_do_padrao;

  // PEV: a diagonal da inversa seletiva das MME no otimo, vezes s2e.
  //
  // As MME estao em unidades de s2e (a penalidade e kron(Cs^-1, K^-1) com Cs = C/s2e),
  // entao [C_s^-1]_ii ja e PEV_i / s2e, a multiplicacao devolve a escala absoluta. Sem
  // esse fator a acuracia sai sistematicamente errada e ainda parece plausivel.
  {
    // no residuo por classe a PEV e os residuos saem das linhas escaladas no ponto final (o
    // R divide residuo e ajustado pelo sqrt(w) da linha, como com os pesos)
    std::vector<double> w_fim;
    Desenho d_fim;
    const bool het_fim = d.modelo.n_residuos() > 1;
    if (het_fim) {
      d_fim = escala_por_classe(d, theta, w_fim);
      cs.mont_dono = nullptr;
    }
    const Desenho& dd = het_fim ? d_fim : d;
    Montado M = monta_mme(dd, theta, &cs);
    if (M.ok) {
      if (!cs.pronto) {
        cs.perm = grau_minimo(M.c);
        cs.sb = simbolica(permuta_sim(M.c, cs.perm));
        cs.pronto = true;
      }
      const std::vector<std::size_t>& perm = cs.perm;
      const Csc& pc = permuta_cache(M.c, cs);
      Csc L;
      if (cholesky(pc, cs.sb, L)) {
        SelInv z = inversa_seletiva(L, 0);
        std::vector<double> diag = z.diagonal();
        R.pev.assign(M.total, std::nan(""));
        for (std::size_t k = 0; k < M.total; k++)
          R.pev[perm[k]] = diag[k] * M.s2e;
        // uma falha aqui (memoria) nao leva junto o ajuste: os residuos ficam NA
        try {
          residuos_no_fator(dd, M, L, z, perm, R);
        } catch (const Interrompido&) {
          throw;
        } catch (const std::exception& e) {
          R.residuo.clear(); R.ajustado.clear(); R.ajustado_fixo.clear(); R.alavanca.clear();
          R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
              "residuals not computed (" + e.what() + "): the fit is kept, residual and "
              "leverage are NA";
        }
      }
    }
  }
  if (cur.ok) {
    const std::vector<std::string> todos = d.modelo.nomes_theta();
    std::vector<std::string> nomes_livres;
    std::vector<double> se_livre;
    for (const std::size_t k : livre) {
      nomes_livres.push_back(todos[k]);
      se_livre.push_back(R.se[k]);
    }
    const std::string av = aviso_informacao(ai_livre, nomes_livres, &se_livre);
    for (std::size_t a = 0; a < livre.size(); a++) R.se[livre[a]] = se_livre[a];
    if (!av.empty()) R.mensagem += (R.mensagem.empty() ? "" : "; ") + av;
    // O COMPONENTE QUASE SEM INFORMACAO fora do conjunto SINGULAR (lista de correcoes, item 20):
    // uma variancia com EP 10 vezes a variancia fenotipica (a soma das var() no ponto) ou mais e
    // falta de informacao do mesmo jeito (medido num materno: var(maternal) 0,020 com EP
    // 12647), e o veredito, que le os autovalores da AI, deixava ela de fora. A regua e a escala
    // fenotipica, e nao a estimativa: uma variancia pequena com EP acima dela cobre o zero, e
    // o EP ali e a informacao que os dados tem
    double vp = 0.0;
    for (std::size_t k = 0; k < nth; k++)
      if (todos[k].compare(0, 4, "var(") == 0 && theta[k] > 0.0) vp += theta[k];
    std::string quase;
    for (std::size_t a = 0; vp > 0.0 && a < livre.size(); a++) {
      const std::size_t k = livre[a];
      if (!(R.se[k] > 0.0) || !std::isfinite(R.se[k])) continue;
      if (todos[k].compare(0, 4, "var(") != 0 || R.se[k] < 10.0 * vp) continue;
      char buf[200];
      std::snprintf(buf, sizeof(buf), "%s (SE %.3g times the phenotypic variance)",
                    todos[k].c_str(), R.se[k] / vp);
      quase += (quase.empty() ? "" : ", ") + std::string(buf);
    }
    if (!quase.empty())
      R.mensagem += (R.mensagem.empty() ? "" : "; ") + std::string("the data hold almost no "
          "information on ") + quase + ": read it as not estimated, as a component of the "
          "SINGULAR verdict";
  }
  if (!nota_partida.empty())
    R.mensagem = R.mensagem.empty() ? nota_partida : R.mensagem + "; " + nota_partida;
  return R;
}

}  // namespace br
