// Multi-caracteristica, parte 2: layout, desenho e o estimador.

#include "mme.h"

#include <cstdio>

// so para o TRACE opcional e a interrupcao: nenhuma conta usa o R
#include <R_ext/Print.h>
#include <R_ext/Utils.h>

namespace br {

// Reindexa o layout do modelo para t caracteristicas: cada grupo passa a dim * t, e o
// residual passa a vech(R0) com t(t+1)/2 parametros.
static void expande_layout_mt(Modelo& m, std::size_t t) {
  std::size_t off = 0;
  for (Grupo& g : m.grupos) {
    g.dim *= t;
    g.offset = off;
    g.nparam = g.dim * (g.dim + 1) / 2;
    off += g.nparam;
  }
  m.offset_residual = off;
  m.ntheta = off + t * (t + 1) / 2;
}

// Nomes de theta no caso multi: coeficiente rotulado "termo[c]@trait".
static std::vector<std::string> nomes_theta_mt(const Modelo& m,
                                               const std::vector<std::string>& alvos) {
  const std::size_t t = alvos.size();
  std::vector<std::string> out(m.ntheta);
  for (const Grupo& g : m.grupos) {
    std::vector<std::string> rot;
    for (std::size_t tau = 0; tau < t; tau++)
      for (std::size_t tm : g.termos) {
        const Termo& te = m.termos[tm];
        for (std::size_t c = 0; c < te.n_coef(); c++) {
          std::string b = te.nome;
          if (te.n_coef() > 1) b += "[" + std::to_string(c) + "]";
          rot.push_back(b + "@" + alvos[tau]);
        }
      }
    for (std::size_t j = 0; j < g.dim; j++)
      for (std::size_t i = j; i < g.dim; i++) {
        const std::size_t k = g.offset + (j * (2 * g.dim - j + 1)) / 2 + (i - j);
        out[k] = (i == j) ? "var(" + rot[i] + ")" : "cov(" + rot[i] + "," + rot[j] + ")";
      }
  }
  for (std::size_t j = 0; j < t; j++)
    for (std::size_t i = j; i < t; i++) {
      const std::size_t k = m.offset_residual + (j * (2 * t - j + 1)) / 2 + (i - j);
      out[k] = (i == j) ? "var(res@" + alvos[i] + ")"
                        : "cov(res@" + alvos[i] + ",res@" + alvos[j] + ")";
    }
  return out;
}

DesenhoMT monta_desenho_mt(Modelo m, const std::vector<std::string>& alvos,
                           const Tabela& tab, const Pedigree* ped,
                            const std::vector<KernelDecl>* kernels) {
  DesenhoMT d;
  d.alvos = alvos;
  d.t = alvos.size();
  if (d.t < 2) throw Erro("multi-trait needs at least two; for one, use model()");
  d.nlin = tab.nlin;

  // y n x t. Caracteristica ausente NAO derruba o registro: marca obs = 0, e o registro
  // participa com a submatriz de R0 do seu padrao. y ausente vira zero — com a inversa
  // embutida zerada naquela linha, o valor nunca e lido, e zero e mais seguro que NaN,
  // que contamina qualquer soma que um defeito futuro deixar passar.
  d.y = Densa(d.nlin, d.t);
  d.usa.assign(d.nlin, 1);
  d.obs.assign(d.nlin * d.t, 1);
  for (std::size_t tau = 0; tau < d.t; tau++) {
    std::vector<double> col = tab.numerico(alvos[tau]);
    for (std::size_t i = 0; i < d.nlin; i++) {
      const bool falta = !std::isfinite(col[i]) ||
          (m.tem_ausente && std::fabs(col[i] - m.codigo_ausente) < 1e-9);
      d.obs[i * d.t + tau] = falta ? 0 : 1;
      d.y.at(i, tau) = falta ? 0.0 : col[i];
    }
  }
  // registro sem NENHUMA caracteristica observada sai
  for (std::size_t i = 0; i < d.nlin; i++) {
    bool alguma = false;
    for (std::size_t tau = 0; tau < d.t; tau++)
      if (d.obs[i * d.t + tau]) { alguma = true; break; }
    if (!alguma) d.usa[i] = 0;
  }

  const bool precisa_ped = [&]{
    for (const Termo& tm : m.termos)
      if (tm.estrutura == Estrutura::Parentesco) return true;
    return false;
  }();
  if (precisa_ped && !ped) throw Erro("there is a term with relationship and no pedigree was given");

  // a K declarada passa pela MESMA reducao do univariado antes de qualquer conta
  std::vector<KernelDecl> kern_red;
  std::vector<std::unordered_set<std::string> > kern_nulos;
  reduz_kernels(m, kernels, kern_red, kern_nulos);

  Csc ainv;
  double ld_ainv = 0.0;
  std::vector<std::string> niveis_ped;
  if (precisa_ped) {
    std::vector<double> f = endogamia(*ped);
    ainv = a_inversa(*ped, f);
    niveis_ped = ped->ids;
    std::vector<std::size_t> perm = grau_minimo(ainv);
    Csc pa = permuta_sim(ainv, perm);
    Simbolica sb = simbolica(pa);
    Csc L;
    if (!cholesky(pa, sb, L)) throw Erro("the pedigree A^-1 is not positive-definite");
    ld_ainv = logdet(L);
  }
  for (const Grupo& g : m.grupos) {
    if (g.estrutura == Estrutura::Parentesco) {
      d.kinv.push_back(ainv);
      d.kinv_logdet.push_back(ld_ainv);
    } else if (g.estrutura == Estrutura::Declarada) {
      kinv_declarada(m, g, kernels, d.kinv, d.kinv_logdet);
    } else {
      d.kinv.push_back(Csc());
      d.kinv_logdet.push_back(0.0);
    }
  }

  // X por registro (a expansao por caracteristica e na montagem)
  std::vector<std::pair<std::string, std::vector<double>>> cols;
  cols.push_back({"intercept", std::vector<double>(d.nlin, 1.0)});
  for (std::size_t k = 0; k < m.termos.size(); k++) {
    if (m.termos[k].aleatorio()) continue;
    DesenhoTermo dt = monta_termo(m, k, tab, nullptr);
    Densa zd = dt.z.densa();
    for (std::size_t j = 0; j < dt.z.ncol; j++) {
      std::vector<double> v(d.nlin);
      for (std::size_t i = 0; i < d.nlin; i++) v[i] = zd.at(i, j);
      cols.push_back({dt.nome + "=" + dt.niveis[j % std::max<std::size_t>(dt.n_niveis, 1)],
                      std::move(v)});
    }
  }
  Densa xfull(d.nlin, cols.size());
  for (std::size_t j = 0; j < cols.size(); j++)
    for (std::size_t i = 0; i < d.nlin; i++) xfull.at(i, j) = cols[j].second[i];
  // POSTO DE X SOBRE AS LINHAS QUE ENTRAM, nao sobre a tabela inteira. E o mesmo conserto
  // que o univariado ja tem (mme.cpp, mesma marca). Um nivel fixo cujos registros TODOS
  // sairam continua com coluna nao-nula na tabela, entao um posto medido sobre tudo o
  // mantinha; na montagem, que so anda nas linhas usadas, essa coluna nao recebe nada, a
  // matriz de coeficientes ganha coluna VAZIA, e a Cholesky morre com "nao
  // positiva-definida" tres passos longe da causa. Medindo o posto onde o modelo de fato
  // vive, a coluna sai como dependente e e REPORTADA em dropped_x.
  std::vector<std::size_t> usadas;
  usadas.reserve(d.n_usadas());
  for (std::size_t i = 0; i < d.nlin; i++) if (d.usa[i]) usadas.push_back(i);
  Densa xusadas(usadas.size(), cols.size());
  for (std::size_t j = 0; j < cols.size(); j++)
    for (std::size_t r = 0; r < usadas.size(); r++)
      xusadas.at(r, j) = cols[j].second[usadas[r]];
  std::vector<std::size_t> fica, sai;
  posto_completo(xusadas, 1e-9, fica, sai);
  d.x = Densa(d.nlin, fica.size());
  for (std::size_t jj = 0; jj < fica.size(); jj++) {
    d.nomes_x.push_back(cols[fica[jj]].first);
    for (std::size_t i = 0; i < d.nlin; i++) d.x.at(i, jj) = xfull.at(i, fica[jj]);
  }
  for (std::size_t j : sai) d.saiu_x.push_back(cols[j].first);

  for (std::size_t k = 0; k < m.termos.size(); k++) {
    if (!m.termos[k].aleatorio()) continue;
    // niveis: do pedigree quando ha parentesco, DA K quando declarada. Todo nivel da
    // estrutura ganha equacao, com ou sem registro, exatamente como em model().
    const std::vector<std::string>* nf = nullptr;
    if (m.termos[k].estrutura == Estrutura::Parentesco) nf = &niveis_ped;
    else if (m.termos[k].estrutura == Estrutura::Declarada) nf = &(*kernels)[k].ids;
    d.aleatorios.push_back(monta_termo(m, k, tab, nf));
  }
  {
    std::vector<DesenhoTermo*> pa;
    for (DesenhoTermo& a : d.aleatorios) pa.push_back(&a);
    casa_niveis_nulos(m, pa, kern_nulos, tab, d.nlin);
  }
  for (const DesenhoTermo& a : d.aleatorios)
    for (std::size_t i = 0; i < d.nlin; i++)
      if (!a.casou[i]) d.usa[i] = 0;
  if (d.n_usadas() == 0) throw Erro("no row enters the analysis");

  expande_layout_mt(m, d.t);
  d.modelo = std::move(m);
  return d;
}

std::vector<std::string> nomes_theta_do_mt(const DesenhoMT& d) {
  return nomes_theta_mt(d.modelo, d.alvos);
}

static std::vector<double> partida_mt(const DesenhoMT& d) {
  const std::size_t t = d.t;
  std::vector<double> var_t(t, 1.0);
  for (std::size_t tau = 0; tau < t; tau++) {
    double soma = 0.0, soma2 = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < d.nlin; i++)
      if (d.usa[i] && d.obs[i * t + tau]) {
        soma += d.y.at(i, tau);
        soma2 += d.y.at(i, tau) * d.y.at(i, tau);
        n++;
      }
    const double media = soma / std::max<std::size_t>(n, 1);
    double v = (soma2 - n * media * media) / std::max(1.0, static_cast<double>(n - 1));
    if (!(v > 0.0) || !std::isfinite(v)) v = 1.0;
    var_t[tau] = v;
  }
  std::vector<double> theta(d.modelo.ntheta, 0.0);
  const double frac = 0.5 / static_cast<double>(d.modelo.grupos.size());
  std::size_t ig = 0;
  for (const Grupo& g : d.modelo.grupos) {
    // A escala da K entra na PARTIDA. Uma K declarada nao tem escala canonica: o usuario
    // pode passar D ou 4D para o mesmo modelo, e a verossimilhanca e identica desde que o
    // componente ande junto. Partir de frac*var(y) nos dois casos e partir de dois pontos
    // DIFERENTES da mesma superficie, e a diferenca sobrevive ate o fim (medido: 204.18
    // contra 215.64 de -2logL no mesmo ajuste, com converged TRUE nos dois). A media
    // geometrica dos autovalores de K sai de graca do log-determinante que o desenho ja
    // guarda, |K^-1|, e dividir por ela deixa a partida EQUIVARIANTE: com cK ela cai por c
    // e o caminho e o mesmo. So a K declarada precisa disso; A tem diagonal ~1 e a
    // diagonal e I, e mexer nelas so deslocaria ajustes que ja estao certos.
    double esc = 1.0;
    if (g.estrutura == Estrutura::Declarada && ig < d.kinv.size() && d.kinv[ig].ncol > 0) {
      const double n = static_cast<double>(d.kinv[ig].ncol);
      const double gm = std::exp(-d.kinv_logdet[ig] / n);
      if (std::isfinite(gm) && gm > 0.0) esc = gm;
    }
    ig++;
    // diagonal por caracteristica: os primeiros n_coef slots sao da trait 0, e assim por diante
    const std::size_t por_t = g.dim / d.t;
    for (std::size_t i = 0; i < g.dim; i++) {
      const std::size_t tau = i / por_t;
      const std::size_t k = g.offset + (i * (2 * g.dim - i + 1)) / 2;   // vech (i,i)
      theta[k] = frac * var_t[tau] / esc;
    }
  }
  for (std::size_t tau = 0; tau < d.t; tau++) {
    const std::size_t k = d.modelo.offset_residual + (tau * (2 * d.t - tau + 1)) / 2;
    theta[k] = 0.5 * var_t[tau];
  }
  return theta;
}

AjusteMT ajusta_mt(const DesenhoMT& d, const std::vector<double>* theta0, std::size_t maxiter,
                   double tol, bool verboso) {
  AjusteMT R;
  // start= do usuario, ou a partida automatica. Vale a mesma razao do univariado: a
  // maneira de conferir que o otimo nao depende de onde a busca comecou, e a saida que a
  // mensagem de parada recomenda.
  std::vector<double> theta = theta0 ? *theta0 : partida_mt(d);
  CacheSimbolica cs;
  AvaliacaoMT cur = avalia_mt(d, theta, &cs);
  if (!cur.ok) {
    R.mensagem = "o theta inicial e INADMISSIVEL";
    return R;
  }

  // O PASSO anda em log-Cholesky, nao em theta cru, e a troca nao e cosmetica. O teste que
  // a mede e a EQUIVARIANCIA: -2logL e invariante sob (K, C) -> (cK, C/c), logo o otimo
  // AJUSTADO tem de ser tambem, e isso e propriedade do caminhante, nao da verossimilhanca.
  // Andando em theta cru, com a mesma K escalada por 25, o ajuste terminava 146 unidades de
  // -2logL longe do de c = 1 (453.55 contra 307.20), com o componente errado por duas
  // ordens de grandeza e converged TRUE nos dois. Em z todo ponto e admissivel, a fronteira
  // det(C_g) = 0 vai para o infinito e a escala de cada direcao acompanha a do proprio
  // componente. Aqui o residuo e um BLOCO de t x t, e nao o escalar do univariado.
  MapaZ mz;
  for (const Grupo& gr : d.modelo.grupos) mz.blocos.push_back(std::make_pair(gr.offset, gr.dim));
  mz.blocos.push_back(std::make_pair(d.modelo.offset_residual, d.t));

  // O certificado final, espelhado do ajustador univariado (o racional completo esta em
  // aireml.cpp): passo relativo pequeno nao prova otimo — AI amortecida pode aceitar um
  // passo minusculo com o gradiente longe de zero. converged exige tambem o decremento
  // de Newton, g' AI^-1 g, ~2x o gap em -2logL perto do otimo. Este laco anda em theta
  // cru, sem pisos e sem log-Cholesky (o limite declarado da etapa 1B), entao as
  // fronteiras que o univariado congela via conjunto ativo aparecem aqui como score que
  // nunca zera. O certificado exclui o que este laco consegue reconhecer: uma variancia
  // diagonal em zero numerico (<= 1e-6 da escala residual) com o score empurrando para
  // baixo, e TODO componente de um grupo cuja C_g esta a 1e-3 (razao das diagonais do
  // Cholesky, o mesmo piso relativo do univariado) da singularidade — este laco nao
  // anda sobre a parede, entao a certificacao e CONDICIONAL a ela e a mensagem diz
  // isso. AI que nao inverte no bloco livre nao certifica nada: conta como recusa.
  const double tol_dec = 2e-4;
  // Quantos componentes o certificado excluiu na ultima avaliacao, para a mensagem. O
  // conjunto sai de PassoZ::na_parede, nas mesmas coordenadas e sobre os mesmos pisos
  // que o passo usa; o teste em theta que havia aqui (lmin < 1e-3 lmax) era a exata
  // desigualdade que o grampo do passo torna falsa, e por isso nunca disparava.
  std::size_t na_fronteira = 0;
  // a escala residual e RECALCULADA a cada ponto, como no univariado: e ela que ancora o
  // piso absoluto, e ela se move durante o ajuste
  auto pecas = [&](const AvaliacaoMT& av) {
    double s = 0.0;
    for (std::size_t tau = 0; tau < d.t; tau++)
      s += theta[d.modelo.offset_residual + (tau * (2 * d.t - tau + 1)) / 2];
    return passo_z(mz, d.modelo.ntheta, theta, av.score, av.ai,
                   std::max(1.0, s / static_cast<double>(d.t)));
  };
  auto conta_parede = [&](const PassoZ& P) {
    std::size_t n = 0;
    for (char c : P.na_parede) n += c;
    return n;
  };

  double lambda = 1e-2;   // multi comeca mais amortecido: as covariancias partem de zero
  int parado = 0;         // iteracoes consecutivas aceitas SEM progresso real
  for (std::size_t it = 1; it <= maxiter; it++) {
    R.iters = it;
    R_CheckUserInterrupt();
    bool aceitou = false;
    const std::size_t ntz = d.modelo.ntheta;
    // z, score e AI em z, pisos e os dois conjuntos ativos, numa chamada compartilhada
    // pelos tres ajustadores. O certificado la embaixo le as MESMAS pecas.
    const PassoZ P = pecas(cur);
    for (int tent = 0; tent < 30 && !aceitou; tent++) {
      Densa m2 = P.ok ? P.az : cur.ai;
      for (std::size_t i = 0; i < m2.nlin; i++) {
        const double di = m2.at(i, i);
        m2.at(i, i) = (di == 0.0) ? lambda : di * (1.0 + lambda);
      }
      if (P.ok)
        for (std::size_t i = 0; i < ntz; i++)
          if (P.congelado[i]) {
            for (std::size_t j = 0; j < ntz; j++) { m2.at(i, j) = 0.0; m2.at(j, i) = 0.0; }
            m2.at(i, i) = 1.0;
          }
      Densa minv;
      try { minv = inv_geral(m2); } catch (const Erro&) { lambda *= 10.0; continue; }
      // a DIRECAO, computada uma vez; o comprimento e que e buscado
      std::vector<double> passo(ntz, 0.0);
      for (std::size_t i = 0; i < ntz; i++)
        for (std::size_t j = 0; j < ntz; j++)
          if (!P.ok || !P.congelado[j]) passo[i] += minv.at(i, j) * (P.ok ? P.sz[j] : cur.score[j]);
      // BUSCA DE COMPRIMENTO na mesma direcao, antes de mexer no amortecimento. Sem ela
      // toda reducao de passo passava por multiplicar lambda por dez, o que encurta E
      // GIRA a direcao: o laco alternava entre dois estados e a verossimilhanca ficava
      // parada (medido: 980 iteracoes com -2logL identico ate a ultima casa).
      for (const double alpha : {1.0, 0.5, 0.25}) {
        std::vector<double> cand;
        if (P.ok) {
          std::vector<double> zn = P.zc;
          for (std::size_t i = 0; i < ntz; i++)
            if (!P.congelado[i]) zn[i] -= alpha * passo[i];
          // PISO na diagonal de cada bloco. Sem ele o passo em z vai ate a fronteira
          // exata (autovalor de G0 medido em 1.5e-14 num ajuste sem sinal genetico
          // nenhum), e ali a inversa do bloco que o MME precisa nao existe mais: os
          // efeitos fixos passam a divergir da mesma GLS montada densa na quinta casa.
          for (std::size_t i = 0; i < ntz; i++)
            if (zn[i] < P.piso[i]) zn[i] = P.piso[i];
          cand = theta_de_z(mz, ntz, zn);
        } else {
          // theta fora do cone (uma C_g ja singular): o passo cru ainda serve de resgate
          cand = theta;
          for (std::size_t i = 0; i < ntz; i++) cand[i] -= alpha * passo[i];
        }
        AvaliacaoMT prox = avalia_mt(d, cand, &cs);
        if (!prox.ok || prox.neg2logl > cur.neg2logl + 1e-9) continue;
        double num = 0.0, den = 0.0;
        for (std::size_t i = 0; i < cand.size(); i++) {
          const double dlt = cand[i] - theta[i];
          num += dlt * dlt;
          den += cand[i] * cand[i];
        }
        // PROGRESSO nao e o mesmo que aceitacao: o teste de aceitacao tem folga de 1e-9
        // e por isso aceita um passo que nao anda nada. Contar esses separadamente e o
        // que impede o laco de virar ponto fixo ate o maxiter.
        const double ganho = cur.neg2logl - prox.neg2logl;
        parado = (ganho > 1e-8 * std::max(1.0, std::fabs(cur.neg2logl))) ? 0 : parado + 1;
        R.reldelta = std::sqrt(num / std::max(den, 1e-300));
        theta = cand;
        cur = std::move(prox);
        lambda = std::max(lambda / 10.0, 1e-10);
        aceitou = true;
      }
      if (!aceitou) lambda *= 10.0;
    }
    if (!aceitou) { R.mensagem = "nenhum passo amortecido melhorou a verossimilhanca"; break; }
    if (verboso) {
      Rprintf("iter %3d  -2logL %.6f  relDelta %.3e\n",
              (int) it, cur.neg2logl, R.reldelta);
      imprime_theta(theta, nomes_theta_do_mt(d), d.modelo);
    }
    if (R.reldelta < tol || parado >= 2) {
      // O CERTIFICADO, nas mesmas coordenadas e sobre os mesmos pisos que o passo. Ele
      // ficou para tras quando o passo foi portado para log-Cholesky, e a assimetria era
      // um defeito de verdade: o teste que havia aqui, lmin < 1e-3 lmax em theta, e a
      // exata desigualdade que o grampo do passo torna FALSA por construcao. No ponto em
      // que o laco para, o grampo deixa lmin = 1e-3 lmax e o `<` estrito nunca dispara,
      // de modo que a direcao grampeada ficava congelada no passo e cobrada integralmente
      // no certificado. Medido: decremento parado em 3.25e+03 com -2logL identico por 980
      // iteracoes, num ponto que uma descida por coordenada melhorava em 0.200 unidade.
      //
      // converged exige as duas coisas, como no univariado: passo pequeno E decremento de
      // Newton dos componentes LIVRES dentro da tolerancia.
      const PassoZ Pc = pecas(cur);
      na_fronteira = conta_parede(Pc);
      R.decremento = decremento_z(Pc);
      if (R.decremento < tol_dec) {
        R.convergiu = true;
        break;
      }
      // A recusa nao encerra o ajuste: o passo pode voltar a andar. Mas ele so volta se o
      // amortecimento for solto — depois de uma sequencia de aceites lambda esta no piso
      // e a direcao e sempre a mesma. `parado` conta os aceites sem ganho real e desiste
      // quando nem soltar o amortecimento adianta.
      if (parado >= 6) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.3g", R.decremento);
        R.mensagem = std::string("did NOT converge: the step stopped making progress, ") +
            "but the Newton decrement of the off-boundary components is " + buf +
            " (tolerance 2e-4 on the -2logL scale), so the point is a stall and not a "
            "certified optimum. Try a different start=, and profile any component the "
            "message reports at a boundary";
        break;
      }
      lambda = 1e-2;
    }
  }
  if (std::isnan(R.decremento)) {
    const PassoZ Pf = pecas(cur);
    na_fronteira = conta_parede(Pf);
    R.decremento = decremento_z(Pf);
  }
  if (na_fronteira > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(na_fronteira) + " component(s) resting at a covariance boundary "
        "excluded from the convergence certificate: a component pinned at a boundary "
        "points out of the cone by construction, so its gradient never vanishes and the "
        "certificate is CONDITIONAL on the pinning; read those components with care";
  R.se.assign(d.modelo.ntheta, std::nan(""));
  try {
    Densa inv = inv_geral(cur.ai);
    for (std::size_t k = 0; k < d.modelo.ntheta; k++) {
      const double v = 2.0 * inv.at(k, k);
      if (v > 0.0) R.se[k] = std::sqrt(v);
    }
  } catch (const Erro&) {}

  if (!R.convergiu && R.mensagem.empty()) {
    char buf[128];
    std::snprintf(buf, sizeof(buf), " (relDelta %.3g, Newton decrement %.3g against the "
                  "2e-4 tolerance)", R.reldelta, R.decremento);
    R.mensagem = "stopped at " + std::to_string(maxiter) + " iteration(s) without a "
        "certified optimum" + buf + ": this model asks for more iterations. Raise maxiter=";
  }
  if (cur.fora_do_padrao > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(cur.fora_do_padrao) + " leitura(s) fora do padrao do fator";

  R.theta = theta;
  R.neg2logl = cur.neg2logl;
  R.solucao = cur.solucao;
  R.pev = cur.pev;
  R.fora_do_padrao = cur.fora_do_padrao;
  return R;
}

}  // namespace br
