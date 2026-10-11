// AR(1)/CAR(1) (Wade e Quaas, 1993), parte 2: desenho e estimador, para t caracteristicas (t = 1 e o caminho
// classico, e os gates antigos sao a regressao disto).

#include "mme.h"

// so para o TRACE opcional e a interrupcao: nenhuma conta usa o R
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <algorithm>

namespace br {

DesenhoAR monta_desenho_ar1(Modelo m, const std::vector<std::string>& alvos,
                            const Tabela& tab, const Pedigree* ped,
                            const std::string& col_sujeito, const std::string& col_tempo,
                            const std::vector<KernelDecl>* kernels) {
  DesenhoAR d;
  d.nlin = tab.nlin;
  d.alvos = alvos;
  d.t = alvos.size();
  const std::size_t t = d.t;
  if (t == 0) throw Erro("the model needs a trait");
  if (t > 32) throw Erro("more than 32 traits in one AR(1) fit");

  // y n x t. REGISTRO COMPLETO por enquanto: com AR(1) a ausencia parcial quebra a
  // separabilidade Gamma (x) R0 (a condicional de um padrao nao e mais um kron), e
  // entrar nisso sem gates proprios seria fingir capacidade. Registro com qualquer
  // caracteristica ausente sai inteiro.
  d.y = Densa(d.nlin, t);
  d.usa.assign(d.nlin, 1);
  d.obs.assign(d.nlin * t, 1);
  for (std::size_t tau = 0; tau < t; tau++) {
    std::vector<double> col = tab.numerico(alvos[tau]);
    for (std::size_t i = 0; i < d.nlin; i++) {
      const bool falta = !std::isfinite(col[i]) ||
          (m.tem_ausente && std::fabs(col[i] - m.codigo_ausente) < 1e-9);
      if (falta) d.obs[i * t + tau] = 0;
      d.y.at(i, tau) = falta ? 0.0 : col[i];
    }
  }
  // o registro so sai quando NENHUMA caracteristica foi observada; celula ausente com outra
  // observada no mesmo registro fica e vira coluna mv (ver DesenhoAR::obs)
  for (std::size_t i = 0; i < d.nlin; i++) {
    bool alguma = false;
    for (std::size_t tau = 0; tau < t; tau++) if (d.obs[i * t + tau]) { alguma = true; break; }
    if (!alguma) { d.usa[i] = 0; d.n_incompletos++; }
  }

  d.tempo = tab.numerico(col_tempo);
  // TEMPO NAO FINITO = registro sem posicao na serie. Antes a linha entrava, o comparador do
  // std::sort recebia NaN (comportamento indefinido) e o ajuste terminava em "theta inicial
  // inadmissivel". Agora sai, e a mensagem conta quantas.
  for (std::size_t i = 0; i < d.nlin; i++)
    if (d.usa[i] && !std::isfinite(d.tempo[i])) { d.usa[i] = 0; d.n_sem_tempo++; }
  // SUJEITO AUSENTE e ERRO, pela regra da coluna de nivel (rotulos_ausentes, em modelo.cpp),
  // em toda linha. O sujeito e o rotulo que junta os registros numa serie, e o NA numerico
  // virava o sujeito "nan" (o texto vazio, o sujeito ""): as linhas sem sujeito formavam UMA
  // serie, com o AR(1) correndo entre registros de animais diferentes. Medido em f2c0d69: 10
  // sujeitos com 2 linhas NA davam n_subjects 11, sem aviso.
  {
    const std::string falta = rotulos_ausentes(tab, col_sujeito);
    if (!falta.empty())
      throw Erro("subject: " + falta + ": a record without a subject belongs to no known "
                 "series, and those rows would form one series; drop those rows or fill "
                 "them");
  }
  std::vector<std::string> suj = tab.rotulos(col_sujeito);

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
      d.kinv_escala.push_back(0.0);
    } else if (g.estrutura == Estrutura::Declarada) {
      kinv_declarada(m, g, kernels, d.kinv, d.kinv_logdet, d.kinv_escala);
    } else {
      d.kinv.push_back(Csc());
      d.kinv_logdet.push_back(0.0);
      d.kinv_escala.push_back(0.0);
    }
  }

  // X (esparso, fixos.cpp) e aleatorios como no caminho comum
  FixosCandidatos fx = colunas_fixas(m, tab);
  // AS LINHAS SEM PAR CAEM ANTES DO POSTO DE X, como em model() (mme.cpp). Antes caiam
  // depois: um CG cujos animais nao estao no pedigree ficava com coluna em X, o posto nao
  // a via como dependente, e a coluna chegava vazia na montagem. model() derrubava o nivel
  // e ajustava; model_ar1() devolvia "did NOT converge" com decremento 0,46 e model_mt()
  // "SINGULAR", dois diagnosticos que apontavam para o lugar errado.
  // niveis: do pedigree quando ha parentesco, DA K quando declarada, da uniao das colunas
  // num grupo sem estrutura de varios termos. Todo nivel do conjunto ganha equacao, com ou
  // sem registro, exatamente como em model() (monta_aleatorios, em mme.cpp).
  d.aleatorios = monta_aleatorios(m, tab, niveis_ped, kernels);
  {
    std::vector<DesenhoTermo*> pa;
    for (DesenhoTermo& a : d.aleatorios) pa.push_back(&a);
    casa_niveis_nulos(m, pa, kern_nulos, tab, d.nlin);
  }
  for (const DesenhoTermo& a : d.aleatorios)
    for (std::size_t i = 0; i < d.nlin; i++)
      if (!a.casou[i]) d.usa[i] = 0;
  // POSTO DE X SOBRE AS LINHAS QUE ENTRAM, nao sobre a tabela inteira. E o mesmo conserto
  // que o univariado ja tem (mme.cpp, mesma marca). Um nivel fixo cujos registros TODOS
  // sairam continua com coluna nao-nula na tabela, entao um posto medido sobre tudo o
  // mantinha; na montagem, que so anda nas linhas usadas, essa coluna nao recebe nada, a
  // matriz de coeficientes ganha coluna VAZIA, e a Cholesky morre com "nao
  // positiva-definida" tres passos longe da causa. Medindo o posto onde o modelo de fato
  // vive, a coluna sai como dependente e e REPORTADA em dropped_x.
  // a referencia global entre os niveis com registro de cada caracteristica, como no
  // model_mt() (o racional esta em monta_desenho_mt): aqui um par sem registro e recusado,
  // e a recusa passa a nomear o nivel que de fato nao tem registro, e nao um vizinho dele
  // que virou a referencia daquela caracteristica. Com t = 1, ou sem nivel incompleto, a
  // ordem e a de sempre. O posto e o esparso (posto_referencia_completa, fixos.cpp).
  // As covariaveis padronizadas antes do posto (padroniza_covariaveis, fixos.cpp), e a volta de
  // b, a mesma em toda caracteristica (aqui toda coluna tem equacao em todas). Uma aninhada
  // cujo centro nao volta pela classe e refeita so escalada, com o posto de novo.
  std::vector<std::size_t> fica, sai, sem_centro;
  for (;;) {
    padroniza_covariaveis(fx, d.usa, sem_centro);
    {
      std::vector<std::string> motivo;
      posto_referencia_completa(fx, d.usa, d.obs, t, 1e-9, fica, sai, motivo);
      d.saiu_motivo = std::move(motivo);
    }
    std::vector<char> mantida(fx.x.ncol, 0);
    for (std::size_t j : fica) mantida[j] = 1;
    const std::vector<std::size_t> falha = aninhadas_sem_volta(fx, mantida, d.usa);
    if (falha.empty()) {
      d.volta_x.assign(t, monta_volta(fx, fica, mantida, d.usa));
      break;
    }
    sem_centro.insert(sem_centro.end(), falha.begin(), falha.end());
    fx = colunas_fixas(m, tab);
  }
  d.x = fixas_por_linha(fx.x, fica);
  for (std::size_t j : fica) d.nomes_x.push_back(fx.nomes[j]);
  for (std::size_t j : sai) d.saiu_x.push_back(fx.nomes[j]);

  // sujeitos: registros usaveis agrupados e ordenados pelo tempo
  {
    std::unordered_map<std::string, std::vector<std::size_t>> mapa;
    for (std::size_t i = 0; i < d.nlin; i++)
      if (d.usa[i]) mapa[suj[i]].push_back(i);
    for (auto& [nome, regs] : mapa) {
      std::sort(regs.begin(), regs.end(),
                [&](std::size_t a, std::size_t b) { return d.tempo[a] < d.tempo[b]; });
      for (std::size_t k = 0; k + 1 < regs.size(); k++)
        if (d.tempo[regs[k + 1]] - d.tempo[regs[k]] <= 0.0)
          throw Erro("subject '" + nome + "' has two records at the SAME time. With AR(1) "
                     "time identifies the record: a tie gives correlation 1 and a singular "
                     "Gamma. Either the data is wrong, or simultaneous repetition calls for "
                     "permanent environment, not AR(1)");
      d.sujeitos.push_back(regs);
    }
  }
  if (d.n_usadas() == 0) throw Erro("no row enters the analysis");
  // o termo do -2logL da padronizacao: toda coluna tem equacao em toda caracteristica aqui
  d.logdet_escala_x = static_cast<double>(t) * logdet_escala(d.volta_x[0].escala);

  // POSTO DE X POR CARACTERISTICA sobre as celulas OBSERVADAS. Com celula ausente, um nivel
  // fixo pode ficar sem nenhuma observacao numa caracteristica, e a coluna dele naquela
  // caracteristica so encontra celulas mv, que a absorvem inteira: a MME fica singular.
  // Nesta etapa isso e RECUSADO nomeando o par; derrubar o par (como o model_mt() faz)
  // pede a numeracao por caracteristica no AR(1), que segue aberta.
  if (t > 1) {
    const Csc xc = fixas_por_coluna(d.x);
    for (std::size_t tau = 0; tau < t; tau++) {
      std::vector<char> lin_tau(d.nlin, 0);
      for (std::size_t i = 0; i < d.nlin; i++) lin_tau[i] = d.usa[i] && d.obs[i * t + tau];
      std::vector<std::size_t> fica_t, sai_t;
      posto_esparso(xc, lin_tau, 1e-9, fica_t, sai_t, false, nullptr, &d.nomes_x);
      if (!sai_t.empty())
        throw Erro("the fixed column '" + d.nomes_x[sai_t[0]] + "' is not estimable for "
                   "trait '" + d.alvos[tau] + "' from the records where that trait was "
                   "observed (" + std::to_string(sai_t.size()) + " column(s) in all): with "
                   "missing cells a level can end up with no observation of one trait. "
                   "Merge that level, or fit that trait's records with model_mt()");
    }
  }

  // layout: grupos com dim expandida por t (como na multi) + vech(R0) + rho
  std::size_t off = 0;
  for (Grupo& g : m.grupos) {
    g.dim *= t;
    g.offset = off;
    g.nparam = g.dim * (g.dim + 1) / 2;
    off += g.nparam;
  }
  m.offset_residual = off;
  m.ntheta = off + t * (t + 1) / 2 + 1;
  d.offset_s2e = off;
  d.offset_rho = off + t * (t + 1) / 2;
  // grade inteira? (ver o campo em mme.h). Basta olhar os intervalos CONSECUTIVOS: se
  // todos sao inteiros, toda soma deles tambem e.
  d.tempo_inteiro = true;
  for (const std::vector<std::size_t>& s : d.sujeitos)
    for (std::size_t k = 0; k + 1 < s.size(); k++) {
      const double g = std::fabs(d.tempo[s[k + 1]] - d.tempo[s[k]]);
      if (std::fabs(g - std::floor(g + 0.5)) > 1e-9) d.tempo_inteiro = false;
    }

  d.modelo = std::move(m);

  // esqueleto: linhas de W em colunas globais COM caracteristica, e colunas por sujeito.
  // Convencao identica a da multi: X coluna j da caracteristica tau em j*t + tau; termo
  // com bloco em colbase, coeficiente (tau, ct), nivel nv em
  // colbase + (tau*n_coef + ct)*nl + nv, e colbase avanca z.ncol*t por termo.
  {
    std::vector<std::size_t> off_grupo(d.modelo.grupos.size());
    std::size_t acc = d.x.ncol * t;
    for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
      off_grupo[g] = acc;
      acc += d.largura(g);
    }
    d.lw.assign(d.nlin, {});
    for (std::size_t r = 0; r < d.nlin; r++)
      for (std::size_t ex = d.x.linptr[r]; ex < d.x.linptr[r + 1]; ex++) {
        const std::size_t j = d.x.coluna[ex];
        const double v = d.x.valor[ex];
        if (v == 0.0) continue;
        for (std::size_t tau = 0; tau < t; tau++)
          d.lw[r].push_back({static_cast<std::uint32_t>(j * t + tau),
                             static_cast<std::uint32_t>(tau), v});
      }
    for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
      std::size_t colbase = off_grupo[g];
      for (std::size_t tm : d.modelo.grupos[g].termos)
        for (const DesenhoTermo& a : d.aleatorios)
          if (a.termo == tm) {
            for (std::size_t c = 0; c < a.z.ncol; c++) {
              const std::size_t ct = c / a.n_niveis, nv = c % a.n_niveis;
              for (std::size_t p = a.z.colptr[c]; p < a.z.colptr[c + 1]; p++)
                for (std::size_t tau = 0; tau < t; tau++)
                  d.lw[a.z.linha[p]].push_back(
                      {static_cast<std::uint32_t>(
                           colbase + (tau * a.n_coef + ct) * a.n_niveis + nv),
                       static_cast<std::uint32_t>(tau), a.z.valor[p]});
            }
            colbase += a.z.ncol * t;
          }
    }
    // uma coluna mv por celula ausente de registro usado, DEPOIS dos grupos, para que os
    // offsets de b, EBV e PEV nao se movam
    d.off_mv = acc;
    d.n_mv = 0;
    if (t > 1)
      for (std::size_t r = 0; r < d.nlin; r++) {
        if (!d.usa[r]) continue;
        for (std::size_t tau = 0; tau < t; tau++)
          if (!d.obs[r * t + tau]) {
            d.lw[r].push_back({static_cast<std::uint32_t>(d.off_mv + d.n_mv),
                               static_cast<std::uint32_t>(tau), 1.0});
            d.n_mv++;
          }
      }
    d.cols_suj.assign(d.sujeitos.size(), {});
    for (std::size_t s = 0; s < d.sujeitos.size(); s++) {
      std::vector<std::uint32_t>& cs = d.cols_suj[s];
      for (std::size_t r : d.sujeitos[s])
        for (const auto& e : d.lw[r]) cs.push_back(e.col);
      std::sort(cs.begin(), cs.end());
      cs.erase(std::unique(cs.begin(), cs.end()), cs.end());
    }
  }
  return d;
}

std::vector<std::string> nomes_theta_ar1(const DesenhoAR& d) {
  const std::size_t t = d.t;
  std::vector<std::string> out(d.modelo.ntheta);
  if (t == 1) {
    std::vector<std::string> base = d.modelo.nomes_theta();
    base.resize(d.modelo.ntheta);
    out = base;
    out[d.offset_s2e] = "var(residual)";
  } else {
    // como na multi: coeficiente rotulado "termo[c]@alvo". A ORDEM segue o layout do Z,
    // termo por fora e caracteristica por dentro (colbase avanca z.ncol*t por termo). Andava
    // ao contrario, e com dois termos no mesmo grupo os componentes saiam com o nome de
    // outro; ver nomes_theta_mt() em multitrait2.cpp, onde isto foi medido.
    for (const Grupo& g : d.modelo.grupos) {
      std::vector<std::string> rot;
      for (std::size_t tm : g.termos)
        for (std::size_t tau = 0; tau < t; tau++) {
          const Termo& te = d.modelo.termos[tm];
          for (std::size_t c = 0; c < te.n_coef(); c++) {
            std::string b = te.nome;
            if (te.n_coef() > 1) b += "[" + std::to_string(c) + "]";
            rot.push_back(b + "@" + d.alvos[tau]);
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
        const std::size_t k = d.offset_s2e + (j * (2 * t - j + 1)) / 2 + (i - j);
        out[k] = (i == j) ? "var(res@" + d.alvos[i] + ")"
                          : "cov(res@" + d.alvos[i] + ",res@" + d.alvos[j] + ")";
      }
  }
  out[d.offset_rho] = "rho(residual)";
  return out;
}

// A partida automatica: variancias de var(y) por caracteristica, R0 diagonal, e o rho
// deslocado de zero (o racional esta no corpo).
static std::vector<double> partida_ar1(const DesenhoAR& d) {
  const std::size_t t = d.t;
  std::vector<double> theta(d.modelo.ntheta, 0.0);
  {
    std::vector<double> var_t(t, 1.0);
    for (std::size_t tau = 0; tau < t; tau++) {
      double soma = 0.0, soma2 = 0.0;
      std::size_t n = 0;
      for (std::size_t i = 0; i < d.nlin; i++)
        if (d.usa[i] && d.obs[i * t + tau]) {
          soma += d.y.at(i, tau); soma2 += d.y.at(i, tau) * d.y.at(i, tau); n++;
        }
      const double media = soma / std::max<std::size_t>(n, 1);
      double v = (soma2 - n * media * media) / std::max(1.0, static_cast<double>(n - 1));
      if (!(v > 0.0) || !std::isfinite(v)) v = 1.0;
      var_t[tau] = v;
    }
    const double frac = 0.5 / static_cast<double>(d.modelo.grupos.size());
    std::size_t ig = 0;
    for (const Grupo& g : d.modelo.grupos) {
      // a mesma equivariancia da partida do multicaracter; o racional esta la
      double esc = 1.0;
      if (g.estrutura == Estrutura::Declarada && ig < d.kinv.size() && d.kinv[ig].ncol > 0) {
        const double gm = escala_kernel(d.kinv_escala, d.kinv_logdet, d.kinv, ig);
        if (std::isfinite(gm) && gm > 0.0) esc = gm;
      }
      ig++;
      // na ordem do layout, termo por fora (o mesmo conserto de partida_mt)
      const std::vector<std::size_t> tr = traco_das_coordenadas(d.modelo, g, t);
      for (std::size_t i = 0; i < g.dim; i++)
        theta[g.offset + (i * (2 * g.dim - i + 1)) / 2] =
            frac * var_t[i < tr.size() ? tr[i] : 0] / esc;
    }
    for (std::size_t tau = 0; tau < t; tau++) {
      const std::size_t k = d.offset_s2e + (tau * (2 * t - tau + 1)) / 2;
      theta[k] = 0.5 * var_t[tau];
    }
    // rho NAO parte de zero. Gamma(dt) = rho^dt tem dGamma/drho = dt rho^(dt-1), que em
    // rho = 0 vale zero para todo dt > 1: se nenhum par de tempos dentro de um sujeito
    // difere de exatamente 1, a linha inteira da AI e o score do rho nascem nulos, o
    // caminhante nunca sai de zero e o ajuste devolve rho = 0 com converged = TRUE e
    // SE = NaN. Medido na mesma serie simulada com rho = 0.6: espacamento 1 recupera
    // 0.563; espacamentos 2 e 7 devolvem 0.000000 com converged TRUE. Medir em dias ou em
    // semanas mudava a resposta, o que e um artefato da unidade e nao um resultado. Em
    // dt fracionario o problema e o oposto e igualmente fatal: a derivada diverge em
    // rho = 0. Uma partida deslocada resolve os dois, e quem quiser outra usa start=.
    // A partida do rho e EQUIVARIANTE na unidade do tempo. Gamma(dt) = rho^dt, entao a
    // mesma serie medida em dias ou em semanas e o mesmo modelo com rho reparametrizado, e
    // uma partida fixa nao e o mesmo ponto nos dois casos. Fixa-se a CORRELACAO no
    // espacamento tipico, e nao o rho: rho0 = 0.3^(1/dt), com dt a mediana dos intervalos
    // entre registros consecutivos do mesmo sujeito.
    double rho0 = 0.3;
    {
      std::vector<double> gaps;
      for (const std::vector<std::size_t>& s : d.sujeitos)
        for (std::size_t k = 0; k + 1 < s.size(); k++) {
          const double g = std::fabs(d.tempo[s[k + 1]] - d.tempo[s[k]]);
          if (g > 0.0 && std::isfinite(g)) gaps.push_back(g);
        }
      if (!gaps.empty()) {
        std::nth_element(gaps.begin(), gaps.begin() + gaps.size() / 2, gaps.end());
        const double dtm = gaps[gaps.size() / 2];
        if (dtm > 0.0) rho0 = std::pow(0.3, 1.0 / dtm);
      }
    }
    theta[d.offset_rho] = rho0;
  }
  return theta;
}

// Os otimos univariados AR(1) de cada caracteristica (layout do model_ar1() de uma:
// grupos, var(residual), rho) viram a partida do multicaracter: grupos e R0 na diagonal por
// caracteristica, covariancias entre caracteristicas zero, e o rho medio. A mesma regra de
// fronteira do multicaracter (partida_univariados_mt): a variancia que o univariado deixou
// abaixo de 1e-4 da partida automatica parte de 1% dela.
std::vector<double> partida_univariados_ar1(const DesenhoAR& d,
                                            const std::vector<std::vector<double> >& uni) {
  const std::size_t t = d.t;
  std::size_t n_u = 2;
  for (const Grupo& g : d.modelo.grupos) {
    const std::size_t du = g.dim / t;
    n_u += du * (du + 1) / 2;
  }
  for (const std::vector<double>& u : uni)
    if (u.size() != n_u) return {};
  std::vector<double> theta(d.modelo.ntheta, 0.0);
  if (!grupos_de_univariados(d.modelo, t, uni, theta)) return {};
  double soma_rho = 0.0;
  for (std::size_t tau = 0; tau < t; tau++) {
    const double v = uni[tau][n_u - 2], r = uni[tau][n_u - 1];
    if (!std::isfinite(v) || !(std::fabs(r) < 1.0)) return {};
    theta[d.offset_s2e + (tau * (2 * t - tau + 1)) / 2] = v;
    soma_rho += r;
  }
  theta[d.offset_rho] = soma_rho / static_cast<double>(t);
  if (!d.tempo_inteiro && theta[d.offset_rho] < 0.0) return {};
  const std::vector<double> auto_ = partida_ar1(d);
  for (std::size_t k = 0; k < theta.size(); k++)
    if (k != d.offset_rho && auto_[k] > 0.0 && theta[k] < 1e-4 * auto_[k])
      theta[k] = 1e-2 * auto_[k];
  return theta;
}

// A CAMINHADA do AR(1), de uma partida ate um veredito; ajusta_ar1 (no fim do arquivo) a
// chama e, com um bloco na parede, explora de dentro do cone (explora_fronteira, aireml.cpp).
// O gancho BREEDINGR_TESTE_SEM_PASSO e o do multicaracter (o racional esta em
// multitrait2.cpp): recusa todo passo AI e todo EM, so para o portao do D2.
// `it0`: as iteracoes ja gastas por uma caminhada anterior do mesmo ajuste (caminha_mt)
static AjusteMT caminha_ar1(DesenhoAR& d, const std::vector<double>* theta0,
                            std::size_t maxiter, double tol, bool verboso, std::size_t it0) {
  AjusteMT R;
  R.iters = it0;
  const std::size_t t = d.t;
  const bool sem_passo = std::getenv("BREEDINGR_TESTE_SEM_PASSO") != nullptr;
  // PADRONIZACAO INTERNA POR CARACTERISTICA, a mesma do multicaracter (o racional e as
  // medidas estao em ajusta_mt, src/multitrait2.cpp). Aqui ela pesa ate com t = 1: o rho nao
  // tem unidade e as variancias tem, entao o relDelta, que soma os dois, e o corte da AI
  // dependiam da unidade de y mesmo num ajuste de uma caracteristica.
  std::vector<double> esc_exato;
  const std::vector<double> esc = desvios_tracos(d.y, d.usa, d.obs, t, &esc_exato);
  const std::vector<double> fator = fatores_escala(d.modelo, t, d.offset_s2e, esc);
  // as unidades do corte de singularidade, com o desvio padrao exato (aviso_informacao)
  const std::vector<double> fator_corte = fatores_escala(d.modelo, t, d.offset_s2e, esc_exato);
  // -2logL original = padronizado + jac, jac = 2 soma_tau (n_tau - p) log s_tau com n_tau
  // as celulas observadas e p as colunas de X (as colunas mv das celulas ausentes cancelam
  // as proprias celulas): so para o TRACE
  double jac = 0.0;
  for (std::size_t tau = 0; tau < t; tau++) {
    double n_tau = 0.0;
    for (std::size_t i = 0; i < d.nlin; i++) if (d.usa[i] && d.obs[i * t + tau]) n_tau++;
    jac += 2.0 * (n_tau - static_cast<double>(d.x.ncol)) * std::log(esc[tau]);
  }
  auto bruto = [&](const std::vector<double>& th) {
    std::vector<double> o = th;
    for (std::size_t k = 0; k < o.size() && k < fator.size(); k++) o[k] *= fator[k];
    return o;
  };
  PadronizaY guarda(d.y, d.usa, d.obs, esc, t);
  // partida: start= do usuario (em unidades originais, levado para as padronizadas), ou a
  // automatica, que ja nasce nas padronizadas
  std::vector<double> theta;
  if (theta0) {
    theta = *theta0;
    for (std::size_t k = 0; k < theta.size() && k < fator.size(); k++) theta[k] /= fator[k];
  } else {
    theta = partida_ar1(d);
  }
  R.partida = bruto(theta);

  CacheSimbolica cs;
  AvaliacaoAR cur = avalia_ar1(d, theta, &cs);
  if (cs.pronto) { R.bloco_denso = bloco_denso_simbolico(cs.sb); R.colunas_fator = cs.sb.n; }
  if (!cur.ok) { R.mensagem = "the starting theta is INADMISSIBLE: some covariance matrix is not positive-definite, or rho is outside (-1, 1)"; return R; }

  // O PASSO anda em log-Cholesky por bloco, como o univariado e o multicaracter. Aqui o
  // residuo e o bloco R0 de t x t e o rho anda em atanh, o que troca a parede |rho| = 1
  // por um infinito: um passo que antes era rejeitado por sair de (-1, 1), com o
  // amortecimento subindo uma ordem de grandeza a cada tentativa, passa a ser um passo
  // grande e admissivel na coordenada certa.
  MapaZ mz;
  for (const Grupo& gr : d.modelo.grupos) mz.blocos.push_back(std::make_pair(gr.offset, gr.dim));
  mz.blocos.push_back(std::make_pair(d.offset_s2e, t));
  mz.correlacoes.push_back(d.offset_rho);

  // O certificado final, espelhado do ajustador univariado (o racional completo esta em
  // aireml.cpp): passo relativo pequeno nao prova otimo. converged exige tambem o
  // decremento de Newton, g' AI^-1 g, ~2x o gap em -2logL perto do otimo. Ele e computado
  // nas MESMAS coordenadas z do passo e sobre os MESMOS pisos, por passo_z/decremento_z:
  // manter dois conjuntos ativos, um por ajustador e em espacos diferentes, foi um defeito
  // medido. Um bloco descansando no piso sai do certificado inteiro, e a certificacao fica
  // CONDICIONAL a esse grampo, dito na mensagem. O rho anda em atanh e nunca esta preso.
  const double tol_dec = 2e-4;
  // Quantos componentes o certificado excluiu, para a mensagem. O conjunto sai de
  // PassoZ::na_parede, nas mesmas coordenadas e sobre os mesmos pisos que o passo usa.
  std::size_t na_fronteira = 0;
  auto pecas = [&](const AvaliacaoAR& av) {
    double s = 0.0;
    for (std::size_t tau = 0; tau < t; tau++)
      s += theta[d.offset_s2e + (tau * (2 * t - tau + 1)) / 2];
    return passo_z(mz, d.modelo.ntheta, theta, av.score, av.ai,
                   std::max(1.0, s / static_cast<double>(t)));
  };
  // RESGATE EM, o que faltava para os espelhos terem a maquinaria completa do univariado.
  // O passo EM e multiplicativo, fica no cone e nao encolhe na fronteira, entao ele anda
  // exatamente onde o passo AI amortecido trava. Aqui ele move so os GRUPOS: o passo M do
  // R0 multivariado, e o do rho no AR(1), nao existem em forma fechada, e por isso a
  // proposta deixa esses componentes onde estao. Isso e um EM GENERALIZADO e continua
  // valido, porque cada candidato so entra se BAIXAR a verossimilhanca.
  //
  // Cada candidato passa pelo grampo em z antes de ser avaliado: um EM que encoste
  // numericamente na singularidade daria um -2logL sem sentido, com det(C_g) em underflow.
  // Ate 15 passos, contra os 50 do univariado, porque aqui cada avaliacao custa mais.
  auto resgate_em = [&](AvaliacaoAR& atual) -> bool {
    bool ganhou = false;
    if (sem_passo) return false;
    for (int e = 0; e < 15; e++) {
      if (atual.em_theta.size() != d.modelo.ntheta) break;
      std::vector<double> ze;
      if (!z_de_theta(mz, d.modelo.ntheta, atual.em_theta, ze)) break;
      const PassoZ Pe = pecas(atual);
      if (Pe.ok)
        for (std::size_t i = 0; i < d.modelo.ntheta; i++)
          if (ze[i] < Pe.piso[i]) ze[i] = Pe.piso[i];
      std::vector<double> cand = theta_de_z(mz, d.modelo.ntheta, ze);
      AvaliacaoAR prox = avalia_ar1(d, cand, &cs);
      // ganho RELEVANTE, e nao qualquer migalha. Este EM move so os grupos: quando o
      // que falta esta no R0 ele nao alcanca, e aceitar um ganho de 1e-8 so desvia um
      // caminho saudavel para um ponto de onde o passo AI nao sai (medido: uma celula
      // que certificava com decremento 1.3e-04 passou a parar em 8.4e-04). Com um piso
      // relativo, o EM entra onde ele de fato resolve, a fronteira, onde o passo
      // amortecido trava, e sai da frente onde nao resolve.
      const double ganho_em = atual.neg2logl - prox.neg2logl;
      if (!prox.ok || ganho_em < 1e-6 * std::max(1.0, std::fabs(atual.neg2logl))) break;
      double num = 0.0, den = 0.0;
      for (std::size_t i = 0; i < cand.size(); i++) {
        const double dlt = cand[i] - theta[i];
        num += dlt * dlt;
        den += cand[i] * cand[i];
      }
      R.reldelta = std::sqrt(num / std::max(den, 1e-300));
      theta = cand;
      atual = std::move(prox);
      ganhou = true;
    }
    return ganhou;
  };
  // o conjunto excluido do ultimo certificado, para a mensagem e para a exploracao
  std::vector<char> parede_final;
  auto conta_parede = [&](const PassoZ& P) {
    parede_final = P.na_parede;
    std::size_t n = 0;
    for (char c : P.na_parede) n += c;
    return n;
  };
  // A SONDA DE KKT pela verossimilhanca, a mesma do multicaracter (o racional esta em
  // candidatos_reentrada, src/aireml.cpp, e o uso em ajusta_mt, src/multitrait2.cpp)
  int n_reentradas = 0;
  double ganho_kkt = 0.0;   // quanto o melhor ponto da ultima sonda recusada baixava
  // A SONDA DA FACE (sonda_face, aireml.cpp), depois da do cone: o certificado tira da conta
  // o bloco inteiro na parede, e um ponto da face que nao e o otimo dela passava junto. Os
  // movimentos sobre a face tem limite proprio, 20, porque cada um e um passo de Newton das
  // diferencas finitas e o caminhante volta a andar depois dele.
  int n_face = 0;
  // a avaliacao em z das sondas da face: o -2logL e o score levado a z pelo Jacobiano
  auto avalia_z = [&](const std::vector<double>& zc, double& v, std::vector<double>& gz) {
    const std::size_t nt = d.modelo.ntheta;
    const AvaliacaoAR av = avalia_ar1(d, theta_de_z(mz, nt, zc), &cs);
    if (!av.ok) return false;
    v = av.neg2logl;
    const Densa J = jacobiano_z(mz, nt, zc);
    gz.assign(nt, 0.0);
    for (std::size_t m = 0; m < nt; m++)
      for (std::size_t k = 0; k < nt; k++) gz[m] += J.at(k, m) * av.score[k];
    return true;
  };
  auto sonda_da_face = [&](const PassoZ& P, bool pode_mover) -> Kkt {
    const std::size_t nt = d.modelo.ntheta;
    std::vector<double> zn;
    double fn = 0.0;
    if (!sonda_face(mz, nt, P, cur.neg2logl, avalia_z, zn, fn)) return Kkt::otimo;
    if (!(fn < cur.neg2logl - tol_reentrada)) return Kkt::otimo;
    ganho_kkt = cur.neg2logl - fn;
    if (!pode_mover || n_face >= 20) return Kkt::limite;
    const std::vector<double> cand = theta_de_z(mz, nt, zn);
    AvaliacaoAR av = avalia_ar1(d, cand, &cs);
    if (!av.ok) return Kkt::otimo;
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < cand.size(); i++) {
      num += (cand[i] - theta[i]) * (cand[i] - theta[i]);
      den += cand[i] * cand[i];
    }
    R.reldelta = std::sqrt(num / std::max(den, 1e-300));
    if (verboso)
      Rprintf("      along the boundary face: -2logL %.6f -> %.6f (the face was not at its "
              "optimum)\n", cur.neg2logl + jac, av.neg2logl + jac);
    theta = cand;
    cur = std::move(av);
    n_face++;
    return Kkt::moveu;
  };
  // A LAPIDACAO DA FACE (o racional esta em tol_lapida, mme.h), so DEPOIS de um certificado
  // aceito com algum bloco na parede: a cada rodada, primeiro o ponto com as covariancias de
  // uma caracteristica no zero zeradas (descorrelada_na_parede, aireml.cpp), depois o passo de
  // Newton da sonda da face; cada um entra se baixa o -2logL por mais que tol_lapida E o ponto
  // novo ainda certifica (decremento dos livres abaixo de tol_dec, nas pecas do proprio
  // ponto). O primeiro ponto que nao certifica fica de fora e o ajuste fica no ultimo que
  // certificou: a lapidacao nunca troca um TRUE por um FALSE e nunca sobe o -2logL. Tentado
  // antes e desfeito: deixar o caminhante andar depois de cada movimento de lapidacao. Na
  // semente 10 (r12 = -0,25) da varredura da revisao, a partida de var(y) foi de TRUE 3,86e-4
  // acima do melhor para FALSE 1,75e-5 acima: o passo AI levava a diagonal na parede de volta
  // para cima do piso (L_22 / L_11 de 1e-3 a 5e-3), o bloco deixava de estar na parede, o
  // score dele perto da singularidade entrava no decremento (0,008) e o caminhante estagnava.
  auto lapida = [&](PassoZ& Pc) {
    const std::size_t nt = d.modelo.ntheta;
    // um candidato vale se baixa por mais que tol_lapida e o ponto novo ainda certifica
    auto aceita = [&](const std::vector<double>& zn, const char* como) -> bool {
      const std::vector<double> cand = theta_de_z(mz, nt, zn);
      AvaliacaoAR av = avalia_ar1(d, cand, &cs);
      if (!av.ok || !(av.neg2logl < cur.neg2logl - tol_lapida)) return false;
      const std::vector<double> antes = theta;
      theta = cand;   // pecas le o theta corrente
      PassoZ Pn = pecas(av);
      const double dn = Pn.ok ? decremento_z(Pn) : tol_dec;
      if (!(dn < tol_dec)) {
        theta = antes;
        return false;
      }
      if (verboso)
        Rprintf("      polishing at the boundary, %s: -2logL %.6f -> %.6f (inside the 1e-4 "
                "slack; the new point still certifies)\n", como, cur.neg2logl + jac,
                av.neg2logl + jac);
      cur = std::move(av);
      Pc = std::move(Pn);
      R.decremento = dn;
      return true;
    };
    for (int k = 0; k < max_lapidacoes; k++) {
      const std::vector<double> zd = descorrelada_na_parede(mz, nt, Pc);
      if (!zd.empty() && aceita(zd, "a trait at zero decorrelated")) continue;
      std::vector<double> zn;
      double fn = 0.0;
      if (!sonda_face(mz, nt, Pc, cur.neg2logl, avalia_z, zn, fn)) return;
      if (!(fn < cur.neg2logl - tol_lapida)) return;
      if (!aceita(zn, "along the face")) return;
    }
  };
  auto sonda_kkt = [&](const PassoZ& P, bool pode_mover) -> Kkt {
    if (!P.ok) return Kkt::otimo;
    bool alguma = false;
    for (char c : P.na_parede) if (c) alguma = true;
    if (!alguma) return Kkt::otimo;
    const std::size_t nt = d.modelo.ntheta;
    double esc_res = 0.0;
    for (std::size_t tau = 0; tau < t; tau++)
      esc_res += theta[d.offset_s2e + (tau * (2 * t - tau + 1)) / 2];
    const CandidatosKKT cands = candidatos_reentrada(
        mz, nt, P.zc, P.na_parede, std::max(1.0, esc_res / static_cast<double>(t)),
        cur.score, cur.ai);
    std::vector<double> z_melhor;
    AvaliacaoAR av_melhor;
    bool achou = false;
    for (const std::vector<std::vector<double> >* lista : {&cands.direcionais, &cands.coordenadas}) {
      for (const std::vector<double>& zc : *lista) {
        AvaliacaoAR av = avalia_ar1(d, theta_de_z(mz, nt, zc), &cs);
        if (!av.ok) continue;
        if (!achou || av.neg2logl < av_melhor.neg2logl) {
          achou = true;
          z_melhor = zc;
          av_melhor = std::move(av);
        }
      }
      if (achou && av_melhor.neg2logl < cur.neg2logl - tol_reentrada) break;
    }
    if (!achou || !(av_melhor.neg2logl < cur.neg2logl - tol_reentrada))
      return sonda_da_face(P, pode_mover);
    ganho_kkt = cur.neg2logl - av_melhor.neg2logl;
    if (!pode_mover || n_reentradas >= max_reentradas) return Kkt::limite;
    const std::vector<double> cand = theta_de_z(mz, nt, z_melhor);
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < cand.size(); i++) {
      num += (cand[i] - theta[i]) * (cand[i] - theta[i]);
      den += cand[i] * cand[i];
    }
    R.reldelta = std::sqrt(num / std::max(den, 1e-300));
    if (verboso)
      Rprintf("      re-entry into the cone: -2logL %.6f -> %.6f (a component at the "
              "boundary was not at a constrained optimum)\n", cur.neg2logl + jac,
              av_melhor.neg2logl + jac);
    theta = cand;
    cur = std::move(av_melhor);
    n_reentradas++;
    return Kkt::moveu;
  };
  auto msg_limite = [&](const std::string& onde) {
    char buf[96];
    std::snprintf(buf, sizeof(buf), "%.3g", ganho_kkt);
    return std::string("did NOT converge: ") + onde + " a component resting at a covariance "
        "boundary is not at a constrained optimum (a point inside the cone or along the "
        "boundary face lowers the -2logL by " + buf + "), and the walk had already moved "
        "off it the most times allowed (" + std::to_string(max_reentradas) + " re-entries "
        "into the cone, 20 moves along the face). Try a start= inside the cone, for example "
        "the univariate fits with the covariances at zero";
  };

  // a partida alem dos pisos do passo vai para eles (o racional e a medida estao em
  // ajusta_mt)
  if (maxiter > 0) {
    const PassoZ P0 = pecas(cur);
    if (P0.ok) {
      std::vector<double> z0 = P0.zc;
      bool mudou = false;
      for (std::size_t i = 0; i < z0.size(); i++)
        if (z0[i] < P0.piso[i]) { z0[i] = P0.piso[i]; mudou = true; }
      if (mudou) {
        const std::vector<double> th0 = theta_de_z(mz, d.modelo.ntheta, z0);
        AvaliacaoAR a0 = avalia_ar1(d, th0, &cs);
        if (a0.ok) {
          if (verboso)
            Rprintf("start beyond the step floors, raised to them: -2logL %.6f -> %.6f\n",
                    cur.neg2logl + jac, a0.neg2logl + jac);
          theta = th0;
          cur = std::move(a0);
          R.partida = bruto(theta);
        }
      }
    }
  }

  double lambda = 1e-2;
  int parado = 0;         // iteracoes consecutivas aceitas SEM progresso real
  bool saiu_do_laco = false;   // o laco parou por um veredito, e nao por esgotar o maxiter
  for (std::size_t it = it0 + 1; it <= maxiter; it++) {
    R.iters = it;
    R_CheckUserInterrupt();
    bool aceitou = false;
    const std::size_t ntz = d.modelo.ntheta;
    PassoZ P = pecas(cur);
    // o ponto corrente abaixo dos proprios pisos vai a eles antes do passo (o racional e a
    // medida estao em caminha_mt, src/multitrait2.cpp)
    if (P.ok) {
      std::vector<double> z0 = P.zc;
      bool abaixo = false;
      for (std::size_t i = 0; i < z0.size(); i++)
        if (z0[i] < P.piso[i]) { z0[i] = P.piso[i]; abaixo = true; }
      if (abaixo) {
        const std::vector<double> th0 = theta_de_z(mz, ntz, z0);
        AvaliacaoAR a0 = avalia_ar1(d, th0, &cs);
        if (a0.ok) {
          theta = th0;
          cur = std::move(a0);
          P = pecas(cur);
        }
      }
    }
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
      std::vector<double> passo(ntz, 0.0);
      for (std::size_t i = 0; i < ntz; i++)
        for (std::size_t j = 0; j < ntz; j++)
          if (!P.ok || !P.congelado[j]) passo[i] += minv.at(i, j) * (P.ok ? P.sz[j] : cur.score[j]);
      // busca de COMPRIMENTO na mesma direcao antes de mexer no amortecimento: subir
      // lambda encurta e GIRA o passo, e o laco alternava entre dois estados sem andar
      for (const double alpha : {1.0, 0.5, 0.25}) {
        std::vector<double> cand;
        if (P.ok) {
          std::vector<double> zn = P.zc;
          for (std::size_t i = 0; i < ntz; i++)
            if (!P.congelado[i]) zn[i] -= alpha * passo[i];
          for (std::size_t i = 0; i < ntz; i++)
            if (zn[i] < P.piso[i]) zn[i] = P.piso[i];
          cand = theta_de_z(mz, ntz, zn);
        } else {
          cand = theta;
          for (std::size_t i = 0; i < ntz; i++) cand[i] -= alpha * passo[i];
        }
        // rho fora de (-1, 1) ou R0 nao positiva-definida: passo rejeitado. Em z o rho
        // anda em atanh e nao ha como sair do intervalo, mas o caminho de resgate cru
        // ainda pode propor um theta inadmissivel.
        if (sem_passo) continue;
        AvaliacaoAR prox = avalia_ar1(d, cand, &cs);
        if (!prox.ok || prox.neg2logl > cur.neg2logl + 1e-9) continue;
        double num = 0.0, den = 0.0;
        for (std::size_t i = 0; i < cand.size(); i++) {
          const double dlt = cand[i] - theta[i];
          num += dlt * dlt;
          den += cand[i] * cand[i];
        }
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
    if (!aceitou) {
      // a causa de verdade, em ingles e com a iteracao (o racional esta em ajusta_mt): antes
      // a string interna em portugues vazava e o bloco do maxiter dizia "stopped at 300
      // iteration(s)" num ajuste parado na 9a
      if (resgate_em(cur)) {
        lambda = 1e-2;
        parado = 0;
        if (verboso)
          Rprintf("      EM rescue: -2logL %.6f (no damped step improved)\n", cur.neg2logl + jac);
        continue;
      }
      // a estagnacao nao certifica (o racional e a medida estao em ajusta_mt)
      const PassoZ Pq = pecas(cur);
      const Kkt kq = sonda_kkt(Pq, true);
      if (kq == Kkt::moveu) {
        lambda = 1e-2;
        parado = 0;
        continue;
      }
      saiu_do_laco = true;
      na_fronteira = conta_parede(Pq);
      R.decremento = decremento_z(Pq);
      if (kq == Kkt::limite) {
        R.mensagem = msg_limite("at iteration " + std::to_string(it) + " no damped step "
                                "improved the likelihood, and");
      } else {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.3g", R.decremento);
        R.mensagem = std::string("did NOT converge: at iteration ") + std::to_string(it) +
            " no damped AI step (30 damping levels, three step lengths each) and no EM "
            "step improved the likelihood, and the Newton decrement of the off-boundary "
            "components is " + buf + (R.decremento < tol_dec
              ? std::string(" (inside the 2e-4 tolerance, but a point no step can "
                            "leave is not certified as an optimum). Try a different "
                            "start=, and profile any component the message reports at a "
                            "boundary")
              : std::string(" (tolerance 2e-4 on the -2logL scale), so the point is a "
                            "stall and not a certified optimum. Try a different start=, "
                            "and profile any component the message reports at a "
                            "boundary"));
      }
      break;
    }
    if (verboso) {
      Rprintf("iter %3d  -2logL %.6f  relDelta %.3e\n",
              (int) it, cur.neg2logl + jac, R.reldelta);
      imprime_theta(bruto(theta), nomes_theta_ar1(d), d.modelo);
    }
    if (R.reldelta < tol || parado >= 2) {
      // O CERTIFICADO, nas mesmas coordenadas e sobre os mesmos pisos que o passo. Ele
      // ficou em theta quando o passo foi para log-Cholesky, e a assimetria era um defeito:
      // o teste que havia aqui, lmin < 1e-3 lmax, e a exata desigualdade que o grampo do
      // passo torna FALSA por construcao, de modo que a direcao grampeada ficava congelada
      // no passo e cobrada integralmente no certificado. converged exige as duas coisas,
      // como no univariado: passo pequeno E decremento de Newton dos LIVRES na tolerancia.
      // o EM primeiro: se ele ainda anda, nao havia convergencia nenhuma, e o
      // amortecimento reinicia para o passo AI voltar a tentar do ponto novo
      if (resgate_em(cur)) {
        lambda = 1e-2;
        parado = 0;
        if (verboso)
          Rprintf("      EM rescue: -2logL %.6f (the AI step had stalled)\n", cur.neg2logl + jac);
        continue;
      }
      PassoZ Pc = pecas(cur);
      na_fronteira = conta_parede(Pc);
      R.decremento = decremento_z(Pc);
      if (R.decremento < tol_dec) {
        // o certificado so vale se nenhum bloco na parede melhora entrando no cone
        const Kkt kc = sonda_kkt(Pc, true);
        if (kc == Kkt::moveu) {
          lambda = 1e-2;
          parado = 0;
          R.decremento = std::numeric_limits<double>::quiet_NaN();
          continue;
        }
        saiu_do_laco = true;
        if (kc == Kkt::limite) {
          R.mensagem = msg_limite("the step had converged, but");
          break;
        }
        lapida(Pc);
        na_fronteira = conta_parede(Pc);
        R.convergiu = true;
        break;
      }
      // O corte de estagnacao ficou mais paciente depois que o resgate EM entrou, porque
      // "nada anda" passou a ser uma afirmacao mais forte: o passo AI travou E o EM
      // tambem nao move. Um EM bem sucedido reinicia o amortecimento, e o passo AI precisa
      // de espaco para retomar do ponto novo antes de o laco declarar estagnacao.
      const Kkt kp = parado >= 25 ? sonda_kkt(Pc, true) : Kkt::otimo;
      if (kp == Kkt::moveu) {
        lambda = 1e-2;
        parado = 0;
        continue;
      }
      if (parado >= 25) {
        saiu_do_laco = true;
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
  // a crista plana no fim do orcamento tambem passa pela sonda de KKT, sem mover (o
  // racional esta em ajusta_mt)
  const bool crista = !R.convergiu && !saiu_do_laco && maxiter > 0 && R.iters >= maxiter &&
                      R.decremento < tol_dec;
  const Kkt k_crista = crista ? sonda_kkt(pecas(cur), false) : Kkt::otimo;
  // de volta as unidades originais: y restaurado bit a bit e a avaliacao final no theta
  // reescalado (o racional esta em ajusta_mt)
  guarda.restaura();
  theta = bruto(theta);
  bool final_ok = true;
  {
    AvaliacaoAR fin = avalia_ar1(d, theta, &cs);
    if (fin.ok) cur = std::move(fin);
    else final_ok = false;
  }
  // EXCLUSAO POR LISTA, dita em voz alta. Um registro com qualquer caracteristica ausente
  // sai inteiro desta rota, e junto pode sair um ponto do meio da serie de um sujeito. O
  // multicaracter comum nao faz isso: la o registro fica e e ajustado contra a submatriz de
  // R0 do padrao dele. A diferenca tem de aparecer, senao o unico sinal e um n_used menor
  // que nrow(data).
  if (d.n_sem_tempo > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(d.n_sem_tempo) + " record(s) left out for a missing or non-finite time";
  if (d.n_incompletos > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(d.n_incompletos) + " record(s) with no trait observed left out";
  if (d.n_mv > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(d.n_mv) + " missing cell(s) kept in the series: the likelihood is "
        "the marginal one of the observed cells, each missing cell carried as its own "
        "fixed effect (the mv device of ASReml), which keeps the residual Gamma (x) R0 whole";
  if (na_fronteira > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(na_fronteira) + " component(s) resting at a covariance boundary "
        "excluded from the convergence certificate: a component pinned at a boundary "
        "points out of the cone by construction, so its gradient never vanishes and the "
        "certificate is CONDITIONAL on the pinning; read those components with care";
  R.se.assign(d.modelo.ntheta, std::nan(""));
  if (final_ok)
    try {
      Densa inv = inv_geral(cur.ai);
      for (std::size_t k = 0; k < d.modelo.ntheta; k++) {
        const double v = 2.0 * inv.at(k, k);
        if (v > 0.0) R.se[k] = std::sqrt(v);
      }
    } catch (const Erro&) {}

  // o mesmo veredicto do model() (aireml.cpp): esgotado o maxiter com o decremento de
  // Newton ja abaixo da tolerancia, o que se tem e uma crista plana, e nao um ajuste
  // que pede mais iteracoes. Aqui ele ainda dizia "Raise maxiter". So quando o laco andou
  // ate o fim (com maxiter = 0 uma avaliacao no start= nao e convergencia, o model() ja dizia
  // FALSE) e a sonda de KKT nao achou ponto melhor dentro do cone.
  if (crista && k_crista == Kkt::limite) {
    char buf[320];
    std::snprintf(buf, sizeof(buf), " (relDelta %.3g, Newton decrement %.3g): a component "
                  "resting at a covariance boundary is not at a constrained optimum, a point "
                  "inside the cone or along the boundary face lowers the -2logL by %.3g", R.reldelta, R.decremento,
                  ganho_kkt);
    R.mensagem = "stopped at " + std::to_string(maxiter) + " iteration(s) without a "
        "certified optimum" + buf + ". Raise maxiter=" +
        (R.mensagem.empty() ? std::string() : "; " + R.mensagem);
  }
  if (crista && k_crista == Kkt::otimo) {
    R.convergiu = true;
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.3g", R.reldelta);
    R.mensagem = std::string(R.mensagem.empty() ? "" : R.mensagem + "; ") +
        "converged by the Newton decrement while the components were still moving "
        "(relDelta " + buf + "): the likelihood is flat along some direction, so the "
        "values along it are about equally supported. Read the standard errors and the "
        "correlation between the components before trusting any one of them";
  }

  // o diagnostico de nao convergencia sai sempre que nenhum outro motivo de falha foi
  // escrito, e vai NA FRENTE das notas informativas (exclusao por lista, fronteira, tempo
  // ausente). Antes so saia com a mensagem vazia, e qualquer nota o escondia.
  if (!R.convergiu && maxiter == 0) {
    // uma avaliacao, nao um ajuste: dito como tal
    char buf[160];
    std::snprintf(buf, sizeof(buf), "maxiter = 0: no iteration was run, so the components "
                  "are the start evaluated, not estimated; the Newton decrement there is "
                  "%.3g", R.decremento);
    R.mensagem = std::string(buf) + (R.mensagem.empty() ? std::string() : "; " + R.mensagem);
  } else if (!R.convergiu && R.iters >= maxiter &&
             R.mensagem.find("did NOT converge") == std::string::npos &&
             R.mensagem.find("stopped at ") != 0) {
    // o teto de verdade: so aqui o ajuste parou POR causa do maxiter
    char buf[128];
    std::snprintf(buf, sizeof(buf), " (relDelta %.3g, Newton decrement %.3g against the "
                  "2e-4 tolerance)", R.reldelta, R.decremento);
    R.mensagem = "stopped at " + std::to_string(maxiter) + " iteration(s) without a "
        "certified optimum" + buf + ": this model asks for more iterations. Raise maxiter=" +
        (R.mensagem.empty() ? std::string() : "; " + R.mensagem);
  }
  if (cur.fora_do_padrao > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(cur.fora_do_padrao) + " value(s) read outside the sparsity pattern of the factor: the result is NOT reliable";

  if (verboso)
    Rprintf("%s at iter %d, -2logL %.6f\n", R.convergiu ? "converged" : "STOPPED",
            (int) R.iters, final_ok ? cur.neg2logl : cur.neg2logl + jac);
  R.theta = theta;
  R.fora_do_padrao = cur.fora_do_padrao;
  R.na_parede = parede_final;
  R.parou = saiu_do_laco && !R.convergiu;
  if (final_ok) {
    R.neg2logl = cur.neg2logl;
    R.solucao = cur.solucao;
    R.pev = cur.pev;
    const std::string av = aviso_informacao(cur.ai, nomes_theta_ar1(d), &R.se, &fator_corte);
    if (!av.empty()) R.mensagem += (R.mensagem.empty() ? "" : "; ") + av;
  } else {
    R.neg2logl = cur.neg2logl + jac;
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        "the final evaluation in the original units failed, so the solutions are not "
        "reported";
  }
  return R;
}

AjusteMT ajusta_ar1(DesenhoAR& d, const std::vector<double>* theta0, std::size_t maxiter,
                    double tol, bool verboso) {
  AjusteMT R = caminha_ar1(d, theta0, maxiter, tol, verboso, 0);
  MapaZ mz;
  for (const Grupo& gr : d.modelo.grupos) mz.blocos.push_back(std::make_pair(gr.offset, gr.dim));
  mz.blocos.push_back(std::make_pair(d.offset_s2e, d.t));
  mz.correlacoes.push_back(d.offset_rho);
  // a partida de var(y) em unidades ORIGINAIS: y ja voltou da padronizacao
  return explora_fronteira(std::move(R), mz, partida_ar1(d), maxiter,
                           [&](const std::vector<double>& th0, std::size_t it0) {
                             return caminha_ar1(d, &th0, maxiter, tol, verboso, it0);
                           }, verboso);
}

}  // namespace br
