// AR(1)/CAR(1) (Wade e Quaas, 1993), parte 2: desenho e estimador, para t caracteristicas (t = 1 e o caminho
// classico, e os gates antigos sao a regressao disto).

#include "mme.h"

// so para o TRACE opcional e a interrupcao: nenhuma conta usa o R
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#include <cstdio>
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
  for (std::size_t tau = 0; tau < t; tau++) {
    std::vector<double> col = tab.numerico(alvos[tau]);
    for (std::size_t i = 0; i < d.nlin; i++) {
      const bool falta = !std::isfinite(col[i]) ||
          (m.tem_ausente && std::fabs(col[i] - m.codigo_ausente) < 1e-9);
      if (falta) d.usa[i] = 0;
      d.y.at(i, tau) = falta ? 0.0 : col[i];
    }
  }

  d.tempo = tab.numerico(col_tempo);
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
    } else if (g.estrutura == Estrutura::Declarada) {
      kinv_declarada(m, g, kernels, d.kinv, d.kinv_logdet);
    } else {
      d.kinv.push_back(Csc());
      d.kinv_logdet.push_back(0.0);
    }
  }

  // X e aleatorios como no caminho comum
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
      for (std::size_t j = 0; j < d.x.ncol; j++) {
        const double v = d.x.at(r, j);
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
    // como na multi: coeficiente rotulado "termo[c]@alvo", tau-major sobre os termos
    for (const Grupo& g : d.modelo.grupos) {
      std::vector<std::string> rot;
      for (std::size_t tau = 0; tau < t; tau++)
        for (std::size_t tm : g.termos) {
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

AjusteMT ajusta_ar1(const DesenhoAR& d, const std::vector<double>* theta0, std::size_t maxiter,
                    double tol, bool verboso) {
  AjusteMT R;
  const std::size_t t = d.t;
  // partida: start= do usuario, ou variancias de var(y) por caracteristica, R0 diagonal,
  // rho = 0
  std::vector<double> theta(d.modelo.ntheta, 0.0);
  if (theta0) {
    theta = *theta0;
  } else {
    std::vector<double> var_t(t, 1.0);
    for (std::size_t tau = 0; tau < t; tau++) {
      double soma = 0.0, soma2 = 0.0;
      std::size_t n = 0;
      for (std::size_t i = 0; i < d.nlin; i++)
        if (d.usa[i]) { soma += d.y.at(i, tau); soma2 += d.y.at(i, tau) * d.y.at(i, tau); n++; }
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
        const double gm = std::exp(-d.kinv_logdet[ig] / static_cast<double>(d.kinv[ig].ncol));
        if (std::isfinite(gm) && gm > 0.0) esc = gm;
      }
      ig++;
      const std::size_t por_t = g.dim / t;
      for (std::size_t i = 0; i < g.dim; i++) {
        const std::size_t tau = i / por_t;
        theta[g.offset + (i * (2 * g.dim - i + 1)) / 2] = frac * var_t[tau] / esc;
      }
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

  CacheSimbolica cs;
  AvaliacaoAR cur = avalia_ar1(d, theta, &cs);
  if (!cur.ok) { R.mensagem = "o theta inicial e INADMISSIVEL"; return R; }

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
  auto conta_parede = [&](const PassoZ& P) {
    std::size_t n = 0;
    for (char c : P.na_parede) n += c;
    return n;
  };

  double lambda = 1e-2;
  int parado = 0;         // iteracoes consecutivas aceitas SEM progresso real
  for (std::size_t it = 1; it <= maxiter; it++) {
    R.iters = it;
    R_CheckUserInterrupt();
    bool aceitou = false;
    const std::size_t ntz = d.modelo.ntheta;
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
    if (!aceitou) { R.mensagem = "nenhum passo amortecido melhorou a verossimilhanca"; break; }
    if (verboso) {
      Rprintf("iter %3d  -2logL %.6f  relDelta %.3e\n",
              (int) it, cur.neg2logl, R.reldelta);
      imprime_theta(theta, nomes_theta_ar1(d), d.modelo);
    }
    if (R.reldelta < tol || parado >= 2) {
      // O CERTIFICADO, nas mesmas coordenadas e sobre os mesmos pisos que o passo. Ele
      // ficou em theta quando o passo foi para log-Cholesky, e a assimetria era um defeito:
      // o teste que havia aqui, lmin < 1e-3 lmax, e a exata desigualdade que o grampo do
      // passo torna FALSA por construcao, de modo que a direcao grampeada ficava congelada
      // no passo e cobrada integralmente no certificado. converged exige as duas coisas,
      // como no univariado: passo pequeno E decremento de Newton dos LIVRES na tolerancia.
      const PassoZ Pc = pecas(cur);
      na_fronteira = conta_parede(Pc);
      R.decremento = decremento_z(Pc);
      if (R.decremento < tol_dec) {
        R.convergiu = true;
        break;
      }
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
