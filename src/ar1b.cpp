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
                            const std::string& col_sujeito, const std::string& col_tempo) {
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
      // Cair no ramo diagonal trocaria a K declarada pela identidade EM SILENCIO — o
      // ajuste convergiria para outra coisa. Erro declarado ate este desenho carregar K.
      throw Erro("kernel() is not available in this fitter yet: in this version only "
                 "model() and eval_internal() carry the declared K");
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
  std::vector<std::size_t> fica, sai;
  posto_completo(xfull, 1e-9, fica, sai);
  d.x = Densa(d.nlin, fica.size());
  for (std::size_t jj = 0; jj < fica.size(); jj++) {
    d.nomes_x.push_back(cols[fica[jj]].first);
    for (std::size_t i = 0; i < d.nlin; i++) d.x.at(i, jj) = xfull.at(i, fica[jj]);
  }
  for (std::size_t j : sai) d.saiu_x.push_back(cols[j].first);
  for (std::size_t k = 0; k < m.termos.size(); k++) {
    if (!m.termos[k].aleatorio()) continue;
    const bool com_ped = m.termos[k].estrutura == Estrutura::Parentesco;
    d.aleatorios.push_back(monta_termo(m, k, tab, com_ped ? &niveis_ped : nullptr));
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

AjusteMT ajusta_ar1(const DesenhoAR& d, std::size_t maxiter, double tol, bool verboso) {
  AjusteMT R;
  const std::size_t t = d.t;
  // partida: variancias de var(y) por caracteristica, R0 diagonal, rho = 0
  std::vector<double> theta(d.modelo.ntheta, 0.0);
  {
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
    for (const Grupo& g : d.modelo.grupos) {
      const std::size_t por_t = g.dim / t;
      for (std::size_t i = 0; i < g.dim; i++) {
        const std::size_t tau = i / por_t;
        theta[g.offset + (i * (2 * g.dim - i + 1)) / 2] = frac * var_t[tau];
      }
    }
    for (std::size_t tau = 0; tau < t; tau++) {
      const std::size_t k = d.offset_s2e + (tau * (2 * t - tau + 1)) / 2;
      theta[k] = 0.5 * var_t[tau];
    }
    theta[d.offset_rho] = 0.0;
  }

  CacheSimbolica cs;
  AvaliacaoAR cur = avalia_ar1(d, theta, &cs);
  if (!cur.ok) { R.mensagem = "o theta inicial e INADMISSIVEL"; return R; }

  // O certificado final, espelhado do ajustador univariado (o racional completo esta em
  // aireml.cpp): passo relativo pequeno nao prova otimo. converged exige tambem o
  // decremento de Newton, g' AI^-1 g, ~2x o gap em -2logL perto do otimo. Este laco anda
  // em theta cru, sem pisos e sem log-Cholesky (o limite declarado da etapa 1B), entao o
  // certificado exclui as fronteiras que consegue reconhecer — variancia diagonal em
  // zero numerico (<= 1e-6 da escala residual) com o score empurrando para baixo, e todo
  // componente de um grupo cuja C_g esta a 1e-3 (razao das diagonais do Cholesky) da
  // singularidade — e a certificacao fica CONDICIONAL a fronteira, dito na mensagem.
  // rho fica sempre no certificado. AI que nao inverte no bloco livre conta como recusa.
  const double tol_dec = 2e-4;
  std::size_t na_fronteira = 0;   // exclusoes da ULTIMA chamada, para a mensagem
  auto decremento = [&](const AvaliacaoAR& av) -> double {
    const std::size_t nt = d.modelo.ntheta;
    std::vector<char> fora(nt, 0);
    double s2m = 0.0;
    for (std::size_t tau = 0; tau < t; tau++)
      s2m += theta[d.offset_s2e + (tau * (2 * t - tau + 1)) / 2];
    s2m = std::max(1.0, s2m / static_cast<double>(t));
    for (const Grupo& gr : d.modelo.grupos) {
      Densa cg(gr.dim, gr.dim);
      for (std::size_t j = 0; j < gr.dim; j++)
        for (std::size_t i = j; i < gr.dim; i++) {
          const std::size_t k = gr.offset + (j * (2 * gr.dim - j + 1)) / 2 + (i - j);
          cg.at(i, j) = theta[k];
          cg.at(j, i) = theta[k];
        }
      Densa L;
      bool sing = !chol_pequena(cg, L);
      if (!sing && gr.dim > 1) {
        double lmax = 0.0, lmin = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < gr.dim; i++) {
          lmax = std::max(lmax, L.at(i, i));
          lmin = std::min(lmin, L.at(i, i));
        }
        sing = lmin < 1e-3 * lmax;
      }
      for (std::size_t j = 0; j < gr.dim; j++)
        for (std::size_t i = j; i < gr.dim; i++) {
          const std::size_t k = gr.offset + (j * (2 * gr.dim - j + 1)) / 2 + (i - j);
          if (sing) fora[k] = 1;
          else if (i == j && theta[k] <= 1e-6 * s2m && av.score[k] > 0.0) fora[k] = 1;
        }
    }
    na_fronteira = 0;
    for (char f2 : fora) na_fronteira += f2;
    if (na_fronteira == nt) return 0.0;
    Densa m2 = av.ai;
    std::vector<double> sc = av.score;
    for (std::size_t i = 0; i < nt; i++)
      if (fora[i]) {
        sc[i] = 0.0;
        for (std::size_t j = 0; j < nt; j++) { m2.at(i, j) = 0.0; m2.at(j, i) = 0.0; }
        m2.at(i, i) = 1.0;
      }
    Densa minv;
    try { minv = inv_geral(m2); } catch (const Erro&) {
      return std::numeric_limits<double>::infinity();
    }
    double dec = 0.0;
    for (std::size_t i = 0; i < nt; i++) {
      if (sc[i] == 0.0) continue;
      for (std::size_t j = 0; j < nt; j++) dec += sc[i] * minv.at(i, j) * sc[j];
    }
    return std::fabs(dec);
  };

  double lambda = 1e-2;
  for (std::size_t it = 1; it <= maxiter; it++) {
    R.iters = it;
    R_CheckUserInterrupt();
    bool aceitou = false;
    for (int tent = 0; tent < 30; tent++) {
      Densa m2 = cur.ai;
      for (std::size_t i = 0; i < m2.nlin; i++) {
        const double di = m2.at(i, i);
        m2.at(i, i) = (di == 0.0) ? lambda : di * (1.0 + lambda);
      }
      Densa minv;
      try { minv = inv_geral(m2); } catch (const Erro&) { lambda *= 10.0; continue; }
      std::vector<double> cand = theta;
      for (std::size_t i = 0; i < cand.size(); i++) {
        double passo = 0.0;
        for (std::size_t j = 0; j < cand.size(); j++) passo += minv.at(i, j) * cur.score[j];
        cand[i] -= passo;
      }
      // rho fora de (-1, 1) ou R0 nao positiva-definida: passo rejeitado, amortecimento
      // sobe — a regra de qualquer theta inadmissivel, nunca projetar em silencio
      AvaliacaoAR prox = avalia_ar1(d, cand, &cs);
      if (prox.ok && prox.neg2logl <= cur.neg2logl + 1e-9) {
        double num = 0.0, den = 0.0;
        for (std::size_t i = 0; i < cand.size(); i++) {
          const double dlt = cand[i] - theta[i];
          num += dlt * dlt;
          den += cand[i] * cand[i];
        }
        R.reldelta = std::sqrt(num / std::max(den, 1e-300));
        theta = cand;
        cur = std::move(prox);
        lambda = std::max(lambda / 10.0, 1e-10);
        aceitou = true;
        break;
      }
      lambda *= 10.0;
    }
    if (!aceitou) { R.mensagem = "nenhum passo amortecido melhorou a verossimilhanca"; break; }
    if (verboso) {
      Rprintf("iter %3d  -2logL %.6f  relDelta %.3e\n",
              (int) it, cur.neg2logl, R.reldelta);
      imprime_theta(theta, nomes_theta_ar1(d));
    }
    if (R.reldelta < tol) {
      // Ao contrario do univariado, converged aqui continua sendo o criterio de passo:
      // este laco anda em theta cru e TRAVA INTEIRO quando uma direcao encosta numa
      // fronteira (o defeito declarado da etapa 1B) — reprovar o certificado so faria o
      // laco travado iterar ate maxiter, empurrando a variancia da fronteira para zero
      // exato sem ganhar verossimilhanca (medido: a celula sem sinal genetico de
      // test-fixed-solutions.R termina com var(animal) = 0 exato e a MME de referencia
      // singular). O decremento e computado e REPORTADO, e um certificado reprovado
      // vira aviso na mensagem; o portao duro fica para o porte log-Cholesky.
      R.decremento = decremento(cur);
      R.convergiu = true;
      if (R.decremento >= tol_dec) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.3g", R.decremento);
        R.mensagem = std::string("the step criterion converged, but the Newton ") +
            "decrement of the off-boundary components is " + buf + " against the 2e-4 "
            "tolerance: this fitter can rest short of the optimum near a covariance "
            "boundary, so treat the estimates as approximate there (the univariate "
            "fitter carries the boundary-following walker and the hard certificate)";
      }
      break;
    }
  }
  if (std::isnan(R.decremento)) R.decremento = decremento(cur);
  if (na_fronteira > 0)
    R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") +
        std::to_string(na_fronteira) + " component(s) at a covariance boundary (a "
        "variance at zero, or a group within 1e-3 of singularity) excluded from the "
        "convergence certificate: this fitter walks in raw theta and cannot follow a "
        "singular boundary, so the certificate is conditional on it; read the "
        "components near that boundary with care";

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
