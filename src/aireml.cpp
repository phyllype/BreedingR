// AI-REML: score analitico, informacao media, EM de aquecimento e passo amortecido.
//
// As regras que nao se quebram, todas vindas de defeito medido no projeto de origem:
//
//  1. convergencia RELATIVA nas componentes, sqrt(soma delta^2 / soma theta^2) < tol. Nunca
//     um limiar absoluto na norma do score, que e O(n_registros): a mesma tolerancia que
//     converge com 500 registros declara nao-convergido com 200 mil, com as estimativas
//     paradas no otimo.
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

  // EM de aquecimento: monotono, fica no cone, e de graca
  for (std::size_t e = 0; e < n_em; e++) {
    Avaliacao prox = avalia(d, cur.em_theta, &cs);
    if (!prox.ok) break;
    theta = cur.em_theta;
    cur = std::move(prox);
  }

  if (verboso && n_em > 0)
    Rprintf("EM warm-up (%d round(s)): -2logL %.6f\n", (int) n_em, cur.neg2logl);
  // CONJUNTO ATIVO. Uma variancia tem piso em zero, e o otimo REML pode estar EM cima
  // dele. Sem tratar isso, a componente encosta no piso, o passo passa a apontar para
  // fora do espaco, cada tentativa e rejeitada, o amortecimento explode e o ajuste
  // queima as iteracoes no lugar — terminando num -2logL PIOR que o do submodelo
  // aninhado, o que e impossivel num otimo. Congelando quem esta no piso com o gradiente
  // empurrando para fora, o resto e otimizado CONDICIONALMENTE a isso, que e a solucao
  // correta do problema com restricao.
  std::vector<char> congelado(theta.size(), 0);
  double lambda = 1e-3;
  for (std::size_t it = 1; it <= maxiter; it++) {
    R.iters = it;
    R_CheckUserInterrupt();
    const double piso = 1e-10 * std::max(1.0, theta[d.modelo.offset_residual]);
    for (const Grupo& gr : d.modelo.grupos)
      for (std::size_t i = 0; i < gr.dim; i++) {
        const std::size_t k = gr.theta_idx(i, i);       // so as VARIANCIAS tem piso
        congelado[k] = (theta[k] <= piso && cur.score[k] > 0.0) ? 1 : 0;
      }
    bool aceitou = false;
    for (int tent = 0; tent < 30; tent++) {
      Densa m = cur.ai;
      for (std::size_t i = 0; i < m.nlin; i++) {
        const double di = m.at(i, i);
        m.at(i, i) = (di == 0.0) ? lambda : di * (1.0 + lambda);
      }
      // congelado: linha e coluna viram identidade e o score sai — e o sistema REDUZIDO
      // dos que ainda podem andar, sem que o bloqueado contamine a direcao dos outros
      std::vector<double> sc = cur.score;
      for (std::size_t i = 0; i < m.nlin; i++)
        if (congelado[i]) {
          sc[i] = 0.0;
          for (std::size_t j = 0; j < m.ncol; j++) { m.at(i, j) = 0.0; m.at(j, i) = 0.0; }
          m.at(i, i) = 1.0;
        }
      Densa minv;
      try { minv = inv_geral(m); } catch (const Erro&) { lambda *= 10.0; continue; }
      std::vector<double> cand = theta;
      for (std::size_t i = 0; i < cand.size(); i++) {
        if (congelado[i]) continue;
        double passo = 0.0;
        for (std::size_t j = 0; j < cand.size(); j++) passo += minv.at(i, j) * sc[j];
        cand[i] -= passo;
      }
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
        cur = std::move(prox);
        lambda = std::max(lambda / 10.0, 1e-10);
        aceitou = true;
        break;
      }
      lambda *= 10.0;
    }
    if (!aceitou) { R.mensagem = "nenhum passo amortecido melhorou a verossimilhanca"; break; }
    if (verboso)
      Rprintf("iter %3d  -2logL %.6f  relDelta %.3e\n",
              (int) it, cur.neg2logl, R.reldelta);
    if (R.reldelta < tol) {
      // O PASSO PEQUENO NAO PROVA OTIMO. Com uma componente encostada no zero a AI fica
      // quase singular naquela direcao, o amortecimento cresce, o passo encolhe — e um
      // criterio que olha so o tamanho do passo declara convergencia com o score longe
      // de zero. O sintoma observado foi um modelo MAIOR parando com -2logL pior que o
      // submodelo aninhado dele, o que e impossivel no otimo.
      //
      // Antes de declarar, tenta um passo EM: ele e multiplicativo, fica dentro do cone
      // e NAO encolhe na fronteira, entao consegue andar exatamente onde o passo AI
      // travou. Se a verossimilhanca ainda melhora ali, nao havia convergencia nenhuma.
      Avaliacao pem = avalia(d, cur.em_theta, &cs);
      if (pem.ok && pem.neg2logl < cur.neg2logl - 1e-8) {
        double num = 0.0, den = 0.0;
        for (std::size_t i = 0; i < theta.size(); i++) {
          const double dlt = cur.em_theta[i] - theta[i];
          num += dlt * dlt;
          den += cur.em_theta[i] * cur.em_theta[i];
        }
        R.reldelta = std::sqrt(num / std::max(den, 1e-300));
        theta = cur.em_theta;
        cur = std::move(pem);
        if (verboso)
          Rprintf("      EM rescue: -2logL %.6f (the AI step had stalled)\n", cur.neg2logl);
        continue;
      }
      R.convergiu = true;
      break;
    }
  }
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

  {
    std::size_t nfix = 0;
    for (char c : congelado) nfix += c;
    if (nfix > 0)
      R.mensagem += std::string(R.mensagem.empty() ? "" : "; ") + std::to_string(nfix) +
          " component(s) at the zero boundary, fixed there while the others were "
          "optimized conditional on that";
  }
  if (!R.convergiu && R.mensagem.empty())
    R.mensagem = "parou em " + std::to_string(maxiter) + " iteracoes sem atingir a tolerancia relativa";
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
