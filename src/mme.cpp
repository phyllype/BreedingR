// As equacoes de modelo misto, montadas esparsas, e a verossimilhanca restrita.
//
// A identidade que amarra tudo, e que os testes exigem dos dois lados:
//
//   -2logL = (n - p) log s2e + log|G*| + log|C_s| + y'Py
//          = (n - p) log s2e + log|M| + log|X'M^-1 X| + y'P*y / s2e ,   M = Z G* Z' + I
//
// O lado esquerdo vem das MME montadas aqui; o direito e a forma V, computada DENSA no
// modulo de referencia. Se os dois nao coincidem, um esta errado — e e assim que a montagem
// e validada sem nunca se conferir contra si mesma.
//
// W'W e acumulada um REGISTRO de cada vez. Cada registro contribui um pequeno clique entre
// as colunas que toca, entao o custo e a soma de nnz_linha^2/2 por registro, e nada
// quadratico no numero de colunas. A penalidade kron(C_g^-1, K^-1) entra depois, na mesma
// lista de triplos — e a montagem SOMAR duplicados e o que faz uma posicao tocada por
// registros e pela penalidade acumular em vez de um apagar o outro.

#include "mme.h"

namespace br {

// C_g do grupo, lida de theta.
Densa cov_grupo(const Modelo& m, const std::vector<double>& theta, std::size_t g) {
  const Grupo& gr = m.grupos[g];
  Densa c(gr.dim, gr.dim);
  for (std::size_t j = 0; j < gr.dim; j++)
    for (std::size_t i = j; i < gr.dim; i++) {
      const double v = theta[gr.theta_idx(i, j)];
      c.at(i, j) = v;
      c.at(j, i) = v;
    }
  return c;
}

// Linhas de uma esparsa, para a montagem andar por registro.
static std::vector<std::vector<std::pair<std::uint32_t, double>>> linhas_de(const Csc& m) {
  std::vector<std::vector<std::pair<std::uint32_t, double>>> out(m.nlin);
  for (std::size_t c = 0; c < m.ncol; c++)
    for (std::size_t k = m.colptr[c]; k < m.colptr[c + 1]; k++)
      out[m.linha[k]].push_back({static_cast<std::uint32_t>(c), m.valor[k]});
  return out;
}

Montado monta_mme(const Desenho& d, const std::vector<double>& theta) {
  Montado M;
  M.s2e = theta[d.modelo.offset_residual];
  if (!(M.s2e > 0.0)) return M;          // inadmissivel: nao e resultado

  M.n_fixo = d.x.ncol;
  M.offset_grupo.resize(d.modelo.grupos.size());
  std::size_t acc = M.n_fixo;
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    M.offset_grupo[g] = acc;
    acc += d.largura(g);
  }
  M.total = acc;

  // (indice do aleatorio, primeira coluna global) por grupo, na ordem dos slots
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> slots(d.modelo.grupos.size());
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    std::size_t col = M.offset_grupo[g];
    for (std::size_t t : d.modelo.grupos[g].termos)
      for (std::size_t a = 0; a < d.aleatorios.size(); a++)
        if (d.aleatorios[a].termo == t) {
          slots[g].push_back({a, col});
          col += d.aleatorios[a].z.ncol;
        }
  }

  std::vector<std::vector<std::vector<std::pair<std::uint32_t, double>>>> zl;
  for (const DesenhoTermo& a : d.aleatorios) zl.push_back(linhas_de(a.z));

  std::vector<std::uint32_t> ti, tj;
  std::vector<double> tv;
  M.rhs.assign(M.total, 0.0);
  std::vector<std::pair<std::uint32_t, double>> lin;
  lin.reserve(64);

  for (std::size_t r = 0; r < d.nlin; r++) {
    if (!d.usa[r]) continue;
    lin.clear();
    for (std::size_t j = 0; j < d.x.ncol; j++) {
      const double v = d.x.at(r, j);
      if (v != 0.0) lin.push_back({static_cast<std::uint32_t>(j), v});
    }
    for (std::size_t g = 0; g < slots.size(); g++)
      for (const auto& [a, col0] : slots[g])
        for (const auto& [c, v] : zl[a][r])
          lin.push_back({static_cast<std::uint32_t>(col0 + c), v});

    const double yv = d.y[r];
    for (std::size_t p = 0; p < lin.size(); p++) {
      M.rhs[lin[p].first] += lin[p].second * yv;
      for (std::size_t q = p; q < lin.size(); q++) {
        std::uint32_t i = lin[p].first, j = lin[q].first;
        if (i > j) std::swap(i, j);
        ti.push_back(i);
        tj.push_back(j);
        tv.push_back(lin[p].second * lin[q].second);
      }
    }
  }

  // penalidade por grupo, so o triangulo superior
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    Densa cg = cov_grupo(d.modelo, theta, g);
    // escala por s2e: C_s = W'W + kron(C_g^-1 s2e, K^-1) e o sistema em unidades de s2e
    Densa cgs = cg;
    for (double& v : cgs.dados) v /= M.s2e;
    Densa cinv;
    try { cinv = inv_pd(cgs); } catch (const Erro&) { return M; }
    const double ld = logdet_pd(cgs);
    if (std::isnan(ld)) return M;

    const std::size_t nl = d.aleatorios.empty() ? 0 :
        [&]{ for (const auto& a : d.aleatorios)
               if (a.termo == d.modelo.grupos[g].termos[0]) return a.n_niveis;
             return static_cast<std::size_t>(0); }();
    M.logdet_g += static_cast<double>(nl) * ld - static_cast<double>(d.modelo.grupos[g].dim) * d.kinv_logdet[g];

    const std::size_t off = M.offset_grupo[g];
    const bool com_k = d.kinv[g].ncol > 0;
    for (std::size_t a = 0; a < d.modelo.grupos[g].dim; a++)
      for (std::size_t b = 0; b < d.modelo.grupos[g].dim; b++) {
        const double f = cinv.at(a, b);
        // f == 0 NAO pula mais: a covariancia que comeca em zero ganharia um padrao
        // menor na primeira avaliacao, e o cache da simbolica congelaria esse padrao
        // errado para o ajuste inteiro. Zero explicito ocupa o slot e nao muda numero.
        // EMITE A MATRIZ CHEIA E MANTEM SO gi <= gj, SEM TROCAR.
        //
        // A primeira versao trocava (swap) e espelhava: com o laco percorrendo a e b
        // completos, cada posicao fora da diagonal era emitida DUAS vezes, e como a
        // montagem soma duplicados, a penalidade saia dobrada. Uma penalidade dobrada
        // ainda e simetrica e definida — converge, para o lugar errado. Mantendo so o
        // triangulo sem trocar, a propria varredura completa de (a,b) garante que cada
        // posicao superior da kron cheia e emitida exatamente uma vez.
        auto poe = [&](std::size_t gi, std::size_t gj, double val) {
          if (gi > gj) return;
          ti.push_back(static_cast<std::uint32_t>(gi));
          tj.push_back(static_cast<std::uint32_t>(gj));
          tv.push_back(val);
        };
        if (com_k) {
          const Csc& k = d.kinv[g];
          for (std::size_t col = 0; col < k.ncol; col++)
            for (std::size_t p = k.colptr[col]; p < k.colptr[col + 1]; p++) {
              const std::size_t rk = k.linha[p];
              // K^-1 vem no triangulo inferior: a matriz cheia tem (rk,col) e (col,rk)
              poe(off + a * nl + rk, off + b * nl + col, f * k.valor[p]);
              if (rk != col) poe(off + a * nl + col, off + b * nl + rk, f * k.valor[p]);
            }
        } else {
          for (std::size_t l = 0; l < nl; l++)
            poe(off + a * nl + l, off + b * nl + l, f);
        }
      }
  }

  M.c = de_triplos(M.total, M.total, ti, tj, tv);
  M.ok = true;
  return M;
}

// -2logL pela via esparsa: monta, fatora, resolve.
Neg2LogL neg2logl_esparsa(const Desenho& d, const std::vector<double>& theta) {
  Neg2LogL r;
  Montado M = monta_mme(d, theta);
  if (!M.ok) return r;

  std::vector<std::size_t> perm = grau_minimo(M.c);
  Csc pc = permuta_sim(M.c, perm);
  Simbolica sb = simbolica(pc);
  Csc L;
  if (!cholesky(pc, sb, L)) return r;

  std::vector<double> pb(M.total);
  for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[perm[k]];
  std::vector<double> px = resolve(L, pb);
  r.solucao.assign(M.total, 0.0);
  for (std::size_t k = 0; k < M.total; k++) r.solucao[perm[k]] = px[k];

  double yy = 0.0;
  for (std::size_t i = 0; i < d.nlin; i++)
    if (d.usa[i]) yy += d.y[i] * d.y[i];
  double bry = 0.0;
  for (std::size_t k = 0; k < M.total; k++) bry += r.solucao[k] * M.rhs[k];

  const double n = static_cast<double>(d.n_usadas());
  const double p = static_cast<double>(M.n_fixo);
  r.logdet_c = logdet(L);
  // O expoente de s2e e n - p e nao n: o REML integra os efeitos fixos. E log|G| carrega
  // log|K|, nao log|K^-1| — a diferenca e constante em theta, entao o OTIMO nao se move e so
  // o VALOR sai errado, o que arruina qualquer comparacao com outro programa sem arruinar as
  // estimativas. Invisivel num teste de recuperacao, fatal num de -2logL.
  r.valor = (n - p) * std::log(M.s2e) + M.logdet_g + r.logdet_c + (yy - bry) / M.s2e;
  r.perm = std::move(perm);
  r.L = std::move(L);
  r.ok = true;
  return r;
}

// -2logL pela FORMA V, densa. E a referencia: um caminho que nao compartilha nada com a
// montagem esparsa alem do desenho. So aguenta problemas pequenos, e e para isso que existe.
double neg2logl_densa_V(const Desenho& d, const std::vector<double>& theta) {
  const double s2e = theta[d.modelo.offset_residual];
  if (!(s2e > 0.0)) return std::nan("");

  // linhas usadas
  std::vector<std::size_t> linhas;
  for (std::size_t i = 0; i < d.nlin; i++)
    if (d.usa[i]) linhas.push_back(i);
  const std::size_t n = linhas.size();

  // M = I + soma_g Z_g (C_g/s2e (x) K) Z_g'
  Densa M(n, n);
  for (std::size_t i = 0; i < n; i++) M.at(i, i) = 1.0;

  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    Densa cg = cov_grupo(d.modelo, theta, g);
    for (double& v : cg.dados) v /= s2e;

    // K densa: inversa da K^-1 esparsa, ou identidade
    const std::size_t nl = [&]{ for (const auto& a : d.aleatorios)
        if (a.termo == d.modelo.grupos[g].termos[0]) return a.n_niveis;
      return static_cast<std::size_t>(0); }();
    Densa K;
    if (d.kinv[g].ncol > 0) K = inv_pd(d.kinv[g].densa_simetrica());
    else { K = Densa(nl, nl); for (std::size_t i = 0; i < nl; i++) K.at(i, i) = 1.0; }

    // Z do grupo, denso nas linhas usadas
    const std::size_t larg = d.largura(g);
    Densa Z(n, larg);
    std::size_t col0 = 0;
    for (std::size_t t : d.modelo.grupos[g].termos)
      for (const DesenhoTermo& a : d.aleatorios)
        if (a.termo == t) {
          Densa zd = a.z.densa();
          for (std::size_t r = 0; r < n; r++)
            for (std::size_t c = 0; c < a.z.ncol; c++)
              Z.at(r, col0 + c) = zd.at(linhas[r], c);
          col0 += a.z.ncol;
        }

    // V_g = kron(C_g, K) na convencao coluna = coef * nl + nivel
    const std::size_t dim = d.modelo.grupos[g].dim;
    Densa Vg(larg, larg);
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t b = 0; b < dim; b++)
        for (std::size_t i = 0; i < nl; i++)
          for (std::size_t j = 0; j < nl; j++)
            Vg.at(a * nl + i, b * nl + j) = cg.at(a, b) * K.at(i, j);

    // M += Z Vg Z'
    Densa ZV(n, larg);
    for (std::size_t r = 0; r < n; r++)
      for (std::size_t c = 0; c < larg; c++) {
        double s = 0.0;
        for (std::size_t k = 0; k < larg; k++) s += Z.at(r, k) * Vg.at(k, c);
        ZV.at(r, c) = s;
      }
    for (std::size_t r = 0; r < n; r++)
      for (std::size_t c = 0; c < n; c++) {
        double s = 0.0;
        for (std::size_t k = 0; k < larg; k++) s += ZV.at(r, k) * Z.at(c, k);
        M.at(r, c) += s;
      }
  }

  const double logdet_M = logdet_pd(M);
  if (std::isnan(logdet_M)) return std::nan("");
  Densa Minv = inv_pd(M);

  // X e y nas linhas usadas
  const std::size_t p = d.x.ncol;
  Densa X(n, p);
  std::vector<double> y(n);
  for (std::size_t r = 0; r < n; r++) {
    y[r] = d.y[linhas[r]];
    for (std::size_t c = 0; c < p; c++) X.at(r, c) = d.x.at(linhas[r], c);
  }

  // X' M^-1 X e X' M^-1 y
  Densa XtMi(p, n);
  for (std::size_t a = 0; a < p; a++)
    for (std::size_t r = 0; r < n; r++) {
      double s = 0.0;
      for (std::size_t k = 0; k < n; k++) s += X.at(k, a) * Minv.at(k, r);
      XtMi.at(a, r) = s;
    }
  Densa XtMiX(p, p);
  std::vector<double> XtMiy(p, 0.0);
  for (std::size_t a = 0; a < p; a++) {
    for (std::size_t b = 0; b < p; b++) {
      double s = 0.0;
      for (std::size_t r = 0; r < n; r++) s += XtMi.at(a, r) * X.at(r, b);
      XtMiX.at(a, b) = s;
    }
    for (std::size_t r = 0; r < n; r++) XtMiy[a] += XtMi.at(a, r) * y[r];
  }
  const double logdet_X = logdet_pd(XtMiX);
  if (std::isnan(logdet_X)) return std::nan("");
  Densa XtMiXinv = inv_pd(XtMiX);

  // y'P*y = y'M^-1 y - (X'M^-1 y)' (X'M^-1X)^-1 (X'M^-1 y)
  double yMiy = 0.0;
  for (std::size_t r = 0; r < n; r++) {
    double s = 0.0;
    for (std::size_t k = 0; k < n; k++) s += Minv.at(r, k) * y[k];
    yMiy += y[r] * s;
  }
  double quad = 0.0;
  for (std::size_t a = 0; a < p; a++)
    for (std::size_t b = 0; b < p; b++) quad += XtMiy[a] * XtMiXinv.at(a, b) * XtMiy[b];

  const double nn = static_cast<double>(n), pp = static_cast<double>(p);
  return (nn - pp) * std::log(s2e) + logdet_M + logdet_X + (yMiy - quad) / s2e;
}



// ------------------------------------------------------------------ montagem do desenho
//
// Junta modelo, tabela e pedigree num Desenho pronto para ajustar. Vive aqui e nao na
// travessia: e logica com decisao numerica (posto completo, exclusao de linha) e tem de
// estar debaixo dos mesmos gates que o resto.

Desenho monta_desenho(const Modelo& m, const Tabela& t, const Pedigree* ped) {
  Desenho d;
  d.modelo = m;
  d.nlin = t.nlin;
  d.y = t.numerico(m.alvo);

  const bool precisa_ped = [&]{
    for (const Termo& tm : m.termos)
      if (tm.estrutura == Estrutura::Parentesco) return true;
    return false;
  }();
  if (precisa_ped && !ped)
    throw Erro("there is a term with relationship and no pedigree was given");

  // K^-1 por grupo. O A^-1 e UM so, partilhado pelos grupos com parentesco, e o log|K^-1|
  // vem da fatoracao esparsa dele.
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
    if (!cholesky(pa, sb, L))
      throw Erro("the pedigree A^-1 is not positive-definite");
    ld_ainv = logdet(L);
  }
  for (const Grupo& g : m.grupos) {
    if (g.estrutura == Estrutura::Parentesco) {
      d.kinv.push_back(ainv);
      d.kinv_logdet.push_back(ld_ainv);
    } else {
      d.kinv.push_back(Csc());
      d.kinv_logdet.push_back(0.0);
    }
  }

  // X: intercepto + termos fixos expandidos, depois posto completo
  std::vector<std::pair<std::string, std::vector<double>>> cols;
  cols.push_back({"intercept", std::vector<double>(d.nlin, 1.0)});
  for (std::size_t k = 0; k < m.termos.size(); k++) {
    if (m.termos[k].aleatorio()) continue;
    DesenhoTermo dt = monta_termo(m, k, t, nullptr);
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

  // aleatorios, com niveis do pedigree quando ha parentesco
  for (std::size_t k = 0; k < m.termos.size(); k++) {
    if (!m.termos[k].aleatorio()) continue;
    const bool com_ped = m.termos[k].estrutura == Estrutura::Parentesco;
    d.aleatorios.push_back(monta_termo(m, k, t, com_ped ? &niveis_ped : nullptr));
  }

  // linhas que entram: nem ausente, nem nivel sem casar no parentesco.
  //
  // NA/NaN na observacao E ausente, com ou sem codigo declarado: e o idioma do R, os
  // caminhos multi e AR ja tratavam assim, e o univariado nao — um NA atravessava a
  // marcacao e virava NaN na verossimilhanca inteira, sem erro nenhum, so um -2logL
  // NaN. O valor tambem TEM de ser zerado: a linha sai de `usa`, mas NaN * 0 continua
  // NaN nas somas que varrem o vetor inteiro.
  d.usa.assign(d.nlin, 1);
  for (std::size_t i = 0; i < d.nlin; i++) {
    const bool falta = !std::isfinite(d.y[i]) ||
        (m.tem_ausente && std::fabs(d.y[i] - m.codigo_ausente) < 1e-9);
    if (falta) { d.usa[i] = 0; d.y[i] = 0.0; }
  }
  for (const DesenhoTermo& a : d.aleatorios)
    for (std::size_t i = 0; i < d.nlin; i++)
      if (!a.casou[i]) d.usa[i] = 0;
  if (d.n_usadas() == 0) throw Erro("no row enters the analysis");
  return d;
}

}  // namespace br
