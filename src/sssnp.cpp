// ssSNPBLUP: o passo unico SEM G, os marcadores como equacoes (Liu et al. 2014).
//
// O modelo equivalente, escrito por extenso porque cada bloco abaixo vem de um termo dele:
//
//   u_g = Z g + p_g,   g ~ N(0, I (1-w) s2u / kd),   p_g ~ N(0, w s2u A22),
//   u_n | u_g ~ o condicional do pedigree,           kd = 2 sum p(1-p).
//
// A precisao conjunta de (u_n, u_g, g) sai da fatoracao p(u_n|u_g) p(u_g|g) p(g), e o
// termo A^gn (A^nn)^-1 A^ng = A^gg - A22^-1 (Schur) transforma tudo em blocos de A^-1:
//
//   Q s2u = [ A^nn   A^ng                          0
//             A^gn   A^gg + ((1-w)/w) A22^-1      -(1/w) A22^-1 Z
//             0     -(1/w) Z' A22^-1               (1/w) Z' A22^-1 Z + (kd/(1-w)) I ]
//
// A parte A^nn/A^ng/A^gg ja E a penalidade que monta_mme escreve com o A^-1 esparso; o
// resto sao os acrescimos deste arquivo. G nunca e formada nem invertida: o sistema e
// resolvido por gradientes conjugados precondicionados (Vandenplas et al. 2018, 2019), e
// cada aplicacao de A22^-1 usa a identidade de Masuda et al. (2017),
//
//   A22^-1 v = A^22 v - A^21 (A^11)^-1 A^12 v,
//
// com UMA fatoracao esparsa do bloco nao-genotipado A^11 feita no comeco e uma resolucao
// triangular por aplicacao. w -> 1 desliga os marcadores e colapsa no BLUP de pedigree,
// e esse colapso e um dos gates.
//
// Mais de um termo e grupos nao escalares. Os marcadores entram em TODOS os grupos com
// parentesco e em cada componente deles, como o H^-1 do caminho genotypes= entra em todos
// (aplica_genomica_em). Num grupo de dimensao q com covariancia K0, u = (u_1..u_q) ~
// N(0, K0 x A), e cada componente tem os seus marcadores: u_g,t = Z g_t + p_t, com
// (g_1..g_q) ~ N(0, K0 x I (1-w)/kd) e (p_1..p_q) ~ N(0, K0 x w A22). A covariancia
// conjunta de (u_n, u_g, g) e K0 x S1, S1 a do caso escalar com s2u = 1, e a precisao e
// K0^-1 x Q1: os acrescimos acima com s2e/s2u trocado pela entrada (a, b) de s2e K0^-1.
// Grupos diferentes nao se tocam. Isso cobre o direto-materno, a norma de reacao, o
// indireto e direto e materno em grupos separados.
//
// Limites declarados: theta e DADO (isto e um resolvedor, como o BLUP com componentes
// fixas da pratica; a REML continua nos caminhos exatos), e G* implicita =
// (1-w) Z Z'/kd + w A22, SEM o ajuste afim do caminho genotypes=, em populacoes fora do
// equilibrio os dois caminhos diferem por construcao.

#include "mme.h"

#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#ifndef FCONE
# define FCONE
#endif

namespace br {

namespace {

// Aplicacao simetrica de uma CSC triangular (superior OU inferior): cada entrada guardada
// (i, j, v) contribui v x_j na linha i e, fora da diagonal, v x_i na linha j.
void symv_tri(const Csc& a, const double* x, double* y) {
  for (std::size_t j = 0; j < a.ncol; j++)
    for (std::size_t k = a.colptr[j]; k < a.colptr[j + 1]; k++) {
      const std::size_t i = a.linha[k];
      const double v = a.valor[k];
      y[i] += v * x[j];
      if (i != j) y[j] += v * x[i];
    }
}

}  // namespace

SnpBlup snp_blup(const Desenho& d, Densa& mg, const std::vector<std::string>& geno_ids,
                 const std::vector<double>& theta, double w, std::size_t maxiter,
                 double tol, bool verboso, bool meio) {
  SnpBlup r;
  if (!(w > 0.0 && w < 1.0))
    throw Erro("rpg (the residual polygenic proportion) must be in (0, 1): 0 leaves the "
               "polygenic residual without a distribution, 1 turns the markers off");

  // As fatias: cada componente de cada grupo com parentesco, na ordem dos slots. Todas
  // indexam os mesmos niveis (o pedigree), e o A^-1 e um so.
  struct Fatia {
    std::size_t grupo, comp;
  };
  std::vector<Fatia> fatias;
  const DesenhoTermo* dt = nullptr;
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    const Grupo& gr = d.modelo.grupos[g];
    if (gr.estrutura != Estrutura::Parentesco) continue;
    std::size_t comp = 0;
    for (std::size_t t : gr.termos) {
      const Termo& tm = d.modelo.termos[t];
      for (std::size_t c = 0; c < tm.n_coef(); c++) {
        fatias.push_back({g, comp++});
        r.fatias.push_back(tm.n_coef() == 1 ? tm.nome : tm.nome + "[" + std::to_string(c) + "]");
      }
      for (const DesenhoTermo& a : d.aleatorios)
        if (a.termo == t) {
          if (!dt) dt = &a;
          else if (a.niveis != dt->niveis)
            throw Erro("the relationship terms do not index the same pedigree levels");
        }
    }
  }
  if (fatias.empty()) throw Erro("there is no relationship term to receive the markers");
  if (!dt) throw Erro("relationship group without an incidence");
  const std::size_t nl = dt->n_niveis;
  const std::size_t nf = fatias.size();
  const std::size_t gpar = fatias[0].grupo;

  // genotipados -> posicao de nivel (== posicao no A^-1)
  std::vector<std::size_t> pos;
  pos.reserve(geno_ids.size());
  {
    std::vector<std::pair<std::string, std::size_t>> tmp;
    tmp.reserve(nl);
    for (std::size_t i = 0; i < nl; i++) tmp.push_back({dt->niveis[i], i});
    std::sort(tmp.begin(), tmp.end());
    for (const std::string& g : geno_ids) {
      auto it = std::lower_bound(tmp.begin(), tmp.end(), std::make_pair(g, std::size_t(0)),
          [](const auto& a, const auto& b){ return a.first < b.first; });
      if (it == tmp.end() || it->first != g)
        throw Erro("genotyped animal '" + g + "' is not in the pedigree");
      pos.push_back(it->second);
    }
  }
  const std::size_t ng = pos.size();

  // Z centrada (imputacao pela media, monomorficos fora), coluna-major para o dgemv. Com
  // metafundadores (meio) e a Z da G05 (Legarra et al., 2015): centrada em 0.5, escala m/2,
  // TODOS os marcadores, a mesma G do caminho genotypes= com metafounders=, entao os dois
  // caminhos resolvem o mesmo sistema (a G* implicita daqui ja nao tem o ajuste afim).
  const std::size_t mm = mg.ncol;
  if (mg.nlin != ng) throw Erro("genotype rows and ids with different lengths");
  r.usa_marcador.assign(mm, 1);
  std::vector<double> pfreq(mm, 0.0), media_m(mm, 0.0);
  double kd = 0.0;
  std::vector<std::size_t> usados;
  for (std::size_t j = 0; j < mm; j++) {
    double soma = 0.0;
    std::size_t k = 0;
    for (std::size_t i = 0; i < ng; i++) {
      const double x = mg.at(i, j);
      if (std::isfinite(x)) { soma += x; k++; }
    }
    if (k == 0) { r.usa_marcador[j] = 0; r.n_monomorficos++; continue; }
    const double media = soma / static_cast<double>(k);
    media_m[j] = media;
    pfreq[j] = meio ? 0.5 : media / 2.0;
    if (!meio && (pfreq[j] <= 0.0 || pfreq[j] >= 1.0)) {
      r.usa_marcador[j] = 0; r.n_monomorficos++; continue;
    }
    kd += 2.0 * pfreq[j] * (1.0 - pfreq[j]);
    usados.push_back(j);
  }
  const std::size_t mu = usados.size();
  if (!(kd > 0.0)) throw Erro("2 sum p(1-p) is not positive: the markers do not vary");
  std::vector<double> zc(ng * mu);
  for (std::size_t jj = 0; jj < mu; jj++) {
    const std::size_t j = usados[jj];
    const double dp = 2.0 * pfreq[j];
    for (std::size_t i = 0; i < ng; i++) {
      double x = mg.at(i, j);
      if (!std::isfinite(x)) { x = media_m[j]; r.n_imputados++; }
      zc[i + ng * jj] = x - dp;
    }
  }

  // blocos de A^-1 pela mascara genotipado/nao: A11 esparsa fatorada uma vez, A12 e A22
  // como triplos para as aplicacoes
  const Csc& ainv = d.kinv[gpar];
  if (ainv.ncol != nl) throw Erro("A^-1 and the term with different sizes");
  std::vector<char> eh_geno(nl, 0);
  for (std::size_t i : pos) eh_geno[i] = 1;
  std::vector<std::size_t> loc(nl);
  std::size_t n1 = 0;
  for (std::size_t i = 0; i < nl; i++) if (!eh_geno[i]) loc[i] = n1++;
  for (std::size_t k = 0; k < ng; k++) loc[pos[k]] = k;

  std::vector<std::uint32_t> i11, j11;
  std::vector<double> v11;
  std::vector<std::uint32_t> i12, j12;      // linha = nao-genotipado local, coluna = genotipado local
  std::vector<double> v12;
  std::vector<std::uint32_t> i22, j22;
  std::vector<double> v22;
  std::vector<double> diag22(ng, 0.0);      // diagonal do BLOCO A^22, o substituto do precondicionador
  for (std::size_t c = 0; c < nl; c++)
    for (std::size_t k = ainv.colptr[c]; k < ainv.colptr[c + 1]; k++) {
      const std::size_t rr = ainv.linha[k];
      const double x = ainv.valor[k];
      const bool rg = eh_geno[rr], cg = eh_geno[c];
      if (!rg && !cg) {
        std::size_t i = loc[rr], j = loc[c];
        if (i < j) std::swap(i, j);
        i11.push_back(static_cast<std::uint32_t>(i));
        j11.push_back(static_cast<std::uint32_t>(j));
        v11.push_back(x);
      } else if (rg && cg) {
        i22.push_back(static_cast<std::uint32_t>(loc[rr]));
        j22.push_back(static_cast<std::uint32_t>(loc[c]));
        v22.push_back(x);
        if (rr == c) diag22[loc[rr]] += x;
      } else {
        i12.push_back(static_cast<std::uint32_t>(rg ? loc[c] : loc[rr]));
        j12.push_back(static_cast<std::uint32_t>(rg ? loc[rr] : loc[c]));
        v12.push_back(x);
      }
    }

  Csc a22blk = de_triplos(ng, ng, i22, j22, v22);
  Csc l11;
  std::vector<std::size_t> perm11;
  if (n1 > 0) {
    Csc a11 = de_triplos(n1, n1, i11, j11, v11);
    perm11 = grau_minimo(a11);
    Csc pa = permuta_sim(a11, perm11);
    Simbolica sb = simbolica(pa);
    if (!cholesky(pa, sb, l11))
      throw Erro("the non-genotyped block of A^-1 is not positive-definite");
  }

  // A22^-1 v = A^22 v - A^21 (A^11)^-1 A^12 v, sem nunca formar A22^-1. O rascunho vem
  // de fora: as aplicacoes de uma iteracao correm em paralelo, cada uma com o seu.
  auto a22inv_vezes = [&](const std::vector<double>& x, std::vector<double>& y,
                          std::vector<double>& t1, std::vector<double>& t2) {
    y.assign(ng, 0.0);
    symv_tri(a22blk, x.data(), y.data());
    if (n1 == 0) return;
    std::fill(t1.begin(), t1.end(), 0.0);
    for (std::size_t k = 0; k < v12.size(); k++) t1[i12[k]] += v12[k] * x[j12[k]];
    for (std::size_t i = 0; i < n1; i++) t2[i] = t1[perm11[i]];
    std::vector<double> px = resolve(l11, t2);
    for (std::size_t i = 0; i < n1; i++) t1[perm11[i]] = px[i];
    for (std::size_t k = 0; k < v12.size(); k++) y[j12[k]] -= v12[k] * t1[i12[k]];
  };

  // as equacoes de base, nas unidades de s2e como todo o resto do motor
  Montado M = monta_mme(d, theta);
  if (!M.ok) throw Erro("theta INADMISSIBLE: some covariance is not positive-definite");
  // s2e K0^-1 por grupo com parentesco, o mesmo fator que monta_mme poe na penalidade
  std::vector<Densa> fg(d.modelo.grupos.size());
  for (const Fatia& f : fatias) {
    if (fg[f.grupo].nlin > 0) continue;
    Densa cgs = cov_grupo(d.modelo, theta, f.grupo);
    for (double& v : cgs.dados) v /= M.s2e;
    fg[f.grupo] = inv_pd(cgs);
  }
  auto col_geno = [&](std::size_t s, std::size_t k) {
    return M.offset_grupo[fatias[s].grupo] + fatias[s].comp * nl + pos[k];
  };

  const std::size_t N = M.total + nf * mu;
  std::vector<double> rhs(N, 0.0);
  for (std::size_t k = 0; k < M.total; k++) rhs[k] = M.rhs[k];

  // precondicionador diagonal: diag(C) mais os acrescimos, com diag(A^22) como
  // substituto de diag(A22^-1), e um majorante (o Schur so subtrai), entao o
  // escalonamento fica do lado conservador e a correcao vem das iteracoes
  std::vector<double> prec(N, 0.0);
  for (std::size_t j = 0; j < M.c.ncol; j++)
    for (std::size_t k = M.c.colptr[j]; k < M.c.colptr[j + 1]; k++)
      if (M.c.linha[k] == j) prec[j] += M.c.valor[k];
  std::vector<double> szz(mu, 0.0);
  for (std::size_t jj = 0; jj < mu; jj++)
    for (std::size_t i = 0; i < ng; i++) {
      const double z = zc[i + ng * jj];
      szz[jj] += diag22[i] * z * z;
    }
  for (std::size_t s = 0; s < nf; s++) {
    const double f = fg[fatias[s].grupo].at(fatias[s].comp, fatias[s].comp);
    for (std::size_t k = 0; k < ng; k++) prec[col_geno(s, k)] += f * (1.0 - w) / w * diag22[k];
    for (std::size_t jj = 0; jj < mu; jj++)
      prec[M.total + s * mu + jj] = f * (szz[jj] / w + kd / (1.0 - w));
  }
  for (double& p : prec) if (!(p > 0.0)) p = 1.0;

  const int ngi = static_cast<int>(ng), mui = static_cast<int>(mu), inc1 = 1;
  const double um = 1.0, zero = 0.0;
  // por fatia: v nas colunas genotipadas, Z v nos marcadores, e as duas A22^-1
  std::vector<std::vector<double>> vg(nf, std::vector<double>(ng)),
      zv(nf, std::vector<double>(ng)), q3(nf), q4(nf);
  const std::size_t njob = 2 * nf;
  std::vector<std::vector<double>> rt1(njob, std::vector<double>(n1)),
      rt2(njob, std::vector<double>(n1));
  const int nth = static_cast<int>(
      std::min<std::size_t>(static_cast<std::size_t>(threads()), njob));
  std::vector<double> hs(ng), gsum(ng), gm(mu);
  auto aplica = [&](const std::vector<double>& v, std::vector<double>& y) {
    y.assign(N, 0.0);
    symv_tri(M.c, v.data(), y.data());
    for (std::size_t s = 0; s < nf; s++) {
      for (std::size_t k = 0; k < ng; k++) vg[s][k] = v[col_geno(s, k)];
      std::fill(zv[s].begin(), zv[s].end(), 0.0);
      if (mu > 0)
        F77_CALL(dgemv)("N", &ngi, &mui, &um, zc.data(), &ngi, v.data() + M.total + s * mu,
                        &inc1, &zero, zv[s].data(), &inc1 FCONE);
    }
    // as 2 nf aplicacoes de A22^-1 sao independentes; cada uma escreve so o seu vetor
#ifdef _OPENMP
#pragma omp parallel for num_threads(nth) schedule(static, 1)
#endif
    for (long jb = 0; jb < static_cast<long>(njob); jb++) {
      const std::size_t s = static_cast<std::size_t>(jb) / 2;
      if (jb % 2 == 0) a22inv_vezes(vg[s], q3[s], rt1[jb], rt2[jb]);
      else a22inv_vezes(zv[s], q4[s], rt1[jb], rt2[jb]);
    }
    // fatia a do grupo g: soma sobre as fatias b do MESMO grupo com f = (s2e K0^-1)(a, b)
    for (std::size_t s = 0; s < nf; s++) {
      const Densa& F = fg[fatias[s].grupo];
      std::fill(hs.begin(), hs.end(), 0.0);
      std::fill(gsum.begin(), gsum.end(), 0.0);
      for (std::size_t s2 = 0; s2 < nf; s2++) {
        if (fatias[s2].grupo != fatias[s].grupo) continue;
        const double f = F.at(fatias[s].comp, fatias[s2].comp);
        for (std::size_t k = 0; k < ng; k++) {
          hs[k] += f * ((1.0 - w) * q3[s2][k] - q4[s2][k]);
          gsum[k] += f * (q3[s2][k] - q4[s2][k]);
        }
        for (std::size_t jj = 0; jj < mu; jj++)
          y[M.total + s * mu + jj] += f * kd / (1.0 - w) * v[M.total + s2 * mu + jj];
      }
      for (std::size_t k = 0; k < ng; k++) y[col_geno(s, k)] += hs[k] / w;
      if (mu > 0) {
        F77_CALL(dgemv)("T", &ngi, &mui, &um, zc.data(), &ngi, gsum.data(), &inc1,
                        &zero, gm.data(), &inc1 FCONE);
        for (std::size_t jj = 0; jj < mu; jj++) y[M.total + s * mu + jj] -= gm[jj] / w;
      }
    }
  };

  // PCG classico com o precondicionador diagonal
  std::vector<double> x(N, 0.0), res(rhs), z(N), pdir(N), q(N);
  double nrhs = 0.0;
  for (double v : rhs) nrhs += v * v;
  nrhs = std::sqrt(nrhs);
  if (nrhs == 0.0) { r.ok = true; r.convergiu = true; r.solucao.assign(N, 0.0); return r; }
  for (std::size_t k = 0; k < N; k++) z[k] = res[k] / prec[k];
  double rho = 0.0;
  for (std::size_t k = 0; k < N; k++) rho += res[k] * z[k];
  pdir = z;
  std::size_t it = 0;
  double rel = 1.0;
  for (; it < maxiter; it++) {
    R_CheckUserInterrupt();
    double nr = 0.0;
    for (double v : res) nr += v * v;
    rel = std::sqrt(nr) / nrhs;
    if (verboso && it > 0 && it % 100 == 0)
      Rprintf("PCG iter %d  relative residual %.3e\n", (int) it, rel);
    if (rel < tol) break;
    aplica(pdir, q);
    double pq = 0.0;
    for (std::size_t k = 0; k < N; k++) pq += pdir[k] * q[k];
    if (!(pq > 0.0)) throw Erro("the ssSNPBLUP system lost positive-definiteness in the PCG");
    const double alfa = rho / pq;
    for (std::size_t k = 0; k < N; k++) { x[k] += alfa * pdir[k]; res[k] -= alfa * q[k]; }
    for (std::size_t k = 0; k < N; k++) z[k] = res[k] / prec[k];
    double rho2 = 0.0;
    for (std::size_t k = 0; k < N; k++) rho2 += res[k] * z[k];
    const double beta = rho2 / rho;
    rho = rho2;
    for (std::size_t k = 0; k < N; k++) pdir[k] = z[k] + beta * pdir[k];
  }

  r.ok = true;
  r.convergiu = rel < tol;
  r.iters = it;
  r.residuo = rel;
  r.solucao.assign(x.begin(), x.begin() + M.total);
  r.efeitos.assign(mm * nf, std::nan(""));
  for (std::size_t s = 0; s < nf; s++)
    for (std::size_t jj = 0; jj < mu; jj++)
      r.efeitos[s * mm + usados[jj]] = x[M.total + s * mu + jj];
  r.n_fixo = M.n_fixo;
  r.offset_grupo = M.offset_grupo;
  return r;
}

}  // namespace br
