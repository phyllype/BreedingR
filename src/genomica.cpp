// Passo unico: G de VanRaden (2008), A22^-1 pelo complemento de Schur, e H^-1
// (Aguilar et al. 2010; Christensen e Lund 2010).
//
// A armadilha silenciosa do passo unico esta toda numa confusao de blocos:
//
//   A22^-1  =  B22 - B21 B11^-1 B12        (Schur sobre A^-1 = [[B11,B12],[B21,B22]])
//
// e isso NAO e o bloco 22 de A^-1. As duas matrizes tem a mesma forma, as duas sao
// simetricas e definidas, e so uma esta certa. O gate que separa as duas exige que a
// fixture as distinga, senao o teste e cego.
//
// B11^-1 B12 nunca e formado denso: B11 e o bloco NAO-genotipado inteiro, e a 20 mil
// animais com 2 mil genotipados isso seria uma densa de 18k x 18k para obter uma resposta
// de 2k x 2k. Uma fatoracao esparsa de B11 e n_geno resolucoes triangulares fazem o mesmo.

#ifdef _OPENMP
#include <omp.h>
#endif
#include "mme.h"
#include <unordered_map>

// O BLAS do proprio R para a G: o produto Z Z' e o unico n^2 m do pacote, e o dsyrk
// faz em minutos o que o laco triplo fazia em horas de uma thread so.
#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

namespace br {

// As frequencias, a imputacao pela media e o denominador 2 sum p(1-p), uma vez so, para a G
// densa e para a rota APY que nunca forma G inteira.
struct FreqZ {
  std::vector<double> p;
  std::vector<char> usa;
  double denom = 0.0;
};

static FreqZ frequencias_z(Densa& m, RelatorioG& rel, bool meio) {
  const std::size_t n = m.nlin, nm = m.ncol;
  if (n == 0 || nm == 0) throw Erro("empty genotypes");
  FreqZ f;
  f.p.assign(nm, 0.0);
  f.usa.assign(nm, 1);
  for (std::size_t j = 0; j < nm; j++) {
    double soma = 0.0;
    std::size_t k = 0;
    for (std::size_t i = 0; i < n; i++) {
      const double x = m.at(i, j);
      if (std::isfinite(x)) { soma += x; k++; }
    }
    if (k == 0) { f.usa[j] = 0; rel.n_monomorficos++; continue; }
    const double media = soma / static_cast<double>(k);
    f.p[j] = meio ? 0.5 : media / 2.0;
    if (!meio && (f.p[j] <= 0.0 || f.p[j] >= 1.0)) { f.usa[j] = 0; rel.n_monomorficos++; continue; }
    for (std::size_t i = 0; i < n; i++)
      if (!std::isfinite(m.at(i, j))) { m.at(i, j) = media; rel.n_imputados++; }
  }
  for (std::size_t j = 0; j < nm; j++)
    if (f.usa[j]) f.denom += 2.0 * f.p[j] * (1.0 - f.p[j]);
  if (!(f.denom > 0.0)) throw Erro("2 sum p(1-p) is not positive: the markers do not vary");
  return f;
}

// G de VanRaden: Z Z' / (2 sum p(1-p)), com Z = M - 2p.
//
// Codigo ausente TEM de chegar como NaN e e imputado pela MEDIA do marcador. Nunca por
// zero: zero e um genotipo valido, e a confusao mudaria as frequencias e a G inteira sem
// nenhum erro visivel.
//
// meio = true e a G05 dos metafundadores (Legarra et al., 2015; Garcia-Baccino et al.,
// 2017): Z = M - 1 (todas as frequencias em 0.5) e escala m/2, TODOS os marcadores, o
// monomorfico inclusive (ele soma a mesma constante a todo par, que e parte da base de
// Gamma). E a G na base dos metafundadores, e por isso ela nao passa pelo ajuste a A22.
Densa vanraden_g(Densa& m, RelatorioG& rel, bool meio) {
  const std::size_t n = m.nlin, nm = m.ncol;
  const FreqZ fz = frequencias_z(m, rel, meio);
  const std::vector<double>& p = fz.p;
  const std::vector<char>& usa = fz.usa;
  const double denom = fz.denom;

  // Z Z' em blocos de marcadores: o bloco Z_b (n x B, coluna-major) e montado ja
  // centrado e acumulado em C += Z_b Z_b', pelos ladrilhos em paralelo (padrao) ou pelo
  // dsyrk do BLAS do R (br_threads(lapack = TRUE)). Os dois escrevem o mesmo triangulo:
  // o superior coluna-major, que e o inferior linha-major. Memoria extra: um bloco + o C,
  // nada proporcional a m.
  const std::size_t B = 2048;
  const int ni = static_cast<int>(n);
  std::vector<double> c(n * n, 0.0), zbuf(n * B);
  double beta = 0.0;
  std::size_t bcol = 0;
  auto acumula = [&]() {
    if (denso_lapack()) {
      const int k = static_cast<int>(bcol);
      const double um = 1.0;
      F77_CALL(dsyrk)("U", "N", &ni, &k, &um, zbuf.data(), &ni, &beta, c.data(), &ni
                      FCONE FCONE);
      beta = 1.0;
    } else {
      produto_ladrilhos(zbuf.data(), n, zbuf.data(), n, bcol, n, n, c.data(), n, true,
                        threads());
    }
    bcol = 0;
  };
  for (std::size_t j = 0; j < nm; j++) {
    if (!usa[j]) continue;
    const double dp = 2.0 * p[j];
    for (std::size_t i = 0; i < n; i++) zbuf[i + n * bcol] = m.at(i, j) - dp;
    if (++bcol == B) acumula();
  }
  if (bcol > 0) acumula();
  Densa g(n, n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t k2 = i; k2 < n; k2++) {
      const double v = c[i + n * k2] / denom;
      g.at(i, k2) = v;
      g.at(k2, i) = v;
    }
  return g;
}

// A22^-1 pelo Schur ESPARSO sobre A^-1, sem nunca formar B11^-1 denso, devolvido como o
// triangulo INFERIOR (na numeracao dos genotipados) so com as entradas que existem: o A22^-1
// e quase todo zero EXATO (medido: 97,6% em 400 genotipados, 98,9% em 800), e e isso que
// deixa a APY esparsa ate o H^-1. B12 tem so os pais e filhos nao genotipados de cada
// genotipado, meia duzia por coluna; guardado denso ele custava memoria n1 n2 e o produto
// B21 x (B11^-1 B12) virava n1 n2^2 flops sobre zeros. Uma resolucao por coluna genotipada,
// as colunas em paralelo (cada uma de uma thread, com a mesma ordem de soma).
Csc a22_inversa_esparsa(const Csc& ainv, const std::vector<std::size_t>& geno) {
  const std::size_t n = ainv.ncol;
  const std::size_t n2 = geno.size();
  std::vector<char> eh_geno(n, 0);
  for (std::size_t i : geno) {
    if (i >= n) throw Erro("genotyped index outside A^-1");
    eh_geno[i] = 1;
  }
  // numeracao: nao-genotipados primeiro
  std::vector<std::size_t> pos(n);
  std::size_t n1 = 0;
  for (std::size_t i = 0; i < n; i++) if (!eh_geno[i]) pos[i] = n1++;
  for (std::size_t k = 0; k < n2; k++) pos[geno[k]] = n1 + k;

  // B11 (triangulo inferior esparso), B12 esparso nas duas direcoes, B22 esparso por coluna
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  std::vector<std::vector<std::pair<std::uint32_t, double>>> b12_col(n2), b12_lin, b22_col(n2);
  for (std::size_t c = 0; c < n; c++)
    for (std::size_t k = ainv.colptr[c]; k < ainv.colptr[c + 1]; k++) {
      const std::size_t r = ainv.linha[k];
      const double x = ainv.valor[k];
      const bool rg = eh_geno[r], cg = eh_geno[c];
      const std::size_t pr = pos[r], pc = pos[c];
      if (!rg && !cg) {
        std::size_t i = pr, j = pc;
        if (i < j) std::swap(i, j);
        li.push_back(static_cast<std::uint32_t>(i));
        cj.push_back(static_cast<std::uint32_t>(j));
        v.push_back(x);
      } else if (rg && cg) {
        std::size_t i = pr - n1, j = pc - n1;
        if (i < j) std::swap(i, j);
        b22_col[j].push_back({static_cast<std::uint32_t>(i), x});
      } else {
        const std::size_t inng = rg ? pc : pr;
        const std::size_t ig = rg ? pr - n1 : pc - n1;
        b12_col[ig].push_back({static_cast<std::uint32_t>(inng), x});
      }
    }
  b12_lin.resize(n1);
  for (std::size_t ig = 0; ig < n2; ig++)
    for (const auto& [i, x] : b12_col[ig])
      b12_lin[i].push_back({static_cast<std::uint32_t>(ig), x});

  Csc L;
  std::vector<std::size_t> perm;
  if (n1 > 0) {
    Csc b11 = de_triplos(n1, n1, li, cj, v);
    perm = grau_minimo(b11);
    Csc pb = permuta_sim(b11, perm);
    Simbolica sb = simbolica(pb);
    if (!cholesky(pb, sb, L))
      throw Erro("the non-genotyped block of A^-1 is not positive-definite");
  }

  // coluna jg do triangulo inferior: B22(ig, jg) - B21 B11^-1 B12(ig, jg), ig >= jg. As
  // resolucoes vao em BLOCOS de W colunas: cada passada pelo fator serve as W de uma vez (a
  // leitura de L, que domina, e dividida por W), com as mesmas operacoes na mesma ordem de
  // uma resolucao sozinha, coluna a coluna. Os blocos em paralelo, cada um de uma thread.
  const std::size_t W = 32;
  std::vector<std::size_t> onde_perm(n1);
  for (std::size_t i = 0; i < n1; i++) onde_perm[perm[i]] = i;
  const std::size_t nblocos = (n2 + W - 1) / W;
  std::vector<std::vector<std::pair<std::uint32_t, double>>> saida(n2);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
  {
    std::vector<double> X(n1 * W), acc(n2, 0.0);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long bb = 0; bb < static_cast<long>(nblocos); bb++) {
      const std::size_t j0 = static_cast<std::size_t>(bb) * W, j1 = std::min(n2, j0 + W);
      const std::size_t w = j1 - j0;
      if (n1 > 0) {
        std::fill(X.begin(), X.begin() + static_cast<std::ptrdiff_t>(n1 * w), 0.0);
        for (std::size_t c = 0; c < w; c++)
          for (const auto& [i, x] : b12_col[j0 + c]) X[onde_perm[i] * w + c] = x;
        // L y = b, L' x = y, W colunas por vez (X linha-major, n1 x w, ja permutado)
        for (std::size_t j = 0; j < n1; j++) {
          double* xj = &X[j * w];
          const double d = L.valor[L.colptr[j]];
          for (std::size_t c = 0; c < w; c++) xj[c] /= d;
          for (std::size_t p = L.colptr[j] + 1; p < L.colptr[j + 1]; p++) {
            const double lp = L.valor[p];
            double* xr = &X[static_cast<std::size_t>(L.linha[p]) * w];
            for (std::size_t c = 0; c < w; c++) xr[c] -= lp * xj[c];
          }
        }
        for (std::size_t j = n1; j-- > 0;) {
          double* xj = &X[j * w];
          for (std::size_t p = L.colptr[j] + 1; p < L.colptr[j + 1]; p++) {
            const double lp = L.valor[p];
            const double* xr = &X[static_cast<std::size_t>(L.linha[p]) * w];
            for (std::size_t c = 0; c < w; c++) xj[c] -= lp * xr[c];
          }
          const double d = L.valor[L.colptr[j]];
          for (std::size_t c = 0; c < w; c++) xj[c] /= d;
        }
      }
      for (std::size_t c = 0; c < w; c++) {
        const std::size_t jg = j0 + c;
        // B21 x, so pelas entradas de B12, linha a linha em i crescente
        if (n1 > 0 && !b12_col[jg].empty())
          for (std::size_t i = 0; i < n1; i++) {
            const double ci = X[onde_perm[i] * w + c];
            if (ci == 0.0) continue;
            for (const auto& [ig, x] : b12_lin[i]) acc[ig] -= x * ci;
          }
        for (const auto& [ig, x] : b22_col[jg]) acc[ig] += x;
        auto& out = saida[jg];
        for (std::size_t ig = jg; ig < n2; ig++)
          if (acc[ig] != 0.0) out.push_back({static_cast<std::uint32_t>(ig), acc[ig]});
        std::fill(acc.begin(), acc.end(), 0.0);
      }
    }
  }
  Csc r(n2, n2);
  for (std::size_t jg = 0; jg < n2; jg++) r.colptr[jg + 1] = r.colptr[jg] + saida[jg].size();
  r.linha.reserve(r.colptr[n2]);
  r.valor.reserve(r.colptr[n2]);
  for (std::size_t jg = 0; jg < n2; jg++)
    for (const auto& [ig, x] : saida[jg]) {
      r.linha.push_back(ig);
      r.valor.push_back(x);
    }
  return r;
}

// A mesma inversa, densa e simetrica: a exportada por a22_inverse() e a da rota de
// metafundadores.
Densa a22_inversa(const Csc& ainv, const std::vector<std::size_t>& geno) {
  const Csc r = a22_inversa_esparsa(ainv, geno);
  Densa out(r.ncol, r.ncol);
  for (std::size_t jg = 0; jg < r.ncol; jg++)
    for (std::size_t k = r.colptr[jg]; k < r.colptr[jg + 1]; k++) {
      out.at(r.linha[k], jg) = r.valor[k];
      out.at(jg, r.linha[k]) = r.valor[k];
    }
  return out;
}

// A22 pelo algoritmo de Colleau (2002): A x = T D T' x, com T = (I - P)^-1 e P os pesos dos
// pais (1/2; 1/4 no avo materno de um pedigree pai/MGS), em tres passadas pelo pedigree
// ordenado. Custa O(n) por coluna genotipada, contra uma resolucao inteira no fator do bloco
// nao genotipado da rota de Schur, e e a rota do preGSf90 (Aguilar et al., 2011). Entrega A22,
// que o passo unico precisa para escalar G. Com metafundadores e a mesma conta sobre
// A(Gamma) = T L T' (Legarra et al., 2015): as linhas dos metafundadores nao tem pais e o
// bloco delas em L e a propria Gamma, cheia, aplicada como K (K' v) com a K K' = Gamma que
// o pedigree ja guarda; nas outras linhas L e a variancia mendeliana, que ja le o F dos pais
// na base de Gamma. As colunas sao independentes, cada uma de uma thread, com soma em ordem
// fixa.
struct Colleau {
  std::vector<double> d, w1, w2, f;
  std::vector<std::size_t> mf_linha;   // coluna de Gamma -> linha do pedigree
};

static Colleau prepara_colleau(const Pedigree& p) {
  const std::size_t n = p.ids.size();
  Colleau c;
  c.f = endogamia(p);
  c.d.resize(n); c.w1.resize(n); c.w2.resize(n);
  c.mf_linha.assign(p.n_mf, 0);
  for (std::size_t i = 0; i < n; i++) {
    const bool mf = !p.eh_mf.empty() && p.eh_mf[i];
    c.d[i] = mf ? 0.0 : variancia_mendeliana(p, c.f, i);
    if (mf) c.mf_linha[static_cast<std::size_t>(p.col_mf[i])] = i;
    c.w1[i] = 0.5;
    c.w2[i] = (!p.mgs.empty() && p.mgs[i]) ? 0.25 : 0.5;
  }
  return c;
}

// v <- A v, no lugar: z = T'x dos mais novos para os mais velhos, w = L z, y = T w dos mais
// velhos para os mais novos
static void colleau_aplica(const Pedigree& p, const Colleau& c, std::vector<double>& v) {
  const std::size_t n = v.size(), q = p.n_mf;
  for (std::size_t i = n; i-- > 0;) {
    if (v[i] == 0.0) continue;
    if (p.pai[i] >= 0) v[static_cast<std::size_t>(p.pai[i])] += c.w1[i] * v[i];
    if (p.mae[i] >= 0) v[static_cast<std::size_t>(p.mae[i])] += c.w2[i] * v[i];
  }
  std::vector<double> kv(q, 0.0);
  for (std::size_t col = 0; col < q; col++)
    for (std::size_t a = 0; a < q; a++)
      kv[col] += p.gama_chol[a * q + col] * v[c.mf_linha[a]];
  for (std::size_t i = 0; i < n; i++) v[i] *= c.d[i];
  for (std::size_t a = 0; a < q; a++) {
    double x = 0.0;
    for (std::size_t col = 0; col < q; col++) x += p.gama_chol[a * q + col] * kv[col];
    v[c.mf_linha[a]] = x;
  }
  for (std::size_t i = 0; i < n; i++) {
    if (p.pai[i] >= 0) v[i] += c.w1[i] * v[static_cast<std::size_t>(p.pai[i])];
    if (p.mae[i] >= 0) v[i] += c.w2[i] * v[static_cast<std::size_t>(p.mae[i])];
  }
}

static Densa a22_colleau(const Pedigree& p, const std::vector<std::size_t>& geno) {
  const std::size_t n = p.ids.size(), n2 = geno.size();
  const Colleau cl = prepara_colleau(p);
  Densa out(n2, n2);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
  {
    std::vector<double> v(n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 8)
#endif
    for (long jj = 0; jj < static_cast<long>(n2); jj++) {
      const std::size_t jg = static_cast<std::size_t>(jj);
      std::fill(v.begin(), v.end(), 0.0);
      v[geno[jg]] = 1.0;
      colleau_aplica(p, cl, v);
      for (std::size_t ig = 0; ig < n2; ig++) out.at(ig, jg) = v[geno[ig]];
    }
  }
  // simetriza contra o arredondamento das duas passadas
  for (std::size_t i = 0; i < n2; i++)
    for (std::size_t j = 0; j < i; j++) {
      const double mdi = 0.5 * (out.at(i, j) + out.at(j, i));
      out.at(i, j) = mdi;
      out.at(j, i) = mdi;
    }
  return out;
}

// So a mistura, G* = (1 - w) G + w A22: o caminho dos metafundadores, em que G ja esta na
// base de A(Gamma).
static Densa mistura_sem_ajuste(const Densa& g, const Densa& a22, double mistura) {
  const std::size_t n = g.nlin;
  if (a22.nlin != n) throw Erro("G and A22 with different sizes");
  if (!(mistura >= 0.0 && mistura <= 1.0)) throw Erro("blend outside [0, 1]");
  Densa out(n, n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      out.at(i, j) = (1.0 - mistura) * g.at(i, j) + mistura * a22.at(i, j);
  return out;
}

// Traz G a escala de A22 (ajuste afim nos momentos) e mistura: G* = (1-w) G_ajustada + w A22.
Densa ajusta_g_para_a22(const Densa& g, const Densa& a22, double mistura) {
  const std::size_t n = g.nlin;
  if (a22.nlin != n) throw Erro("G and A22 with different sizes");
  if (!(mistura >= 0.0 && mistura <= 1.0)) throw Erro("blend outside [0, 1]");

  double mg_diag = 0.0, ma_diag = 0.0, mg_off = 0.0, ma_off = 0.0;
  for (std::size_t i = 0; i < n; i++) {
    mg_diag += g.at(i, i);
    ma_diag += a22.at(i, i);
    for (std::size_t j = 0; j < i; j++) { mg_off += g.at(i, j); ma_off += a22.at(i, j); }
  }
  mg_diag /= n; ma_diag /= n;
  const double noff = n > 1 ? static_cast<double>(n) * (n - 1) / 2.0 : 1.0;
  mg_off /= noff; ma_off /= noff;

  // b (mg - mg_off) + a = ma  na diagonal e fora: resolve o afim a + b G
  const double b = (mg_diag - mg_off) != 0.0
      ? (ma_diag - ma_off) / (mg_diag - mg_off) : 1.0;
  const double a = ma_off - b * mg_off;

  Densa out(n, n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      out.at(i, j) = (1.0 - mistura) * (a + b * g.at(i, j)) + mistura * a22.at(i, j);
  return out;
}

// H^-1 = A^-1 + [0 0; 0 G*^-1 - A22^-1], devolvida como triplos do triangulo inferior.
Csc constroi_hinv(const Csc& ainv, const std::vector<std::size_t>& geno,
                  const Densa& gstar_inv, const Densa& a22_inv) {
  const std::size_t n = ainv.ncol;
  const std::size_t n2 = geno.size();
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  li.reserve(ainv.nnz() + n2 * (n2 + 1) / 2);
  cj.reserve(ainv.nnz() + n2 * (n2 + 1) / 2);
  v.reserve(ainv.nnz() + n2 * (n2 + 1) / 2);
  for (std::size_t c = 0; c < n; c++)
    for (std::size_t k = ainv.colptr[c]; k < ainv.colptr[c + 1]; k++) {
      li.push_back(ainv.linha[k]);
      cj.push_back(static_cast<std::uint32_t>(c));
      v.push_back(ainv.valor[k]);
    }
  // ZERO EXATO NAO VIRA TRIPLO. A APY (e a de Vecchia) produzem uma G^-1 com estrutura: o
  // bloco jovem x jovem sai so na diagonal, porque Mnn e diagonal, e esses zeros sao
  // exatos, nunca escritos. O A22^-1 tambem e quase todo zero exato (medido: 97.6% em 400
  // genotipados, 98.9% em 800). Empurrar os n2(n2+1)/2 pares inteiros para os triplos
  // enterrava tudo isso: o H^-1 saia denso no bloco genotipado, a fatoracao do MME ficava
  // O(n_geno^3) por iteracao e a APY comprava so estabilidade numerica, nenhum tempo
  // (medido: 12.10 s/iter denso contra 11.76 s/iter com nucleo de 300, em 2400
  // genotipados, dentro do ruido).
  //
  // O filtro e por zero EXATO e nao por tolerancia, de proposito: nao ha aproximacao
  // nenhuma aqui. A fracao de nao-zeros da diferenca e a mesma contando `!= 0` e contando
  // `|x| > 1e-12` (75.12% e 43.84% nos dois testes), entao o que sai e exatamente o que
  // nao existe. Um valor pequeno mas real continua entrando.
  for (std::size_t a = 0; a < n2; a++)
    for (std::size_t b = 0; b <= a; b++) {
      const double x = gstar_inv.at(a, b) - a22_inv.at(a, b);
      if (x == 0.0) continue;
      std::size_t i = geno[a], j = geno[b];
      if (i < j) std::swap(i, j);
      li.push_back(static_cast<std::uint32_t>(i));
      cj.push_back(static_cast<std::uint32_t>(j));
      v.push_back(x);
    }
  return de_triplos(n, n, li, cj, v);
}

// H^-1 com APY SEM nenhuma matriz n_geno x n_geno. A rota densa forma G, A22, G*, a inversa
// APY e o A22^-1 inteiros: com 50 mil genotipados sao 20 GB por matriz, e a APY existe
// justamente para esse tamanho. Aqui so aparece o que a APY usa (Misztal, Legarra e Aguilar,
// 2014): G nas linhas do NUCLEO (c x n, pelos ladrilhos em blocos de marcadores) e a
// diagonal; as medias do ajuste afim por somas (1'G1 = ||Z'1||^2 / k, e 1'A22 1 por uma
// passada de Colleau); A22 nas colunas do nucleo por Colleau; os blocos da inversa APY
// (nucleo x nucleo, nucleo x jovem, diagonal dos jovens); e o A22^-1 esparso pelo Schur. A
// memoria fica O(c n + nnz(A22^-1)). As mesmas contas da rota densa, na mesma ordem de
// definicao: o portao de exatidao da APY parcial continua valendo.
//
// Com metafundadores (com_mf) G e a G05 e nao ha ajuste afim, como na rota densa: a0 = 0,
// b = 1, e a A(Gamma)22 das linhas do nucleo sai do mesmo Colleau, que ja aplica Gamma.
static Csc h_inversa_apy(const Pedigree& ped, const Csc& ainv, const std::vector<std::size_t>& idx,
                         Densa& m, double mistura, const std::vector<std::size_t>& nuc,
                         bool com_mf, RelatorioG& rel) {
  const std::size_t n = m.nlin, nm = m.ncol, c = nuc.size();
  if (c < 2) throw Erro("the APY core needs at least 2 animals");
  if (!(mistura >= 0.0 && mistura <= 1.0)) throw Erro("blend outside [0, 1]");
  std::vector<char> eh_nuc(n, 0);
  for (std::size_t a : nuc) {
    if (a >= n) throw Erro("APY core index outside the genotyped set");
    eh_nuc[a] = 1;
  }
  std::vector<std::size_t> jov;
  for (std::size_t i = 0; i < n; i++) if (!eh_nuc[i]) jov.push_back(i);
  const std::size_t nj = jov.size();
  const int nth = threads();
  const FreqZ fz = frequencias_z(m, rel, com_mf);

  // 1. G nas linhas do nucleo, a diagonal de G e 1'ZZ'1, por blocos de marcadores
  Densa gs(c, n);
  std::vector<double> gdiag(n, 0.0);
  double soma1 = 0.0;
  const std::size_t Bm = 1024;
  std::vector<double> zb(n * Bm), zcb(c * Bm);
  std::size_t bcol = 0;
  auto acumula = [&]() {
    produto_ladrilhos(zcb.data(), c, zb.data(), n, bcol, c, n, gs.dados.data(), n, false, nth);
    bcol = 0;
  };
  for (std::size_t j = 0; j < nm; j++) {
    if (!fz.usa[j]) continue;
    const double dp = 2.0 * fz.p[j];
    double* zk = &zb[bcol * n];
    double sj = 0.0;
    for (std::size_t i = 0; i < n; i++) {
      const double z = m.at(i, j) - dp;
      zk[i] = z;
      gdiag[i] += z * z;
      sj += z;
    }
    double* zck = &zcb[bcol * c];
    for (std::size_t a = 0; a < c; a++) zck[a] = zk[nuc[a]];
    soma1 += sj * sj;
    if (++bcol == Bm) acumula();
  }
  if (bcol > 0) acumula();
  std::vector<double>().swap(zb);
  std::vector<double>().swap(zcb);
  for (double& x : gs.dados) x /= fz.denom;
  for (double& x : gdiag) x /= fz.denom;
  soma1 /= fz.denom;

  // 2. o ajuste afim G -> A22 pelas medias, sem G nem A22 inteiras (nenhum com Gamma)
  const Colleau cl = prepara_colleau(ped);
  const double nn = static_cast<double>(n);
  double a0 = 0.0, b = 1.0;
  if (!com_mf) {
    double tr_g = 0.0, tr_a = 0.0;
    for (std::size_t i = 0; i < n; i++) { tr_g += gdiag[i]; tr_a += 1.0 + cl.f[idx[i]]; }
    double soma_a = 0.0;
    std::vector<double> um(ped.ids.size(), 0.0);
    for (std::size_t i = 0; i < n; i++) um[idx[i]] = 1.0;
    colleau_aplica(ped, cl, um);
    for (std::size_t i = 0; i < n; i++) soma_a += um[idx[i]];
    const double noff = n > 1 ? nn * (nn - 1.0) : 1.0;
    const double mg_diag = tr_g / nn, ma_diag = tr_a / nn;
    const double mg_off = (soma1 - tr_g) / noff, ma_off = (soma_a - tr_a) / noff;
    b = (mg_diag - mg_off) != 0.0 ? (ma_diag - ma_off) / (mg_diag - mg_off) : 1.0;
    a0 = ma_off - b * mg_off;
  }

  // 3. G* nas linhas do nucleo (A22 delas por Colleau) e na diagonal
#ifdef _OPENMP
#pragma omp parallel num_threads(nth)
#endif
  {
    std::vector<double> v(ped.ids.size());
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif
    for (long aa = 0; aa < static_cast<long>(c); aa++) {
      const std::size_t a = static_cast<std::size_t>(aa);
      std::fill(v.begin(), v.end(), 0.0);
      v[idx[nuc[a]]] = 1.0;
      colleau_aplica(ped, cl, v);
      double* g = gs.linha(a);
      for (std::size_t i = 0; i < n; i++)
        g[i] = (1.0 - mistura) * (a0 + b * g[i]) + mistura * v[idx[i]];
    }
  }
  std::vector<double> gsd(n);
  for (std::size_t i = 0; i < n; i++)
    gsd[i] = (1.0 - mistura) * (a0 + b * gdiag[i]) + mistura * (1.0 + cl.f[idx[i]]);
  rel.diag_gstar = gsd;

  // 4. APY: Gcc^-1, P = Gcc^-1 Gcn, residuos mendelianos, e os blocos da inversa
  Densa gcc(c, c);
  for (std::size_t a = 0; a < c; a++)
    for (std::size_t bb = 0; bb < c; bb++)
      gcc.at(a, bb) = 0.5 * (gs.at(a, nuc[bb]) + gs.at(bb, nuc[a]));
  Densa gcci = inv_pd(gcc);
  Densa pm(c, n);
  produto_ladrilhos(gcci.dados.data(), c, gs.dados.data(), n, c, c, n, pm.dados.data(), n,
                    false, nth);
  double diag_media = 0.0;
  for (double x : gsd) diag_media += x;
  diag_media /= nn;
  const double limiar = 1e-8 * diag_media;
  std::vector<double> mend(nj);
  for (std::size_t k = 0; k < nj; k++) mend[k] = gsd[jov[k]];
  for (std::size_t a = 0; a < c; a++) {
    const double* g = gs.linha(a);
    const double* pa = pm.linha(a);
    for (std::size_t k = 0; k < nj; k++) mend[k] -= g[jov[k]] * pa[jov[k]];
  }
  std::vector<double> minv(nj);
  for (std::size_t k = 0; k < nj; k++) {
    if (mend[k] <= limiar)
      throw Erro("degenerate Mendelian residual in APY: a non-core animal is collinear with the core "
                 "(clone or duplicate). Enlarge the core or the blend.");
    minv[k] = 1.0 / mend[k];
  }
  Densa().dados.swap(gs.dados);
  // nucleo x nucleo: Gcc^-1 + P M^-1 P', pelos ladrilhos com Ps' (nj x c) "por k"
  Densa cc = gcci;
  {
    std::vector<double> pst(nj * c);
    for (std::size_t a = 0; a < c; a++) {
      const double* pa = pm.linha(a);
      for (std::size_t k = 0; k < nj; k++) pst[k * c + a] = pa[jov[k]] * std::sqrt(minv[k]);
    }
    produto_ladrilhos(pst.data(), c, pst.data(), c, nj, c, c, cc.dados.data(), c, true, nth);
  }

  // 5. A22^-1 esparso e os triplos do H^-1
  const Csc a22i = a22_inversa_esparsa(ainv, idx);
  const std::size_t nh = ainv.ncol;
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  const std::size_t total = ainv.nnz() + c * (c + 1) / 2 + c * nj + nj + a22i.nnz();
  li.reserve(total); cj.reserve(total); v.reserve(total);
  auto poe = [&](std::size_t i, std::size_t j, double x) {
    if (x == 0.0) return;
    if (i < j) std::swap(i, j);
    li.push_back(static_cast<std::uint32_t>(i));
    cj.push_back(static_cast<std::uint32_t>(j));
    v.push_back(x);
  };
  for (std::size_t col = 0; col < nh; col++)
    for (std::size_t k = ainv.colptr[col]; k < ainv.colptr[col + 1]; k++)
      poe(ainv.linha[k], col, ainv.valor[k]);
  for (std::size_t a = 0; a < c; a++)
    for (std::size_t bb = 0; bb <= a; bb++) poe(idx[nuc[a]], idx[nuc[bb]], cc.at(a, bb));
  for (std::size_t a = 0; a < c; a++) {
    const double* pa = pm.linha(a);
    for (std::size_t k = 0; k < nj; k++) poe(idx[jov[k]], idx[nuc[a]], -pa[jov[k]] * minv[k]);
  }
  for (std::size_t k = 0; k < nj; k++) poe(idx[jov[k]], idx[jov[k]], minv[k]);
  for (std::size_t jg = 0; jg < a22i.ncol; jg++)
    for (std::size_t k = a22i.colptr[jg]; k < a22i.colptr[jg + 1]; k++)
      poe(idx[a22i.linha[k]], idx[jg], -a22i.valor[k]);
  return de_triplos(nh, nh, li, cj, v);
}

// Substitui o K^-1 dos grupos com parentesco pelo H^-1, e o logdet correspondente.
//
// Genotipado fora do pedigree e ERRO explicito, nao descarte silencioso: um genotipado sem
// linha em A nao tem onde entrar em H^-1, e some-lo mudaria a analise sem aviso.
// A inversa APY (Misztal, Legarra e Aguilar, 2014) de uma G densa ja ajustada: nucleo
// exato, jovens por recursao condicional.
//
//   G_APY^-1 = [ Gcc^-1 + P Mnn^-1 P\'   -P Mnn^-1 ]      P = Gcc^-1 Gcn
//              [ -Mnn^-1 P\'              Mnn^-1   ]      Mnn = diag(g_ii - g_ic P_i)
//
// O residuo mendeliano degenerado tem limiar RELATIVO a escala de G: um clone entre os
// jovens da residuo de +-1e-13, e dividir por ele poria 1e13 dentro da inversa.
static Densa apy_de(const Densa& g, const std::vector<std::size_t>& nucleo) {
  const std::size_t n = g.nlin;
  std::vector<char> eh_nucleo(n, 0);
  for (std::size_t c : nucleo) {
    if (c >= n) throw Erro("APY core index outside the genotyped set");
    eh_nucleo[c] = 1;
  }
  std::vector<std::size_t> jovens;
  for (std::size_t i = 0; i < n; i++) if (!eh_nucleo[i]) jovens.push_back(i);
  const std::size_t nc = nucleo.size(), nj = jovens.size();
  if (nc < 2) throw Erro("the APY core needs at least 2 animals");

  Densa gcc(nc, nc);
  for (std::size_t a = 0; a < nc; a++)
    for (std::size_t b = 0; b < nc; b++) gcc.at(a, b) = g.at(nucleo[a], nucleo[b]);
  Densa gcc_inv = inv_pd(gcc);

  Densa out(n, n);
  if (nj == 0) {
    for (std::size_t a = 0; a < nc; a++)
      for (std::size_t b = 0; b < nc; b++) out.at(nucleo[a], nucleo[b]) = gcc_inv.at(a, b);
    return out;
  }

  // Os dois produtos nc^2 nj vao aos ladrilhos em paralelo, ou ao BLAS do R com
  // br_threads(lapack = TRUE). P = Gcc^-1 Gcn: P(a, b) = soma_k Gcc^-1(k, a) Gcn(k, b), com
  // as duas "por k" como estao guardadas (Gcc^-1 simetrica). Para o BLAS coluna-major cada
  // Densa linha-major e a transposta: P' (nj x nc) = Gcn' Gcc^-1.
  Densa gcn(nc, nj), pmat(nc, nj);
  for (std::size_t a = 0; a < nc; a++)
    for (std::size_t b = 0; b < nj; b++) gcn.at(a, b) = g.at(nucleo[a], jovens[b]);
  if (denso_lapack()) {
    const int m_ = static_cast<int>(nj), n_ = static_cast<int>(nc);
    const double um = 1.0, zero = 0.0;
    F77_CALL(dgemm)("N", "N", &m_, &n_, &n_, &um, gcn.dados.data(), &m_,
                    gcc_inv.dados.data(), &n_, &zero, pmat.dados.data(), &m_ FCONE FCONE);
  } else {
    produto_ladrilhos(gcc_inv.dados.data(), nc, gcn.dados.data(), nj, nc, nc, nj,
                      pmat.dados.data(), nj, false, threads());
  }

  double diag_media = 0.0;
  for (std::size_t i = 0; i < n; i++) diag_media += g.at(i, i);
  diag_media /= static_cast<double>(n);
  const double limiar = 1e-8 * diag_media;

  // residuo mendeliano de cada jovem, com as linhas de Gcn e P lidas em sequencia
  std::vector<double> mendv(nj);
  for (std::size_t b = 0; b < nj; b++) mendv[b] = g.at(jovens[b], jovens[b]);
  for (std::size_t k = 0; k < nc; k++) {
    const double* gk = gcn.linha(k);
    const double* pk = pmat.linha(k);
    for (std::size_t b = 0; b < nj; b++) mendv[b] -= gk[b] * pk[b];
  }
  std::vector<double> minv(nj);
  for (std::size_t b = 0; b < nj; b++) {
    const double mend = mendv[b];
    if (mend <= limiar)
      throw Erro("degenerate Mendelian residual in APY: a non-core animal is collinear with the core "
                 "(clone or duplicate). Enlarge the core or the blend.");
    minv[b] = 1.0 / mend;
  }

  // bloco do nucleo: Gcc^-1 + P Mnn^-1 P' = Gcc^-1 + Ps Ps', Ps = P Mnn^-1/2 (Mnn > 0 pelo
  // limiar acima), por dsyrk: C = A' A com A = Ps' (nj x nc) coluna-major
  {
    Densa ps = pmat;
    for (std::size_t a = 0; a < nc; a++) {
      double* pa = ps.linha(a);
      for (std::size_t k = 0; k < nj; k++) pa[k] *= std::sqrt(minv[k]);
    }
    Densa cc = gcc_inv;
    if (denso_lapack()) {
      const int n_ = static_cast<int>(nc), k_ = static_cast<int>(nj);
      const double um = 1.0;
      F77_CALL(dsyrk)("U", "T", &n_, &k_, &um, ps.dados.data(), &k_, &um, cc.dados.data(), &n_
                      FCONE FCONE);
    } else {
      // os ladrilhos querem "por k" (k = jovem): Ps' linha-major, nj x nc
      Densa pst(nj, nc);
      for (std::size_t a = 0; a < nc; a++)
        for (std::size_t k = 0; k < nj; k++) pst.at(k, a) = ps.at(a, k);
      produto_ladrilhos(pst.dados.data(), nc, pst.dados.data(), nc, nj, nc, nc,
                        cc.dados.data(), nc, true, threads());
    }
    // "U" coluna-major = triangulo inferior da linha-major: espelha a partir dele
    for (std::size_t a = 0; a < nc; a++)
      for (std::size_t b = 0; b <= a; b++) {
        out.at(nucleo[a], nucleo[b]) = cc.at(a, b);
        out.at(nucleo[b], nucleo[a]) = cc.at(a, b);
      }
  }
  for (std::size_t a = 0; a < nc; a++)
    for (std::size_t b = 0; b < nj; b++) {
      const double v = -pmat.at(a, b) * minv[b];
      out.at(nucleo[a], jovens[b]) = v;
      out.at(jovens[b], nucleo[a]) = v;
    }
  for (std::size_t b = 0; b < nj; b++) out.at(jovens[b], jovens[b]) = minv[b];
  return out;
}

// A inversa de Vecchia (1988): cada animal condiciona nos SEUS k vizinhos mais proximos entre
// os anteriores, nao num nucleo global. E a generalizacao da APY, e do proprio A^-1 de
// Henderson (1976), que e exatamente Vecchia com os PAIS como conjunto de condicionamento
// (exato porque o pedigree e markoviano). Schafer, Katzfuss & Owhadi (2021) mostram que,
// dado o padrao, o fator esparso assim construido minimiza a divergencia KL; conjuntos
// maiores nunca pioram.
//
//   coluna i de U:  b = G[c,c]^-1 G[c,i],  d = g_ii - G[c,i]' b   (residuo mendeliano)
//                   U[c,i] = -b/sqrt(d),   U[i,i] = 1/sqrt(d),    G^-1 ~ U U'
//
// A ORDEM importa e e a ordem das linhas dos genotipos: ancestrais primeiro e o
// natural, porque condicionar no passado e o que a recursao significa. A vizinhanca e
// escolhida pelo MAIOR |g_ij| entre os anteriores (parentesco mais forte), com
// desempate deterministico pelo indice.
static Densa vecchia_de(const Densa& g, std::size_t k) {
  const std::size_t n = g.nlin;
  if (k < 1) throw Erro("vecchia_k must be at least 1");
  double diag_media = 0.0;
  for (std::size_t i = 0; i < n; i++) diag_media += g.at(i, i);
  diag_media /= static_cast<double>(n);
  const double limiar = 1e-8 * diag_media;

  Densa out(n, n);
  std::vector<std::size_t> cand;
  std::vector<double> u(n);
  for (std::size_t i = 0; i < n; i++) {
    const std::size_t m2 = std::min(k, i);
    cand.resize(i);
    for (std::size_t j = 0; j < i; j++) cand[j] = j;
    if (m2 < i)
      std::partial_sort(cand.begin(), cand.begin() + m2, cand.end(),
          [&](std::size_t a, std::size_t b) {
            const double fa = std::fabs(g.at(a, i)), fb = std::fabs(g.at(b, i));
            if (fa != fb) return fa > fb;
            return a < b;
          });
    cand.resize(m2);
    std::sort(cand.begin(), cand.end());

    double d = g.at(i, i);
    std::vector<double> b(m2, 0.0);
    if (m2 > 0) {
      Densa s(m2, m2);
      for (std::size_t a = 0; a < m2; a++)
        for (std::size_t c = 0; c < m2; c++) s.at(a, c) = g.at(cand[a], cand[c]);
      std::vector<double> rhs(m2);
      for (std::size_t a = 0; a < m2; a++) rhs[a] = g.at(cand[a], i);
      if (!chol_densa(s))
        throw Erro("a Vecchia conditioning block is not positive-definite; raise the blend");
      // resolve L L' b = rhs, denso e pequeno (m2 x m2)
      for (std::size_t a = 0; a < m2; a++) {
        double acc = rhs[a];
        for (std::size_t c = 0; c < a; c++) acc -= s.at(a, c) * b[c];
        b[a] = acc / s.at(a, a);
      }
      for (std::size_t a = m2; a-- > 0;) {
        double acc = b[a];
        for (std::size_t c = a + 1; c < m2; c++) acc -= s.at(c, a) * b[c];
        b[a] = acc / s.at(a, a);
      }
      for (std::size_t a = 0; a < m2; a++) d -= g.at(cand[a], i) * b[a];
    }
    if (d <= limiar)
      throw Erro("degenerate Mendelian residual in the Vecchia recursion: an animal is "
                 "collinear with its neighborhood (clone or duplicate). Raise the blend "
                 "or remove the duplicate.");
    const double sd = std::sqrt(d);
    // coluna i de U sobre (cand, i); acumula o produto externo em G^-1 ~ U U'
    for (std::size_t a = 0; a < m2; a++) u[a] = -b[a] / sd;
    const double uii = 1.0 / sd;
    for (std::size_t a = 0; a < m2; a++) {
      for (std::size_t c = 0; c <= a; c++) {
        const double v = u[a] * u[c];
        out.at(cand[a], cand[c]) += v;
        if (a != c) out.at(cand[c], cand[a]) += v;
      }
      const double v = u[a] * uii;
      out.at(cand[a], i) += v;
      out.at(i, cand[a]) += v;
    }
    out.at(i, i) += uii * uii;
  }
  return out;
}

// O nucleo, comum aos tres desenhos: tudo o que a genomica toca sao os grupos do modelo
// e os K^-1 com seus log-determinantes, uni, multi e AR(1) carregam exatamente esses
// campos, entao o passo unico e UM caminho, nao tres.
// A H^-1 do passo unico, de um pedigree e seus genotipos: A^-1 + [0 0; 0 G*^-1 - A22^-1],
// com G* a G de VanRaden trazida a escala de A22 e misturada, invertida exata, pela APY ou
// por Vecchia. E o nucleo que os ajustadores usam (aplica_genomica_em) e o que h_inverse()
// exporta, para que os motores em R (limiar, sobrevivencia) tenham o MESMO passo unico.
Csc h_inversa(const Pedigree& ped, const Csc& ainv, const std::vector<std::string>& geno_ids,
              Densa& m, double mistura, const std::vector<std::string>& nucleo_apy,
              std::size_t vecchia_k, RelatorioG& rel) {
  std::vector<std::size_t> idx;
  idx.reserve(geno_ids.size());
  {
    std::size_t n = ped.ids.size();
    // mapa id -> posicao
    std::vector<std::pair<std::string, std::size_t>> tmp;
    tmp.reserve(n);
    for (std::size_t i = 0; i < n; i++) tmp.push_back({ped.ids[i], i});
    std::sort(tmp.begin(), tmp.end());
    for (const std::string& g : geno_ids) {
      auto it = std::lower_bound(tmp.begin(), tmp.end(), std::make_pair(g, std::size_t(0)),
          [](const auto& a, const auto& b){ return a.first < b.first; });
      if (it == tmp.end() || it->first != g)
        throw Erro("genotyped animal '" + g + "' is not in the pedigree");
      idx.push_back(it->second);
    }
  }

  const bool com_mf = !ped.eh_mf.empty() &&
                      std::any_of(ped.eh_mf.begin(), ped.eh_mf.end(), [](char c) { return c != 0; });
  // Com metafundadores (restricao #15): G05 na base de Gamma, A22 e a A(Gamma)22 pela rota
  // de Schur sobre a A(Gamma)^-1, e G* = (1 - w) G05 + w A22 SEM o ajuste afim, que e
  // exatamente a correcao de base que Gamma ja faz (fazer os dois corrige a base duas vezes).
  if (vecchia_k == 0 && !nucleo_apy.empty()) {
    std::vector<std::size_t> nc_idx;
    nc_idx.reserve(nucleo_apy.size());
    std::unordered_map<std::string, std::size_t> onde;
    onde.reserve(geno_ids.size());
    for (std::size_t k = 0; k < geno_ids.size(); k++) onde.emplace(geno_ids[k], k);
    for (const std::string& nid : nucleo_apy) {
      const auto it = onde.find(nid);
      if (it == onde.end())
        throw Erro("APY core animal '" + nid + "' is not among the genotyped");
      nc_idx.push_back(it->second);
    }
    rel.linha_ped = idx;
    return h_inversa_apy(ped, ainv, idx, m, mistura, nc_idx, com_mf, rel);
  }
  Densa g = vanraden_g(m, rel, com_mf);
  Densa a22, a22i;
  if (com_mf) {
    a22i = a22_inversa(ainv, idx);
    a22 = a22_colleau(ped, idx);
  } else {
    a22 = a22_colleau(ped, idx);
    // Com Vecchia a G^-1 e esparsa, e o A22^-1 tem de chegar com os zeros EXATOS do Schur:
    // a inversa numerica da A22 poe ruido de arredondamento onde o verdadeiro e zero, e o
    // bloco do H^-1 enche (medido: 330.529 contra 306.512 nao-zeros, 800 genotipados, k = 20).
    // Na rota exata a G^-1 ja e densa, e a inversa da A22 de Colleau basta.
    a22i = vecchia_k > 0 ? a22_inversa(ainv, idx) : inv_pd(a22);
  }
  Densa gstar = com_mf ? mistura_sem_ajuste(g, a22, mistura) : ajusta_g_para_a22(g, a22, mistura);
  // a priori de cada genotipado, guardada AQUI porque este e o unico ponto em que G*
  // existe formada; accuracy() a usa no lugar de 1 + F do pedigree
  rel.diag_gstar.resize(gstar.nlin);
  for (std::size_t q = 0; q < gstar.nlin; q++) rel.diag_gstar[q] = gstar.at(q, q);
  rel.linha_ped = idx;

  // Com nucleo APY declarado, a inversa de G* e a APY; com vecchia_k, a de Vecchia;
  // sem, a exata. Os gates que sustentam os caminhos aproximados sao os colapsos:
  // nucleo = TODOS e k >= n - 1 tem de reproduzir o ajuste exato identicamente,
  // porque nesses casos a aproximacao E a propria inversa.
  if (!nucleo_apy.empty() && vecchia_k > 0)
    throw Erro("apy_core and vecchia_k are two approximations of the same inverse: "
               "declare one");
  Densa gstar_inv;
  if (vecchia_k > 0) {
    gstar_inv = vecchia_de(gstar, vecchia_k);
  } else if (nucleo_apy.empty()) {
    gstar_inv = inv_pd(gstar);
  } else {
    std::vector<std::size_t> nc_idx;
    nc_idx.reserve(nucleo_apy.size());
    std::unordered_map<std::string, std::size_t> onde;
    onde.reserve(geno_ids.size());
    for (std::size_t k = 0; k < geno_ids.size(); k++) onde.emplace(geno_ids[k], k);
    for (const std::string& nid : nucleo_apy) {
      const auto it = onde.find(nid);
      if (it == onde.end())
        throw Erro("APY core animal '" + nid + "' is not among the genotyped");
      nc_idx.push_back(it->second);
    }
    gstar_inv = apy_de(gstar, nc_idx);
  }
  return constroi_hinv(ainv, idx, gstar_inv, a22i);
}

static RelatorioG aplica_genomica_em(const Modelo& modelo, std::vector<Csc>& kinv,
                                     std::vector<double>& kinv_logdet, const Pedigree& ped,
                                     const std::vector<std::string>& geno_ids, Densa& m,
                                     double mistura,
                                     const std::vector<std::string>& nucleo_apy,
                                     std::size_t vecchia_k) {
  RelatorioG rel;
  // o A^-1 ja esta no desenho (primeiro grupo com parentesco)
  const Csc* ainv = nullptr;
  for (std::size_t g = 0; g < modelo.grupos.size(); g++)
    if (modelo.grupos[g].estrutura == Estrutura::Parentesco) { ainv = &kinv[g]; break; }
  if (!ainv) throw Erro("there is no relationship group to receive the genomics");

  Csc hinv = h_inversa(ped, *ainv, geno_ids, m, mistura, nucleo_apy, vecchia_k, rel);

  // log|H^-1| pela fatoracao esparsa
  std::vector<std::size_t> perm = grau_minimo(hinv);
  Csc ph = permuta_sim(hinv, perm);
  Simbolica sb = simbolica(ph);
  Csc L;
  if (!cholesky(ph, sb, L))
    throw Erro("H^-1 is not positive-definite; raise the blend");
  const double ld = logdet(L);

  for (std::size_t gg = 0; gg < modelo.grupos.size(); gg++)
    if (modelo.grupos[gg].estrutura == Estrutura::Parentesco) {
      kinv[gg] = hinv;
      kinv_logdet[gg] = ld;
    }
  return rel;
}

RelatorioG aplica_genomica(Desenho& d, const Pedigree& ped,
                           const std::vector<std::string>& geno_ids, Densa& m,
                           double mistura, const std::vector<std::string>& nucleo_apy,
                           std::size_t vecchia_k) {
  return aplica_genomica_em(d.modelo, d.kinv, d.kinv_logdet, ped, geno_ids, m, mistura,
                            nucleo_apy, vecchia_k);
}

RelatorioG aplica_genomica(DesenhoMT& d, const Pedigree& ped,
                           const std::vector<std::string>& geno_ids, Densa& m,
                           double mistura, const std::vector<std::string>& nucleo_apy,
                           std::size_t vecchia_k) {
  return aplica_genomica_em(d.modelo, d.kinv, d.kinv_logdet, ped, geno_ids, m, mistura,
                            nucleo_apy, vecchia_k);
}

RelatorioG aplica_genomica(DesenhoAR& d, const Pedigree& ped,
                           const std::vector<std::string>& geno_ids, Densa& m,
                           double mistura, const std::vector<std::string>& nucleo_apy,
                           std::size_t vecchia_k) {
  return aplica_genomica_em(d.modelo, d.kinv, d.kinv_logdet, ped, geno_ids, m, mistura,
                            nucleo_apy, vecchia_k);
}

}  // namespace br
