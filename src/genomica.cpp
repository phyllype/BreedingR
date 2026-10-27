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
#include <cstdio>
#include <cstring>
#include <unordered_map>

// O BLAS do proprio R para a G: o produto Z Z' e o unico n^2 m do pacote, e o dsyrk
// faz em minutos o que o laco triplo fazia em horas de uma thread so.
#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#ifndef FCONE
# define FCONE
#endif

namespace br {

// As frequencias e o denominador 2 sum p(1-p), uma vez so, para a G densa, para a rota APY
// que nunca forma G inteira e para o Lanczos. A MEDIA de cada marcador fica guardada: e ela
// que entra no lugar do ausente quando o marcador e lido (bloco_z), sem imputar no lugar.
struct FreqZ {
  std::vector<double> p, media;
  std::vector<char> usa;
  double denom = 0.0;
};

static FreqZ frequencias_z(const Genotipos& gt, RelatorioG& rel, bool meio) {
  const std::size_t n = gt.n, nm = gt.m;
  if (n == 0 || nm == 0) throw Erro("empty genotypes");
  FreqZ f;
  f.p.assign(nm, 0.0);
  f.media.assign(nm, 0.0);
  f.usa.assign(nm, 1);
  std::vector<std::size_t> faltam(nm, 0);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
  {
    std::vector<double> col(n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
    for (long jj = 0; jj < static_cast<long>(nm); jj++) {
      const std::size_t j = static_cast<std::size_t>(jj);
      gt.coluna(j, col.data());
      double soma = 0.0, mn = 0.0, mx = 0.0;
      std::size_t k = 0;
      for (std::size_t i = 0; i < n; i++)
        if (std::isfinite(col[i])) {
          if (k == 0) mn = mx = col[i];
          mn = std::min(mn, col[i]);
          mx = std::max(mx, col[i]);
          soma += col[i];
          k++;
        }
      if (k == 0) { f.usa[j] = 0; continue; }
      // o MESMO genotipo em todo animal (todos 1: uma falha de chamada, ou um painel imputado)
      // nao informa parentesco nenhum, e com p = 0,5 nao e monomorfico pela frequencia: entrava
      // so no denominador, e um painel inteiro assim dava G = 0. Fica fora, como o monomorfico
      // (a G05 de meio, com todos os marcadores por construcao, segue como esta)
      if (!meio && mn == mx) { f.usa[j] = 0; continue; }
      const double media = soma / static_cast<double>(k);
      f.media[j] = media;
      f.p[j] = meio ? 0.5 : media / 2.0;
      if (!meio && (f.p[j] <= 0.0 || f.p[j] >= 1.0)) { f.usa[j] = 0; continue; }
      faltam[j] = n - k;
    }
  }
  for (std::size_t j = 0; j < nm; j++) {
    if (!f.usa[j]) rel.n_monomorficos++;
    else rel.n_imputados += faltam[j];
  }
  for (std::size_t j = 0; j < nm; j++)
    if (f.usa[j]) f.denom += 2.0 * f.p[j] * (1.0 - f.p[j]);
  if (!(f.denom > 0.0))
    throw Erro("no marker varies among the genotyped animals (each one is monomorphic, or has "
               "the same genotype in every animal, as a panel of failed calls): there is no G");
  return f;
}

// Z nos marcadores cols[0..nb), centrada e com o ausente na media, por marcador:
// zb[c n + a]. Um marcador fora de uso sai com zeros.
static void bloco_z(const Genotipos& gt, const FreqZ& fz, const std::size_t* cols,
                    std::size_t nb, double* zb, int nth) {
  const std::size_t n = gt.n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nth) schedule(dynamic, 4)
#endif
  for (long cc = 0; cc < static_cast<long>(nb); cc++) {
    const std::size_t c = static_cast<std::size_t>(cc), j = cols[c];
    double* z = zb + c * n;
    if (!fz.usa[j]) { std::fill(z, z + n, 0.0); continue; }
    gt.coluna(j, z);
    const double med = fz.media[j], dp = 2.0 * fz.p[j];
    for (std::size_t a = 0; a < n; a++) z[a] = (std::isfinite(z[a]) ? z[a] : med) - dp;
  }
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
static void vanraden_em(const Genotipos& gt, const FreqZ& fz, double* c);

Densa vanraden_g(const Genotipos& gt, RelatorioG& rel, bool meio) {
  const FreqZ fz = frequencias_z(gt, rel, meio);
  Densa g(gt.n, gt.n);
  vanraden_em(gt, fz, g.dados.data());
  return g;
}

// A mesma G escrita no buffer de quem chama (a matriz do R de g_matrix(), sem copia n x n)
void vanraden_g_em(const Genotipos& gt, RelatorioG& rel, bool meio, double* out) {
  const FreqZ fz = frequencias_z(gt, rel, meio);
  vanraden_em(gt, fz, out);
}

// Z Z' (ou W W', a de dominancia) em blocos de marcadores, no buffer c do CHAMADOR (n x n):
// o codificador monta o bloco (n x B, coluna-major, um marcador por coluna) ja centrado e
// imputado, e o bloco e acumulado em C += Z_b Z_b', pelos ladrilhos em paralelo (padrao) ou
// pelo dsyrk do BLAS do R (br_threads(lapack = TRUE)). Os dois escrevem o mesmo triangulo:
// c[i + n k] com i <= k, o superior coluna-major, que e o inferior linha-major. Memoria
// extra: um bloco, nada proporcional a m, e nenhuma copia n x n: quem chama pode passar a
// propria matriz do R. A interrupcao e conferida entre os blocos, fora das regioes paralelas.
template <class Codifica>
static void gram_blocos(std::size_t n, const std::vector<std::size_t>& usados,
                        Codifica codifica, double* c) {
  const std::size_t B = 2048;
  const int ni = static_cast<int>(n);
  std::fill(c, c + n * n, 0.0);
  std::vector<double> zbuf(n * std::min(B, std::max<std::size_t>(usados.size(), 1)));
  double beta = 0.0;
  for (std::size_t u0 = 0; u0 < usados.size(); u0 += B) {
    checa_interrupcao();
    const std::size_t bcol = std::min(B, usados.size() - u0);
    codifica(usados.data() + u0, bcol, zbuf.data());
    if (denso_lapack()) {
      const int k = static_cast<int>(bcol);
      const double um = 1.0;
      F77_CALL(dsyrk)("U", "N", &ni, &k, &um, zbuf.data(), &ni, &beta, c, &ni
                      FCONE FCONE);
      beta = 1.0;
    } else {
      produto_ladrilhos(zbuf.data(), n, zbuf.data(), n, bcol, n, n, c, n, true, threads());
    }
  }
}

// O triangulo c[i + n k] (i <= k) de gram_blocos, dividido pela escala e espelhado no proprio
// buffer: sai a matriz cheia e simetrica (a mesma lida por linhas ou por colunas).
static void espelha_divide(double* c, std::size_t n, double escala) {
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t k2 = i; k2 < n; k2++) {
      const double v = c[i + n * k2] / escala;
      c[i + n * k2] = v;
      c[k2 + n * i] = v;
    }
}

static void vanraden_em(const Genotipos& gt, const FreqZ& fz, double* c) {
  std::vector<std::size_t> usados;
  for (std::size_t j = 0; j < gt.m; j++) if (fz.usa[j]) usados.push_back(j);
  gram_blocos(gt.n, usados,
              [&](const std::size_t* cols, std::size_t nb, double* zb) {
                bloco_z(gt, fz, cols, nb, zb, threads());
              }, c);
  espelha_divide(c, gt.n, fz.denom);
}

// ---- as portas de genotipos: F e D genomicos, Gamma, efeitos de SNP

// Heterozigotos e observados por animal (F genomico = 1 - het/obs): faixas de 1024 animais,
// cada faixa de uma thread, lendo o trecho dela de cada marcador. Contagens inteiras: exatas.
void homozigose(const Genotipos& gt, double* het, double* obs) {
  const std::size_t n = gt.n, m = gt.m, F = 1024;
  const long nf = static_cast<long>((n + F - 1) / F);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
  {
    std::vector<double> buf(F);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long ff = 0; ff < nf; ff++) {
      const std::size_t a0 = static_cast<std::size_t>(ff) * F, a1 = std::min(n, a0 + F);
      std::vector<std::size_t> h(a1 - a0, 0), o(a1 - a0, 0);
      for (std::size_t j = 0; j < m; j++) {
        gt.trecho(j, a0, a1, buf.data());
        for (std::size_t a = 0; a < a1 - a0; a++) {
          const double x = buf[a];
          if (!std::isfinite(x)) continue;
          o[a]++;
          if (x == 1.0) h[a]++;
        }
      }
      for (std::size_t a = a0; a < a1; a++) {
        het[a] = static_cast<double>(h[a - a0]);
        obs[a] = static_cast<double>(o[a - a0]);
      }
    }
  }
}

// As contagens de 0, 1 e 2 de cada marcador nos animais marcados em linha (todos com linha
// nula): o controle de qualidade tira dai a taxa de chamada do marcador, a frequencia, o MAF e
// o qui-quadrado de Hardy-Weinberg, e os ajustadores sabem quais marcadores entram na G. Cada
// marcador e de uma thread so; blocos de marcadores, com a interrupcao conferida entre eles.
// Contagens inteiras: o mesmo resultado com qualquer numero de threads.
void contagens_marcadores(const Genotipos& gt, const char* linha, double* n0, double* n1,
                          double* n2) {
  const std::size_t n = gt.n, m = gt.m, B = 1024;
  for (std::size_t j0 = 0; j0 < m; j0 += B) {
    checa_interrupcao();
    const long jb = static_cast<long>(j0), je = static_cast<long>(std::min(m, j0 + B));
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
    {
      std::vector<double> col(n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
      for (long jj = jb; jj < je; jj++) {
        const std::size_t j = static_cast<std::size_t>(jj);
        gt.coluna(j, col.data());
        std::size_t c0 = 0, c1 = 0, c2 = 0;
        for (std::size_t a = 0; a < n; a++) {
          if (linha && !linha[a]) continue;
          const double x = col[a];
          c0 += x == 0.0;
          c1 += x == 1.0;
          c2 += x == 2.0;
        }
        n0[j] = static_cast<double>(c0);
        n1[j] = static_cast<double>(c1);
        n2[j] = static_cast<double>(c2);
      }
    }
  }
}

// As mesmas contagens so nas linhas dadas (base 0): o que o controle de qualidade subtrai das
// contagens dos marcadores quando tira esses animais, lendo k x m genotipos e nao a matriz
// inteira a cada passada. Cada marcador e de uma thread so, e as linhas sao lidas na ordem
// dada; contagens inteiras, identicas com qualquer numero de threads. Com poucas linhas o
// laco fica numa thread so (abrir a regiao paralela custaria mais que a conta).
void contagens_nas_linhas(const Genotipos& gt, const std::vector<std::size_t>& linhas,
                          double* n0, double* n1, double* n2) {
  const std::size_t m = gt.m, k = linhas.size(), B = 4096;
  for (std::size_t a : linhas)
    if (a >= gt.n) throw Erro("marker counts: row outside the genotypes");
  for (std::size_t j0 = 0; j0 < m; j0 += B) {
    checa_interrupcao();
    const long jb = static_cast<long>(j0), je = static_cast<long>(std::min(m, j0 + B));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads()) if (k >= 64)
#endif
    for (long jj = jb; jj < je; jj++) {
      const std::size_t j = static_cast<std::size_t>(jj);
      std::size_t c0 = 0, c1 = 0, c2 = 0;
      for (std::size_t a : linhas) {
        const double x = gt.valor(j, a);
        c0 += x == 0.0;
        c1 += x == 1.0;
        c2 += x == 2.0;
      }
      n0[j] = static_cast<double>(c0);
      n1[j] = static_cast<double>(c1);
      n2[j] = static_cast<double>(c2);
    }
  }
}

// As chamadas de cada animal nos marcadores marcados em usa (a taxa de chamada por animal,
// preGSf90 callrateAnim e PLINK --mind; quem chama passa os marcadores que entram na G).
// Blocos de marcadores com a interrupcao entre eles; dentro, faixas de 1024 animais, cada
// faixa de uma thread so, dona das suas contagens. Contagens inteiras: identico com qualquer
// numero de threads. Um bloco sem marcador marcado nem abre a regiao paralela: o controle de
// qualidade passa aqui so as poucas colunas que sairam numa passada.
void chamadas_por_animal(const Genotipos& gt, const char* usa, double* chamadas) {
  const std::size_t n = gt.n, m = gt.m, F = 1024, B = 256;
  const long nf = static_cast<long>((n + F - 1) / F);
  std::vector<std::size_t> c(n, 0);
  for (std::size_t j0 = 0; j0 < m; j0 += B) {
    checa_interrupcao();
    const std::size_t j1 = std::min(m, j0 + B);
    if (usa && std::none_of(usa + j0, usa + j1, [](char u) { return u != 0; })) continue;
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
    {
      std::vector<double> buf(F);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (long ff = 0; ff < nf; ff++) {
        const std::size_t a0 = static_cast<std::size_t>(ff) * F, a1 = std::min(n, a0 + F);
        for (std::size_t j = j0; j < j1; j++) {
          if (usa && !usa[j]) continue;
          gt.trecho(j, a0, a1, buf.data());
          for (std::size_t a = 0; a < a1 - a0; a++) c[a0 + a] += std::isfinite(buf[a]);
        }
      }
    }
  }
  for (std::size_t a = 0; a < n; a++) chamadas[a] = static_cast<double>(c[a]);
}

// Homozigotos opostos entre pai e filho genotipados (o teste de exclusao de preGSf90,
// verify_parentage, e de Hayes 2011): num marcador chamado nos dois, filho 0 e pai 2 (ou o
// contrario) e impossivel por Mendel, fora erro de genotipagem. So os marcadores marcados em
// usa entram (todos com usa nulo): o cromossomo sexual fica de fora, como no preGSf90. O laco
// de fora e por blocos de marcadores, com a interrupcao conferida entre eles; dentro, cada
// thread e dona de uma faixa fixa de pares e le a coluna j nas linhas dos seus pares (a
// coluna inteira cabe no cache). As somas sem desvio (o += do teste) evitam o erro de
// previsao do salto quando metade dos pares conflita. Contagens inteiras, cada par de uma
// thread so: identico com qualquer numero de threads.
void conflitos_parentesco(const Genotipos& gt, const std::vector<std::size_t>& filho,
                          const std::vector<std::size_t>& pai, const char* usa, double* ambos,
                          double* opostos) {
  const std::size_t np = filho.size(), m = gt.m, B = 256;
  if (pai.size() != np) throw Erro("parentage: offspring and parent of different lengths");
  for (std::size_t k = 0; k < np; k++)
    if (filho[k] >= gt.n || pai[k] >= gt.n) throw Erro("parentage: row outside the genotypes");
  std::vector<std::size_t> a(np, 0), o(np, 0);
  for (std::size_t j0 = 0; j0 < m; j0 += B) {
    checa_interrupcao();
    const std::size_t j1 = std::min(m, j0 + B);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
    {
      std::size_t nt = 1, t = 0;
#ifdef _OPENMP
      nt = static_cast<std::size_t>(omp_get_num_threads());
      t = static_cast<std::size_t>(omp_get_thread_num());
#endif
      const std::size_t k0 = np * t / nt, k1 = np * (t + 1) / nt;
      for (std::size_t j = j0; j < j1; j++) {
        if (usa && !usa[j]) continue;
        for (std::size_t k = k0; k < k1; k++) {
          const double x = gt.valor(j, filho[k]), y = gt.valor(j, pai[k]);
          const bool ok = std::isfinite(x) && std::isfinite(y);
          a[k] += ok;
          o[k] += ok && std::fabs(x - y) == 2.0;
        }
      }
    }
  }
  for (std::size_t k = 0; k < np; k++) {
    ambos[k] = static_cast<double>(a[k]);
    opostos[k] = static_cast<double>(o[k]);
  }
}

// D de Vitezica et al. (2013): cada genotipo 0, 1, 2 codificado -2p^2, 2pq, -2q^2 e
// D = W W' / soma (2pq)^2, frequencias dos proprios genotipados. O ausente leva a media dos
// codigos observados do marcador, (n0 c0 + n1 c1 + n2 c2) / k. A media do marcador e a do
// motor (frequencias_z): soma/k em double, a mesma p da G; o R somava em long double, entao
// D difere da versao em R no ultimo bit quando k passa de ~2000 (portao por tolerancia).
struct FreqD {
  std::vector<double> c0, c1, c2, imp;
  std::vector<char> usa;
  double denom = 0.0;
};

static FreqD frequencias_d(const Genotipos& gt) {
  const std::size_t n = gt.n, nm = gt.m;
  // sem marcador nenhum cai na recusa de g_dominancia (nenhum polimorfico)
  if (n == 0) throw Erro("empty genotypes");
  FreqD f;
  f.c0.assign(nm, 0.0); f.c1.assign(nm, 0.0); f.c2.assign(nm, 0.0); f.imp.assign(nm, 0.0);
  f.usa.assign(nm, 0);
  std::vector<double> h(nm, 0.0);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
  {
    std::vector<double> col(n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
    for (long jj = 0; jj < static_cast<long>(nm); jj++) {
      const std::size_t j = static_cast<std::size_t>(jj);
      gt.coluna(j, col.data());
      std::size_t n0 = 0, n1 = 0, n2 = 0;
      for (std::size_t a = 0; a < n; a++) {
        const double x = col[a];
        if (!std::isfinite(x)) continue;
        if (x == 0.0) n0++;
        else if (x == 1.0) n1++;
        else n2++;
      }
      const std::size_t k = n0 + n1 + n2;
      if (k == 0) continue;
      // a soma inteira e exata: a mesma media de frequencias_z
      const double p = (static_cast<double>(n1 + 2 * n2) / static_cast<double>(k)) / 2.0;
      if (!(p > 0.0 && p < 1.0)) continue;
      const double q = 1.0 - p;
      f.usa[j] = 1;
      f.c0[j] = -2.0 * (p * p);
      f.c1[j] = (2.0 * p) * q;
      f.c2[j] = -2.0 * (q * q);
      f.imp[j] = (static_cast<double>(n0) * f.c0[j] + static_cast<double>(n1) * f.c1[j] +
                  static_cast<double>(n2) * f.c2[j]) / static_cast<double>(k);
      h[j] = f.c1[j] * f.c1[j];
    }
  }
  for (std::size_t j = 0; j < nm; j++)
    if (f.usa[j]) f.denom += h[j];
  return f;
}

// W nos marcadores cols[0..nb): o codigo de cada genotipo, ou a media dos codigos no ausente
static void bloco_w(const Genotipos& gt, const FreqD& fd, const std::size_t* cols,
                    std::size_t nb, double* wb) {
  const std::size_t n = gt.n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads()) schedule(dynamic, 4)
#endif
  for (long cc = 0; cc < static_cast<long>(nb); cc++) {
    const std::size_t c = static_cast<std::size_t>(cc), j = cols[c];
    double* w = wb + c * n;
    gt.coluna(j, w);
    for (std::size_t a = 0; a < n; a++) {
      const double x = w[a];
      w[a] = !std::isfinite(x) ? fd.imp[j]
           : x == 0.0 ? fd.c0[j] : x == 1.0 ? fd.c1[j] : fd.c2[j];
    }
  }
}

void g_dominancia(const Genotipos& gt, double* out) {
  const FreqD fd = frequencias_d(gt);
  std::vector<std::size_t> usados;
  for (std::size_t j = 0; j < gt.m; j++) if (fd.usa[j]) usados.push_back(j);
  if (usados.empty()) throw Erro("every marker is monomorphic: there is no D to build");
  const std::size_t n = gt.n;
  gram_blocos(n, usados,
              [&](const std::size_t* cols, std::size_t nb, double* wb) {
                bloco_w(gt, fd, cols, nb, wb);
              }, out);
  espelha_divide(out, n, fd.denom);
}

// T = (X - c)' V, m x k: out[j + b m] = soma em a crescente de (x_aj - c_j) v[a + b n], com o
// ausente valendo imp_j, isto e, contribuindo (imp_j - c_j) v (nada quando imp_j = c_j). Cada
// marcador de uma thread e um acumulador so em ordem de animal, a ordem do dgemm de referencia
// ('T', 'N'): o crossprod(Z, V) do R sem formar Z, e o mesmo resultado com qualquer numero
// de threads. Em fatias de 64 marcadores por thread (4 pedacos de 16 por thread, para todas
// trabalharem com qualquer numero delas), com a interrupcao conferida entre as fatias (a fatia
// nao muda conta nenhuma: cada marcador continua de um dono so).
void genotipos_tv(const Genotipos& gt, const double* centro, const double* imp,
                  const double* v, std::size_t k, double* out) {
  const std::size_t n = gt.n, m = gt.m;
  const std::size_t fatia = 64 * static_cast<std::size_t>(std::max(1, threads()));
  for (std::size_t f0 = 0; f0 < m; f0 += fatia) {
    checa_interrupcao();
    const long f1 = static_cast<long>(std::min(m, f0 + fatia));
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
    {
      std::vector<double> col(n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
      for (long jj = static_cast<long>(f0); jj < f1; jj++) {
        const std::size_t j = static_cast<std::size_t>(jj);
        gt.coluna(j, col.data());
        const double cj = centro[j], dz = imp[j] - centro[j];
        for (std::size_t b = 0; b < k; b++) {
          const double* vb = v + b * n;
          double t = 0.0;
          for (std::size_t a = 0; a < n; a++) {
            const double x = col[a];
            if (std::isfinite(x)) t += (x - cj) * vb[a];
            else if (dz != 0.0) t += dz * vb[a];
          }
          out[j + b * m] = t;
        }
      }
    }
  }
}

// Quantos ausentes (NaN, NA_INTEGER, 5 no raw): a recusa de estimate_gamma() sem converter
double conta_ausentes(const Genotipos& gt) {
  double soma = 0.0;
  std::vector<double> col(gt.n);
  for (std::size_t j = 0; j < gt.m; j++) {
    gt.coluna(j, col.data());
    for (double x : col) if (!std::isfinite(x)) soma += 1.0;
  }
  return soma;
}

// A matriz de Gram do lado MENOR de Z (centrada, imputada, so os marcadores em uso, sem dividir
// pela escala): Z Z' (n x n) com n <= m, ou Z'Z (m x m) no outro caso, que e o que o eigen
// exato de apy_core_select() decompoe. Z'Z vai por blocos de ANIMAIS, cada bloco lido do R
// marcador a marcador, sem nunca formar Z inteira.
Densa gram_menor(const Genotipos& gt, RelatorioG& rel, double& escala) {
  const std::size_t n = gt.n;
  const FreqZ fz = frequencias_z(gt, rel, false);
  escala = fz.denom;
  std::vector<std::size_t> usados;
  for (std::size_t j = 0; j < gt.m; j++) if (fz.usa[j]) usados.push_back(j);
  const std::size_t mu = usados.size();
  const int nth = threads();
  const std::size_t d = std::min(n, mu);
  std::vector<double> c(d * d, 0.0);
  if (n <= mu) {
    const std::size_t B = 2048;
    std::vector<double> zb(n * std::min(B, mu));
    for (std::size_t u0 = 0; u0 < mu; u0 += B) {
      const std::size_t nb = std::min(B, mu - u0);
      bloco_z(gt, fz, usados.data() + u0, nb, zb.data(), nth);
      produto_ladrilhos(zb.data(), n, zb.data(), n, nb, n, n, c.data(), n, true, nth);
    }
  } else {
    const std::size_t A = 1024;
    std::vector<double> zr(std::min(A, n) * mu);
    for (std::size_t a0 = 0; a0 < n; a0 += A) {
      const std::size_t na = std::min(A, n - a0);
#ifdef _OPENMP
#pragma omp parallel num_threads(nth)
#endif
      {
        std::vector<double> buf(na);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
        for (long uu = 0; uu < static_cast<long>(mu); uu++) {
          const std::size_t u = static_cast<std::size_t>(uu), j = usados[u];
          gt.trecho(j, a0, a0 + na, buf.data());
          const double med = fz.media[j], dp = 2.0 * fz.p[j];
          for (std::size_t a = 0; a < na; a++)
            zr[a * mu + u] = (std::isfinite(buf[a]) ? buf[a] : med) - dp;
        }
      }
      produto_ladrilhos(zr.data(), mu, zr.data(), mu, na, mu, mu, c.data(), mu, true, nth);
    }
  }
  Densa g(d, d);
  for (std::size_t i = 0; i < d; i++)
    for (std::size_t k = 0; k <= i; k++) {
      g.at(i, k) = c[i * d + k];
      g.at(k, i) = c[i * d + k];
    }
  return g;
}

// Lanczos sobre G = ZZ'/k para a quadratura de Lanczos estocastica (Ubaru, Chen e Saad,
// 2017, SIAM J. Matrix Anal. Appl. 38:1075): com sondas de Rademacher v normalizadas, os nos e
// os pesos de Gauss da tridiagonal de cada sonda estimam a medida espectral de G, e dela sai
// quantos autovalores explicam 98% do traco sem decompor G. Os autovalores nao nulos de ZZ'
// e Z'Z sao os mesmos, entao o operador e o do lado menor. Z nunca existe inteira: cada
// passo le os genotipos do R por blocos de marcadores, centra e imputa o bloco, e os dois
// produtos do bloco vao aos ladrilhos com todas as sondas de uma vez. Sem reortogonalizacao: a
// quadratura continua certa na perda de ortogonalidade (Knizhnerman 1996; e o que o artigo
// usa), e o portao mede isso contra o eigen exato. O erro sistematico que sobra tambem nao vem
// dos fantasmas: e o da regra de Gauss no limiar (os nos de Ritz caem quase no mesmo lugar em
// todas as sondas), que reortogonalizar nao tira; o R o tira com a media das regras das
// ultimas tridiagonais lideres de cada sonda (nos_lanczos, em R/apy.R), que ja estao aqui.
LanczosG lanczos_g(const Genotipos& gt, const Densa& sondas, std::size_t passos,
                   RelatorioG& rel) {
  const std::size_t n = gt.n, nm = gt.m;
  const FreqZ fz = frequencias_z(gt, rel, false);
  const int nth = threads();
  const std::size_t B = std::min<std::size_t>(1024, nm);
  std::vector<std::size_t> todos(nm);
  for (std::size_t j = 0; j < nm; j++) todos[j] = j;
  std::vector<double> zb(n * B);
  double traco = 0.0;
  for (std::size_t j0 = 0; j0 < nm; j0 += B) {
    const std::size_t nb = std::min(B, nm - j0);
    bloco_z(gt, fz, todos.data() + j0, nb, zb.data(), nth);
    for (std::size_t k = 0; k < nb * n; k++) traco += zb[k] * zb[k];
  }
  const bool lado_n = n <= nm;
  const std::size_t dim = lado_n ? n : nm, nv = sondas.ncol;
  if (sondas.nlin != dim) throw Erro("probes with the wrong number of rows");
  if (passos < 2) throw Erro("Lanczos needs at least 2 steps");
  LanczosG r;
  r.dim = dim;
  r.traco = traco / fz.denom;
  r.alfa = Densa(passos, nv);
  r.beta = Densa(passos, nv);
  r.passos.assign(nv, passos);
  std::vector<char> viva(nv, 1);

  // w = A q, A = ZZ'/k (lado n) ou Z'Z/k (lado m), bloco a bloco de marcadores
  std::vector<double> w(dim * nv), wb(B * nv), t(lado_n ? 0 : n * nv);
  auto aplica = [&](const std::vector<double>& q) {
    std::fill(w.begin(), w.end(), 0.0);
    if (lado_n) {
      for (std::size_t j0 = 0; j0 < nm; j0 += B) {
        const std::size_t nb = std::min(B, nm - j0);
        bloco_z(gt, fz, todos.data() + j0, nb, zb.data(), nth);
        std::fill(wb.begin(), wb.begin() + nb * nv, 0.0);
        produto_ladrilhos(zb.data(), n, q.data(), nv, n, nb, nv, wb.data(), nv, false, nth,
                          true);
        produto_ladrilhos(zb.data(), n, wb.data(), nv, nb, n, nv, w.data(), nv, false, nth);
      }
    } else {
      std::fill(t.begin(), t.end(), 0.0);
      for (std::size_t j0 = 0; j0 < nm; j0 += B) {
        const std::size_t nb = std::min(B, nm - j0);
        bloco_z(gt, fz, todos.data() + j0, nb, zb.data(), nth);
        produto_ladrilhos(zb.data(), n, q.data() + j0 * nv, nv, nb, n, nv, t.data(), nv, false,
                          nth);
      }
      for (std::size_t j0 = 0; j0 < nm; j0 += B) {
        const std::size_t nb = std::min(B, nm - j0);
        bloco_z(gt, fz, todos.data() + j0, nb, zb.data(), nth);
        produto_ladrilhos(zb.data(), n, t.data(), nv, n, nb, nv, w.data() + j0 * nv, nv, false,
                          nth, true);
      }
    }
    for (double& x : w) x /= fz.denom;
  };

// q normalizada por coluna; q0 a anterior
  std::vector<double> q(dim * nv), q0(dim * nv, 0.0);
  for (std::size_t b = 0; b < nv; b++) {
    double s2 = 0.0;
    for (std::size_t i = 0; i < dim; i++) s2 += sondas.at(i, b) * sondas.at(i, b);
    if (!(s2 > 0.0)) throw Erro("a probe is zero");
    const double s = std::sqrt(s2);
    for (std::size_t i = 0; i < dim; i++) q[i * nv + b] = sondas.at(i, b) / s;
  }
  for (std::size_t j = 0; j < passos; j++) {
    R_CheckUserInterrupt();
    aplica(q);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nth) schedule(static)
#endif
    for (long bb = 0; bb < static_cast<long>(nv); bb++) {
      const std::size_t b = static_cast<std::size_t>(bb);
      if (!viva[b]) continue;
      double a = 0.0;
      for (std::size_t i = 0; i < dim; i++) a += w[i * nv + b] * q[i * nv + b];
      const double bprev = j > 0 ? r.beta.at(j - 1, b) : 0.0;
      double s2 = 0.0;
      for (std::size_t i = 0; i < dim; i++) {
        const double x = w[i * nv + b] - a * q[i * nv + b] - bprev * q0[i * nv + b];
        w[i * nv + b] = x;
        s2 += x * x;
      }
      r.alfa.at(j, b) = a;
      const double bt = std::sqrt(s2);
      r.beta.at(j, b) = bt;
      // subespaco invariante: a tridiagonal desta sonda para aqui, exata
      if (!(bt > 1e-12 * std::fabs(a)) || j + 1 == passos) {
        r.passos[b] = j + 1;
        viva[b] = 0;
        continue;
      }
      for (std::size_t i = 0; i < dim; i++) {
        q0[i * nv + b] = q[i * nv + b];
        q[i * nv + b] = w[i * nv + b] / bt;
      }
    }
    if (std::none_of(viva.begin(), viva.end(), [](char c) { return c != 0; })) break;
  }
  return r;
}

// A22^-1 pelo Schur ESPARSO sobre A^-1, sem nunca formar B11^-1 denso, devolvido como o
// triangulo INFERIOR (na numeracao dos genotipados) so com as entradas que existem: o A22^-1
// e quase todo zero EXATO (medido: 97,6% em 400 genotipados, 98,9% em 800), e e isso que
// deixa a APY esparsa ate o H^-1. B12 tem so os pais e filhos nao genotipados de cada
// genotipado, meia duzia por coluna; guardado denso ele custava memoria n1 n2 e o produto
// B21 x (B11^-1 B12) virava n1 n2^2 flops sobre zeros. Uma resolucao por coluna genotipada,
// as colunas em paralelo (cada uma de uma thread, com a mesma ordem de soma).
// O preparo do Schur, comum a A22^-1 explicita e as resolucoes sem formar A22^-1 (as de
// estimate_gamma(): A22^-1 X de poucas colunas e a forma quadratica dos marcadores do ML, em
// que A22^-1 pode ser 100% densa quando so a ultima geracao tem genotipo). Numeracao:
// nao-genotipados primeiro; B11 fatorado (grau minimo), B12 esparso nas duas direcoes, B22
// esparso por coluna (triangulo inferior, a diagonal inclusive).
struct SchurA22 {
  std::size_t n1 = 0, n2 = 0;
  Csc L;
  std::vector<std::size_t> onde_perm;   // linha nao genotipada -> posicao no fator
  std::vector<std::vector<std::pair<std::uint32_t, double>>> b12_col, b12_lin, b22_col;
};

static SchurA22 prepara_schur(const Csc& ainv, const std::vector<std::size_t>& geno,
                              CacheSimbolica& (*cache)(const Csc&) = nullptr) {
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

  SchurA22 s;
  s.n1 = n1;
  s.n2 = n2;
  // B11 (triangulo inferior esparso), B12 esparso nas duas direcoes, B22 esparso por coluna
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  s.b12_col.resize(n2);
  s.b22_col.resize(n2);
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
        s.b22_col[j].push_back({static_cast<std::uint32_t>(i), x});
      } else {
        const std::size_t inng = rg ? pc : pr;
        const std::size_t ig = rg ? pr - n1 : pc - n1;
        s.b12_col[ig].push_back({static_cast<std::uint32_t>(inng), x});
      }
    }
  s.b12_lin.resize(n1);
  for (std::size_t ig = 0; ig < n2; ig++)
    for (const auto& [i, x] : s.b12_col[ig])
      s.b12_lin[i].push_back({static_cast<std::uint32_t>(ig), x});

  std::vector<std::size_t> perm;
  if (n1 > 0) {
    Csc b11 = de_triplos(n1, n1, li, cj, v);
    bool ok;
    if (cache) {
      // o mesmo padrao a cada iteracao do pseudo-EM (so Gamma muda os valores): o grau minimo
      // e a simbolica guardados pela assinatura do padrao, o mesmo resultado bit a bit
      CacheSimbolica& cs = cache(b11);
      perm = cs.perm;
      ok = cholesky(permuta_cache(b11, cs), cs.sb, s.L);
    } else {
      perm = grau_minimo(b11);
      Csc pb = permuta_sim(b11, perm);
      Simbolica sb = simbolica(pb);
      ok = cholesky(pb, sb, s.L);
    }
    if (!ok) throw Erro("the non-genotyped block of A^-1 is not positive-definite");
  }
  s.onde_perm.assign(n1, 0);
  for (std::size_t i = 0; i < n1; i++) s.onde_perm[perm[i]] = i;
  return s;
}

// L y = b no lugar, W colunas por vez (X linha-major, n1 x w, ja permutado)
// (interrompe: confere a interrupcao a cada 65536 linhas; so para quem chama fora de regiao
// paralela, porque a excecao nao pode sair de uma)
static void schur_frente(const Csc& L, std::size_t n1, std::size_t w, double* X,
                         bool interrompe = false) {
  for (std::size_t j = 0; j < n1; j++) {
    if (interrompe && (j & 0xFFFF) == 0) checa_interrupcao();
    double* xj = &X[j * w];
    const double d = L.valor[L.colptr[j]];
    for (std::size_t c = 0; c < w; c++) xj[c] /= d;
    for (std::size_t p = L.colptr[j] + 1; p < L.colptr[j + 1]; p++) {
      const double lp = L.valor[p];
      double* xr = &X[static_cast<std::size_t>(L.linha[p]) * w];
      for (std::size_t c = 0; c < w; c++) xr[c] -= lp * xj[c];
    }
  }
}

// L' x = y no lugar, no mesmo arranjo
static void schur_tras(const Csc& L, std::size_t n1, std::size_t w, double* X,
                       bool interrompe = false) {
  for (std::size_t j = n1; j-- > 0;) {
    if (interrompe && (j & 0xFFFF) == 0) checa_interrupcao();
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

// Y = P B12 X para as colunas c0..c0+w de X (n2 x q coluna-major), em Y (n1 x w linha-major,
// na ordem do fator): cada linha nao genotipada soma as suas entradas de B12 em ordem fixa
static void schur_b12(const SchurA22& s, const double* X, std::size_t c0, std::size_t w,
                      double* Y) {
  for (std::size_t i = 0; i < s.n1; i++) {
    double* yi = &Y[s.onde_perm[i] * w];
    for (std::size_t c = 0; c < w; c++) {
      const double* xc = X + (c0 + c) * s.n2;
      double t = 0.0;
      for (const auto& [ig, x] : s.b12_lin[i]) t += x * xc[ig];
      yi[c] = t;
    }
  }
}

// A22^-1 X = B22 X - B21 B11^-1 B12 X, X (n2 x q coluna-major) em out (o mesmo arranjo), sem
// formar A22^-1: uma fatoracao para as q colunas, em blocos de W colunas, cada bloco de uma
// thread e cada coluna de saida de um bloco so, com as somas em ordem fixa.
static void schur_aplica(const SchurA22& s, const double* X, std::size_t q, double* out) {
  const std::size_t W = 32, n1 = s.n1, n2 = s.n2;
  const std::size_t nblocos = (q + W - 1) / W;
  const auto bloco = [&](std::size_t bb, double* Y, bool interrompe) {
    const std::size_t c0 = bb * W, w = std::min(W, q - c0);
    if (n1 > 0) {
      schur_b12(s, X, c0, w, Y);
      schur_frente(s.L, n1, w, Y, interrompe);
      schur_tras(s.L, n1, w, Y, interrompe);
    }
    for (std::size_t c = 0; c < w; c++) {
      const double* xc = X + (c0 + c) * n2;
      double* o = out + (c0 + c) * n2;
      std::fill(o, o + n2, 0.0);
      for (std::size_t jg = 0; jg < n2; jg++)
        for (const auto& [ig, x] : s.b22_col[jg]) {
          o[ig] += x * xc[jg];
          if (ig != jg) o[jg] += x * xc[ig];
        }
      if (n1 > 0)
        for (std::size_t ig = 0; ig < n2; ig++)
          for (const auto& [i, x] : s.b12_col[ig]) o[ig] -= x * Y[s.onde_perm[i] * w + c];
    }
  };
  // um bloco so (as q <= 32 colunas do pseudo-EM, do GLS e do 1'A22^-1 do ML): fora de regiao
  // paralela, que nao teria o que dividir, e com a interrupcao conferida dentro das
  // resolucoes, que num pedigree fundo levam segundos. As mesmas contas na mesma ordem.
  if (nblocos == 1) {
    std::vector<double> Y(n1 * W);
    bloco(0, Y.data(), true);
    return;
  }
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
  {
    std::vector<double> Y(n1 * W);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long bb = 0; bb < static_cast<long>(nblocos); bb++)
      bloco(static_cast<std::size_t>(bb), Y.data(), false);
  }
}

Csc a22_inversa_esparsa(const Csc& ainv, const std::vector<std::size_t>& geno) {
  const SchurA22 s = prepara_schur(ainv, geno);
  const std::size_t n1 = s.n1, n2 = s.n2;
  const Csc& L = s.L;
  const auto& b12_col = s.b12_col;
  const auto& b12_lin = s.b12_lin;
  const auto& b22_col = s.b22_col;
  const auto& onde_perm = s.onde_perm;

  // coluna jg do triangulo inferior: B22(ig, jg) - B21 B11^-1 B12(ig, jg), ig >= jg. As
  // resolucoes vao em BLOCOS de W colunas: cada passada pelo fator serve as W de uma vez (a
  // leitura de L, que domina, e dividida por W), com as mesmas operacoes na mesma ordem de
  // uma resolucao sozinha, coluna a coluna. Os blocos em paralelo, cada um de uma thread.
  const std::size_t W = 32;
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
        schur_frente(L, n1, w, X.data());
        schur_tras(L, n1, w, X.data());
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

// A22^-1 X sem formar A22^-1 (estimate_gamma(): o V do pseudo-EM e do GLS), X n2 x q
// coluna-major na ordem de geno
std::vector<double> a22_inversa_vezes(const Csc& ainv, const std::vector<std::size_t>& geno,
                                      const std::vector<double>& X, std::size_t q,
                                      CacheSimbolica& (*cache)(const Csc&)) {
  if (X.size() != geno.size() * q) throw Erro("A22^-1 X: X with the wrong size");
  const SchurA22 s = prepara_schur(ainv, geno, cache);
  std::vector<double> out(X.size(), 0.0);
  schur_aplica(s, X.data(), q, out.data());
  return out;
}

// Os tres resumos do ML de um metafundador (Legarra et al. 2024b), sem A22^-1 densa e sem
// Z Z': s = 1'A22^-1 1, t = soma_j z_j' A22^-1 z_j e r = soma_j (z_j' A22^-1 1)^2, com
// z_j = m_j - 1 (o ausente nao chega: estimate_gamma recusa antes). O A22^-1 pode ser 100%
// denso (genotipados so na ultima geracao), entao a forma quadratica sai do Schur por bloco
// de W marcadores: z'A22^-1 z = z'B22 z - ||L^-1 P B12 z||^2, so a substituicao para a frente
// no fator de B11. Cada marcador de um bloco, cada bloco de uma thread, e as parciais por
// marcador somadas em ordem de marcador: o mesmo resultado com qualquer numero de threads.
ResumosMl gamma_ml_resumos(const Csc& ainv, const std::vector<std::size_t>& geno,
                           const Genotipos& gt) {
  const SchurA22 s = prepara_schur(ainv, geno);
  const std::size_t n1 = s.n1, n2 = s.n2, m = gt.m;
  if (gt.n != n2) throw Erro("genotypes and genotyped animals with different counts");
  ResumosMl r;
  std::vector<double> um(n2, 1.0), w(n2, 0.0);
  schur_aplica(s, um.data(), 1, w.data());
  for (std::size_t i = 0; i < n2; i++) r.s += w[i];
  std::vector<double> centro(m, 1.0), t(m, 0.0);
  genotipos_tv(gt, centro.data(), centro.data(), w.data(), 1, t.data());
  for (std::size_t j = 0; j < m; j++) r.r += t[j] * t[j];

  // em fatias de 8 blocos por thread, com a interrupcao conferida entre as fatias (fora da
  // regiao paralela): num pedigree fundo cada bloco custa uma passada inteira pelo fator. A
  // fatia nao muda conta nenhuma (cada marcador de um bloco so, somados em ordem no fim).
  const std::size_t W = 32;
  const std::size_t nblocos = (m + W - 1) / W;
  const std::size_t fatia = 8 * static_cast<std::size_t>(std::max(1, threads()));
  std::vector<double> parc(m, 0.0);
  for (std::size_t f0 = 0; f0 < nblocos; f0 += fatia) {
    checa_interrupcao();
    const long f1 = static_cast<long>(std::min(nblocos, f0 + fatia));
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
    {
      std::vector<double> Z(n2 * W), Y(n1 * W);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (long bb = static_cast<long>(f0); bb < f1; bb++) {
        const std::size_t j0 = static_cast<std::size_t>(bb) * W, w2 = std::min(W, m - j0);
        for (std::size_t c = 0; c < w2; c++) {
          double* z = &Z[c * n2];
          gt.coluna(j0 + c, z);
          for (std::size_t a = 0; a < n2; a++) z[a] = std::isfinite(z[a]) ? z[a] - 1.0 : 0.0;
        }
        if (n1 > 0) {
          schur_b12(s, Z.data(), 0, w2, Y.data());
          schur_frente(s.L, n1, w2, Y.data());
        }
        for (std::size_t c = 0; c < w2; c++) {
          const double* z = &Z[c * n2];
          double q22 = 0.0;
          for (std::size_t jg = 0; jg < n2; jg++)
            for (const auto& [ig, x] : s.b22_col[jg])
              q22 += (ig == jg ? 1.0 : 2.0) * x * z[ig] * z[jg];
          double y2 = 0.0;
          for (std::size_t i = 0; i < n1; i++) y2 += Y[i * w2 + c] * Y[i * w2 + c];
          parc[j0 + c] = q22 - y2;
        }
      }
    }
  }
  for (std::size_t j = 0; j < m; j++) r.t += parc[j];
  return r;
}

// Fracoes de genes dos metafundadores em cada animal, linhas em ORDEM TOPOLOGICA (pais antes
// dos filhos): q_mf = e_k na linha do metafundador k, q_i = (q_pai + q_mae) / 2 nas outras. O
// laco por linha do pedigree que estimate_gamma() fazia em R, com a mesma conta: Q n x q por
// linhas.
std::vector<double> fracoes_genes(const std::vector<std::int64_t>& pai,
                                  const std::vector<std::int64_t>& mae,
                                  const std::vector<std::int64_t>& col_mf, std::size_t q) {
  const std::size_t n = pai.size();
  if (mae.size() != n || col_mf.size() != n) throw Erro("gene fractions: columns of different lengths");
  std::vector<double> Q(n * q, 0.0);
  for (std::size_t i = 0; i < n; i++) {
    double* qi = &Q[i * q];
    if (col_mf[i] >= 0) {
      if (static_cast<std::size_t>(col_mf[i]) >= q) throw Erro("metafounder column outside Gamma");
      qi[col_mf[i]] = 1.0;
      continue;
    }
    for (const std::int64_t pp : {pai[i], mae[i]})
      if (pp >= static_cast<std::int64_t>(i))
        throw Erro("gene fractions: the pedigree is not in topological order");
    for (std::size_t k = 0; k < q; k++) {
      const double s = pai[i] >= 0 ? Q[static_cast<std::size_t>(pai[i]) * q + k] : 0.0;
      const double d = mae[i] >= 0 ? Q[static_cast<std::size_t>(mae[i]) * q + k] : 0.0;
      qi[k] = (s + d) / 2.0;
    }
  }
  return Q;
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

Densa a22_colleau(const Pedigree& p, const std::vector<std::size_t>& geno) {
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
// Os coeficientes do afim a + b G que levam a media da diagonal e a de fora dela aos de A22
static void coef_afim(const Densa& g, const Densa& a22, double& a, double& b) {
  const std::size_t n = g.nlin;
  double mg_diag = 0.0, ma_diag = 0.0, mg_off = 0.0, ma_off = 0.0;
  for (std::size_t i = 0; i < n; i++) {
    mg_diag += g.at(i, i);
    ma_diag += a22.at(i, i);
    for (std::size_t j = 0; j < i; j++) { mg_off += g.at(i, j); ma_off += a22.at(i, j); }
  }
  mg_diag /= n; ma_diag /= n;
  const double noff = n > 1 ? static_cast<double>(n) * (n - 1) / 2.0 : 1.0;
  mg_off /= noff; ma_off /= noff;
  b = (mg_diag - mg_off) != 0.0 ? (ma_diag - ma_off) / (mg_diag - mg_off) : 1.0;
  a = ma_off - b * mg_off;
}

Densa ajusta_g_para_a22(const Densa& g, const Densa& a22, double mistura) {
  const std::size_t n = g.nlin;
  if (a22.nlin != n) throw Erro("G and A22 with different sizes");
  if (!(mistura >= 0.0 && mistura <= 1.0)) throw Erro("blend outside [0, 1]");
  // b (mg - mg_off) + a = ma  na diagonal e fora: resolve o afim a + b G
  double a = 0.0, b = 1.0;
  coef_afim(g, a22, a, b);
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
// Os passos 1 a 4 da rota APY (G* nas linhas do nucleo, Gcc^-1, P = Gcc^-1 Gcn e os residuos
// mendelianos dos jovens), comuns ao H^-1 e ao retro-calculo dos efeitos de SNP, que aplica a
// MESMA inversa APY ao vetor de valores geneticos sem formar nada n x n.
struct ApyGstar {
  FreqZ fz;
  Densa gcci, pm;                    // c x c e c x n
  std::vector<double> minv, gsd;     // 1/residuo mendeliano dos jovens; diagonal de G*
  std::vector<std::size_t> jov;
  double a0 = 0.0, b = 1.0;
};

static ApyGstar apy_gstar(const Pedigree& ped, const std::vector<std::size_t>& idx,
                          const Genotipos& gt, double mistura,
                          const std::vector<std::size_t>& nuc, bool com_mf, RelatorioG& rel) {
  const std::size_t n = gt.n, nm = gt.m, c = nuc.size();
  if (c < 2) throw Erro("the APY core needs at least 2 animals");
  if (!(mistura >= 0.0 && mistura <= 1.0)) throw Erro("blend outside [0, 1]");
  std::vector<char> eh_nuc(n, 0);
  for (std::size_t a : nuc) {
    if (a >= n) throw Erro("APY core index outside the genotyped set");
    eh_nuc[a] = 1;
  }
  ApyGstar ap;
  std::vector<std::size_t>& jov = ap.jov;
  for (std::size_t i = 0; i < n; i++) if (!eh_nuc[i]) jov.push_back(i);
  const std::size_t nj = jov.size();
  const int nth = threads();
  ap.fz = frequencias_z(gt, rel, com_mf);
  const FreqZ& fz = ap.fz;

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
    const double dp = 2.0 * fz.p[j], med = fz.media[j];
    double* zk = &zb[bcol * n];
    gt.coluna(j, zk);
    double sj = 0.0;
    for (std::size_t i = 0; i < n; i++) {
      const double z = (std::isfinite(zk[i]) ? zk[i] : med) - dp;
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
  double& a0 = ap.a0;
  double& b = ap.b;
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
  std::vector<double>& gsd = ap.gsd;
  gsd.assign(n, 0.0);
  for (std::size_t i = 0; i < n; i++)
    gsd[i] = (1.0 - mistura) * (a0 + b * gdiag[i]) + mistura * (1.0 + cl.f[idx[i]]);
  rel.diag_gstar = gsd;

  // 4. APY: Gcc^-1, P = Gcc^-1 Gcn, residuos mendelianos, e os blocos da inversa
  Densa gcc(c, c);
  for (std::size_t a = 0; a < c; a++)
    for (std::size_t bb = 0; bb < c; bb++)
      gcc.at(a, bb) = 0.5 * (gs.at(a, nuc[bb]) + gs.at(bb, nuc[a]));
  ap.gcci = inv_pd(gcc);
  const Densa& gcci = ap.gcci;
  ap.pm = Densa(c, n);
  Densa& pm = ap.pm;
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
  std::vector<double>& minv = ap.minv;
  minv.assign(nj, 0.0);
  for (std::size_t k = 0; k < nj; k++) {
    if (mend[k] <= limiar)
      throw Erro("degenerate Mendelian residual in APY: a non-core animal is collinear with the core "
                 "(clone or duplicate). Enlarge the core or the blend.");
    minv[k] = 1.0 / mend[k];
  }
  return ap;
}

// G*^-1 u pela inversa APY, sem forma-la: com t = M^-1 (u_n - P' u_c), o bloco do nucleo e
// Gcc^-1 u_c - P t e o dos jovens e t, que e [Gcc^-1 + P M^-1 P', -P M^-1; -M^-1 P', M^-1] u.
static std::vector<double> apy_aplica(const ApyGstar& ap, const std::vector<std::size_t>& nuc,
                                      const std::vector<double>& u) {
  const std::size_t c = nuc.size(), nj = ap.jov.size();
  std::vector<double> uc(c), t(nj), out(u.size(), 0.0);
  for (std::size_t a = 0; a < c; a++) uc[a] = u[nuc[a]];
  for (std::size_t k = 0; k < nj; k++) {
    const std::size_t i = ap.jov[k];
    double s = u[i];
    for (std::size_t a = 0; a < c; a++) s -= ap.pm.at(a, i) * uc[a];
    t[k] = s * ap.minv[k];
    out[i] = t[k];
  }
  for (std::size_t a = 0; a < c; a++) {
    const double* ga = ap.gcci.linha(a);
    const double* pa = ap.pm.linha(a);
    double s = 0.0;
    for (std::size_t b2 = 0; b2 < c; b2++) s += ga[b2] * uc[b2];
    for (std::size_t k = 0; k < nj; k++) s -= pa[ap.jov[k]] * t[k];
    out[nuc[a]] = s;
  }
  return out;
}

static Csc h_inversa_apy(const Pedigree& ped, const Csc& ainv, const std::vector<std::size_t>& idx,
                         const Genotipos& gt, double mistura, const std::vector<std::size_t>& nuc,
                         bool com_mf, RelatorioG& rel) {
  const ApyGstar ap = apy_gstar(ped, idx, gt, mistura, nuc, com_mf, rel);
  const std::size_t c = nuc.size(), nj = ap.jov.size();
  const int nth = threads();
  const Densa& pm = ap.pm;
  const std::vector<double>& minv = ap.minv;
  const std::vector<std::size_t>& jov = ap.jov;
  // nucleo x nucleo: Gcc^-1 + P M^-1 P', pelos ladrilhos com Ps' (nj x c) "por k"
  Densa cc = ap.gcci;
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
              const Genotipos& gt, double mistura, const std::vector<std::string>& nucleo_apy,
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
        throw Erro("genotyped animal '" + g + "' is not in the pedigree" +
                   dica_cientifica(g, ped.ids));
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
    return h_inversa_apy(ped, ainv, idx, gt, mistura, nc_idx, com_mf, rel);
  }
  Densa g = vanraden_g(gt, rel, com_mf);
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

// A linha do pedigree de cada genotipado (a mesma busca e a mesma mensagem de h_inversa)
std::vector<std::size_t> linhas_genotipados(const Pedigree& ped,
                                            const std::vector<std::string>& geno_ids) {
  std::unordered_map<std::string, std::size_t> onde;
  onde.reserve(ped.ids.size());
  for (std::size_t i = 0; i < ped.ids.size(); i++) onde.emplace(ped.ids[i], i);
  std::vector<std::size_t> idx;
  idx.reserve(geno_ids.size());
  for (const std::string& g : geno_ids) {
    const auto it = onde.find(g);
    if (it == onde.end())
      throw Erro("genotyped animal '" + g + "' is not in the pedigree" +
                 dica_cientifica(g, ped.ids));
    idx.push_back(it->second);
  }
  return idx;
}

// A diagonal de G sem formar G: G_ii = soma_j z_ij^2 / k, por faixas de 1024 animais (cada
// faixa de uma thread, os marcadores em ordem), com o ausente na media do marcador.
static std::vector<double> diag_g_faixas(const Genotipos& gt, const FreqZ& fz) {
  const std::size_t n = gt.n, m = gt.m, F = 1024;
  const long nf = static_cast<long>((n + F - 1) / F);
  std::vector<double> d(n, 0.0);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
  {
    std::vector<double> buf(F);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long ff = 0; ff < nf; ff++) {
      const std::size_t a0 = static_cast<std::size_t>(ff) * F, a1 = std::min(n, a0 + F);
      double* da = d.data() + a0;
      for (std::size_t j = 0; j < m; j++) {
        if (!fz.usa[j]) continue;
        gt.trecho(j, a0, a1, buf.data());
        const double med = fz.media[j], dp = 2.0 * fz.p[j];
        for (std::size_t a = 0; a < a1 - a0; a++) {
          const double z = (std::isfinite(buf[a]) ? buf[a] : med) - dp;
          da[a] += z * z;
        }
      }
    }
  }
  for (std::size_t i = 0; i < n; i++) d[i] /= fz.denom;
  return d;
}

// A diagonal de G* sem formar G nem A22, O(n m + n_ped): G_ii por faixas, 1'G 1 = soma_j
// (soma_i z_ij)^2 / k por genotipos_tv, 1 + F e 1'A22 1 por uma passada de Colleau, e o mesmo
// afim das rotas (sem afim com metafundadores). E a conferencia contra a diagonal que o ajuste
// guardou ANTES do trabalho n^2 m e n^3 das rotas exata e de Vecchia; as somas daqui e as das
// rotas diferem na ordem (ulps), o que nao pesa numa conferencia a 1e-8.
static std::vector<double> diag_gstar_previa(const Pedigree& ped, const Colleau& cl,
                                             const std::vector<std::size_t>& idx,
                                             const Genotipos& gt, const FreqZ& fz,
                                             double mistura, bool com_mf) {
  const std::size_t n = gt.n, m = gt.m;
  const double nn = static_cast<double>(n);
  std::vector<double> dg = diag_g_faixas(gt, fz);
  double a0 = 0.0, b = 1.0;
  if (!com_mf) {
    std::vector<double> centro(m), imp(m), um(n, 1.0), t(m);
    for (std::size_t j = 0; j < m; j++) { centro[j] = 2.0 * fz.p[j]; imp[j] = fz.media[j]; }
    genotipos_tv(gt, centro.data(), imp.data(), um.data(), 1, t.data());
    double soma_g = 0.0, tr_g = 0.0, tr_a = 0.0, soma_a = 0.0;
    for (std::size_t j = 0; j < m; j++) if (fz.usa[j]) soma_g += t[j] * t[j];
    soma_g /= fz.denom;
    std::vector<double> v(ped.ids.size(), 0.0);
    for (std::size_t i = 0; i < n; i++) {
      tr_g += dg[i];
      tr_a += 1.0 + cl.f[idx[i]];
      v[idx[i]] = 1.0;
    }
    colleau_aplica(ped, cl, v);
    for (std::size_t i = 0; i < n; i++) soma_a += v[idx[i]];
    const double pares = n > 1 ? nn * (nn - 1.0) : 1.0;
    const double mg_diag = tr_g / nn, mg_off = (soma_g - tr_g) / pares;
    const double ma_diag = tr_a / nn, ma_off = (soma_a - tr_a) / pares;
    b = (mg_diag - mg_off) != 0.0 ? (ma_diag - ma_off) / (mg_diag - mg_off) : 1.0;
    a0 = ma_off - b * mg_off;
  }
  for (std::size_t i = 0; i < n; i++)
    dg[i] = (1.0 - mistura) * (a0 + b * dg[i]) + mistura * (1.0 + cl.f[idx[i]]);
  return dg;
}

// A G* refeita contra a do ajuste, pela diagonal que ele guardou (fit$h_prior, na ordem dos
// genotipados; vazia = sem conferencia): um pedigree, uma matriz de genotipos ou uma mistura
// diferentes dariam outra G*, e o retro-calculo seria de outro modelo.
static void confere_diag_ajuste(const std::vector<double>& dg, const std::vector<double>& hp) {
  if (hp.empty()) return;
  if (hp.size() != dg.size())
    throw Erro("the stored diagonal of G* and the genotyped animals differ in number");
  double dif = 0.0;
  for (std::size_t i = 0; i < dg.size(); i++) {
    const double d = std::fabs(dg[i] - hp[i]) / std::max(std::fabs(hp[i]), 1.0);
    if (std::isnan(d) || d > dif) dif = d;
  }
  if (!(dif <= 1e-8)) {
    char txt[64];
    std::snprintf(txt, sizeof txt, "%.3g", dif);
    throw Erro(std::string("the G* rebuilt from this pedigree and these genotypes is not the ") +
               "one the fit used (its diagonal differs by up to " + txt + "): give the " +
               "same pedigree and genotypes as in the fit");
  }
}

// G*^-1 u pela rota exata, com UMA matriz n x n: G no buffer (vanraden_em, a mesma do
// ajuste), os momentos do ajuste afim (os de G pelo triangulo, como coef_afim; os de A22 por
// somas, 1 + F na diagonal e 1'A22 1 por uma passada de Colleau, como a rota APY), G* escrita
// por cima de G no triangulo inferior, coluna a coluna de A22 por Colleau (cada coluna de uma
// thread, em fatias com a interrupcao conferida entre elas). O triangulo inferior por LINHAS,
// compactado e com a ordem invertida, e o triangulo inferior compactado por COLUNAS de P G* P
// (P a reversao): fatora-se esse, no lugar, pela Cholesky em ladrilhos, e resolve-se com P u.
// O pico e n^2 doubles, o da propria G.
static std::vector<double> gstar_resolve_exata(const Pedigree& ped, const Colleau& cl,
                                               const std::vector<std::size_t>& idx,
                                               const Genotipos& gt, const FreqZ& fz,
                                               double mistura, bool com_mf,
                                               const std::vector<double>& u, double& a0,
                                               double& b, std::vector<double>& dg) {
  const std::size_t n = gt.n, np = ped.ids.size();
  const double nn = static_cast<double>(n);
  std::vector<double> c(n * n);
  vanraden_em(gt, fz, c.data());
  checa_interrupcao();
  a0 = 0.0;
  b = 1.0;
  if (!com_mf) {
    double mg_diag = 0.0, mg_off = 0.0;
    for (std::size_t i = 0; i < n; i++) {
      mg_diag += c[i * n + i];
      for (std::size_t j = 0; j < i; j++) mg_off += c[i * n + j];
    }
    mg_diag /= nn;
    mg_off /= n > 1 ? nn * (nn - 1.0) / 2.0 : 1.0;
    double tr_a = 0.0, soma_a = 0.0;
    std::vector<double> um(np, 0.0);
    for (std::size_t i = 0; i < n; i++) { tr_a += 1.0 + cl.f[idx[i]]; um[idx[i]] = 1.0; }
    colleau_aplica(ped, cl, um);
    for (std::size_t i = 0; i < n; i++) soma_a += um[idx[i]];
    const double ma_diag = tr_a / nn;
    const double ma_off = (soma_a - tr_a) / (n > 1 ? nn * (nn - 1.0) : 1.0);
    b = (mg_diag - mg_off) != 0.0 ? (ma_diag - ma_off) / (mg_diag - mg_off) : 1.0;
    a0 = ma_off - b * mg_off;
  }
  const std::size_t fatia = 64 * static_cast<std::size_t>(std::max(1, threads()));
  for (std::size_t f0 = 0; f0 < n; f0 += fatia) {
    checa_interrupcao();
    const long f1 = static_cast<long>(std::min(n, f0 + fatia));
#ifdef _OPENMP
#pragma omp parallel num_threads(threads())
#endif
    {
      std::vector<double> v(np);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif
      for (long jj = static_cast<long>(f0); jj < f1; jj++) {
        const std::size_t jg = static_cast<std::size_t>(jj);
        std::fill(v.begin(), v.end(), 0.0);
        v[idx[jg]] = 1.0;
        colleau_aplica(ped, cl, v);
        for (std::size_t ig = jg; ig < n; ig++) {
          double& x = c[ig * n + jg];
          x = (1.0 - mistura) * (a0 + b * x) + mistura * v[idx[ig]];
        }
      }
    }
  }
  dg.assign(n, 0.0);
  for (std::size_t i = 0; i < n; i++) dg[i] = c[i * n + i];
  checa_interrupcao();

  const std::size_t N = n * (n + 1) / 2;
  for (std::size_t i = 1; i < n; i++)
    std::memmove(&c[i * (i + 1) / 2], &c[i * n], (i + 1) * sizeof(double));
  c.resize(N);
  std::reverse(c.begin(), c.end());
  if (!cholesky_empacotada(c, n, threads()))
    throw Erro("G* rebuilt from the pedigree and the genotypes is not positive-definite");
  // L L' y = P u, para a frente e para tras na coluna compactada; alfa = P y
  std::vector<double> y(n);
  for (std::size_t i = 0; i < n; i++) y[i] = u[n - 1 - i];
  std::size_t off = 0;
  for (std::size_t j = 0; j < n; j++) {
    const double* lj = &c[off];
    y[j] /= lj[0];
    const double yj = y[j];
    for (std::size_t i = j + 1; i < n; i++) y[i] -= lj[i - j] * yj;
    off += n - j;
  }
  for (std::size_t j = n; j-- > 0;) {
    off -= n - j;
    const double* lj = &c[off];
    double s = y[j];
    for (std::size_t i = j + 1; i < n; i++) s -= lj[i - j] * y[i];
    y[j] = s / lj[0];
  }
  std::vector<double> alfa(n);
  for (std::size_t i = 0; i < n; i++) alfa[i] = y[n - 1 - i];
  return alfa;
}

// Efeitos de marcador retro-calculados de um ajuste de passo unico (VanRaden 2008; Wang et al.
// 2012): a = b (1 - w) Z' G*^-1 u / k, com a MESMA G* do ajuste. Z centrada em 2p e o ausente
// na media do marcador (frequencias_z), k = 2 soma p(1-p), a + b G o ajuste afim a A22 e w a
// mistura. Com metafundadores G e a G05 (Z = M - 1, k = m/2, todos os marcadores) e nao ha
// ajuste afim (b = 1): o modelo de marcadores equivalente e u = (M - 1) a com a de variancia
// sigma2 (1 - w) / (m/2), e a mesma conta da a esperanca condicional de a dado u. O ausente
// entra como a G05 o leu, (media - 1), e nao como zero.
//
// G*^-1 u pela rota que o ajuste usou: exata (Cholesky de G*), APY (a inversa APY aplicada
// ao vetor, O(c n), nada n x n) ou Vecchia (o mesmo vecchia_de do ajuste sobre a mesma G*
// densa, entao a mesma vizinhanca e a mesma inversa).
//
// hp e a diagonal de G* que o ajuste guardou, na ordem de geno_ids (vazia = sem conferencia):
// nas rotas exata e de Vecchia ela e conferida pela diagonal previa, antes de formar G; na APY
// pela diagonal que apy_gstar entrega, antes de aplicar a inversa.
EfeitosSnp efeitos_snp(const Pedigree& ped, const std::vector<std::string>& geno_ids,
                       const Genotipos& gt, const std::vector<double>& u, double mistura,
                       const std::vector<std::string>& nucleo_apy, std::size_t vecchia_k,
                       const std::vector<double>& hp) {
  const std::size_t n = gt.n, m = gt.m;
  if (geno_ids.size() != n) throw Erro("genotyped ids and genotype rows differ in number");
  if (u.size() != n) throw Erro("one genetic value per genotyped animal is needed");
  if (!(mistura >= 0.0 && mistura <= 1.0)) throw Erro("blend outside [0, 1]");
  if (!nucleo_apy.empty() && vecchia_k > 0)
    throw Erro("apy_core and vecchia_k are two approximations of the same inverse: "
               "declare one");
  const std::vector<std::size_t> idx = linhas_genotipados(ped, geno_ids);
  const bool com_mf = !ped.eh_mf.empty() &&
                      std::any_of(ped.eh_mf.begin(), ped.eh_mf.end(), [](char c) { return c != 0; });
  RelatorioG rel;
  EfeitosSnp r;
  FreqZ fz;
  std::vector<double> alfa;
  if (!nucleo_apy.empty()) {
    std::unordered_map<std::string, std::size_t> onde;
    onde.reserve(n);
    for (std::size_t k = 0; k < n; k++) onde.emplace(geno_ids[k], k);
    std::vector<std::size_t> nc_idx;
    nc_idx.reserve(nucleo_apy.size());
    for (const std::string& nid : nucleo_apy) {
      const auto it = onde.find(nid);
      if (it == onde.end())
        throw Erro("APY core animal '" + nid + "' is not among the genotyped");
      nc_idx.push_back(it->second);
    }
    ApyGstar ap = apy_gstar(ped, idx, gt, mistura, nc_idx, com_mf, rel);
    confere_diag_ajuste(ap.gsd, hp);
    alfa = apy_aplica(ap, nc_idx, u);
    r.a = ap.a0;
    r.b = ap.b;
    r.diag_gstar = ap.gsd;
    fz = std::move(ap.fz);
    r.rota = 2;
  } else if (vecchia_k > 0) {
    fz = frequencias_z(gt, rel, com_mf);
    if (!hp.empty())
      confere_diag_ajuste(diag_gstar_previa(ped, prepara_colleau(ped), idx, gt, fz, mistura,
                                            com_mf), hp);
    Densa gstar;
    {
      Densa g(n, n);
      vanraden_em(gt, fz, g.dados.data());
      const Densa a22 = a22_colleau(ped, idx);
      if (com_mf) {
        gstar = mistura_sem_ajuste(g, a22, mistura);
      } else {
        coef_afim(g, a22, r.a, r.b);
        gstar = ajusta_g_para_a22(g, a22, mistura);
      }
    }
    r.diag_gstar.assign(n, 0.0);
    for (std::size_t i = 0; i < n; i++) r.diag_gstar[i] = gstar.at(i, i);
    const Densa gi = vecchia_de(gstar, vecchia_k);
    alfa.assign(n, 0.0);
    for (std::size_t i = 0; i < n; i++) {
      const double* gl = gi.linha(i);
      double s = 0.0;
      for (std::size_t j = 0; j < n; j++) s += gl[j] * u[j];
      alfa[i] = s;
    }
    r.rota = 3;
  } else {
    fz = frequencias_z(gt, rel, com_mf);
    const Colleau cl = prepara_colleau(ped);
    if (!hp.empty())
      confere_diag_ajuste(diag_gstar_previa(ped, cl, idx, gt, fz, mistura, com_mf), hp);
    alfa = gstar_resolve_exata(ped, cl, idx, gt, fz, mistura, com_mf, u, r.a, r.b,
                               r.diag_gstar);
    r.rota = 1;
  }
  checa_interrupcao();

  std::vector<double> centro(m), imp(m), t(m);
  for (std::size_t j = 0; j < m; j++) { centro[j] = 2.0 * fz.p[j]; imp[j] = fz.media[j]; }
  genotipos_tv(gt, centro.data(), imp.data(), alfa.data(), 1, t.data());
  r.denom = fz.denom;
  for (std::size_t j = 0; j < m; j++) {
    if (!fz.usa[j]) continue;
    r.usados.push_back(j);
    r.freq.push_back(fz.p[j]);
    r.efeito.push_back(t[j] / fz.denom * r.b * (1.0 - mistura));
  }
  return r;
}

static RelatorioG aplica_genomica_em(const Modelo& modelo, std::vector<Csc>& kinv,
                                     std::vector<double>& kinv_logdet, const Pedigree& ped,
                                     const std::vector<std::string>& geno_ids, const Genotipos& gt,
                                     double mistura,
                                     const std::vector<std::string>& nucleo_apy,
                                     std::size_t vecchia_k) {
  RelatorioG rel;
  // o A^-1 ja esta no desenho (primeiro grupo com parentesco)
  const Csc* ainv = nullptr;
  for (std::size_t g = 0; g < modelo.grupos.size(); g++)
    if (modelo.grupos[g].estrutura == Estrutura::Parentesco) { ainv = &kinv[g]; break; }
  if (!ainv) throw Erro("there is no relationship group to receive the genomics");

  Csc hinv = h_inversa(ped, *ainv, geno_ids, gt, mistura, nucleo_apy, vecchia_k, rel);

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
                           const std::vector<std::string>& geno_ids, const Genotipos& gt,
                           double mistura, const std::vector<std::string>& nucleo_apy,
                           std::size_t vecchia_k) {
  return aplica_genomica_em(d.modelo, d.kinv, d.kinv_logdet, ped, geno_ids, gt, mistura,
                            nucleo_apy, vecchia_k);
}

RelatorioG aplica_genomica(DesenhoMT& d, const Pedigree& ped,
                           const std::vector<std::string>& geno_ids, const Genotipos& gt,
                           double mistura, const std::vector<std::string>& nucleo_apy,
                           std::size_t vecchia_k) {
  return aplica_genomica_em(d.modelo, d.kinv, d.kinv_logdet, ped, geno_ids, gt, mistura,
                            nucleo_apy, vecchia_k);
}

RelatorioG aplica_genomica(DesenhoAR& d, const Pedigree& ped,
                           const std::vector<std::string>& geno_ids, const Genotipos& gt,
                           double mistura, const std::vector<std::string>& nucleo_apy,
                           std::size_t vecchia_k) {
  return aplica_genomica_em(d.modelo, d.kinv, d.kinv_logdet, ped, geno_ids, gt, mistura,
                            nucleo_apy, vecchia_k);
}

}  // namespace br
