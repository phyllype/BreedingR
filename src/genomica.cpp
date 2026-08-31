// Passo unico: G de VanRaden, A22^-1 pelo complemento de Schur, e H^-1.
//
// A armadilha silenciosa do passo unico esta toda numa confusao de blocos:
//
//   A22^-1  =  B22 - B21 B11^-1 B12        (Schur sobre A^-1 = [[B11,B12],[B21,B22]])
//
// e isso NAO e o bloco 22 de A^-1. As duas matrizes tem a mesma forma, as duas sao
// simetricas e definidas, e so uma esta certa. O gate que separa as duas exige que a
// fixture as distinga — senao o teste e cego.
//
// B11^-1 B12 nunca e formado denso: B11 e o bloco NAO-genotipado inteiro, e a 20 mil
// animais com 2 mil genotipados isso seria uma densa de 18k x 18k para obter uma resposta
// de 2k x 2k. Uma fatoracao esparsa de B11 e n_geno resolucoes triangulares fazem o mesmo.

#include "mme.h"

// O BLAS do proprio R para a G: o produto Z Z' e o unico n^2 m do pacote, e o dsyrk
// faz em minutos o que o laco triplo fazia em horas de uma thread so.
#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

namespace br {

// G de VanRaden: Z Z' / (2 sum p(1-p)), com Z = M - 2p.
//
// Codigo ausente TEM de chegar como NaN e e imputado pela MEDIA do marcador. Nunca por
// zero: zero e um genotipo valido, e a confusao mudaria as frequencias e a G inteira sem
// nenhum erro visivel.
Densa vanraden_g(Densa& m, RelatorioG& rel) {
  const std::size_t n = m.nlin, nm = m.ncol;
  if (n == 0 || nm == 0) throw Erro("empty genotypes");

  std::vector<double> p(nm);
  std::vector<char> usa(nm, 1);
  for (std::size_t j = 0; j < nm; j++) {
    double soma = 0.0;
    std::size_t k = 0;
    for (std::size_t i = 0; i < n; i++) {
      const double x = m.at(i, j);
      if (std::isfinite(x)) { soma += x; k++; }
    }
    if (k == 0) { usa[j] = 0; rel.n_monomorficos++; continue; }
    const double media = soma / static_cast<double>(k);
    p[j] = media / 2.0;
    if (p[j] <= 0.0 || p[j] >= 1.0) { usa[j] = 0; rel.n_monomorficos++; continue; }
    for (std::size_t i = 0; i < n; i++)
      if (!std::isfinite(m.at(i, j))) { m.at(i, j) = media; rel.n_imputados++; }
  }

  double denom = 0.0;
  for (std::size_t j = 0; j < nm; j++)
    if (usa[j]) denom += 2.0 * p[j] * (1.0 - p[j]);
  if (!(denom > 0.0)) throw Erro("2 sum p(1-p) is not positive: the markers do not vary");

  // Z Z' por dsyrk em blocos de marcadores: o bloco Z_b (n x B, coluna-major) e montado
  // ja centrado, e o dsyrk acumula C += Z_b Z_b'. Memoria extra: um bloco + o C
  // coluna-major, nada proporcional a m.
  const std::size_t B = 2048;
  const int ni = static_cast<int>(n);
  std::vector<double> c(n * n, 0.0), zbuf(n * B);
  double beta = 0.0;
  std::size_t bcol = 0;
  for (std::size_t j = 0; j < nm; j++) {
    if (!usa[j]) continue;
    const double dp = 2.0 * p[j];
    for (std::size_t i = 0; i < n; i++) zbuf[i + n * bcol] = m.at(i, j) - dp;
    if (++bcol == B) {
      const int k = static_cast<int>(bcol);
      const double um = 1.0;
      F77_CALL(dsyrk)("U", "N", &ni, &k, &um, zbuf.data(), &ni, &beta, c.data(), &ni
                      FCONE FCONE);
      beta = 1.0;
      bcol = 0;
    }
  }
  if (bcol > 0) {
    const int k = static_cast<int>(bcol);
    const double um = 1.0;
    F77_CALL(dsyrk)("U", "N", &ni, &k, &um, zbuf.data(), &ni, &beta, c.data(), &ni
                    FCONE FCONE);
  }
  Densa g(n, n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t k2 = i; k2 < n; k2++) {
      const double v = c[i + n * k2] / denom;
      g.at(i, k2) = v;
      g.at(k2, i) = v;
    }
  return g;
}

// A22^-1 pelo Schur ESPARSO sobre A^-1, sem nunca formar B11^-1 denso.
Densa a22_inversa(const Csc& ainv, const std::vector<std::size_t>& geno) {
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

  // B11 (triangulo inferior esparso), B12 denso por coluna genotipada, B22 denso
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  Densa b12(n1, n2), b22(n2, n2);
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
        b22.at(pr - n1, pc - n1) = x;
        if (r != c) b22.at(pc - n1, pr - n1) = x;
      } else {
        // um em cada bloco: entra em B12 (nao-genotipado nas linhas)
        const std::size_t inng = rg ? pc : pr;
        const std::size_t ig = rg ? pr - n1 : pc - n1;
        b12.at(inng, ig) = x;
      }
    }

  if (n1 == 0) return b22;   // todos genotipados: A22^-1 = A^-1 inteiro

  Csc b11 = de_triplos(n1, n1, li, cj, v);
  std::vector<std::size_t> perm = grau_minimo(b11);
  Csc pb = permuta_sim(b11, perm);
  Simbolica sb = simbolica(pb);
  Csc L;
  if (!cholesky(pb, sb, L))
    throw Erro("the non-genotyped block of A^-1 is not positive-definite");

  // A22^-1 = B22 - B21 (B11^-1 B12): uma resolucao por coluna genotipada
  Densa out = b22;
  std::vector<double> col(n1), pcol(n1);
  for (std::size_t jg = 0; jg < n2; jg++) {
    for (std::size_t i = 0; i < n1; i++) col[i] = b12.at(i, jg);
    for (std::size_t i = 0; i < n1; i++) pcol[i] = col[perm[i]];
    std::vector<double> px = resolve(L, pcol);
    for (std::size_t i = 0; i < n1; i++) col[perm[i]] = px[i];
    for (std::size_t ig = 0; ig < n2; ig++) {
      double s = 0.0;
      for (std::size_t i = 0; i < n1; i++) s += b12.at(i, ig) * col[i];
      out.at(ig, jg) -= s;
    }
  }
  // simetriza contra o arredondamento das resolucoes
  for (std::size_t i = 0; i < n2; i++)
    for (std::size_t j = 0; j < i; j++) {
      const double mdi = 0.5 * (out.at(i, j) + out.at(j, i));
      out.at(i, j) = mdi;
      out.at(j, i) = mdi;
    }
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
  for (std::size_t a = 0; a < n2; a++)
    for (std::size_t b = 0; b <= a; b++) {
      std::size_t i = geno[a], j = geno[b];
      double x = gstar_inv.at(a, b) - a22_inv.at(a, b);
      if (i < j) std::swap(i, j);
      li.push_back(static_cast<std::uint32_t>(i));
      cj.push_back(static_cast<std::uint32_t>(j));
      v.push_back(x);
    }
  return de_triplos(n, n, li, cj, v);
}

// Substitui o K^-1 dos grupos com parentesco pelo H^-1, e o logdet correspondente.
//
// Genotipado fora do pedigree e ERRO explicito, nao descarte silencioso: um genotipado sem
// linha em A nao tem onde entrar em H^-1, e some-lo mudaria a analise sem aviso.
// A inversa APY de uma G densa ja ajustada: nucleo exato, jovens por recursao condicional.
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

  Densa gcn(nc, nj), pmat(nc, nj);
  for (std::size_t a = 0; a < nc; a++)
    for (std::size_t b = 0; b < nj; b++) gcn.at(a, b) = g.at(nucleo[a], jovens[b]);
  for (std::size_t a = 0; a < nc; a++)
    for (std::size_t b = 0; b < nj; b++) {
      double sacc = 0.0;
      for (std::size_t k = 0; k < nc; k++) sacc += gcc_inv.at(a, k) * gcn.at(k, b);
      pmat.at(a, b) = sacc;
    }

  double diag_media = 0.0;
  for (std::size_t i = 0; i < n; i++) diag_media += g.at(i, i);
  diag_media /= static_cast<double>(n);
  const double limiar = 1e-8 * diag_media;

  std::vector<double> minv(nj);
  for (std::size_t b = 0; b < nj; b++) {
    double mend = g.at(jovens[b], jovens[b]);
    for (std::size_t k = 0; k < nc; k++) mend -= gcn.at(k, b) * pmat.at(k, b);
    if (mend <= limiar)
      throw Erro("degenerate Mendelian residual in APY: a non-core animal is collinear with the core "
                 "(clone or duplicate). Enlarge the core or the blend.");
    minv[b] = 1.0 / mend;
  }

  for (std::size_t a = 0; a < nc; a++)
    for (std::size_t b = 0; b < nc; b++) {
      double sacc = gcc_inv.at(a, b);
      for (std::size_t k = 0; k < nj; k++) sacc += pmat.at(a, k) * minv[k] * pmat.at(b, k);
      out.at(nucleo[a], nucleo[b]) = sacc;
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

// A inversa de Vecchia: cada animal condiciona nos SEUS k vizinhos mais proximos entre
// os anteriores, nao num nucleo global. E a generalizacao da APY — e do proprio A^-1 de
// Henderson, que e exatamente Vecchia com os PAIS como conjunto de condicionamento
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
// e os K^-1 com seus log-determinantes — uni, multi e AR(1) carregam exatamente esses
// campos, entao o passo unico e UM caminho, nao tres.
static RelatorioG aplica_genomica_em(const Modelo& modelo, std::vector<Csc>& kinv,
                                     std::vector<double>& kinv_logdet, const Pedigree& ped,
                                     const std::vector<std::string>& geno_ids, Densa& m,
                                     double mistura,
                                     const std::vector<std::string>& nucleo_apy,
                                     std::size_t vecchia_k) {
  RelatorioG rel;
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

  // o A^-1 ja esta no desenho (primeiro grupo com parentesco)
  const Csc* ainv = nullptr;
  for (std::size_t g = 0; g < modelo.grupos.size(); g++)
    if (modelo.grupos[g].estrutura == Estrutura::Parentesco) { ainv = &kinv[g]; break; }
  if (!ainv) throw Erro("there is no relationship group to receive the genomics");

  Densa g = vanraden_g(m, rel);
  Densa a22i = a22_inversa(*ainv, idx);
  Densa a22 = inv_pd(a22i);
  Densa gstar = ajusta_g_para_a22(g, a22, mistura);

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
    for (const std::string& nid : nucleo_apy) {
      bool achou = false;
      for (std::size_t k = 0; k < geno_ids.size(); k++)
        if (geno_ids[k] == nid) { nc_idx.push_back(k); achou = true; break; }
      if (!achou)
        throw Erro("APY core animal '" + nid + "' is not among the genotyped");
    }
    gstar_inv = apy_de(gstar, nc_idx);
  }
  Csc hinv = constroi_hinv(*ainv, idx, gstar_inv, a22i);

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
