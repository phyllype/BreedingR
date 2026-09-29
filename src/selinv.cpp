// Inversa seletiva de Takahashi (Takahashi, Fagan e Chin, 1973): os elementos de C^-1 no padrao do fator.
//
// E o que a PEV, a acuracia e os tracos do score precisam, e e a unica parte de C^-1 que se
// pode pagar. A inversa cheia da matriz de coeficientes e densa: a 20.400 colunas sao 3,3 GB,
// e cada posicao fora do padrao e uma que ninguem le.
//
// ## Por que o padrao de L basta, e fecha
//
// Com C = L D L' e L unitaria inferior, a recorrencia de Takahashi calcula, de baixo para
// cima,
//
//   Z[i,j] = - soma_{k > j, L[k,j] != 0}  Z[i,k] . L[k,j]      para i > j
//   Z[j,j] = 1/D[j] - soma_{k > j, L[k,j] != 0}  L[k,j] . Z[k,j]
//
// Todo Z que ela le do lado direito esta numa posicao (i,k) com i e k ambos no padrao da
// coluna j. A propriedade de preenchimento do fator de Cholesky diz que esse par tambem
// esta no padrao. Logo a recorrencia nunca pede um elemento que nao calculou, e nunca
// precisa de um fora de L. Esse fechamento e a razao de tudo isto funcionar.
//
// Aqui o fator e L L', entao D[j] = L[j,j]^2 e a unitaria e L~[k,j] = L[k,j] / L[j,j].

#include "mme.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace br {

// Quantas colunas finais do fator sao completamente densas.
//
// Num modelo de passo unico o bloco genomico e um clique denso, e o produto de Kronecker com
// a covariancia do grupo o propaga: as ultimas colunas de L ficam cheias. Detectar esse bloco
// e o que permite a forma fechada, e ele so existe se a ordenacao deixou o clique no FIM —
// que e por que o pacote ordena por grau minimo e nao por Cuthill-McKee reverso.
std::size_t bloco_denso_simbolico(const Simbolica& sb) {
  const std::size_t n = sb.n;
  std::size_t t = n;
  while (t > 0) {
    const std::size_t j = t - 1;
    if (sb.colptr[j + 1] - sb.colptr[j] != n - j) break;
    t--;
  }
  return n - t;
}

std::size_t bloco_denso_final(const Csc& L) {
  const std::size_t n = L.ncol;
  std::size_t t = n;
  while (t > 0) {
    const std::size_t j = t - 1;
    if (L.colptr[j + 1] - L.colptr[j] != n - j) break;
    t--;
  }
  return n - t;
}

// Metodos de SelInv (declarada em mme.h), fora da classe para nao virarem inline de uma
// definicao que os outros modulos nao veem.
bool SelInv::get(std::size_t i, std::size_t j, double& out) const {
  std::size_t r = i, c = j;
  if (r < c) std::swap(r, c);
  if (c >= n) return false;
  const std::size_t lo = colptr[c], hi = colptr[c + 1];
  auto it = std::lower_bound(linha.begin() + lo, linha.begin() + hi,
                             static_cast<std::uint32_t>(r));
  if (it != linha.begin() + hi && *it == static_cast<std::uint32_t>(r)) {
    out = valor[static_cast<std::size_t>(it - linha.begin())];
    return true;
  }
  return false;
}

std::vector<double> SelInv::diagonal() const {
  std::vector<double> d(n);
  for (std::size_t j = 0; j < n; j++) d[j] = valor[colptr[j]];
  return d;
}

// `bloco` de 0 significa detectar; passe um valor explicito so para forcar a recorrencia pura.
SelInv inversa_seletiva(const Csc& L, std::size_t bloco = 0) {
  const std::size_t n = L.ncol;
  SelInv z;
  z.n = n;
  z.colptr = L.colptr;
  z.linha = L.linha;
  z.valor.assign(L.nnz(), 0.0);

  const std::size_t t = bloco == 0 ? bloco_denso_final(L) : std::min(bloco, n);
  const std::size_t inicio = n - std::min(t, n);

  // Forma fechada no bloco denso final.
  //
  // Para um bloco final denso T o complemento de Schur e exatamente L[T,T] L[T,T]', logo
  // Z[T,T] = inv(L[T,T] L[T,T]'). Como L[T,T] JA E o fator desse produto, formar o produto
  // e fatorar de novo seria fazer a mesma fatoracao duas vezes. O bloco vai ao nucleo
  // compactado em ladrilhos (inversa_empacotada), em paralelo: copiado das entradas de L
  // (zero onde o padrao nao tem entrada, quando `bloco` pede mais que a cauda densa) e
  // devolvido so nas posicoes do padrao.
  const int nth = threads();
  if (t >= 2) {
    std::vector<std::size_t> off(t + 1, 0);
    for (std::size_t j = 0; j < t; j++) off[j + 1] = off[j] + (t - j);
    std::vector<double> lp(off[t], 0.0), zp(off[t]);
    for (std::size_t j = inicio; j < n; j++)
      for (std::size_t p = L.colptr[j]; p < L.colptr[j + 1]; p++) {
        const std::size_t i = L.linha[p];
        if (i >= inicio) lp[off[j - inicio] + (i - j)] = L.valor[p];
      }
    inversa_empacotada(lp.data(), t, zp.data(), nth);
    for (std::size_t j = inicio; j < n; j++)
      for (std::size_t p = L.colptr[j]; p < L.colptr[j + 1]; p++) {
        const std::size_t i = L.linha[p];
        if (i >= inicio) z.valor[p] = zp[off[j - inicio] + (i - j)];
      }
  }

  // Recorrencia para cima nas colunas restantes.
  //
  // A leitura ingenua da identidade pede, para cada par ORDENADO (i,k) de linhas abaixo da
  // diagonal, o elemento Z[max,min], e acha-lo por busca binaria custa m^2 log sondagens
  // aleatorias numa coluna com m elementos. Nao precisa: a coluna c guarda Z[r,c] para r >= c
  // em ordem crescente, e pela simetria esse UNICO elemento serve aos dois pares ordenados.
  // Percorrer as colunas guardadas com um vetor de selo dizendo quais linhas estao no padrao
  // faz a mesma soma dupla sem busca nenhuma e com leitura sequencial.
  //
  // Um par fora do padrao continua contribuindo zero, que e a mesma semantica que a busca
  // binaria tinha ao nao encontrar.
  //
  // Em paralelo POR NIVEL da arvore de eliminacao: a coluna j le so as colunas das linhas do
  // seu padrao, que sao ancestrais dela na arvore, e escreve so a propria coluna. Colunas do
  // mesmo nivel (mesma distancia ate a cauda ou a raiz) sao independentes; cada uma e feita
  // inteira por uma thread, com a mesma ordem de soma do laco serial, entao o resultado e bit
  // a bit o mesmo com qualquer numero de threads.
  const std::size_t primeiro = (t >= 2) ? inicio : n;
  std::vector<std::uint32_t> nivel(primeiro, 0);
  std::uint32_t maior = 0;
  for (std::size_t j = primeiro; j-- > 0;) {
    const std::size_t lo = L.colptr[j], hi = L.colptr[j + 1];
    if (hi - lo > 1) {
      const std::size_t pai = L.linha[lo + 1];
      nivel[j] = pai >= primeiro ? 0 : nivel[pai] + 1;
    }
    maior = std::max(maior, nivel[j]);
  }
  std::vector<std::size_t> ini_nivel(static_cast<std::size_t>(maior) + 2, 0);
  for (std::size_t j = 0; j < primeiro; j++) ini_nivel[nivel[j] + 1]++;
  for (std::size_t v = 0; v + 1 < ini_nivel.size(); v++) ini_nivel[v + 1] += ini_nivel[v];
  std::vector<std::size_t> por_nivel(primeiro);
  {
    std::vector<std::size_t> pos = ini_nivel;
    for (std::size_t j = primeiro; j-- > 0;) por_nivel[pos[nivel[j]]++] = j;
  }

  struct Rascunho {
    std::vector<std::uint32_t> marca, ondice, linhas;
    std::vector<double> acc, vals;
    std::uint32_t selo = 0;
  };
  std::vector<Rascunho> rasc(static_cast<std::size_t>(std::max(1, nth)));
  auto coluna = [&](std::size_t j, Rascunho& r) {
    if (r.marca.empty()) { r.marca.assign(n, 0); r.ondice.assign(n, 0); r.acc.assign(n, 0.0); }
    const std::size_t lo = L.colptr[j], hi = L.colptr[j + 1];
    const double djj = L.valor[lo];
    const double d = djj * djj;
    r.linhas.clear();
    r.vals.clear();
    for (std::size_t p = lo + 1; p < hi; p++) {
      r.linhas.push_back(L.linha[p]);
      r.vals.push_back(L.valor[p] / djj);
    }
    const std::size_t m = r.linhas.size();
    if (m == 0) { z.valor[lo] = 1.0 / d; return; }
    r.selo++;
    for (std::size_t idx = 0; idx < m; idx++) {
      r.marca[r.linhas[idx]] = r.selo;
      r.ondice[r.linhas[idx]] = static_cast<std::uint32_t>(idx);
      r.acc[idx] = 0.0;
    }
    const std::uint32_t ultima = r.linhas.back();
    for (std::size_t ic = 0; ic < m; ic++) {
      const std::size_t c = r.linhas[ic];
      const double lc = r.vals[ic];
      for (std::size_t p = z.colptr[c]; p < z.colptr[c + 1]; p++) {
        const std::uint32_t rr = z.linha[p];
        if (rr > ultima) break;                // a coluna esta ordenada
        if (r.marca[rr] != r.selo) continue;
        const std::size_t ir = r.ondice[rr];
        const double zv = z.valor[p];
        if (ir == ic) r.acc[ic] += zv * lc;
        else { r.acc[ir] += zv * lc; r.acc[ic] += zv * r.vals[ir]; }
      }
    }
    for (std::size_t idx = 0; idx < m; idx++) z.valor[lo + 1 + idx] = -r.acc[idx];
    double dg = 0.0;
    for (std::size_t idx = 0; idx < m; idx++) dg += r.vals[idx] * z.valor[lo + 1 + idx];
    z.valor[lo] = 1.0 / d - dg;
  };
  for (std::size_t v = 0; v + 1 < ini_nivel.size(); v++) {
    const std::size_t a = ini_nivel[v], b = ini_nivel[v + 1];
    if (b - a < 64 || nth <= 1) {
      for (std::size_t q = a; q < b; q++) coluna(por_nivel[q], rasc[0]);
      continue;
    }
#ifdef _OPENMP
#pragma omp parallel num_threads(nth)
#endif
    {
#ifdef _OPENMP
      Rascunho& r = rasc[static_cast<std::size_t>(omp_get_thread_num())];
#pragma omp for schedule(dynamic, 16)
#else
      Rascunho& r = rasc[0];
#endif
      for (long q = static_cast<long>(a); q < static_cast<long>(b); q++)
        coluna(por_nivel[static_cast<std::size_t>(q)], r);
    }
  }

  return z;
}

}  // namespace br
