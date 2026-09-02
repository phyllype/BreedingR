// Cholesky esparsa: arvore de eliminacao, fatoracao simbolica, fatoracao numerica
// olhando-para-cima, resolucao triangular, permutacao simetrica e grau minimo.
//
// E o que tira o pacote do caminho denso. A Cholesky densa de estruturas.cpp fica como
// REFERENCIA: tudo aqui e conferido contra ela, nunca contra si mesmo.
//
// Duas ideias decidem o formato do codigo (Davis, 2006, Direct Methods for Sparse Linear Systems):
//
//  - a ARVORE DE ELIMINACAO diz, para cada coluna, a que coluna anterior o primeiro
//    elemento fora da diagonal pertence. Subir a arvore da o padrao de nao-zeros de uma
//    LINHA de L sem formar nada.
//  - a fatoracao OLHANDO-PARA-CIMA calcula L uma linha de cada vez, resolvendo um sistema
//    triangular contra a parte ja pronta. Ela precisa do padrao da linha primeiro, que e
//    exatamente o que a arvore da, e custa O(nnz(L)) em vez de algo quadratico.
//
// A entrada e o triangulo SUPERIOR de uma simetrica em CSC. Para uma simetrica isso e a
// transposta do triangulo inferior.

#include "mme.h"

namespace br {

// Triangulo superior, incluindo a diagonal.
Csc triu(const Csc& a) {
  Csc out(a.nlin, a.ncol);
  out.linha.reserve(a.nnz() / 2 + a.ncol);
  out.valor.reserve(a.nnz() / 2 + a.ncol);
  for (std::size_t c = 0; c < a.ncol; c++) {
    for (std::size_t k = a.colptr[c]; k < a.colptr[c + 1]; k++)
      if (a.linha[k] <= c) { out.linha.push_back(a.linha[k]); out.valor.push_back(a.valor[k]); }
    out.colptr[c + 1] = out.linha.size();
  }
  return out;
}

// Arvore de eliminacao, com compressao de caminho.
std::vector<std::int64_t> arvore(const Csc& au) {
  const std::size_t n = au.ncol;
  std::vector<std::int64_t> pai(n, -1), ancestral(n, -1);
  for (std::size_t k = 0; k < n; k++)
    for (std::size_t p = au.colptr[k]; p < au.colptr[k + 1]; p++) {
      std::size_t i = au.linha[p];
      if (i >= k) continue;
      while (true) {
        const std::int64_t prox = ancestral[i];
        ancestral[i] = static_cast<std::int64_t>(k);
        if (prox < 0) { pai[i] = static_cast<std::int64_t>(k); break; }
        if (static_cast<std::size_t>(prox) == k) break;
        i = static_cast<std::size_t>(prox);
      }
    }
  return pai;
}

// Padrao de nao-zeros da linha k de L, em ordem topologica, devolvido em s[topo..n).
static std::size_t alcance(const Csc& au, std::size_t k, const std::vector<std::int64_t>& pai,
                           std::vector<std::uint32_t>& s, std::vector<std::uint32_t>& marca) {
  const std::size_t n = au.ncol;
  std::size_t topo = n;
  marca[k] = static_cast<std::uint32_t>(k) + 1;
  for (std::size_t p = au.colptr[k]; p < au.colptr[k + 1]; p++) {
    std::size_t i = au.linha[p];
    if (i > k) continue;
    std::size_t len = 0;
    while (marca[i] != static_cast<std::uint32_t>(k) + 1) {
      s[len++] = static_cast<std::uint32_t>(i);
      marca[i] = static_cast<std::uint32_t>(k) + 1;
      if (pai[i] < 0) break;
      i = static_cast<std::size_t>(pai[i]);
    }
    while (len > 0) s[--topo] = s[--len];
  }
  return topo;
}

Simbolica simbolica(const Csc& au) {
  const std::size_t n = au.ncol;
  Simbolica sb;
  sb.n = n;
  sb.pai = arvore(au);
  std::vector<std::size_t> conta(n, 0);
  std::vector<std::uint32_t> s(n), marca(n, 0);
  for (std::size_t k = 0; k < n; k++) {
    const std::size_t topo = alcance(au, k, sb.pai, s, marca);
    for (std::size_t t = topo; t < n; t++) conta[s[t]]++;
    conta[k]++;                      // a diagonal
  }
  sb.colptr.assign(n + 1, 0);
  for (std::size_t j = 0; j < n; j++) sb.colptr[j + 1] = sb.colptr[j] + conta[j];
  return sb;
}

// Fatoracao numerica olhando-para-cima. Devolve false se a matriz nao for positiva-definida
// — o que, neste engine, e informacao sobre theta e nao uma falha a reportar como tal.
bool cholesky(const Csc& au, const Simbolica& sb, Csc& L) {
  const std::size_t n = sb.n;
  const std::size_t nz = sb.nnz();
  std::vector<std::size_t> colptr = sb.colptr;
  std::vector<std::uint32_t> linha(nz, 0);
  std::vector<double> valor(nz, 0.0);
  std::vector<std::size_t> prox(colptr.begin(), colptr.begin() + n);

  std::vector<double> x(n, 0.0);
  std::vector<std::uint32_t> s(n), marca(n, 0);

  for (std::size_t k = 0; k < n; k++) {
    const std::size_t topo = alcance(au, k, sb.pai, s, marca);

    double d = 0.0;
    for (std::size_t p = au.colptr[k]; p < au.colptr[k + 1]; p++) {
      const std::size_t i = au.linha[p];
      if (i < k) x[i] = au.valor[p];
      else if (i == k) d = au.valor[p];
    }

    for (std::size_t t = topo; t < n; t++) {
      const std::size_t j = s[t];
      const double lkj = x[j] / valor[colptr[j]];   // L[j,j] e o primeiro da coluna
      x[j] = 0.0;
      for (std::size_t p = colptr[j] + 1; p < prox[j]; p++) x[linha[p]] -= valor[p] * lkj;
      d -= lkj * lkj;
      // Uma simbolica menor que o padrao real transbordaria a coluna j para a j+1 e
      // corromperia L sem erro nenhum. E invariante de quem chama, entao para aqui.
      if (prox[j] >= colptr[j + 1])
        throw Erro("symbolic factorization too small: the matrix pattern grew afterwards");
      linha[prox[j]] = static_cast<std::uint32_t>(k);
      valor[prox[j]] = lkj;
      prox[j]++;
    }

    if (!(d > 0.0) || !std::isfinite(d)) return false;
    linha[prox[k]] = static_cast<std::uint32_t>(k);
    valor[prox[k]] = std::sqrt(d);   // a diagonal vai PRIMEIRO, e o laco acima conta com isso
    prox[k]++;
  }

  L = Csc(n, n);
  L.linha.reserve(nz);
  L.valor.reserve(nz);
  for (std::size_t j = 0; j < n; j++) {
    for (std::size_t p = colptr[j]; p < prox[j]; p++) {
      L.linha.push_back(linha[p]);
      L.valor.push_back(valor[p]);
    }
    L.colptr[j + 1] = L.linha.size();
  }
  return true;
}

double logdet(const Csc& L) {
  double s = 0.0;
  for (std::size_t j = 0; j < L.ncol; j++) s += std::log(L.valor[L.colptr[j]]);
  return 2.0 * s;
}

// Resolve L L' x = b.
std::vector<double> resolve(const Csc& L, const std::vector<double>& b) {
  const std::size_t n = L.ncol;
  std::vector<double> x = b;
  for (std::size_t j = 0; j < n; j++) {          // para a frente
    x[j] /= L.valor[L.colptr[j]];
    for (std::size_t p = L.colptr[j] + 1; p < L.colptr[j + 1]; p++)
      x[L.linha[p]] -= L.valor[p] * x[j];
  }
  for (std::size_t j = n; j-- > 0;) {            // para tras
    for (std::size_t p = L.colptr[j] + 1; p < L.colptr[j + 1]; p++)
      x[j] -= L.valor[p] * x[L.linha[p]];
    x[j] /= L.valor[L.colptr[j]];
  }
  return x;
}

// Permuta simetricamente e devolve o triangulo SUPERIOR de P A P'.
//
// A ENTRADA TEM DE TER CADA PAR GUARDADO UMA VEZ SO — um triangulo, qualquer um deles.
// Duas armadilhas ja apareceram aqui, em direcoes opostas:
//
//  - com uma matriz CHEIA, (i,j) e (j,i) caem na mesma posicao permutada e a montagem por
//    triplos SOMA as duas: o resultado fica com o dobro fora da diagonal e a fatoracao passa
//    a dizer que a matriz nao e positiva-definida.
//  - filtrando por um triangulo fixo, uma entrada guardada no OUTRO triangulo e descartada
//    inteira. Foi o que aconteceu na primeira versao deste arquivo: a entrada vinha no
//    inferior, o filtro mantinha o superior, e sobrava so a diagonal. L saia diagonal, o
//    bloco denso nao aparecia, e a resolucao dava numero errado sem erro nenhum.
//
// Por isso nao ha filtro: cada entrada e mapeada uma vez e colocada no triangulo superior do
// indice novo, venha ela de onde vier.
Csc permuta_sim(const Csc& a, const std::vector<std::size_t>& perm) {
  const std::size_t n = a.ncol;
  if (perm.size() != n) throw Erro("permutation of the wrong size");
  std::vector<std::size_t> inv(n);
  for (std::size_t k = 0; k < n; k++) {
    if (perm[k] >= n) throw Erro("permutation out of range");
    inv[perm[k]] = k;
  }
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  li.reserve(a.nnz()); cj.reserve(a.nnz()); v.reserve(a.nnz());
  for (std::size_t c = 0; c < n; c++)
    for (std::size_t k = a.colptr[c]; k < a.colptr[c + 1]; k++) {
      const std::size_t r = a.linha[k];
      std::size_t i = inv[r], j = inv[c];
      if (i > j) std::swap(i, j);                // triangulo superior no novo indice
      li.push_back(static_cast<std::uint32_t>(i));
      cj.push_back(static_cast<std::uint32_t>(j));
      v.push_back(a.valor[k]);
    }
  return de_triplos(n, n, li, cj, v);
}

// Grau minimo (Tinney e Walker, 1967) sobre o grafo quociente.
//
// Medido, nao assumido. Num pedigree de 3.000 animais mais um clique genomico o preenchimento
// foi: natural 1.332.217, Cuthill-McKee reverso 746.017, Cuthill-McKee simples 2.506.473,
// grau minimo 201.487. O RCM (Cuthill e McKee, 1969) reduz pela metade mas poe os nos de grau alto PRIMEIRO, que e
// exatamente onde a forma fechada do bloco denso final nao os enxerga.
//
// O grafo quociente e o que torna isto pagavel: uma variavel eliminada vira um ELEMENTO que
// guarda o clique que ela criou, entao o preenchimento nunca e materializado.
std::vector<std::size_t> grau_minimo(const Csc& a) {
  const std::size_t n = a.ncol;
  std::vector<std::size_t> ordem;
  if (n == 0) return ordem;
  ordem.reserve(n);

  std::vector<std::vector<std::uint32_t>> av(n), ev(n), le(n);
  for (std::size_t c = 0; c < n; c++)
    for (std::size_t k = a.colptr[c]; k < a.colptr[c + 1]; k++) {
      const std::size_t r = a.linha[k];
      if (r == c) continue;
      av[c].push_back(static_cast<std::uint32_t>(r));
      av[r].push_back(static_cast<std::uint32_t>(c));
    }
  for (auto& l : av) { std::sort(l.begin(), l.end()); l.erase(std::unique(l.begin(), l.end()), l.end()); }

  std::vector<char> vivo(n, 1);
  std::vector<std::uint32_t> marca(n, 0), no_lp(n, 0);
  std::uint32_t selo = 0, selo_lp = 0;

  // vizinhanca de uma variavel viva: vizinhos-variavel mais os membros de cada elemento
  auto viz = [&](std::size_t i) {
    selo++;
    std::vector<std::uint32_t> out;
    marca[i] = selo;
    for (std::uint32_t j : av[i]) if (vivo[j] && marca[j] != selo) { marca[j] = selo; out.push_back(j); }
    for (std::uint32_t e : ev[i])
      for (std::uint32_t j : le[e]) if (vivo[j] && marca[j] != selo) { marca[j] = selo; out.push_back(j); }
    return out;
  };

  std::vector<std::size_t> grau(n);
  for (std::size_t i = 0; i < n; i++) grau[i] = av[i].size();

  std::vector<std::vector<std::uint32_t>> baldes(n + 1);
  for (std::size_t i = 0; i < n; i++) baldes[std::min(grau[i], n)].push_back(static_cast<std::uint32_t>(i));
  std::size_t lo = 0;

  while (ordem.size() < n) {
    std::size_t p = static_cast<std::size_t>(-1);
    while (lo <= n) {
      while (!baldes[lo].empty()) {
        const std::uint32_t cand = baldes[lo].back();
        baldes[lo].pop_back();
        if (vivo[cand] && std::min(grau[cand], n) == lo) { p = cand; break; }
      }
      if (p != static_cast<std::size_t>(-1)) break;
      lo++;
    }
    if (p == static_cast<std::size_t>(-1)) {          // o que sobrou esta isolado
      for (std::size_t i = 0; i < n; i++) if (vivo[i]) { vivo[i] = 0; ordem.push_back(i); }
      break;
    }

    std::vector<std::uint32_t> lp = viz(p);
    vivo[p] = 0;
    ordem.push_back(p);
    le[p] = lp;

    // pertencer a L_p por SELO, e nao varrendo L_p. A versao com varredura era quadratica
    // dentro de um laco que ja e sobre L_p, e a ordenacao virava a parte mais cara do ajuste.
    selo_lp++;
    for (std::uint32_t x : lp) no_lp[x] = selo_lp;

    for (std::uint32_t iu : lp) {
      const std::size_t i = iu;
      auto& a_i = av[i];
      a_i.erase(std::remove_if(a_i.begin(), a_i.end(),
                               [&](std::uint32_t x) { return x == p || !vivo[x]; }), a_i.end());
      auto& e_i = ev[i];
      e_i.erase(std::remove(e_i.begin(), e_i.end(), static_cast<std::uint32_t>(p)), e_i.end());
      e_i.push_back(static_cast<std::uint32_t>(p));
      // absorcao: um elemento cujos membros vivos ja estao dentro de L_p esta morto
      e_i.erase(std::remove_if(e_i.begin(), e_i.end(), [&](std::uint32_t e) {
        if (e == p) return false;
        for (std::uint32_t x : le[e]) if (vivo[x] && no_lp[x] != selo_lp) return false;
        return true;
      }), e_i.end());
    }
    for (std::uint32_t iu : lp) {
      const std::size_t i = iu;
      const std::size_t d = viz(i).size();
      grau[i] = d;
      const std::size_t b = std::min(d, n);
      baldes[b].push_back(iu);
      if (b < lo) lo = b;
    }
  }
  return ordem;
}

}  // namespace br
