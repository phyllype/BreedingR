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

  // no ADIADO (linha densa, ver abaixo) sai do grafo e nao so da fila: fica vivo para ser
  // posto no fim da ordem, mas nao conta em grau, nao entra em L_p e nao segura absorcao.
  // Medido na APY (3.000 animais, 1.500 genotipados, nucleo de 414): com o nucleo adiado mas
  // ainda nas listas, cada jovem eliminado criava um elemento com os 414 dentro, e todo
  // viz() e toda checagem de absorcao varriam os 414 de novo. Foram 1,3e10 passos de viz()
  // e 5e9 de absorcao contra 7e7 no caminho exato, 17 s por ordenacao em vez de 0,3 s. E o
  // nucleo contado no grau e nos vivos arrastava 341 jovens para o bloco denso junto com
  // ele (752 colunas densas em vez de nucleo + 1). O AMD tira a linha densa do grafo
  // (Amestoy, Davis e Duff, 1996); aqui ela so saia da fila.
  std::vector<char> adiado(n, 0);
  std::size_t n_adiados = 0;
  auto ativo = [&](std::uint32_t j) { return vivo[j] && !adiado[j]; };

  // a lista de um elemento so perde membros (quem morre ou e adiado nao volta), entao ela
  // e compactada ao ser lida: cada entrada morta custa uma leitura so, e nao uma a cada
  // viz() e a cada checagem de absorcao que passar por ela
  auto compacta = [&](std::uint32_t e) {
    auto& l = le[e];
    l.erase(std::remove_if(l.begin(), l.end(), [&](std::uint32_t x) { return !ativo(x); }),
            l.end());
  };

  // vizinhanca de uma variavel viva: vizinhos-variavel mais os membros de cada elemento
  auto viz = [&](std::size_t i) {
    selo++;
    std::vector<std::uint32_t> out;
    marca[i] = selo;
    for (std::uint32_t j : av[i]) if (ativo(j) && marca[j] != selo) { marca[j] = selo; out.push_back(j); }
    for (std::uint32_t e : ev[i]) {
      compacta(e);
      for (std::uint32_t j : le[e]) if (marca[j] != selo) { marca[j] = selo; out.push_back(j); }
    }
    return out;
  };

  std::vector<std::size_t> grau(n);
  for (std::size_t i = 0; i < n; i++) grau[i] = av[i].size();

  // LINHA DENSA SAI DO JOGO DE GRAUS (o "dense row handling" do AMD).
  //
  // Medido em dado real, 5.930 animais com 5.924 genotipados: o H^-1 sai 99,7% denso,
  // porque com genotipagem quase completa o A22^-1 e denso, e grau minimo num clique e o
  // PIOR caso desta estrutura de dados. Cada eliminacao chama viz() para cada membro de
  // L_p, e viz() custa O(k); com |L_p| = k dentro de k eliminacoes o custo e k^3 e nao nnz.
  // Foram 10 minutos parados aqui num ajuste que nem tinha chegado a primeira iteracao do
  // AI-REML, e a proxima chamada (a das MME, em aireml.cpp) e sobre uma matriz 2x maior.
  //
  // O conserto e nao deixar esses nos entrarem no jogo: quem passa do limiar vai direto
  // para o fim da ordem. Num clique todos passam, e o que sobra e exatamente a ordem que a
  // medicao de enchimento ja tinha apontado como otima aqui: elimina a parte esparsa por
  // grau minimo e empilha o clique no fim (dense_block = nucleo+1, enchimento ZERO).
  //
  // ORDENACAO NAO MUDA RESULTADO, so velocidade e enchimento: e semelhanca por permutacao.
  // O pior caso de errar o limiar e um ajuste lento, nunca um numero diferente. Isso e o
  // que torna esta troca barata de gatilhar e segura de errar.
  //
  // O LIMIAR E RELATIVO AO QUE AINDA ESTA VIVO, e nao a sqrt(n). Medido: com 10*sqrt(n), o
  // classico do AMD, a ordenacao ficou 8x a 13x mais rapida E PIOR, porque na APY o animal
  // JOVEM tem grau igual ao tamanho do nucleo (ele se liga a todo o nucleo) e passava do
  // limiar junto com o clique. O bloco denso saltou do nucleo para n_genotipados e o
  // enchimento subiu 7% a 12%. Isso troca custo de ordenacao, que e uma vez por ajuste, por
  // custo de fatoracao, que e TODA iteracao: com nucleo de 2.533 em 5.924 genotipados seriam
  // 12,8x mais flops por iteracao. Troca ruim.
  //
  // A assinatura da patologia nao e "grau alto", e "grau igual ao que sobrou": um clique que
  // abrange quase todos os vivos. Um no com grau 600 entre 1.200 vivos e denso e util; o
  // mesmo grau 600 entre 700 vivos e um clique e nao ha ordem que o salve.
  const double fracao = 0.8;
  auto denso_demais = [&](std::size_t d, std::size_t vivos) {
    return vivos > 64 && static_cast<double>(d) > fracao * static_cast<double>(vivos);
  };

  // NA PARTIDA, o corte relativo aos vivos chega tarde quando so parte dos animais e
  // genotipada. Medido na APY com 1.500 genotipados em 3.000: o nucleo (grau ~1.500, metade
  // dos vivos) so passava de 0,8 no passo 1.019, e ate la cada jovem eliminado atualizava
  // os 414 do nucleo, 2,7 s de ordenacao contra 0,13 s no caminho exato. O 10*sqrt(n) do
  // AMD pega o nucleo, mas sozinho tambem pega o JOVEM quando o nucleo e grande (o grau
  // dele e o tamanho do nucleo; ver acima). Os dois cortes juntos separam os dois casos:
  // o nucleo tem o grau do topo do grafo (liga-se a todo genotipado), o jovem tem o
  // tamanho do nucleo, e nos dois regimes medidos (nucleo 414 de 1.500, nucleo 2.533 de
  // 5.924) o jovem fica abaixo de 0,8 do topo. Num grafo esparso comum o corte absoluto
  // do AMD nunca e atingido e nada muda.
  //
  // O TOPO E O PERCENTIL 99 DO GRAU, e nao o maximo. Medido nas MME do mesmo caso: o
  // intercepto liga todo animal com registro (grau 2.700), virava o maximo, e 0,8 dele
  // passava do grau do nucleo; a ordenacao das MME ficava nos 2,5 s. Um punhado de
  // equacoes de efeito fixo nao pode definir o que e denso para os animais.
  std::size_t grau_topo = 0;
  {
    std::vector<std::size_t> g2(grau);
    const std::size_t pos = n - 1 - std::min(n - 1, n / 100);
    std::nth_element(g2.begin(), g2.begin() + pos, g2.end());
    grau_topo = g2[pos];
  }
  const double corte_abs = std::max(16.0, 10.0 * std::sqrt(static_cast<double>(n)));
  auto denso_na_partida = [&](std::size_t d) {
    return static_cast<double>(d) > corte_abs &&
           static_cast<double>(d) > fracao * static_cast<double>(grau_topo);
  };

  // A REGRA DE PARTIDA SO VALE QUANDO O CONJUNTO DENSO E PEQUENO. Tirar a linha densa do
  // grafo e o que o AMD faz, e e inocuo quando elas sao poucas: na APY o nucleo e 14% dos
  // nos e se liga a todo jovem, entao ignora-lo nao muda a ordem dos jovens. No passo unico
  // EXATO o clique genotipado e metade do grafo, e sem ele o grau minimo ordena os nao
  // genotipados sem saber quantos vizinhos genotipados cada um tem: medido, o ajuste exato
  // ficou ~1,5x mais lento por iteracao. Uma tentativa de devolver essa informacao contando
  // os vizinhos adiados na chave piorou a APY (bloco denso 1.012 -> 1.211) sem consertar o
  // exato, e foi descartada. Acima de 25% o conjunto fica com o corte relativo aos vivos,
  // que ja pega o clique quando os nao genotipados saem.
  {
    std::size_t na_partida = 0;
    for (std::size_t i = 0; i < n; i++) if (denso_na_partida(grau[i])) na_partida++;
    const bool vale_partida = na_partida <= n / 4;
    for (std::size_t i = 0; i < n; i++)
      if (denso_demais(grau[i], n) || (vale_partida && denso_na_partida(grau[i]))) {
        adiado[i] = 1; n_adiados++;
      }
  }
  std::vector<std::vector<std::uint32_t>> baldes(n + 1);
  for (std::size_t i = 0; i < n; i++)
    if (!adiado[i]) baldes[std::min(grau[i], n)].push_back(static_cast<std::uint32_t>(i));
  std::size_t lo = 0;

  // GRAU PREGUICOSO PARA OS HUBS. Medido no multicaracter com pedigree real de estrutura
  // comum, poucos pais para muitos filhos: a ordenacao era 98% do tempo de uma avaliacao e
  // crescia ~n^1.8 (0.06 s com 2 mil animais, 2.36 s com 16 mil), enquanto a Cholesky
  // numerica levava 0.014 s. A causa e o hub: um touro com centenas de filhos, ou a
  // equacao de um grupo contemporaneo ligada a todos os animais dele. A cada filho
  // eliminado, o touro estava em L_p, e o laco limpava a lista dele e recalculava o grau
  // dele com viz(), as duas coisas O(filhos). Sao n eliminacoes vezes O(n / touros).
  //
  // E o grau de um hub nao importa ate o fim: ele nunca e o minimo enquanto tem filhos.
  // Entao o no de lista grande so e MARCADO sujo em L_p, e a limpeza e o viz() acontecem
  // uma vez, quando ele e retirado do balde como candidato. No de lista pequena continua
  // sendo atualizado na hora, exatamente como antes, e onde nao ha hub a ordem sai igual.
  // Ordenacao nao muda resultado, so enchimento; o que se confere e que o enchimento nao
  // piorou.
  const std::size_t lista_grande = 16;
  // e so quando o elemento NOVO e pequeno. A preguica compensa porque cada filho do touro
  // traz um elemento de dois ou tres membros, e adiar a absorcao deles e barato. No nucleo
  // da APY cada jovem eliminado traz um elemento do tamanho do nucleo inteiro; adiado, isso
  // incha a lista ate o proximo recalculo, e medido, deixou a fatoracao da APY 2x MAIS
  // LENTA com enchimento identico. Com L_p grande vale o caminho antigo, exato.
  const std::size_t lp_pequeno = 32;
  std::vector<char> sujo(n, 0);
  auto limpa = [&](std::size_t i) {
    auto& a_i = av[i];
    a_i.erase(std::remove_if(a_i.begin(), a_i.end(),
                             [&](std::uint32_t x) { return !ativo(x); }), a_i.end());
    auto& e_i = ev[i];
    std::sort(e_i.begin(), e_i.end());
    e_i.erase(std::unique(e_i.begin(), e_i.end()), e_i.end());
  };

  while (ordem.size() < n) {
    std::size_t p = static_cast<std::size_t>(-1);
    while (lo <= n) {
      if (baldes[lo].empty()) { lo++; continue; }
      const std::uint32_t cand = baldes[lo].back();
      baldes[lo].pop_back();
      if (!vivo[cand] || adiado[cand] || std::min(grau[cand], n) != lo) continue;
      if (sujo[cand]) {
        // o hub chegou a frente com um grau velho: agora sim ele e limpo e medido, e volta
        // ao balde certo. Se o grau verdadeiro for menor que lo, a busca recomeca dali
        limpa(cand);
        sujo[cand] = 0;
        const std::size_t d = viz(cand).size();
        grau[cand] = d;
        if (denso_demais(d, n - ordem.size() - n_adiados)) { adiado[cand] = 1; n_adiados++; continue; }
        const std::size_t b = std::min(d, n);
        if (b != lo) {
          baldes[b].push_back(cand);
          if (b < lo) lo = b;
          continue;
        }
      }
      p = cand;
      break;
    }
    if (p == static_cast<std::size_t>(-1)) {          // o que sobrou esta isolado ou adiado
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
      // O MESMO DESVIO TEM DE ESTAR NOS DOIS LACOS, e a primeira versao so o pos no de
      // baixo. Medido: com a guarda so la, o nucleo era adiado certo na inicializacao, mas
      // cada eliminacao de um NAO-nucleo ainda arrastava os k do nucleo por aqui, e av[i] de
      // um no do clique tem k entradas. Da O(k^2) por eliminacao, 2.533 x 5.923 num caso
      // real, e o ganho caiu para ~2x em vez dos ~1000x que a conta pedia.
      //
      // Pular isto e seguro por um invariante: av[i] e ev[i] so sao LIDOS por viz(i), e
      // viz(i) nunca roda para um no adiado, porque ele nao entra em balde (nao vira pivo) e
      // o laco de baixo o desvia. Manter a adjacencia de quem nunca mais sera consultado e
      // trabalho jogado fora. A absorcao le no_lp[x] para x em le[e], e um le[e] antigo
      // ainda pode trazer um no adiado depois dele; por isso ela filtra por ativo(), e
      // o no adiado nao segura a morte de elemento nenhum.
      if (adiado[i]) continue;
      if (lp.size() <= lp_pequeno && (sujo[i] || av[i].size() + ev[i].size() > lista_grande)) {
        // hub: so registra o elemento novo; limpeza e grau ficam para quando ele sair do
        // balde. viz() ja filtra os mortos, entao a lista suja nao muda o grau medido
        ev[i].push_back(static_cast<std::uint32_t>(p));
        sujo[i] = 1;
        continue;
      }
      if (sujo[i]) { limpa(i); sujo[i] = 0; }
      auto& a_i = av[i];
      a_i.erase(std::remove_if(a_i.begin(), a_i.end(),
                               [&](std::uint32_t x) { return x == p || !ativo(x); }), a_i.end());
      auto& e_i = ev[i];
      e_i.erase(std::remove(e_i.begin(), e_i.end(), static_cast<std::uint32_t>(p)), e_i.end());
      e_i.push_back(static_cast<std::uint32_t>(p));
      // absorcao: um elemento cujos membros vivos ja estao dentro de L_p esta morto
      e_i.erase(std::remove_if(e_i.begin(), e_i.end(), [&](std::uint32_t e) {
        if (e == p) return false;
        compacta(e);
        for (std::uint32_t x : le[e]) if (no_lp[x] != selo_lp) return false;
        return true;
      }), e_i.end());
    }
    for (std::uint32_t iu : lp) {
      const std::size_t i = iu;
      if (adiado[i] || sujo[i]) continue;
      const std::size_t d = viz(i).size();
      grau[i] = d;
      // um no VIRA clique durante a eliminacao: o grau dele nao cresce, o numero de vivos e
      // que encolhe ate a razao passar do corte. E assim que o nucleo da APY e pego, depois
      // que os jovens ja sairam, e nao antes deles
      if (denso_demais(d, n - ordem.size() - n_adiados)) { adiado[i] = 1; n_adiados++; continue; }
      const std::size_t b = std::min(d, n);
      baldes[b].push_back(iu);
      if (b < lo) lo = b;
    }
  }
  return ordem;
}

}  // namespace br
