// As equacoes de modelo misto (Henderson, 1950), montadas esparsas, e a
// verossimilhanca restrita (Patterson e Thompson, 1971).
//
// A identidade que amarra tudo, e que os testes exigem dos dois lados:
//
//   -2logL = (n - p) log s2e + log|G*| + log|C_s| + y'Py
//          = (n - p) log s2e + log|M| + log|X'M^-1 X| + y'P*y / s2e ,   M = Z G* Z' + I
//
// O lado esquerdo vem das MME montadas aqui; o direito e a forma V, computada DENSA no
// modulo de referencia. Se os dois nao coincidem, um esta errado, e e assim que a montagem
// e validada sem nunca se conferir contra si mesma.
//
// W'W e acumulada um REGISTRO de cada vez. Cada registro contribui um pequeno clique entre
// as colunas que toca, entao o custo e a soma de nnz_linha^2/2 por registro, e nada
// quadratico no numero de colunas. A penalidade kron(C_g^-1, K^-1) entra depois, na mesma
// lista de triplos, e a montagem SOMAR duplicados e o que faz uma posicao tocada por
// registros e pela penalidade acumular em vez de um apagar o outro.

#include "mme.h"

#include <R_ext/Print.h>
#include <algorithm>
#include <cstring>
#include <cstdio>

#include <unordered_set>

namespace br {

// LINHA ESTRUTURALMENTE NULA numa K declarada: a linha inteira zero, diagonal inclusive.
// E o padrao da inversa generalizada das matrizes parciais multirraca (Mrode & Pocrnic,
// 4a ed., p.243-244): o nivel NAO contribui para este termo, nao ganha equacao, e os
// registros dele ficam na analise com incidencia zero AQUI, em vez de sair. A distincao
// com o id AUSENTE da K e deliberada: a linha zero e uma declaracao ("este animal nao
// carrega genes desta raca"), a ausencia e uma lacuna, e lacuna continua excluindo a linha
// como sempre. A K reduzida (so os ids nao nulos) e o que segue para a fatoracao e para os
// NIVEIS do termo, que e a parte que os espelhos nao faziam: sem ela o termo tirava os
// niveis da tabela, e uma K de 15 ids contra 12 niveis de dados virava "triplet outside
// the matrix" na montagem. Uma implementacao so, chamada pelos tres.
void reduz_kernels(const Modelo& m, const std::vector<KernelDecl>*& kernels,
                   std::vector<KernelDecl>& kern_red,
                   std::vector<std::unordered_set<std::string> >& kern_nulos) {
  if (!kernels) return;
  kern_red = *kernels;
  kern_nulos.assign(kern_red.size(), std::unordered_set<std::string>());
  for (std::size_t k = 0; k < kern_red.size(); k++) {
    KernelDecl& kd = kern_red[k];
    if (kd.vazia()) continue;
    const std::size_t nk = kd.ids.size();
    if (kd.k.nlin != nk || kd.k.ncol != nk) continue;  // shape falha adiante, com erro
    std::vector<char> nulo(nk, 0);
    std::size_t nz = 0;
    for (std::size_t i = 0; i < nk; i++) {
      bool zera = true;
      for (std::size_t j = 0; j < nk && zera; j++)
        if (kd.k.at(i, j) != 0.0) zera = false;
      if (zera) { nulo[i] = 1; nz++; }
    }
    for (std::size_t i = 0; i < nk; i++)
      if (!nulo[i] && kd.k.at(i, i) == 0.0)
        throw Erro("K of term '" + m.termos[k].nome + "': the row of '" + kd.ids[i] +
                   "' has a zero diagonal with nonzero covariances. A null "
                   "contribution is a whole row of zeros; anything else is not a "
                   "covariance matrix");
    if (nz == 0) continue;
    if (nz == nk)
      throw Erro("the K of term '" + m.termos[k].nome + "' is entirely zero: there is "
                 "no covariance to declare. Drop the term instead");
    std::vector<std::string> ids2;
    ids2.reserve(nk - nz);
    for (std::size_t i = 0; i < nk; i++) {
      if (nulo[i]) kern_nulos[k].insert(kd.ids[i]);
      else ids2.push_back(kd.ids[i]);
    }
    Densa k2(ids2.size(), ids2.size());
    for (std::size_t j = 0, jj = 0; j < nk; j++) {
      if (nulo[j]) continue;
      for (std::size_t i = 0, ii = 0; i < nk; i++) {
        if (nulo[i]) continue;
        k2.at(ii, jj) = kd.k.at(i, j);
        ii++;
      }
      jj++;
    }
    kd.ids = std::move(ids2);
    kd.k = std::move(k2);
  }
  kernels = &kern_red;
}

// O MESMO NUMERO ESCRITO DE DOIS JEITOS nao e um animal fora do pedigree. O R escreve o
// double 100000 como "1e+05" (as.character(), factor(), rownames<-), e o motor rotula a
// coluna numerica dos dados como "100000". Quando um lado passou por essa conversao e o
// outro nao, o nivel nao casava com o pedigree (ou com a K), e os registros daquele animal
// saiam da analise em silencio, so com o n_used menor. Aqui isso vira erro declarado, que diz
// como escrever o id, nos dois sentidos: dado em notacao cientifica com o inteiro entre os
// niveis, e nivel em notacao cientifica com o inteiro nos dados. So olha as linhas que nao
// casaram, entao nao custa nada quando tudo casa.
//
// Os ids conferidos sao os niveis do termo E, num kernel(), os das linhas NULAS da K, que
// reduz_kernels() tira dos niveis e guarda em kern_nulos: um registro de nivel nulo fica na
// analise com incidencia zero (casa_niveis_nulos), e um id nulo escrito "1e+05" na K contra
// 100000 nos dados nao casava com nada e saia calado, porque a conferencia so via os niveis.
static void confere_rotulos_cientificos(
    const Modelo& m, const std::vector<DesenhoTermo*>& aleatorios,
    const std::vector<std::unordered_set<std::string> >& kern_nulos, const Tabela& t,
    std::size_t nlin) {
  for (const DesenhoTermo* a : aleatorios) {
    const Termo& tm = m.termos[a->termo];
    if (tm.estrutura != Estrutura::Parentesco && tm.estrutura != Estrutura::Declarada) continue;
    bool algum = false;
    for (std::size_t i = 0; i < nlin && !algum; i++) algum = !a->casou[i];
    if (!algum) continue;
    std::vector<std::string> ids = a->niveis;
    if (tm.estrutura == Estrutura::Declarada && a->termo < kern_nulos.size())
      ids.insert(ids.end(), kern_nulos[a->termo].begin(), kern_nulos[a->termo].end());
    // a coluna que da o nivel da linha: a classe de aninhamento numa covariavel aninhada; no
    // termo social e no pai / avo materno o aninhado e a baia e o avo, e o nivel e a coluna
    const std::string& col =
        (!tm.aninhado.empty() && !tm.social && !tm.mgs) ? tm.aninhado : tm.coluna;
    const std::vector<std::string> rot = t.rotulos(col);
    const std::string onde = tm.estrutura == Estrutura::Parentesco ? "pedigree" : "K";
    auto recusa = [&](const std::string& nos_dados, const std::string& nos_niveis) {
      throw Erro("term '" + tm.nome + "': the data has the id '" + nos_dados + "' and the " +
                 onde + " has '" + nos_niveis + "', the same number written two ways ('" +
                 (nos_dados.find("e+") != std::string::npos ? nos_dados : nos_niveis) +
                 "' is how as.character() and factor() write a round number). Its records "
                 "would leave the analysis: write the ids the same way on both sides, as "
                 "numbers or with format(x, scientific = FALSE, trim = TRUE)");
    };
    // nas linhas, na ordem dos dados, para a mensagem apontar sempre o mesmo id
    std::unordered_set<std::string> fora, niveis;
    for (std::size_t i = 0; i < nlin; i++) {
      if (a->casou[i]) continue;
      fora.insert(rot[i]);
      const std::string c = inteiro_de_cientifico(rot[i]);
      if (c.empty()) continue;
      if (niveis.empty()) niveis.insert(ids.begin(), ids.end());
      if (niveis.count(c)) recusa(rot[i], c);
    }
    for (const std::string& nv : ids) {
      const std::string c = inteiro_de_cientifico(nv);
      if (!c.empty() && fora.count(c)) recusa(c, nv);
    }
  }
}

// Registro de um nivel declarado NULO na K: fica na analise, com incidencia zero neste
// termo. A linha zero da K diz que o efeito e exatamente zero, entao nao ha nada a somar e
// nada a excluir. Vale para os tres ajustadores, e e o ultimo passo do casamento de niveis:
// o que sobra sem casar passa pela conferencia dos ids em notacao cientifica.
void casa_niveis_nulos(const Modelo& m, const std::vector<DesenhoTermo*>& aleatorios,
                       const std::vector<std::unordered_set<std::string> >& kern_nulos,
                       const Tabela& t, std::size_t nlin) {
  for (DesenhoTermo* a : aleatorios) {
    if (m.termos[a->termo].estrutura != Estrutura::Declarada) continue;
    if (a->termo >= kern_nulos.size() || kern_nulos[a->termo].empty()) continue;
    const std::vector<std::string> rot = t.rotulos(m.termos[a->termo].coluna);
    for (std::size_t i = 0; i < nlin; i++)
      if (!a->casou[i] && kern_nulos[a->termo].count(rot[i])) a->casou[i] = 1;
  }
  confere_rotulos_cientificos(m, aleatorios, kern_nulos, t, nlin);
}

// Ordem dos niveis de um grupo sem estrutura de varios termos: os rotulos inteiros primeiro,
// em ordem NUMERICA ("9" antes de "10", como o dado ja renumerado espera), depois o resto em
// ordem de bytes. E uma ordem total sobre rotulos distintos: sinal, magnitude sem zeros a
// esquerda e, no empate ("007" e "7"), os proprios bytes.
static bool rotulo_inteiro(const std::string& s) {
  std::size_t i = (!s.empty() && s[0] == '-') ? 1 : 0;
  if (i == s.size()) return false;
  for (; i < s.size(); i++)
    if (s[i] < '0' || s[i] > '9') return false;
  return true;
}

static bool rotulo_antes(const std::string& a, const std::string& b) {
  const bool ia = rotulo_inteiro(a), ib = rotulo_inteiro(b);
  if (ia != ib) return ia;
  if (!ia) return a < b;
  const bool na = a[0] == '-', nb = b[0] == '-';
  if (na != nb) return na;
  auto magnitude = [](const std::string& s) {
    std::size_t i = (s[0] == '-') ? 1 : 0;
    while (i + 1 < s.size() && s[i] == '0') i++;
    return s.substr(i);
  };
  const std::string ma = magnitude(a), mb = magnitude(b);
  if (ma != mb) {
    const bool menor = ma.size() != mb.size() ? ma.size() < mb.size() : ma < mb;
    return na ? !menor : menor;
  }
  return a < b;
}

// Os termos ALEATORIOS do desenho, com o conjunto de niveis de cada um. Um lugar so para os
// tres montadores (model(), model_mt(), model_ar1(); gibbs() e snp_blup() passam por
// monta_desenho).
//
// Todo termo de um grupo tem de indexar O MESMO conjunto de niveis, na mesma ordem: a
// penalidade e kron(C_g^-1, K^-1) sobre a coluna coef * n_niveis + nivel, entao o nivel l de
// um termo so covaria com o nivel l do outro, e montagem, traco da AI, EM, Gibbs e os nomes
// dos EBV leem n_niveis do PRIMEIRO termo. No parentesco o conjunto e o pedigree; no kernel()
// sao os ids da K, e kinv_declarada ja exige a mesma K em todo termo do grupo. No grupo SEM
// estrutura (random(), pe()) cada termo tirava os niveis da propria coluna, na ordem de
// aparicao, e o nivel i de um pareava com o nivel i do outro, que e outro nivel: medido no
// caso de 20 + 20 niveis de tests/testthat/test-grupo-iid-niveis.R, o mesmo dado e o mesmo
// theta davam -2logL 340.573 na ordem original das linhas e 340.078 na permutacao daquele
// portao, contra 338.860 da MME densa pareada pelo nome; com 8 touros e 16 vacas a kron
// tinha a dimensao do primeiro termo e os EBV das vacas saiam com rotulos de touro. Agora o
// grupo usa a UNIAO dos rotulos dos seus termos, em rotulo_antes(), que nao depende da ordem
// das linhas (um termo aleatorio aninhado e recusado em monta_modelo, entao o rotulo e o da
// propria coluna). Um nivel que so aparece na coluna de um termo continua sendo efeito do
// outro, sem registro ali, so com a priori (e com a informacao que a covariancia do grupo
// traz do nivel irmao). O grupo de UM termo segue com a ordem de aparicao de sempre: sem par
// nao ha pareamento a errar.
std::vector<DesenhoTermo> monta_aleatorios(const Modelo& m, const Tabela& t,
                                           const std::vector<std::string>& niveis_ped,
                                           const std::vector<KernelDecl>* kernels) {
  std::vector<std::vector<std::string>> uniao(m.termos.size());
  std::vector<char> tem_uniao(m.termos.size(), 0);
  for (const Grupo& g : m.grupos) {
    if (g.estrutura != Estrutura::Diagonal || g.termos.size() < 2) continue;
    std::vector<std::string> niveis;
    std::unordered_set<std::string> visto;
    for (std::size_t tk : g.termos)
      for (const std::string& r : t.rotulos(m.termos[tk].coluna))
        if (visto.insert(r).second) niveis.push_back(r);
    std::sort(niveis.begin(), niveis.end(), rotulo_antes);
    for (std::size_t tk : g.termos) {
      uniao[tk] = niveis;
      tem_uniao[tk] = 1;
    }
  }

  std::vector<DesenhoTermo> out;
  for (std::size_t k = 0; k < m.termos.size(); k++) {
    if (!m.termos[k].aleatorio()) continue;
    const std::vector<std::string>* nf = nullptr;
    if (m.termos[k].estrutura == Estrutura::Parentesco) nf = &niveis_ped;
    else if (m.termos[k].estrutura == Estrutura::Declarada) nf = &(*kernels)[k].ids;
    else if (tem_uniao[k]) nf = &uniao[k];
    out.push_back(monta_termo(m, k, t, nf));
  }

  // O invariante que a montagem supoe, conferido aqui: quebrado, ele nao da erro em lugar
  // nenhum adiante, so estima outra coisa.
  for (const Grupo& g : m.grupos) {
    const DesenhoTermo* ref = nullptr;
    for (std::size_t tk : g.termos)
      for (const DesenhoTermo& a : out) {
        if (a.termo != tk) continue;
        if (!ref) ref = &a;
        else if (a.niveis != ref->niveis)
          throw Erro("group '" + g.nome + "': terms '" + ref->nome + "' (" +
                     std::to_string(ref->n_niveis) + " levels) and '" + a.nome + "' (" +
                     std::to_string(a.n_niveis) + " levels) do not index the same levels, "
                     "and the group covariance pairs them level by level");
      }
  }
  return out;
}

// A COVARIANCIA DE DOIS TERMOS IID SO EXISTE NA VEROSSIMILHANCA ATRAVES DOS NIVEIS COMUNS. Num
// grupo sem estrutura, Var(y) leva Z_a (C_ab I) Z_b', cuja entrada (i, j) e C_ab quando o
// nivel de a na linha i e o nivel de b na linha j sao o mesmo rotulo, e zero nos outros
// casos. Sem nenhum nivel com registro nos dois termos (touro e vaca com ids que nunca se
// repetem entre os sexos), C_ab some da verossimilhanca: ela e plana nessa direcao, e o
// REML andava as 300 iteracoes, saia com converged = FALSE e erro padrao NaN em todos os
// componentes. A recusa vem antes do ajuste, com o par de termos. Avaliar num theta dado
// (maxiter = 0, theta_fixed = do gibbs()) continua valendo: ali C_ab e um valor conhecido, e
// o BLUP de um nivel que so tem registro num termo usa essa covariancia. So conta linha que
// entra na analise (`usa`), com incidencia nao nula.
void confere_covariancias_iid(const Modelo& m, const std::vector<DesenhoTermo>& aleatorios,
                              const std::vector<char>& usa) {
  for (const Grupo& g : m.grupos) {
    if (g.estrutura != Estrutura::Diagonal || g.termos.size() < 2) continue;
    // com_registro[k][l]: o nivel l tem ao menos uma linha usada no k-esimo termo do grupo
    std::vector<std::vector<char>> com_registro;
    std::vector<const DesenhoTermo*> dt;
    for (std::size_t tk : g.termos)
      for (const DesenhoTermo& a : aleatorios) {
        if (a.termo != tk) continue;
        std::vector<char> tem(a.n_niveis, 0);
        for (std::size_t c = 0; c < a.z.ncol; c++)
          for (std::size_t q = a.z.colptr[c]; q < a.z.colptr[c + 1]; q++)
            if (a.z.valor[q] != 0.0 && a.z.linha[q] < usa.size() && usa[a.z.linha[q]])
              tem[c % a.n_niveis] = 1;
        com_registro.push_back(std::move(tem));
        dt.push_back(&a);
      }
    for (std::size_t p = 0; p < dt.size(); p++)
      for (std::size_t q = p + 1; q < dt.size(); q++) {
        bool comum = false;
        for (std::size_t l = 0; l < com_registro[p].size() && !comum; l++)
          comum = com_registro[p][l] && com_registro[q][l];
        if (comum) continue;
        const Termo& tp = m.termos[dt[p]->termo];
        const Termo& tq = m.termos[dt[q]->termo];
        throw Erro("group '" + g.nome + "': no level has records in both '" + tp.nome +
                   "' (column '" + tp.coluna + "') and '" + tq.nome + "' (column '" +
                   tq.coluna + "'). Their covariance pairs equal levels, and with no label "
                   "in common it does not enter the likelihood: it cannot be estimated, and "
                   "REML would end without convergence or standard errors. Put the two "
                   "terms in separate groups, or hold the components at known values "
                   "(start = with maxiter = 0 and n_em = 0 in model(), theta_fixed = in "
                   "gibbs())");
      }
  }
}


// Imprime o ESTADO DAS ESTIMATIVAS a cada iteracao, e nao so o tamanho do passo.
//
// POR QUE: um verbose que so diz "-2logL caiu, o passo encolheu" informa que o ajuste
// ANDA, nao PARA ONDE. As tres coisas que decidem se vale esperar ou matar a rodada sao
// uma variancia caminhando para zero, uma correlacao subindo para +/-1 e um componente
// explodindo, e as tres so aparecem no vetor. Sem imprimi-lo, elas so viram visiveis no
// fim, quando o tempo ja foi gasto.
//
// PARA QUE: em ajuste longo (multicaracter com muitos tracos, AR(1) em serie comprida,
// passo unico) a decisao de interromper e trocar start= vale horas.
//
// Quebra em varias linhas com recuo, para nao arruinar o terminal quando o modelo tem
// dezenas de componentes. %.4g porque a leitura aqui e de ORDEM DE GRANDEZA e de
// direcao de caminhada; o valor exato sai no fim, com erro-padrao.
// A K DECLARADA de um grupo, em UM lugar so.
//
// Estava escrita dentro de monta_desenho() e os dois espelhos (multicaracter e AR(1))
// recusavam kernel() por nao a terem. Duplicar o bloco resolveria os dois e criaria o
// problema seguinte: tres copias da mesma inversao densa, do mesmo teste de
// positividade e da mesma convencao de triangulo, divergindo com o tempo.
void kinv_declarada(const Modelo& mo, const Grupo& g,
                    const std::vector<KernelDecl>* kernels,
                    std::vector<Csc>& kinv, std::vector<double>& kinv_logdet) {
        // K DECLARADA (kernel): a matriz veio pronta do usuario, D de dominancia, G_AA de
        // epistasia, uma parcial por raca. A inversao aqui e DENSA de proposito: a K
        // declarada tem o tamanho do problema que o usuario montou, e a rota esparsa para D
        // de pedigree grande (Hoeschele & VanRaden 1991) segue por fazer.
        if (!kernels)
          throw Erro("a kernel() term declares its own covariance matrix, and this fitting "
                     "route does not carry it: kernel() reaches model(), model_mt() and "
                     "model_ar1() with their eval_internal() companions and gibbs(), but not snp_blup(), "
                     "whose marker equations have no place for a declared K");
        // O nome do termo so serve para a mensagem, e le-lo de um Modelo que o chamador
        // ainda nao preencheu foi um acesso fora de faixa que so aparecia no caminho de
        // ERRO: com K boa a linha nunca era executada. Por isso o indice e conferido.
        auto nome_do = [&](std::size_t tk) {
          return tk < mo.termos.size() ? mo.termos[tk].nome : std::string("kernel");
        };
        if (g.termos.empty())
          throw Erro("group '" + g.nome + "' is declared kernel and carries no term");
        const KernelDecl* kd = nullptr;
        for (std::size_t tk : g.termos) {
          if (tk >= kernels->size() || (*kernels)[tk].vazia())
            throw Erro("term '" + nome_do(tk) +
                       "' is declared kernel and no K reached the engine");
          const KernelDecl& c = (*kernels)[tk];
          if (!kd) kd = &c;
          else if (kd->ids != c.ids || kd->k.dados != c.k.dados)
            throw Erro("group '" + g.nome + "' has kernel() terms with different K: the "
                       "penalty of a group is a single kron(C, K); give the terms the same "
                       "K or separate groups");
        }
        const std::size_t nk = kd->ids.size();
        if (kd->k.nlin != nk || kd->k.ncol != nk)
          throw Erro("kernel K of " + std::to_string(kd->k.nlin) + " x " +
                     std::to_string(kd->k.ncol) + " for " + std::to_string(nk) + " ids");
        // Positiva-definida com FOLGA, nao so por passar na Cholesky: uma K singular (a G crua
        // de mais animais que marcadores independentes) tem o menor pivo em zero a menos do
        // arredondamento, e o sinal desse arredondamento muda com a plataforma e com a ordem
        // de soma de quem montou a K (medido: a mesma G de posto 14 em 15 recusada no Windows
        // e aceita no Linux). O corte e relativo a escala: pivo^2 < 1e-12 da maior diagonal e
        // singular para todos os efeitos, a inversa perderia doze digitos.
        Densa lk = kd->k;
        const bool pd = chol_densa(lk);
        double ldk = 0.0, piv_min = std::numeric_limits<double>::infinity(), dmax = 0.0;
        for (std::size_t i = 0; pd && i < nk; i++) {
          ldk += 2.0 * std::log(lk.at(i, i));
          piv_min = std::min(piv_min, lk.at(i, i) * lk.at(i, i));
          dmax = std::max(dmax, kd->k.at(i, i));
        }
        if (!pd || !(piv_min > 1e-12 * dmax))
          throw Erro("the K of term '" + nome_do(g.termos[0]) + "' is not "
                     "positive-definite" + std::string(pd ? " (numerically singular)" : "") +
                     ": if it comes from markers, add a small ridge to "
                     "the diagonal (K + 0.01 I) before the call");
        Densa ki = inv_pd(kd->k);
        // para a Csc do triangulo inferior, a MESMA convencao do A^-1 do pedigree
        std::vector<std::uint32_t> li, cj;
        std::vector<double> v;
        for (std::size_t j = 0; j < nk; j++)
          for (std::size_t i = j; i < nk; i++)
            if (ki.at(i, j) != 0.0) {
              li.push_back(static_cast<std::uint32_t>(i));
              cj.push_back(static_cast<std::uint32_t>(j));
              v.push_back(ki.at(i, j));
            }
        kinv.push_back(de_triplos(nk, nk, li, cj, v));
        // o campo guarda log|K^-1|, e |K^-1| = 1/|K|
        kinv_logdet.push_back(-ldk);
}

void imprime_theta(const std::vector<double>& th, const std::vector<std::string>& nomes,
                   const Modelo& m) {
  // UMA LINHA POR TERMO, e nao uma fila corrida de componentes. Num grupo de covariancia
  // os parametros sao o vech de UMA matriz (direto-materno, ou o bloco aditivo de um
  // multicaracter), e imprimi-los emendados com os de outro termo obriga quem le a
  // reconstruir a mao quais numeros formam qual matriz. Agrupando por termo, a matriz do
  // aditivo, a do ambiente permanente e a do residual sao lidas separadamente, que e a
  // leitura que decide se o ajuste vai bem: a diagonal diz onde cada variancia esta, e o
  // fora da diagonal diz se uma correlacao esta caminhando para +/-1.
  auto emite = [&](const std::string& rotulo, std::size_t de, std::size_t ate) {
    if (de >= ate || de >= th.size()) return;
    std::string linha = "    " + rotulo;
    while (linha.size() < 15) linha += " ";
    const std::size_t recuo = 15;
    for (std::size_t k = de; k < ate && k < th.size(); k++) {
      char buf[96];
      const char* nome = (k < nomes.size()) ? nomes[k].c_str() : "?";
      std::snprintf(buf, sizeof(buf), "%s=%.4g", nome, th[k]);
      if (linha.size() > recuo && linha.size() + std::strlen(buf) + 2 > 78) {
        Rprintf("%s\n", linha.c_str());
        linha.assign(recuo, ' ');
      }
      if (linha.size() > recuo) linha += "  ";
      linha += buf;
    }
    if (linha.size() > recuo) Rprintf("%s\n", linha.c_str());
  };

  for (const Grupo& g : m.grupos) {
    std::string rot = g.nome;
    if (rot.empty() && !g.termos.empty() && g.termos[0] < m.termos.size())
      rot = m.termos[g.termos[0]].nome;
    if (rot.empty()) rot = "grupo";
    emite("[" + rot + "]", g.offset, g.offset + g.nparam);
  }
  // Do offset do residual ate o fim: no multicaracter o residual e uma MATRIZ entre
  // caracteristicas, e no AR(1) o rho mora depois dela. Os dois pertencem a mesma leitura.
  emite("[residual]", m.offset_residual, th.size());
}

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

// Montagem das MME. W'W nao depende de theta nem de y e, com um cache (o ajuste AI-REML e a
// cadeia de Gibbs passam o seu), e montada UMA vez com zeros explicitos nas posicoes da
// penalidade, junto com a posicao de cada entrada da penalidade no padrao final; nas
// avaliacoes seguintes so os valores da penalidade sao somados no lugar. W'y e refeito a cada
// chamada, coluna a coluna: na cadeia probit y e a liabilidade, sorteada de novo a cada
// iteracao NO MESMO desenho (guardar W'y deixava a cadeia com o y da primeira iteracao, e o
// portao do probit pegou). O resultado e bit a bit o da montagem por triplos: cada posicao
// superior recebe UMA entrada de penalidade, (soma dos dados + 0) + penalidade e a mesma soma
// na mesma ordem, e cada entrada de W'y soma os registros na mesma ordem crescente.
Montado monta_mme(const Desenho& d, const std::vector<double>& theta, CacheSimbolica* cache) {
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
  const bool pronto = cache && cache->mont_dono == &d && cache->mont_base.ncol == M.total;

  std::vector<std::uint32_t> ti, tj;
  std::vector<double> tv;
  if (pronto) {
    M.rhs.assign(M.total, 0.0);
    for (std::size_t j = 0; j < d.x.ncol; j++) {
      double s = 0.0;
      for (std::size_t r = 0; r < d.nlin; r++) {
        if (!d.usa[r]) continue;
        const double v = d.x.at(r, j);
        if (v != 0.0) s += v * d.y[r];
      }
      M.rhs[j] = s;
    }
    for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
      std::size_t col = M.offset_grupo[g];
      for (std::size_t t : d.modelo.grupos[g].termos)
        for (const DesenhoTermo& a : d.aleatorios)
          if (a.termo == t) {
            for (std::size_t c = 0; c < a.z.ncol; c++) {
              double s = 0.0;
              for (std::size_t p = a.z.colptr[c]; p < a.z.colptr[c + 1]; p++)
                if (d.usa[a.z.linha[p]]) s += a.z.valor[p] * d.y[a.z.linha[p]];
              M.rhs[col + c] = s;
            }
            col += a.z.ncol;
          }
    }
  } else {
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
  }

  // os fatores da penalidade por grupo, e o inadmissivel antes de emitir qualquer coisa
  std::vector<Densa> cinvs(d.modelo.grupos.size());
  std::vector<std::size_t> nls(d.modelo.grupos.size(), 0);
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    Densa cg = cov_grupo(d.modelo, theta, g);
    // escala por s2e: C_s = W'W + kron(C_g^-1 s2e, K^-1) e o sistema em unidades de s2e
    Densa cgs = cg;
    for (double& v : cgs.dados) v /= M.s2e;
    try { cinvs[g] = inv_pd(cgs); } catch (const Erro&) { return M; }
    const double ld = logdet_pd(cgs);
    if (std::isnan(ld)) return M;

    const std::size_t nl = d.aleatorios.empty() ? 0 :
        [&]{ for (const auto& a : d.aleatorios)
               if (a.termo == d.modelo.grupos[g].termos[0]) return a.n_niveis;
             return static_cast<std::size_t>(0); }();
    nls[g] = nl;
    M.logdet_g += static_cast<double>(nl) * ld - static_cast<double>(d.modelo.grupos[g].dim) * d.kinv_logdet[g];
  }

  // A penalidade, emitida sempre na MESMA ordem: f(i, j, valor) para cada posicao do
  // triangulo superior. f == 0 NAO pula: a covariancia que comeca em zero ganharia um padrao
  // menor na primeira avaliacao, e o cache da simbolica congelaria esse padrao errado para o
  // ajuste inteiro. Zero explicito ocupa o slot e nao muda numero. EMITE A MATRIZ CHEIA E
  // MANTEM SO gi <= gj, SEM TROCAR: a primeira versao trocava (swap) e espelhava, e com o
  // laco percorrendo a e b completos cada posicao fora da diagonal era emitida DUAS vezes; a
  // montagem soma duplicados, e a penalidade saia dobrada (ainda simetrica e definida,
  // converge, para o lugar errado). Mantendo so o triangulo sem trocar, a varredura completa
  // de (a,b) emite cada posicao superior da kron cheia exatamente uma vez.
  auto emite = [&](auto&& f) {
    for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
      const Densa& cinv = cinvs[g];
      const std::size_t nl = nls[g], off = M.offset_grupo[g];
      const bool com_k = d.kinv[g].ncol > 0;
      for (std::size_t a = 0; a < d.modelo.grupos[g].dim; a++)
        for (std::size_t b = 0; b < d.modelo.grupos[g].dim; b++) {
          const double fab = cinv.at(a, b);
          auto poe = [&](std::size_t gi, std::size_t gj, double val) {
            if (gi > gj) return;
            f(gi, gj, val);
          };
          if (com_k) {
            const Csc& k = d.kinv[g];
            for (std::size_t col = 0; col < k.ncol; col++)
              for (std::size_t p = k.colptr[col]; p < k.colptr[col + 1]; p++) {
                const std::size_t rk = k.linha[p];
                // K^-1 vem no triangulo inferior: a matriz cheia tem (rk,col) e (col,rk)
                poe(off + a * nl + rk, off + b * nl + col, fab * k.valor[p]);
                if (rk != col) poe(off + a * nl + col, off + b * nl + rk, fab * k.valor[p]);
              }
          } else {
            for (std::size_t l = 0; l < nl; l++)
              poe(off + a * nl + l, off + b * nl + l, fab);
          }
        }
    }
  };

  if (!cache) {
    emite([&](std::size_t i, std::size_t j, double v) {
      ti.push_back(static_cast<std::uint32_t>(i));
      tj.push_back(static_cast<std::uint32_t>(j));
      tv.push_back(v);
    });
    M.c = de_triplos(M.total, M.total, ti, tj, tv);
    M.ok = true;
    return M;
  }
  if (!pronto) {
    emite([&](std::size_t i, std::size_t j, double) {
      ti.push_back(static_cast<std::uint32_t>(i));
      tj.push_back(static_cast<std::uint32_t>(j));
      tv.push_back(0.0);
    });
    cache->mont_base = de_triplos(M.total, M.total, ti, tj, tv);
    std::vector<std::uint32_t>().swap(ti);
    std::vector<std::uint32_t>().swap(tj);
    std::vector<double>().swap(tv);
    const Csc& b = cache->mont_base;
    cache->mont_pos.clear();
    emite([&](std::size_t i, std::size_t j, double) {
      const auto ini = b.linha.begin() + static_cast<std::ptrdiff_t>(b.colptr[j]);
      const auto fim = b.linha.begin() + static_cast<std::ptrdiff_t>(b.colptr[j + 1]);
      const auto it = std::lower_bound(ini, fim, static_cast<std::uint32_t>(i));
      cache->mont_pos.push_back(static_cast<std::size_t>(it - b.linha.begin()));
    });
    cache->mont_dono = &d;
  }
  M.c = cache->mont_base;
  std::size_t idx = 0;
  emite([&](std::size_t, std::size_t, double v) { M.c.valor[cache->mont_pos[idx++]] += v; });
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
  // log|K|, nao log|K^-1|, a diferenca e constante em theta, entao o OTIMO nao se move e so
  // o VALOR sai errado, o que arruina qualquer comparacao com outro programa sem arruinar as
  // estimativas. Invisivel num teste de recuperacao, fatal num de -2logL.
  r.valor = (n - p) * std::log(M.s2e) + M.logdet_g + r.logdet_c + (yy - bry) / M.s2e
            - d.logdet_peso;
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

Desenho monta_desenho(const Modelo& m, const Tabela& t, const Pedigree* ped,
                      const std::vector<double>* peso,
                      const std::vector<KernelDecl>* kernels) {
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

  // a reducao da K por linha nula e os niveis que ela define: ver reduz_kernels()
  std::vector<KernelDecl> kern_red;
  std::vector<std::unordered_set<std::string> > kern_nulos;
  reduz_kernels(m, kernels, kern_red, kern_nulos);

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
    } else if (g.estrutura == Estrutura::Declarada) {
      kinv_declarada(m, g, kernels, d.kinv, d.kinv_logdet);
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

  // aleatorios, com niveis do pedigree quando ha parentesco, da K quando declarada e da
  // uniao das colunas num grupo sem estrutura de varios termos: todo nivel do conjunto
  // ganha equacao, com ou sem registro. Ver monta_aleatorios().
  d.aleatorios = monta_aleatorios(m, t, niveis_ped, kernels);

  {
    std::vector<DesenhoTermo*> pa;
    for (DesenhoTermo& a : d.aleatorios) pa.push_back(&a);
    casa_niveis_nulos(m, pa, kern_nulos, t, d.nlin);
  }

  // linhas que entram: nem ausente, nem nivel sem casar no parentesco.
  //
  // NA/NaN na observacao E ausente, com ou sem codigo declarado: e o idioma do R, os
  // caminhos multi e AR ja tratavam assim, e o univariado nao, um NA atravessava a
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

  // PESOS, como escala de linha por sqrt(w). Um registro de peso w tem residual s2e/w;
  // multiplicar a linha inteira (y, X e cada Z) por sqrt(w) transforma o modelo
  // ponderado no modelo homocedastico das mesmas equacoes normais. A verossimilhanca do
  // modelo ORIGINAL difere da escalada pelo jacobiano soma(log w), constante em theta:
  // nao move o otimo, mas e guardada para o -2logL sair no valor certo.
  if (peso) {
    if (peso->size() != d.nlin) throw Erro("weights with the wrong length");
    for (std::size_t i = 0; i < d.nlin; i++) {
      if (!d.usa[i]) continue;
      const double w = (*peso)[i];
      if (!(w > 0.0) || !std::isfinite(w))
        throw Erro("weight that is not finite and positive at row " + std::to_string(i + 1));
      d.logdet_peso += std::log(w);
    }
    for (std::size_t i = 0; i < d.nlin; i++) {
      const double r = d.usa[i] ? std::sqrt((*peso)[i]) : 1.0;
      if (r == 1.0) continue;
      d.y[i] *= r;
      for (std::size_t j = 0; j < xfull.ncol; j++) xfull.at(i, j) *= r;
      for (DesenhoTermo& a : d.aleatorios)
        for (std::size_t c = 0; c < a.z.ncol; c++)
          for (std::size_t k = a.z.colptr[c]; k < a.z.colptr[c + 1]; k++)
            if (a.z.linha[k] == i) a.z.valor[k] *= r;
    }
  }

  // POSTO DE X SOBRE AS LINHAS QUE ENTRAM, nao sobre a tabela inteira.
  //
  // Um nivel fixo cujos registros TODOS sairam (ausentes, ou nivel sem par no pedigree)
  // continua com coluna nao-nula na tabela, entao um posto medido sobre tudo o mantinha;
  // na montagem, que so anda nas linhas usadas, essa coluna nao recebe nada e a matriz
  // de coeficientes ganha uma coluna VAZIA. Dali o efeito era em cascata: a Cholesky
  // acha diagonal zero e devolve "nao positiva-definida", o ajuste para na primeira
  // avaliacao sem theta, e a travessia ainda tentava montar as MME com esse theta vazio.
  // Medindo o posto onde o modelo de fato vive, a coluna sai como dependente e e
  // REPORTADA em dropped_x, que e o comportamento correto e visivel.
  std::vector<std::size_t> usadas;
  usadas.reserve(d.n_usadas());
  for (std::size_t i = 0; i < d.nlin; i++) if (d.usa[i]) usadas.push_back(i);
  Densa xusadas(usadas.size(), cols.size());
  for (std::size_t r = 0; r < usadas.size(); r++)
    for (std::size_t j = 0; j < cols.size(); j++) xusadas.at(r, j) = xfull.at(usadas[r], j);
  std::vector<std::size_t> fica, sai;
  posto_completo(xusadas, 1e-9, fica, sai);
  d.x = Densa(d.nlin, fica.size());
  for (std::size_t jj = 0; jj < fica.size(); jj++) {
    d.nomes_x.push_back(cols[fica[jj]].first);
    for (std::size_t i = 0; i < d.nlin; i++) d.x.at(i, jj) = xfull.at(i, fica[jj]);
  }
  for (std::size_t j : sai) d.saiu_x.push_back(cols[j].first);
  return d;
}

}  // namespace br
