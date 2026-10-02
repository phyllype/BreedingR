// Inversa esparsa da matriz de dominancia por subclasses pai x mae: Hoeschele e VanRaden
// (1991), J Dairy Sci 74:557-569, na forma geral, com endogamia, geracoes sobrepostas, pais
// desconhecidos e casais repetidos.
//
// O MODELO AUMENTADO. O desvio de dominancia do animal i com pai S e mae D conhecidos e
// d_i = h_SD + delta_i: h_SD e o efeito da subclasse S x D (a dominancia media de infinitos
// irmaos completos) e delta_i o desvio dentro dela, independente. Var(h) = (F / 4) s2d, com
// F_{SD,XY} = a_SX a_DY + a_SY a_DX, e Var(delta_i) = Delta_i s2d,
//   Delta_i = 1 - F_{SD,SD} / 4 = 1 - [(1 + F_S)(1 + F_D) + 4 F_i^2] / 4   (a_SD = 2 F_i).
// Sem os dois pais nao ha subclasse e Delta_i = 1 (Cockerham zera a linha). A identidade
// D = W (F/4) W' + diag(Delta) vale para a D de Cockerham de dominance_matrix(), com
// endogamia inclusive; e a MESMA D, nao uma correcao dela (a D de Cockerham com diagonal 1
// continua sendo a aproximacao classica sob endogamia).
//
// A PRECISAO de [d; h] e
//   Q = [[Delta^-1, -Delta^-1 W], [-W' Delta^-1, 4 F^-1 + W' Delta^-1 W]],
// log|K| = sum log Delta + log|F| - n_h log 4, e Var([d; h]) = s2d Q^-1: o componente e o
// mesmo s2d de kernel(id, K = D), e -2logL, score e AI sao identicos aos daquela rota.
//
// QUANDO NAO EXISTE: Delta_i <= 0. Com autofecundacao isso acontece ja na segunda geracao
// (Delta = -0.125) e com irmao x irma na sexta (pais com F = 0,59), enquanto a D de Cockerham
// do mesmo pedigree ainda pode ser positiva-definida. A D existe; a representacao por
// subclasse nao. A rota recusa, e quando o animal nao tem registro basta deixa-lo fora de
// Q (animals = os animais com registro): o Delta dele nao entra em lugar nenhum.
//
// DUAS ROTAS para F^-1, ambas exatas para esta D:
//  - densa: F das subclasses cheias pela A entre os pais (Colleau), invertida densa. Custo
//    n_sub^3 / 3 e memoria n_sub^2. Ganha com poucas subclasses de muitos filhos (leitegada).
//  - esparsa: recorrencia de pares. Com f_xy = (g_xy + g_yx) / sqrt(2), g ~ N(0, A (x) A), a
//    linha (x, y) de (I - P) (x) (I - P) e a eq. [5] do artigo,
//      f_xy = .5 sum_{p in pais(x)} f_py + .5 sum_{q in pais(y)} f_xq - .25 sum_{p,q} f_pq + e_xy,
//    com Var(e_xy) = d_x d_y (x != y) e 2 d_x^2 no par proprio (x = y), d_x o coeficiente
//    mendeliano da A ja com a endogamia dos pais. Termos repetidos se acumulam (no par proprio,
//    +1 em {p, x}, -.5 em {s_x, d_x}, -.25 em {s_x, s_x} e {d_x, d_x}). O fecho ancestral dos
//    pares cheios e uma marginal exata (conjunto ancestralmente fechado de uma rede gaussiana
//    em DAG). A poda tambem e exata: um par sem animal e com no maximo UM filho e integrado
//    na linha do filho (q_ck += q_cu q_uk, m_c += q_cu^2 m_u). Ganha quando cada subclasse
//    tem poucos filhos (estrutura leiteira). A regra substitui a poda heuristica do artigo.
//
// A ESCOLHA AUTOMATICA olha a fatoracao simbolica das MME CONJUNTAS de um modelo animal +
// dominancia com um registro por animal de Q (A^-1, Q e o acoplamento a_i - d_i), e nao o
// bloco de dominancia isolado: soma de c_j^2 (flops) de cada rota, com a parte ESPARSA do
// fator pesando PESO_ESPARSO vezes a cauda densa. O peso e uma constante da regra: a
// cholesky() do pacote fatora a cauda densa final em ladrilhos com bloco 4 x 4 em
// registradores e o resto coluna a coluna, escalar e com acesso indireto, varias vezes mais
// devagar por operacao. Contar flops sem peso escolheu a esparsa num pedigree leiteiro de
// 7.216 animais (3,3e10 contra 5,6e10) que avaliou em 45 s contra 18 s da densa. O peso nao
// depende do numero de threads, de proposito: a rota muda os niveis do termo, e o resultado
// nao pode depender das threads.
//
// O CUSTO DA PROPRIA DECISAO e limitado pela rota densa, que e avaliada primeiro e e barata
// de avaliar (o bloco h e um clique de contagem fechada). Tres saidas antecipadas levam a
// densa sem ordenar a esparsa inteira:
//  - o fecho passa das n_sub (n_sub + 1) / 2 entradas do bloco denso ou da memoria de pico da
//    densa (BYTES_POR_PAR por par);
//  - o bloco de pares da Q esparsa seria montado com mais triplos que o triangulo inteiro do
//    bloco denso. E a mesma contagem dos dois lados (as duas rotas montam Q em triplos, e a
//    esparsa os cria com repeticao), so do bloco h: a A^-1, a parte d e o acoplamento entram
//    iguais nas duas candidatas e ficam fora da comparacao. Comparar as MME inteiras da
//    esparsa (com a A^-1) contra o pico de montar a Q densa (sem ela) empurrava para a densa
//    todo pedigree grande com poucas subclasses;
//  - a ORDENACAO da candidata esparsa passa de um teto de TRABALHO (entradas de lista lidas
//    pelo grau minimo, uma contagem deterministica): o trabalho da ordenacao da candidata
//    densa mais o equivalente ao tempo de UMA fatoracao da densa. Medido em pedigrees de
//    leitegadas, o grau minimo exato numa candidata de muito enchimento andava a ~3e7 c_j^2
//    por segundo, 200 vezes mais devagar que a fatoracao em ladrilhos: com o teto so em c_j^2
//    (o custo da densa) a decisao levava 40 a 350 s para escolher uma densa que se monta em 1
//    a 3 s. Uma candidata cuja ordenacao custa mais que a da densa e uma fatoracao da densa
//    perde a comparacao por ela, porque o ajuste paga essa mesma ordenacao antes da primeira
//    iteracao; a decisao custa entao no maximo cerca de uma fatoracao da densa alem das
//    duas ordenacoes, mais o fecho. O trabalho da densa entra no teto pelo mesmo motivo da
//    saida anterior: a ordenacao da A^-1 e paga pelas duas rotas e so o excesso conta. A
//    regra e uma heuristica: num pedigree de 60.100 animais com 800 subclasses (cinco
//    geracoes de um filho por casal) a ordenacao da esparsa passou do teto e a densa
//    escolhida custou 1,9 vezes a esparsa num ajuste (validation/dominance_hv91_routes.R);
//    um teto de varias fatoracoes talvez acertasse ali (nao testado), e deixaria a decisao
//    mais lenta nas leitegadas.
// A ordenacao tambem para quando a soma parcial de c_j^2 passa do custo da densa (custo >=
// flops, a escolha nao muda).
//
// Tudo aqui e sequencial salvo a A entre os pais (Colleau, uma coluna por thread, dono
// unico) e a Cholesky/inversa em ladrilhos: o resultado nao depende do numero de threads.

#include "mme.h"
#include <unordered_map>
#include <queue>
#include <cmath>
#include <algorithm>
#include <limits>
#include <chrono>
#include <cstdio>

namespace br {

namespace {

// chave de um par nao ordenado de linhas do pedigree topologico (base 0): maior * n + menor.
// Os pares pais de um par tem sempre chave MENOR, entao a ordem das chaves e topologica.
// int64: com n acima de 46.340 animais a chave passa de 2^31.
inline std::int64_t chave(std::int64_t a, std::int64_t b, std::int64_t n) {
  return a > b ? a * n + b : b * n + a;
}

struct NoPar {
  std::int64_t x = 0, y = 0;                         // x <= y
  bool cheio = false, vivo = true, feito = false;
  double m = 0.0;                                    // Var(e) do par, em unidades de F
  std::vector<std::pair<std::uint32_t, double> > q;  // pais: (no, coeficiente) em f = sum q f + e
  std::vector<std::uint32_t> filhos;
};

struct Fecho {
  std::vector<NoPar> nos;
  std::unordered_map<std::int64_t, std::uint32_t> indice;
  std::size_t processados = 0;
  bool abortado = false;
};

void tira(std::vector<std::uint32_t>& v, std::uint32_t x) {
  for (std::size_t k = 0; k < v.size(); k++)
    if (v[k] == x) { v[k] = v.back(); v.pop_back(); return; }
}

// O fecho de pares com eliminacao exata, processado por chave decrescente: quando um par
// sai da fila, todos os filhos dele ja foram processados (tem chave maior), entao o
// conjunto de filhos e final e a regra "no maximo um filho" e decidida com certeza.
Fecho fecho_pares(const Pedigree& p, const std::vector<double>& dx,
                  const std::vector<std::pair<std::int64_t, std::int64_t> >& cheios,
                  std::size_t limite) {
  const std::int64_t n = static_cast<std::int64_t>(p.ids.size());
  Fecho fe;
  std::priority_queue<std::int64_t> fila;
  auto pega = [&](std::int64_t a, std::int64_t b) -> std::uint32_t {
    const std::int64_t k = chave(a, b, n);
    auto it = fe.indice.find(k);
    if (it != fe.indice.end()) return it->second;
    if (fe.nos.size() >= static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
      throw Erro("the pair closure outgrew 32-bit node indices");
    const std::uint32_t j = static_cast<std::uint32_t>(fe.nos.size());
    NoPar no;
    no.x = std::min(a, b);
    no.y = std::max(a, b);
    fe.nos.push_back(std::move(no));
    fe.indice.emplace(k, j);
    fila.push(k);
    return j;
  };
  for (const auto& c : cheios) fe.nos[pega(c.first, c.second)].cheio = true;

  std::vector<std::pair<std::int64_t, double> > qk;
  while (!fila.empty()) {
    const std::int64_t k = fila.top();
    fila.pop();
    const std::uint32_t u = fe.indice.find(k)->second;
    if (fe.nos[u].feito) continue;
    fe.nos[u].feito = true;
    fe.processados++;
    if ((fe.processados & 0xFFFF) == 0) checa_interrupcao();
    const std::int64_t x = fe.nos[u].x, y = fe.nos[u].y;
    const std::int64_t px[2] = {p.pai[x], p.mae[x]}, py[2] = {p.pai[y], p.mae[y]};
    // os coeficientes do par, ACUMULADOS por chave: no par proprio e na autofecundacao
    // varios termos caem no mesmo par e se somam
    qk.clear();
    auto soma = [&](std::int64_t a, std::int64_t b, double c) {
      const std::int64_t kk = chave(a, b, n);
      for (auto& e : qk)
        if (e.first == kk) { e.second += c; return; }
      qk.push_back({kk, c});
    };
    for (int a = 0; a < 2; a++) if (px[a] >= 0) soma(px[a], y, 0.5);
    for (int b = 0; b < 2; b++) if (py[b] >= 0) soma(x, py[b], 0.5);
    for (int a = 0; a < 2; a++)
      for (int b = 0; b < 2; b++)
        if (px[a] >= 0 && py[b] >= 0) soma(px[a], py[b], -0.25);
    // ordem fixa de criacao dos nos novos: a numeracao interna nao depende de nada externo
    std::sort(qk.begin(), qk.end());
    for (const auto& e : qk) {
      if (e.second == 0.0) continue;
      const std::uint32_t j = pega(e.first / n, e.first % n);
      fe.nos[u].q.push_back({j, e.second});
      fe.nos[j].filhos.push_back(u);
    }
    fe.nos[u].m = (x == y) ? 2.0 * dx[x] * dx[x] : dx[x] * dx[y];
    if (fe.nos.size() > limite) {
      fe.abortado = true;
      return fe;
    }
    // ELIMINACAO EXATA de um par sem animal com no maximo um filho
    if (fe.nos[u].cheio || fe.nos[u].filhos.size() > 1) continue;
    fe.nos[u].vivo = false;
    for (const auto& pk : fe.nos[u].q) tira(fe.nos[pk.first].filhos, u);
    if (!fe.nos[u].filhos.empty()) {
      const std::uint32_t c = fe.nos[u].filhos[0];
      NoPar& nc = fe.nos[c];
      double b = 0.0;
      for (std::size_t t = 0; t < nc.q.size(); t++)
        if (nc.q[t].first == u) {
          b = nc.q[t].second;
          nc.q.erase(nc.q.begin() + static_cast<std::ptrdiff_t>(t));
          break;
        }
      for (const auto& pk : fe.nos[u].q) {
        bool achou = false;
        for (std::size_t t = 0; t < nc.q.size(); t++) {
          if (nc.q[t].first != pk.first) continue;
          achou = true;
          nc.q[t].second += b * pk.second;
          // coeficientes diadicos cancelam a zero exato; o zero sai do padrao
          if (nc.q[t].second == 0.0) {
            nc.q.erase(nc.q.begin() + static_cast<std::ptrdiff_t>(t));
            tira(fe.nos[pk.first].filhos, c);
          }
          break;
        }
        if (!achou) {
          nc.q.push_back({pk.first, b * pk.second});
          fe.nos[pk.first].filhos.push_back(c);
        }
      }
      nc.m += b * b * fe.nos[u].m;
    }
    std::vector<std::pair<std::uint32_t, double> >().swap(fe.nos[u].q);
    std::vector<std::uint32_t>().swap(fe.nos[u].filhos);
  }
  return fe;
}

// a_xy da A tabular, pela recursao de Henderson memorizada por par e com pilha explicita (um
// pedigree profundo estouraria a recursao): para x < y, a_xy = (a_{x,s_y} + a_{x,d_y}) / 2 e
// a_xx = 1 + F_x. Visita no maximo o fecho completo dos pares pedidos.
class ParentescoPares {
 public:
  ParentescoPares(const Pedigree& p, const std::vector<double>& f)
      : p_(p), f_(f), n_(static_cast<std::int64_t>(p.ids.size())) {}
  double operator()(std::int64_t a, std::int64_t b) {
    std::vector<std::pair<std::int64_t, std::int64_t> > pilha;
    pilha.push_back({std::min(a, b), std::max(a, b)});
    while (!pilha.empty()) {
      if ((++passos_ & 0xFFFFF) == 0) checa_interrupcao();
      const std::int64_t x = pilha.back().first, y = pilha.back().second;
      if (x == y || memo_.count(chave(x, y, n_))) { pilha.pop_back(); continue; }
      const std::int64_t s = p_.pai[y], d = p_.mae[y];
      double v[2] = {0.0, 0.0};
      bool falta = false;
      const std::int64_t pp[2] = {s, d};
      for (int t = 0; t < 2; t++) {
        if (pp[t] < 0) continue;
        if (!pega(x, pp[t], v[t])) {
          pilha.push_back({std::min(x, pp[t]), std::max(x, pp[t])});
          falta = true;
        }
      }
      if (falta) continue;
      memo_.emplace(chave(x, y, n_), 0.5 * (v[0] + v[1]));
      pilha.pop_back();
    }
    double out = 0.0;
    pega(a, b, out);
    return out;
  }

 private:
  bool pega(std::int64_t a, std::int64_t b, double& v) const {
    if (a == b) { v = 1.0 + f_[static_cast<std::size_t>(a)]; return true; }
    auto it = memo_.find(chave(a, b, n_));
    if (it == memo_.end()) return false;
    v = it->second;
    return true;
  }
  const Pedigree& p_;
  const std::vector<double>& f_;
  std::int64_t n_;
  std::unordered_map<std::int64_t, double> memo_;
  std::uint64_t passos_ = 0;
};

const double PESO_ESPARSO = 8.0;
// memoria de um par vivo do fecho: o no (~64), os vetores de pais e filhos (~100 com a
// reserva), a entrada do hash (~60) e a da fila; uma ESTIMATIVA pelo leiaute, nao medida
const double BYTES_POR_PAR = 320.0;
// entradas de lista lidas pelo grau minimo no tempo de uma operacao c_j^2 da cauda densa em
// ladrilhos: converte o custo da densa no trabalho que a ordenacao da esparsa pode ler alem
// do que a ordenacao da densa leu. Ordem de grandeza
// medida em pedigrees de leitegadas (corridas unicas): o grau minimo le ~1e8 a 2e8 entradas
// por segundo nessas candidatas, e com este teto a decisao levou de 0,1 a 1 s onde a rota
// densa se monta em 0,2 a 2 s (validation/dominance_hv91_routes.R tem as replicas). E uma
// constante, e nao uma medida da maquina: a rota (e os niveis do termo) nao pode depender
// do relogio.
const double RAZAO_TRABALHO = 0.05;

// por que a rota automatica escolheu o que escolheu (o R le pela posicao)
enum Motivo {
  M_DADA = 0,          // a rota foi pedida (route = "dense" / "sparse")
  M_COMPARADA = 1,     // as duas foram custeadas, ficou a mais barata
  M_FECHO = 2,         // o fecho passou das entradas ou da memoria da densa
  M_MEMORIA = 3,       // o bloco de pares da esparsa teria mais triplos que o triangulo denso
  M_TRABALHO = 4,      // a ordenacao da esparsa passou do teto de trabalho
  M_CUSTO = 5,         // a ordenacao da esparsa passou do custo da densa em c_j^2
  M_DENSA_GRANDE = 6,  // a densa passaria de 4 GB
  M_SEM_SUB = 7        // nenhum animal de Q tem os dois pais
};

struct CustoMME {
  double flops = 0.0;     // soma de c_j^2 de todas as colunas do fator
  double custo = 0.0;     // a parte esparsa pesando PESO_ESPARSO vezes a cauda densa
  double trabalho = 0.0;  // entradas de lista lidas pela ordenacao
  int corte = 0;          // 0 ordenou tudo, 1 teto de c_j^2, 2 teto de trabalho
};

// Fatoracao simbolica das MME conjuntas do modelo animal + dominancia de referencia: A^-1
// sobre o pedigree, a precisao Q (triangulo inferior, niveis 0..n_q-1 animais e n_q.. os
// pares), e o acoplamento de cada animal de Q com o proprio efeito aditivo (um registro por
// animal). Com h_densa o bloco h-h e um clique: o padrao recebido traz so a diagonal dele, os
// niveis h vao para o fim da ordem e a parte densa entra pela contagem fechada. A cauda e a
// do motor (bloco_denso_simbolico): as colunas finais completamente cheias, a partir de 128.
// teto > 0: o grau minimo desiste quando a soma parcial de c_j^2 passa dele, e o custo volta
// infinito. Como custo >= flops (PESO_ESPARSO >= 1), desistir acima do custo da outra rota
// nao muda a escolha: so poupa uma ordenacao que custaria tanto quanto a fatoracao perdedora.
// teto_trabalho > 0: desiste tambem quando a ordenacao le mais entradas de lista que isso.
CustoMME flops_conjuntas(const Csc& ainv, const std::vector<std::size_t>& linha_ped,
                         const std::vector<std::uint32_t>& li, const std::vector<std::uint32_t>& cj,
                         std::size_t n_q, std::size_t n_h, bool h_densa, double teto,
                         double teto_trabalho) {
  const std::size_t na = ainv.ncol, nt = na + n_q + n_h;
  std::vector<std::uint32_t> ri, rc;
  std::vector<double> rv;
  ri.reserve(ainv.nnz() + li.size() + n_q);
  rc.reserve(ri.capacity());
  for (std::size_t c = 0; c < na; c++)
    for (std::size_t k = ainv.colptr[c]; k < ainv.colptr[c + 1]; k++) {
      ri.push_back(ainv.linha[k]);
      rc.push_back(static_cast<std::uint32_t>(c));
    }
  for (std::size_t t = 0; t < li.size(); t++) {
    ri.push_back(static_cast<std::uint32_t>(na + li[t]));
    rc.push_back(static_cast<std::uint32_t>(na + cj[t]));
  }
  for (std::size_t a = 0; a < n_q; a++) {
    ri.push_back(static_cast<std::uint32_t>(na + a));
    rc.push_back(static_cast<std::uint32_t>(linha_ped[a]));
  }
  rv.assign(ri.size(), 1.0);
  const Csc jt = de_triplos(nt, nt, ri, rc, rv);
  std::vector<std::uint32_t>().swap(ri);
  std::vector<std::uint32_t>().swap(rc);
  std::vector<double>().swap(rv);
  double trabalho = 0.0;
  std::vector<std::size_t> perm = grau_minimo(jt, teto, teto_trabalho, &trabalho);
  if (perm.empty() && nt > 0) {
    CustoMME inf;
    inf.flops = inf.custo = std::numeric_limits<double>::infinity();
    inf.trabalho = trabalho;
    inf.corte = (teto_trabalho > 0.0 && trabalho > teto_trabalho) ? 2 : 1;
    return inf;
  }
  if (h_densa)
    std::stable_partition(perm.begin(), perm.end(),
                          [&](std::size_t v) { return v < na + n_q; });
  const Simbolica sb = simbolica(permuta_sim(jt, perm));
  std::vector<std::size_t> cnt(nt);
  for (std::size_t j = 0; j < nt; j++) cnt[j] = sb.colptr[j + 1] - sb.colptr[j];
  if (h_densa)
    for (std::size_t t = 0; t < n_h; t++) cnt[nt - n_h + t] = n_h - t;
  std::size_t T = 0;
  while (T < nt && cnt[nt - 1 - T] == T + 1) T++;
  if (T < 128) T = 0;
  CustoMME c;
  c.trabalho = trabalho;
  double cauda = 0.0;
  for (std::size_t j = 0; j < nt; j++) {
    const double q = static_cast<double>(cnt[j]);
    c.flops += q * q;
    if (j >= nt - T) cauda += q * q;
  }
  c.custo = PESO_ESPARSO * (c.flops - cauda) + cauda;
  return c;
}

// log|M| de uma simetrica positiva-definida esparsa (triangulo inferior), pela mesma
// Cholesky das MME; NaN se nao for positiva-definida
double logdet_esparsa(const Csc& m) {
  const std::vector<std::size_t> perm = grau_minimo(m);
  const Csc pa = permuta_sim(m, perm);
  const Simbolica sb = simbolica(pa);
  Csc L;
  if (!cholesky(pa, sb, L)) return std::nan("");
  return logdet(L);
}

}  // namespace

DominanciaInversa dominancia_inversa(const Pedigree& p, const std::vector<char>& em_q,
                                     int rota, double max_pares) {
  const std::size_t n = p.ids.size();
  if (em_q.size() != n) throw Erro("dominance inverse: one flag per pedigree animal");
  // o fecho indexa os pares em 32 bits: um teto acima disso nao limita nada, e a conversao de
  // um double acima de 2^64 para inteiro nem e definida (no MinGW virava 0)
  if (!(max_pares >= 1.0)) throw Erro("max_pairs must be at least 1");
  max_pares = std::min(max_pares, 4294967295.0);
  if (!p.mgs.empty())
    for (char c : p.mgs)
      if (c) throw Erro("the dominance inverse needs sire and dam; a sire / maternal-grandsire "
                        "pedigree has no sire x dam subclass");
  // defensivo: R_dominancia_inversa monta o pedigree sem metafundadores, e um codigo de pai
  // sem linha propria ja e recusado por constroi_pedigree
  if (p.n_mf > 0)
    throw Erro("the dominance inverse does not take metafounders: the Cockerham D is "
               "defined on a base of unrelated, non-inbred founders");
  const std::vector<double> f = endogamia(p);
  std::vector<double> dx(n);
  for (std::size_t i = 0; i < n; i++) dx[i] = variancia_mendeliana(p, f, i);

  DominanciaInversa r;
  // OS ANIMAIS DE Q, em ordem topologica, e o Delta de cada um
  for (std::size_t i = 0; i < n; i++)
    if (em_q[i]) r.animal.push_back(i);
  const std::size_t n_q = r.animal.size();
  if (n_q == 0) throw Erro("dominance inverse: no animal to put in the precision");
  std::vector<double> delta(n_q, 1.0);
  std::vector<std::int64_t> chave_sub(n_q, -1);
  const std::int64_t nn = static_cast<std::int64_t>(n);
  std::vector<std::pair<std::int64_t, std::int64_t> > cheios;
  {
    // quantos animais de Q caem em cada subclasse: com dois ou mais, Delta <= 0 torna a
    // propria D indefinida (dois irmaos completos tem D_ij = F_cc / 4 = 1 - Delta, e o bloco
    // 2 x 2 deles tem autovalor Delta); com um so, a D ainda pode ser positiva-definida
    std::unordered_map<std::int64_t, std::size_t> conta;
    std::size_t ruim = n_q;
    for (std::size_t a = 0; a < n_q; a++) {
      const std::size_t i = r.animal[a];
      const std::int64_t s = p.pai[i], d = p.mae[i];
      if (s < 0 || d < 0) continue;
      const double fs = f[static_cast<std::size_t>(s)], fd = f[static_cast<std::size_t>(d)];
      delta[a] = 1.0 - 0.25 * ((1.0 + fs) * (1.0 + fd) + 4.0 * f[i] * f[i]);
      if (!(delta[a] > 1e-10) && ruim == n_q) ruim = a;
      chave_sub[a] = chave(s, d, nn);
      if (conta[chave_sub[a]]++ == 0) cheios.push_back({std::min(s, d), std::max(s, d)});
    }
    if (ruim < n_q) {
      const std::size_t i = r.animal[ruim];
      const std::size_t s = static_cast<std::size_t>(p.pai[i]);
      const std::size_t d = static_cast<std::size_t>(p.mae[i]);
      const std::size_t irmaos = conta[chave_sub[ruim]];
      char buf[64];
      std::snprintf(buf, sizeof buf, "%.6g", delta[ruim]);
      std::string msg = "animal '" + p.ids[i] + "' (sire '" + p.ids[s] + "', dam '" + p.ids[d] +
                        "') has Delta = " + buf + " = 1 - [(1 + F_s)(1 + F_d) + 4 F^2] / 4 <= 0. ";
      if (irmaos >= 2)
        msg += "Its subclass has " + std::to_string(irmaos) + " animals in Q, and full sibs have "
               "D_ij = 1 - Delta: D of these animals is itself not positive-definite, no route "
               "fits it. Leave animals without records out (animals =) or drop inbred records";
      else
        msg += "D may still be positive-definite, but it has no subclass representation. If the "
               "animal has no record use animals = <the animals with records>; otherwise "
               "kernel(id, K = dominance_matrix(ped))";
      throw Erro(msg);
    }
  }
  std::sort(cheios.begin(), cheios.end(),
            [&](const std::pair<std::int64_t, std::int64_t>& a,
                const std::pair<std::int64_t, std::int64_t>& b) {
              return chave(a.first, a.second, nn) < chave(b.first, b.second, nn);
            });
  const std::size_t n_sub = cheios.size();
  r.n_sub = n_sub;
  double soma_log_delta = 0.0;
  for (double v : delta) soma_log_delta += std::log(v);
  // os pais das subclasses cheias, linhas da A que a rota densa pede ao Colleau
  std::vector<std::size_t> pais;
  for (const auto& c : cheios) {
    pais.push_back(static_cast<std::size_t>(c.first));
    pais.push_back(static_cast<std::size_t>(c.second));
  }
  std::sort(pais.begin(), pais.end());
  pais.erase(std::unique(pais.begin(), pais.end()), pais.end());

  // A ROTA. O pico de memoria da densa, contado e nao estimado por baixo: a A entre os pais
  // (8 n_pais^2 bytes, solta logo depois de formar F); F, F^-1 e a copia de trabalho da
  // inversa em ladrilhos (tres triangulos de 4 n_sub^2 bytes); e, na saida, os triplos do
  // triangulo de Q (16 bytes por entrada) mais a Csc e os indices de de_triplos (10 bytes por
  // entrada), 18 n_sub^2 no total. Acima de 4 GB ela nao e oferecida e a esparsa e a unica.
  const double ns2 = static_cast<double>(n_sub) * static_cast<double>(n_sub);
  const double np2 = static_cast<double>(pais.size()) * static_cast<double>(pais.size());
  const double densa_bytes = 8.0 * np2 + 18.0 * ns2;
  const bool densa_cabe = densa_bytes <= 4.0e9;
  const double dense_lim = 0.5 * static_cast<double>(n_sub) * static_cast<double>(n_sub + 1);
  using relogio = std::chrono::steady_clock;
  const auto seg = [](relogio::time_point a, relogio::time_point b) {
    return std::chrono::duration<double>(b - a).count();
  };
  const relogio::time_point t_ini = relogio::now();
  Fecho fe;
  bool tem_fecho = false;
  if (rota == 2 || (rota == 0 && n_sub > 0)) {
    // na rota automatica o fecho desiste para a densa quando passa das entradas do bloco
    // denso OU da memoria da densa: um par vivo custa ~BYTES_POR_PAR (o no, os vetores de
    // pais e filhos e a entrada do hash), 40 vezes uma entrada densa. Num pedigree de
    // leitegadas de 48.660 animais e 4.800 subclasses, so com o primeiro limite (11,5 milhoes
    // de pares) o processo passou de 1,8 GB ainda no fecho e foi interrompido, quando a rota
    // densa inteira pede ~0,6 GB pela conta de densa_bytes.
    const double teto = rota == 0 && densa_cabe
        ? std::min(max_pares, std::max(std::min(dense_lim, densa_bytes / BYTES_POR_PAR), 1e5))
        : max_pares;
    fe = fecho_pares(p, dx, cheios, static_cast<std::size_t>(teto));
    r.n_criados = fe.nos.size();
    r.n_processados = fe.processados;
    r.fecho_abortado = fe.abortado;
    if (fe.abortado && (rota == 2 || !densa_cabe))
      throw Erro("the pair closure passed max_pairs = " +
                 std::to_string(static_cast<unsigned long long>(teto)) +
                 " pairs" + std::string(densa_cabe ? "" : ", and the dense subclass route "
                 "would need more than 4 GB") + ". Restrict the precision to the animals with "
                 "records (animals =), or raise max_pairs if the memory allows");
    tem_fecho = !fe.abortado;
  }
  const relogio::time_point t_fecho = relogio::now();
  r.seg_fecho = seg(t_ini, t_fecho);
  if (rota == 1 && !densa_cabe)
    throw Erro("the dense subclass route needs about " +
               std::to_string((long long) std::ceil(densa_bytes / 1e9)) + " GB for " +
               std::to_string(n_sub) + " subclasses; use route = \"sparse\" or animals =");

  // OS PARES VIVOS da rota esparsa, em ordem de chave, e o padrao de Q
  std::vector<std::uint32_t> vivos;
  std::vector<std::int64_t> posicao;
  if (tem_fecho) {
    for (std::size_t j = 0; j < fe.nos.size(); j++)
      if (fe.nos[j].vivo) vivos.push_back(static_cast<std::uint32_t>(j));
    std::sort(vivos.begin(), vivos.end(), [&](std::uint32_t a, std::uint32_t b) {
      return chave(fe.nos[a].x, fe.nos[a].y, nn) < chave(fe.nos[b].x, fe.nos[b].y, nn);
    });
    posicao.assign(fe.nos.size(), -1);
    for (std::size_t t = 0; t < vivos.size(); t++) posicao[vivos[t]] = static_cast<std::int64_t>(t);
  }
  // nivel h da subclasse de cada animal, nas duas rotas
  auto nivel_sub_esparsa = [&](std::size_t a) -> std::int64_t {
    if (chave_sub[a] < 0) return -1;
    return posicao[fe.indice.find(chave_sub[a])->second];
  };
  std::unordered_map<std::int64_t, std::size_t> pos_cheio;
  for (std::size_t c = 0; c < n_sub; c++)
    pos_cheio.emplace(chave(cheios[c].first, cheios[c].second, nn), c);

  // o triangulo inferior da parte d e do acoplamento, comum as duas rotas (com o nivel h)
  auto parte_d = [&](std::size_t base_h, bool esparsa,
                     std::vector<std::uint32_t>& li, std::vector<std::uint32_t>& cj,
                     std::vector<double>& v, std::vector<double>& s_h) {
    for (std::size_t a = 0; a < n_q; a++) {
      li.push_back(static_cast<std::uint32_t>(a));
      cj.push_back(static_cast<std::uint32_t>(a));
      v.push_back(1.0 / delta[a]);
      if (chave_sub[a] < 0) continue;
      const std::size_t h = esparsa ? static_cast<std::size_t>(nivel_sub_esparsa(a))
                                    : pos_cheio.find(chave_sub[a])->second;
      li.push_back(static_cast<std::uint32_t>(base_h + h));
      cj.push_back(static_cast<std::uint32_t>(a));
      v.push_back(-1.0 / delta[a]);
      s_h[h] += 1.0 / delta[a];
    }
  };

  // A DECISAO da rota automatica: flops da fatoracao das MME conjuntas de cada uma, com as
  // saidas antecipadas descritas no cabecalho
  int escolhida = rota;
  r.motivo = M_DADA;
  r.bytes_densa = densa_bytes;
  r.entradas_densa = dense_lim;
  if (n_sub == 0) {
    escolhida = 1;
    r.motivo = M_SEM_SUB;
  } else if (rota == 0) {
    if (!densa_cabe) {
      // sem a densa nao ha o que comparar: a analise simbolica da esparsa custaria uma
      // ordenacao inteira das MME (40 s num pedigree leiteiro de 50 mil animais) para nada
      escolhida = 2;
      r.motivo = M_DENSA_GRANDE;
    } else if (!tem_fecho) {
      escolhida = 1;
      r.motivo = M_FECHO;
    } else {
      const Csc ainv = a_inversa(p, f);
      // a densa primeiro: o bloco h dela e um clique de contagem fechada e a ordenacao sai
      // barata. O custo dela vira os dois tetos da ordenacao da esparsa, que e a cara.
      {
        std::vector<std::uint32_t> ld, cd;
        std::vector<double> vd, sd(n_sub, 0.0);
        parte_d(n_q, false, ld, cd, vd, sd);
        for (std::size_t c = 0; c < n_sub; c++) {
          ld.push_back(static_cast<std::uint32_t>(n_q + c));
          cd.push_back(static_cast<std::uint32_t>(n_q + c));
        }
        const CustoMME cd2 = flops_conjuntas(ainv, r.animal, ld, cd, n_q, n_sub, true, 0.0, 0.0);
        r.flops_densa = cd2.flops;
        r.custo_densa = cd2.custo;
        r.trabalho_densa = cd2.trabalho;
      }
      // os triplos do bloco de pares da Q esparsa, CONTADOS antes de serem montados, como a
      // montagem os cria: por par com k pais, o triangulo de (1 + k) x (1 + k) da linha de
      // (I - Q), com repeticao, mais a diagonal de W' Delta^-1 W. A densa monta o triangulo
      // inteiro, dense_lim; a A^-1 e a parte d sao comuns as duas e nao entram
      double n_ent_h = static_cast<double>(vivos.size());
      for (std::size_t t = 0; t < vivos.size(); t++) {
        const double k = static_cast<double>(fe.nos[vivos[t]].q.size());
        n_ent_h += 0.5 * (1.0 + k) * (2.0 + k);
      }
      r.entradas_esparsa = n_ent_h;
      if (n_ent_h > dense_lim) {
        escolhida = 1;
        r.motivo = M_MEMORIA;
      } else {
        std::vector<std::uint32_t> li, cj;
        std::vector<double> v, s_h(vivos.size(), 0.0);
        parte_d(n_q, true, li, cj, v, s_h);
        std::vector<double>().swap(v);
        for (std::size_t t = 0; t < vivos.size(); t++) {
          const NoPar& no = fe.nos[vivos[t]];
          const std::uint32_t ic = static_cast<std::uint32_t>(n_q + t);
          li.push_back(ic); cj.push_back(ic);
          for (std::size_t a = 0; a < no.q.size(); a++) {
            const std::uint32_t ia = static_cast<std::uint32_t>(n_q + posicao[no.q[a].first]);
            li.push_back(std::max(ic, ia)); cj.push_back(std::min(ic, ia));
            for (std::size_t b = 0; b < a; b++) {
              const std::uint32_t ib = static_cast<std::uint32_t>(n_q + posicao[no.q[b].first]);
              li.push_back(std::max(ia, ib)); cj.push_back(std::min(ia, ib));
            }
            li.push_back(ia); cj.push_back(ia);
          }
        }
        r.teto_trabalho = r.trabalho_densa + RAZAO_TRABALHO * r.custo_densa;
        const CustoMME ce = flops_conjuntas(ainv, r.animal, li, cj, n_q, vivos.size(), false,
                                            r.custo_densa, r.teto_trabalho);
        r.flops_esparsa = ce.flops;
        r.custo_esparsa = ce.custo;
        r.trabalho_ordem = ce.trabalho;
        if (ce.corte == 2) {
          escolhida = 1;
          r.motivo = M_TRABALHO;
        } else if (ce.corte == 1) {
          escolhida = 1;
          r.motivo = M_CUSTO;
        } else {
          escolhida = r.custo_densa < r.custo_esparsa ? 1 : 2;
          r.motivo = M_COMPARADA;
        }
      }
    }
  }
  r.rota = escolhida;
  // escolhida a densa, o fecho nao serve mais: sai antes de a densa alocar o dela, e o pico
  // da rota automatica fica o maior dos dois, nao a soma
  if (escolhida == 1 && tem_fecho) {
    Fecho().nos.swap(fe.nos);
    std::unordered_map<std::int64_t, std::uint32_t>().swap(fe.indice);
    std::vector<std::uint32_t>().swap(vivos);
    std::vector<std::int64_t>().swap(posicao);
    tem_fecho = false;
  }
  const relogio::time_point t_decisao = relogio::now();
  r.seg_decisao = seg(t_fecho, t_decisao);

  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  std::vector<double> logdet_qhh_v;
  double log_f = 0.0;
  std::size_t n_h = 0;
  if (escolhida == 2 && n_sub > 0) {
    // ROTA ESPARSA: Q a partir do fecho
    n_h = vivos.size();
    std::vector<double> s_h(n_h, 0.0);
    parte_d(n_q, true, li, cj, v, s_h);
    const std::size_t ini_h = li.size();
    for (std::size_t t = 0; t < n_h; t++) {
      const NoPar& no = fe.nos[vivos[t]];
      if (!(no.m > 0.0)) throw Erro("dominance inverse: a pair with a non-positive residual");
      log_f += std::log(no.m);
      // linha de (I - Q): 1 no proprio par, -q nos pais; contribui c c' / m, vezes 4
      std::vector<std::pair<std::uint32_t, double> > lin;
      lin.push_back({static_cast<std::uint32_t>(t), 1.0});
      for (const auto& pk : no.q) {
        if (posicao[pk.first] < 0) throw Erro("dominance inverse: an eliminated pair is still referenced");
        lin.push_back({static_cast<std::uint32_t>(posicao[pk.first]), -pk.second});
      }
      const double w = 4.0 / no.m;
      for (std::size_t a = 0; a < lin.size(); a++)
        for (std::size_t b = 0; b < lin.size(); b++) {
          if (lin[a].first < lin[b].first) continue;
          li.push_back(static_cast<std::uint32_t>(n_q + lin[a].first));
          cj.push_back(static_cast<std::uint32_t>(n_q + lin[b].first));
          v.push_back(w * lin[a].second * lin[b].second);
        }
    }
    for (std::size_t t = 0; t < n_h; t++) {
      li.push_back(static_cast<std::uint32_t>(n_q + t));
      cj.push_back(static_cast<std::uint32_t>(n_q + t));
      v.push_back(s_h[t]);
    }
    // Q_hh = 4 F^-1 + W' Delta^-1 W, para log|D| = log|K| + log|Q_hh|
    {
      std::vector<std::uint32_t> hi, hj;
      std::vector<double> hv;
      for (std::size_t t = ini_h; t < li.size(); t++) {
        hi.push_back(li[t] - static_cast<std::uint32_t>(n_q));
        hj.push_back(cj[t] - static_cast<std::uint32_t>(n_q));
        hv.push_back(v[t]);
      }
      r.logdet_qhh = logdet_esparsa(de_triplos(n_h, n_h, hi, hj, hv));
    }
    // os pares e a priori de cada um, (a_xx a_yy + a_xy^2) / 4
    ParentescoPares axy(p, f);
    for (std::size_t t = 0; t < n_h; t++) {
      const NoPar& no = fe.nos[vivos[t]];
      r.par_x.push_back(no.x);
      r.par_y.push_back(no.y);
      r.par_cheio.push_back(no.cheio ? 1 : 0);
      const double axx = 1.0 + f[static_cast<std::size_t>(no.x)];
      const double ayy = 1.0 + f[static_cast<std::size_t>(no.y)];
      const double a = axy(no.x, no.y);
      r.priori_h.push_back(0.25 * (axx * ayy + a * a));
    }
  } else if (n_sub > 0) {
    // ROTA DENSA: F das subclasses cheias pela A entre os pais (Colleau)
    n_h = n_sub;
    std::vector<std::size_t> off(n_sub + 1, 0);
    for (std::size_t j = 0; j < n_sub; j++) off[j + 1] = off[j] + (n_sub - j);
    std::vector<double> lf(off[n_sub]);
    {
      // a A entre os pais so vive ate F estar formada
      std::unordered_map<std::size_t, std::size_t> pp;
      for (std::size_t k = 0; k < pais.size(); k++) pp.emplace(pais[k], k);
      const Densa ap = a22_colleau(p, pais);
      std::vector<std::size_t> sx(n_sub), sy(n_sub);
      for (std::size_t c = 0; c < n_sub; c++) {
        sx[c] = pp.find(static_cast<std::size_t>(cheios[c].first))->second;
        sy[c] = pp.find(static_cast<std::size_t>(cheios[c].second))->second;
      }
      for (std::size_t j = 0; j < n_sub; j++)
        for (std::size_t i = j; i < n_sub; i++)
          lf[off[j] + (i - j)] = ap.at(sx[i], sx[j]) * ap.at(sy[i], sy[j]) +
                                 ap.at(sx[i], sy[j]) * ap.at(sy[i], sx[j]);
    }
    for (std::size_t c = 0; c < n_sub; c++) {
      r.par_x.push_back(cheios[c].first);
      r.par_y.push_back(cheios[c].second);
      r.par_cheio.push_back(1);
      r.priori_h.push_back(0.25 * lf[off[c]]);
    }
    const int nth = threads();
    if (!cholesky_empacotada(lf, n_sub, nth))
      throw Erro("dominance inverse: the subclass relationship F is not positive-definite");
    for (std::size_t j = 0; j < n_sub; j++) log_f += 2.0 * std::log(lf[off[j]]);
    std::vector<double> zf(off[n_sub]);
    inversa_empacotada(lf.data(), n_sub, zf.data(), nth);
    std::vector<double>().swap(lf);
    std::vector<double> s_h(n_sub, 0.0);
    // reservados de uma vez: o crescimento por dobra passaria o pico contado acima
    li.reserve(off[n_sub] + 2 * n_q);
    cj.reserve(off[n_sub] + 2 * n_q);
    v.reserve(off[n_sub] + 2 * n_q);
    parte_d(n_q, false, li, cj, v, s_h);
    for (std::size_t j = 0; j < n_sub; j++) {
      zf[off[j]] = 4.0 * zf[off[j]] + s_h[j];
      for (std::size_t i = j + 1; i < n_sub; i++) zf[off[j] + (i - j)] *= 4.0;
    }
    for (std::size_t j = 0; j < n_sub; j++)
      for (std::size_t i = j; i < n_sub; i++) {
        const double x = zf[off[j] + (i - j)];
        if (x == 0.0) continue;
        li.push_back(static_cast<std::uint32_t>(n_q + i));
        cj.push_back(static_cast<std::uint32_t>(n_q + j));
        v.push_back(x);
      }
    if (!cholesky_empacotada(zf, n_sub, nth))
      throw Erro("dominance inverse: the subclass block of the precision is not positive-definite");
    double lq = 0.0;
    for (std::size_t j = 0; j < n_sub; j++) lq += 2.0 * std::log(zf[off[j]]);
    r.logdet_qhh = lq;
    std::vector<double>().swap(zf);
  } else {
    // nenhum animal de Q com os dois pais: Q = I
    for (std::size_t a = 0; a < n_q; a++) {
      li.push_back(static_cast<std::uint32_t>(a));
      cj.push_back(static_cast<std::uint32_t>(a));
      v.push_back(1.0);
    }
    r.logdet_qhh = 0.0;
  }
  r.q = de_triplos(n_q + n_h, n_q + n_h, li, cj, v);
  r.logdet_k = soma_log_delta + log_f - static_cast<double>(n_h) * std::log(4.0);
  // log|D| pelo lema do determinante: |Q| = |Q_hh| / |D|
  r.logdet_d = r.logdet_k + r.logdet_qhh;
  r.seg_montagem = seg(t_decisao, relogio::now());
  return r;
}

}  // namespace br
