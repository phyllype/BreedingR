// Declaracoes compartilhadas entre a montagem (mme.cpp), o estimador (aireml.cpp) e a
// travessia (entrada.cpp).
//
// Um unico lugar, de proposito: a primeira versao redefinia `Desenho` em dois .cpp, o que e
// violacao de ODR — os dois corpos tem de ser identicos token a token, e qualquer edicao
// futura num deles quebraria o outro em silencio, no vinculo ou pior, em tempo de execucao.

#ifndef BREEDINGR_MME_H
#define BREEDINGR_MME_H

#include <unordered_set>
#include "modelo.h"

#include <limits>

namespace br {

// ---- pedigree.cpp
struct Pedigree {
  std::vector<std::string> ids;
  std::vector<std::int64_t> pai, mae;
  // metafundadores (Legarra et al. 2015). eh_mf marca as linhas-base virtuais, prefixadas
  // ao pedigree, e gama traz a DIAGONAL de Gamma por linha. O truque que faz a variancia
  // mendeliana generalizar sozinha: F(mf) = gamma_ii - 1, NEGATIVO quando gamma < 1 e
  // legitimamente negativo, entao nao pode ser grampeado em zero em lugar nenhum.
  //
  // ATENCAO ao pseudo-codigo IMPRESSO no artigo, que traz "F(i) = 1 - gamma(i)": o sinal
  // esta trocado la, e com ele o proprio Exemplo 1 do artigo nao fecha. O texto corrido diz
  // o certo: "a self-relationship of a11 = gamma and an individual inbreeding coefficient
  // of Fi = a11 - 1 = gamma - 1".
  //
  // Gamma pode ser CHEIA. O gamma_jk fora da diagonal e a ancestralidade compartilhada
  // entre duas bases, que e o parametro que interessa em analise multirraca, e ele entra em
  // TRES lugares, dos quais so o primeiro e obvio:
  //   1. o bloco Gamma^-1 de A^-1 (gama_inv), SOMADO por cima das contribuicoes que os
  //      filhos ja jogam nas linhas dos metafundadores;
  //   2. a endogamia, pelo termo l_mf' Gamma l_mf, computado como ||K' l_mf||^2 com K a
  //      Cholesky inferior (gama_chol), que e a forma publicada e a numericamente estavel;
  //   3. INDIRETAMENTE na variancia mendeliana, sempre que um pai for mestico, porque o F
  //      desse pai ja carrega gamma_jk dentro.
  // Com Gamma diagonal os tres colapsam no comportamento anterior por construcao. Nao ha
  // ramo separado para o caso diagonal de proposito: um ramo duplicado e onde os dois
  // comportamentos divergiriam com o tempo.
  //
  // col_mf existe porque a ordenacao topologica REORDENA as linhas: o metafundador k de
  // Gamma nao fica necessariamente na linha k do pedigree ordenado. col_mf[i] devolve a
  // coluna de Gamma da linha i, ou -1. Sem isso, um desalinhamento com Gamma diagonal so
  // trocaria gamma_ii de lugar; com Gamma cheia ele corrompe todos os cruzados de um jeito
  // dificil de perceber.
  std::vector<char> eh_mf;
  std::vector<double> gama;
  std::vector<std::int64_t> col_mf;
  std::size_t n_mf = 0;
  // K com K K' = Gamma, de decomposicao ESPECTRAL e portanto NAO triangular: e o que
  // permite Gamma singular (gamma = 0, o limite de grupo de pais desconhecidos, e dois
  // metafundadores para a mesma populacao). Quem soma ||K' l||^2 tem de varrer a coluna
  // INTEIRA, nao so o triangulo.
  std::vector<double> gama_chol;   // n_mf x n_mf, por linhas
  std::vector<double> gama_inv;    // Gamma^-1 (pseudo-inversa se singular), por linhas
};
Pedigree constroi_pedigree(const std::vector<std::string>&, const std::vector<std::string>&,
                           const std::vector<std::string>&,
                           const std::vector<std::string>& = {},
                           const std::vector<double>& = {});
std::vector<double> endogamia(const Pedigree&);
// estado das ESTIMATIVAS no verbose, compartilhado pelos tres ajustadores iterativos
void imprime_theta(const std::vector<double>&, const std::vector<std::string>&,
                   const Modelo&);
Csc a_inversa(const Pedigree&, const std::vector<double>&);

// ---- linalg.cpp
struct Simbolica {
  std::vector<std::int64_t> pai;
  std::vector<std::size_t> colptr;
  std::size_t n = 0;
  std::size_t nnz() const { return colptr.empty() ? 0 : colptr.back(); }
};
Csc triu(const Csc&);
Simbolica simbolica(const Csc&);
bool cholesky(const Csc&, const Simbolica&, Csc&);

// Cache da fatoracao simbolica para os lacos de ajuste: o PADRAO de C nao depende de
// theta (a montagem empurra sempre os mesmos slots, zeros explicitos inclusive), entao
// o grau minimo e a analise simbolica valem para o ajuste inteiro e sao refeitos hoje a
// cada avaliacao por nada. Ponteiro nulo preserva o comportamento avulso dos gates.
struct CacheSimbolica {
  bool pronto = false;
  std::vector<std::size_t> perm;
  Simbolica sb;
};
double logdet(const Csc&);
std::vector<double> resolve(const Csc&, const std::vector<double>&);
Csc permuta_sim(const Csc&, const std::vector<std::size_t>&);
std::vector<std::size_t> grau_minimo(const Csc&);

// ---- selinv.cpp
struct SelInv {
  std::size_t n = 0;
  std::vector<std::size_t> colptr;
  std::vector<std::uint32_t> linha;
  std::vector<double> valor;
  bool get(std::size_t, std::size_t, double&) const;
  std::vector<double> diagonal() const;
};
std::size_t bloco_denso_final(const Csc&);
SelInv inversa_seletiva(const Csc&, std::size_t);

// Uma matriz de covariancia DECLARADA pelo usuario (kernel(id, K=)): os ids que nomeiam
// linhas e colunas, e a K densa. O vetor de declaracoes anda PARALELO aos termos do
// modelo — entrada vazia para termo que nao e kernel — para que o desenho ache a K do
// termo sem um mapa a parte.
struct KernelDecl {
  std::vector<std::string> ids;
  Densa k;
  bool vazia() const { return ids.empty(); }
};

// ---- o desenho completo
struct Desenho {
  Modelo modelo;
  Densa x;
  std::vector<std::string> nomes_x;
  std::vector<std::string> saiu_x;
  std::vector<DesenhoTermo> aleatorios;
  std::vector<Csc> kinv;
  std::vector<double> kinv_logdet;
  std::vector<double> y;
  std::vector<char> usa;
  std::size_t nlin = 0;
  // Pesos: um registro de peso w tem residual s2e/w. Entram como ESCALA DE LINHA por
  // sqrt(w) em y, X e Z — as equacoes normais da tabela escalada SAO as equacoes
  // ponderadas, entao score, AI e inversa seletiva continuam valendo sem uma linha de
  // mudanca. So a verossimilhanca precisa do jacobiano, que e esta constante.
  double logdet_peso = 0.0;

  std::size_t n_usadas() const {
    std::size_t s = 0;
    for (char u : usa) s += u;
    return s;
  }
  std::size_t largura(std::size_t g) const {
    std::size_t s = 0;
    for (std::size_t t : modelo.grupos[g].termos)
      for (const DesenhoTermo& d : aleatorios)
        if (d.termo == t) s += d.z.ncol;
    return s;
  }
  std::size_t total_colunas() const {
    std::size_t s = x.ncol;
    for (std::size_t g = 0; g < modelo.grupos.size(); g++) s += largura(g);
    return s;
  }
};

// ---- mme.cpp
struct Montado {
  Csc c;
  std::vector<double> rhs;
  double logdet_g = 0.0;
  double s2e = 0.0;
  std::size_t n_fixo = 0;
  std::size_t total = 0;
  std::vector<std::size_t> offset_grupo;
  bool ok = false;
};
// a K declarada de um grupo, compartilhada pelos tres ajustadores
void reduz_kernels(const Modelo& m, const std::vector<KernelDecl>*& kernels,
                   std::vector<KernelDecl>& kern_red,
                   std::vector<std::unordered_set<std::string> >& kern_nulos);
void casa_niveis_nulos(const Modelo& m, const std::vector<DesenhoTermo*>& aleatorios,
                       const std::vector<std::unordered_set<std::string> >& kern_nulos,
                       const Tabela& t, std::size_t nlin);
void kinv_declarada(const Modelo&, const Grupo&, const std::vector<KernelDecl>*,
                    std::vector<Csc>&, std::vector<double>&);
Montado monta_mme(const Desenho&, const std::vector<double>&);
Densa cov_grupo(const Modelo&, const std::vector<double>&, std::size_t);

// ---- aireml.cpp: Cholesky of a small dense covariance group; false when not PD.
bool chol_pequena(const Densa&, Densa&);

struct Neg2LogL {
  bool ok = false;
  double valor = 0.0;
  std::vector<double> solucao;
  double logdet_c = 0.0;
  std::vector<std::size_t> perm;
  Csc L;
};
Neg2LogL neg2logl_esparsa(const Desenho&, const std::vector<double>&);
double neg2logl_densa_V(const Desenho&, const std::vector<double>&);

// ---- aireml.cpp
struct Avaliacao {
  bool ok = false;
  double neg2logl = 0.0;
  std::vector<double> solucao;
  std::vector<double> score;
  Densa ai;
  std::vector<double> em_theta;
  std::size_t fora_do_padrao = 0;
};
Avaliacao avalia(const Desenho&, const std::vector<double>&, CacheSimbolica* = nullptr);

struct Ajuste {
  bool convergiu = false;
  std::size_t iters = 0;
  double reldelta = 0.0;
  // Newton decrement of the FREE components at the final point, g' AI^-1 g with the
  // active set excluded: the convergence certificate (~ twice the -2logL gap to the
  // optimum near it). NaN when the fit died before a first evaluation.
  double decremento = std::numeric_limits<double>::quiet_NaN();
  double neg2logl = 0.0;
  std::vector<double> theta, se, solucao;
  std::vector<double> pev;   // diagonal de [C_s^-1] vezes s2e, na numeracao das colunas
  std::vector<double> score; // score no ponto final: a evidencia de otimo
  std::vector<double> vcov;  // 2 AI^-1, a covariancia amostral das componentes
  std::string mensagem;
  std::size_t fora_do_padrao = 0;
};
Desenho monta_desenho(const Modelo&, const Tabela&, const Pedigree*,
                      const std::vector<double>* = nullptr,
                      const std::vector<KernelDecl>* = nullptr);

// ---- genomica.cpp
struct RelatorioG {
  std::size_t n_imputados = 0;
  std::size_t n_monomorficos = 0;
  // Diagonal de G* e a LINHA do pedigree de cada genotipado, para que accuracy() possa
  // dividir pela variancia a priori CERTA. Sob passo unico a priori de um genotipado e a
  // diagonal de H, que naquele bloco e a de G*, e nao 1 + F do pedigree: medido numa
  // populacao simulada de 510 animais todos genotipados, as duas diferem ate 0,18 e movem
  // uma acuracia individual ate 0,067. A media nao muda (0,6966 contra 0,6964); o que muda
  // e o INDIVIDUO, e com ele o ranqueamento por acuracia.
  std::vector<double> diag_gstar;
  std::vector<std::size_t> linha_ped;
};
Densa vanraden_g(Densa&, RelatorioG&);
Densa a22_inversa(const Csc&, const std::vector<std::size_t>&);
Densa ajusta_g_para_a22(const Densa&, const Densa&, double);
Csc constroi_hinv(const Csc&, const std::vector<std::size_t>&, const Densa&, const Densa&);

// Substitui o K^-1 dos grupos com parentesco pelo H^-1 do passo unico.
RelatorioG aplica_genomica(Desenho&, const Pedigree&, const std::vector<std::string>&,
                           Densa&, double, const std::vector<std::string>&,
                           std::size_t = 0);
std::vector<double> partida(const Desenho&);

// ---- sssnp.cpp: o passo unico sem G — marcadores como equacoes, PCG, A22^-1 livre de
// matriz (Liu et al. 2014; Masuda et al. 2017; Vandenplas et al. 2018, 2019). theta e
// dado: e um resolvedor.
struct SnpBlup {
  bool ok = false, convergiu = false;
  std::size_t iters = 0;
  double residuo = 0.0;              // ||r|| / ||rhs|| ao parar
  std::vector<double> solucao;       // as colunas do desenho de base (fixos + grupos)
  std::vector<double> efeitos;       // m efeitos de marcador, NaN nos monomorficos
  std::vector<char> usa_marcador;
  std::size_t n_imputados = 0, n_monomorficos = 0;
  std::size_t n_fixo = 0;
  std::vector<std::size_t> offset_grupo;
};
SnpBlup snp_blup(const Desenho&, Densa&, const std::vector<std::string>&,
                 const std::vector<double>&, double, std::size_t, double, bool = false);

// ---- gibbs.cpp: o lado bayesiano, bloco unico de localizacao + condicionais conjugadas
struct GibbsSaida {
  bool ok = false;
  std::size_t n_amostras = 0, ntheta = 0;
  std::vector<double> amostras;    // n_amostras x ntheta, por linha
  std::vector<double> media_loc, var_loc;
  std::string mensagem;
};
GibbsSaida gibbs(const Desenho&, std::size_t n_iter, std::size_t burnin, std::size_t thin,
                 bool loc_fixa, const std::vector<double>* theta_fixo, bool verboso = false);
Ajuste ajusta(const Desenho&, const std::vector<double>*, std::size_t, std::size_t, double,
              bool = false);


// ---- multitrait.cpp / multitrait2.cpp
struct DesenhoMT {
  Modelo modelo;
  std::vector<std::string> alvos;
  std::size_t t = 0;
  Densa x;
  std::vector<std::string> nomes_x;
  std::vector<std::string> saiu_x;
  std::vector<DesenhoTermo> aleatorios;
  std::vector<Csc> kinv;
  std::vector<double> kinv_logdet;
  Densa y;
  std::vector<char> usa;
  // obs[r*t + tau]: a caracteristica tau foi observada no registro r. Um registro pode
  // participar com um subconjunto — e a R0 dele e a submatriz do padrao.
  std::vector<char> obs;
  std::size_t nlin = 0;
  std::size_t n_usadas() const { std::size_t s = 0; for (char u : usa) s += u; return s; }
  std::size_t largura(std::size_t g) const {
    std::size_t s = 0;
    for (std::size_t tm : modelo.grupos[g].termos)
      for (const DesenhoTermo& d : aleatorios)
        if (d.termo == tm) s += d.z.ncol * t;
    return s;
  }
  std::size_t total_colunas() const {
    std::size_t s = x.ncol * t;
    for (std::size_t g = 0; g < modelo.grupos.size(); g++) s += largura(g);
    return s;
  }
};
struct AvaliacaoMT {
  bool ok = false;
  double neg2logl = 0.0;
  std::vector<double> solucao;
  std::vector<double> pev;   // diagonal de C^-1 na numeracao original (unidades absolutas)
  std::vector<double> score;
  Densa ai;
  std::size_t fora_do_padrao = 0;
};
AvaliacaoMT avalia_mt(const DesenhoMT&, const std::vector<double>&, CacheSimbolica* = nullptr);
double neg2logl_densa_V_mt(const DesenhoMT&, const std::vector<double>&);
DesenhoMT monta_desenho_mt(Modelo, const std::vector<std::string>&, const Tabela&,
                           const Pedigree*,
                           const std::vector<KernelDecl>* = nullptr);
std::vector<std::string> nomes_theta_do_mt(const DesenhoMT&);
struct AjusteMT {
  bool convergiu = false;
  std::size_t iters = 0;
  double reldelta = 0.0;
  // Newton decrement at the final point (these walkers have no active set, so it is
  // g' AI^-1 g whole); NaN when the fit died before a first evaluation.
  double decremento = std::numeric_limits<double>::quiet_NaN();
  double neg2logl = 0.0;
  std::vector<double> theta, se, solucao;
  std::vector<double> pev;   // diagonal de C^-1 no otimo, ja em unidades absolutas
  std::string mensagem;
  std::size_t fora_do_padrao = 0;
};
AjusteMT ajusta_mt(const DesenhoMT&, std::size_t, double, bool = false);
// passo unico na multi: o mesmo nucleo do uni, sobre os mesmos campos
RelatorioG aplica_genomica(DesenhoMT&, const Pedigree&, const std::vector<std::string>&,
                           Densa&, double, const std::vector<std::string>&,
                           std::size_t = 0);


// ---- ar1.cpp / ar1b.cpp
// Entrada de linha do esqueleto AR: coluna global, caracteristica e valor.
struct EntAR {
  std::uint32_t col;
  std::uint32_t trait;
  double val;
};
struct DesenhoAR {
  Modelo modelo;
  std::vector<std::string> alvos;   // t = alvos.size(); t = 1 e o caminho classico
  std::size_t t = 1;
  Densa x;
  std::vector<std::string> nomes_x;
  std::vector<std::string> saiu_x;
  std::vector<DesenhoTermo> aleatorios;
  std::vector<Csc> kinv;
  std::vector<double> kinv_logdet;
  Densa y;                          // nlin x t
  std::vector<char> usa;
  std::size_t nlin = 0;
  std::vector<std::vector<std::size_t>> sujeitos;
  std::vector<double> tempo;
  // Esqueleto pre-computado no desenho (nada disto depende de theta): as linhas de W ja
  // em colunas globais COM a caracteristica, e as colunas distintas de cada sujeito.
  std::vector<std::vector<EntAR>> lw;
  std::vector<std::vector<std::uint32_t>> cols_suj;
  // offset_s2e = inicio do vech(R0) (t = 1: o proprio s2e); offset_rho = logo depois
  std::size_t offset_s2e = 0, offset_rho = 0;
  std::size_t n_usadas() const { std::size_t s = 0; for (char u : usa) s += u; return s; }
  std::size_t largura(std::size_t g) const {
    std::size_t s = 0;
    for (std::size_t tm : modelo.grupos[g].termos)
      for (const DesenhoTermo& d : aleatorios)
        if (d.termo == tm) s += d.z.ncol * t;
    return s;
  }
  std::size_t total_colunas() const {
    std::size_t s = x.ncol * t;
    for (std::size_t g = 0; g < modelo.grupos.size(); g++) s += largura(g);
    return s;
  }
};
struct AvaliacaoAR {
  bool ok = false;
  double neg2logl = 0.0;
  std::vector<double> solucao, score;
  std::vector<double> pev;   // diagonal de C^-1 na numeracao original (unidades absolutas)
  Densa ai;
  std::size_t fora_do_padrao = 0;
};
AvaliacaoAR avalia_ar1(const DesenhoAR&, const std::vector<double>&, CacheSimbolica* = nullptr);
double neg2logl_densa_V_ar1(const DesenhoAR&, const std::vector<double>&);
DesenhoAR monta_desenho_ar1(Modelo, const std::vector<std::string>&, const Tabela&,
                            const Pedigree*, const std::string&, const std::string&,
                            const std::vector<KernelDecl>* = nullptr);
std::vector<std::string> nomes_theta_ar1(const DesenhoAR&);
AjusteMT ajusta_ar1(const DesenhoAR&, std::size_t, double, bool = false);
// passo unico no AR(1): o mesmo nucleo do uni, sobre os mesmos campos
RelatorioG aplica_genomica(DesenhoAR&, const Pedigree&, const std::vector<std::string>&,
                           Densa&, double, const std::vector<std::string>&,
                           std::size_t = 0);

}  // namespace br

#endif
