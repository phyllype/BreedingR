// Declaracoes compartilhadas entre a montagem (mme.cpp), o estimador (aireml.cpp) e a
// travessia (entrada.cpp).
//
// Um unico lugar, de proposito: a primeira versao redefinia `Desenho` em dois .cpp, o que e
// violacao de ODR — os dois corpos tem de ser identicos token a token, e qualquer edicao
// futura num deles quebraria o outro em silencio, no vinculo ou pior, em tempo de execucao.

#ifndef BREEDINGR_MME_H
#define BREEDINGR_MME_H

#include "modelo.h"

namespace br {

// ---- pedigree.cpp
struct Pedigree {
  std::vector<std::string> ids;
  std::vector<std::int64_t> pai, mae;
  // metafundadores (Legarra et al. 2015), Gamma DIAGONAL nesta versao: eh_mf marca as
  // linhas-base virtuais e gama traz o gamma de cada uma (0 < gamma < 2). O truque que
  // faz tudo generalizar: F(mf) = gamma - 1, e as formulas de variancia mendeliana
  // existentes ja leem o resto sozinhas.
  std::vector<char> eh_mf;
  std::vector<double> gama;
};
Pedigree constroi_pedigree(const std::vector<std::string>&, const std::vector<std::string>&,
                           const std::vector<std::string>&,
                           const std::vector<std::string>& = {},
                           const std::vector<double>& = {});
std::vector<double> endogamia(const Pedigree&);
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
Montado monta_mme(const Desenho&, const std::vector<double>&);
Densa cov_grupo(const Modelo&, const std::vector<double>&, std::size_t);

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
  double neg2logl = 0.0;
  std::vector<double> theta, se, solucao;
  std::vector<double> pev;   // diagonal de [C_s^-1] vezes s2e, na numeracao das colunas
  std::string mensagem;
  std::size_t fora_do_padrao = 0;
};
Desenho monta_desenho(const Modelo&, const Tabela&, const Pedigree*);

// ---- genomica.cpp
struct RelatorioG {
  std::size_t n_imputados = 0;
  std::size_t n_monomorficos = 0;
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
// matriz (Liu et al. 2014; Vandenplas et al. 2019). theta e dado: e um resolvedor.
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
                           const Pedigree*);
std::vector<std::string> nomes_theta_do_mt(const DesenhoMT&);
struct AjusteMT {
  bool convergiu = false;
  std::size_t iters = 0;
  double reldelta = 0.0;
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
                            const Pedigree*, const std::string&, const std::string&);
std::vector<std::string> nomes_theta_ar1(const DesenhoAR&);
AjusteMT ajusta_ar1(const DesenhoAR&, std::size_t, double, bool = false);
// passo unico no AR(1): o mesmo nucleo do uni, sobre os mesmos campos
RelatorioG aplica_genomica(DesenhoAR&, const Pedigree&, const std::vector<std::string>&,
                           Densa&, double, const std::vector<std::string>&,
                           std::size_t = 0);

}  // namespace br

#endif
