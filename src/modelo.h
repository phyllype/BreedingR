// O modelo como DADO, e o desenho que sai dele.
//
// A unidade de layout e o GRUPO DE COVARIANCIA, nao o termo. Essa e a decisao que separa
// este desenho do anterior: com um bloco por termo, dois efeitos aleatorios nunca podem ser
// correlacionados e o modelo direto-materno e inexprimivel. Com o grupo como unidade,
//
//   Var(u_g) = C_g (x) K
//
// uma regressao aleatoria vira o caso de UM termo com m coeficientes, e o direto-materno vira
// o caso de DOIS termos escalares no mesmo grupo. A mesma penalidade kron(C_g^-1, K^-1), sem
// caso especial nenhum.
//
// A convencao de coluna que faz a penalidade sair exatamente assim:
//
//   coluna de Z = coeficiente * n_niveis + nivel
//
// com o NIVEL variando mais rapido. Trocar a ordem nao quebra nada visivelmente: o modelo
// monta, converge, e estima outra coisa.

#ifndef BREEDINGR_MODELO_H
#define BREEDINGR_MODELO_H

#include "estruturas.h"

namespace br {

enum class Efeito { Classe, Covariavel };
// Declarada: o termo traz a PROPRIA matriz de covariancia (kernel(id, K=)) — D de
// dominancia, G_AA de epistasia, uma parcial por raca. O K^-1 do grupo vem da K declarada
// em vez de A ou H, e todo o resto (penalidade kron, score, AI) nem sabe a diferenca.
enum class Estrutura { Fixo, Diagonal, Parentesco, Declarada };

struct Termo {
  std::string nome;
  std::string coluna;
  Efeito efeito = Efeito::Classe;
  Estrutura estrutura = Estrutura::Fixo;
  // colunas que formam a base; vazio = escalar, um coeficiente por nivel
  std::vector<std::string> base;
  // covariavel estimada DENTRO de cada nivel desta classe; num termo social, a coluna da
  // BAIA (o grupo de convivencia)
  std::string aninhado;
  // Efeito genetico INDIRETO (modelo associativo: Griffing 1967; Muir e Schinckel 2002;
  // Bijma et al. 2007): a incidencia da linha i
  // marca os COMPANHEIROS de baia do animal i, nao o proprio animal. O fenotipo de i carrega
  // o efeito social de cada colega; o efeito direto continua no termo animal comum, e os
  // dois dividem um grupo de covariancia com a correlacao direto-social estimada.
  bool social = false;
  // Diluicao do efeito indireto com o tamanho da baia (Bijma 2010, Genetics 186:1013-1028):
  // a entrada de Z_S de cada companheiro vale (n_i - 1)^(-diluicao), com n_i o numero de
  // animais DISTINTOS na baia do registro i. diluicao = 0 reproduz a soma do livro
  // (coeficiente 1 por companheiro, Mrode e Pocrnic 2023 cap. 9, n fixo); diluicao = 1 e a
  // media dos companheiros. So tem sentido num termo social; validado em monta_modelo.
  double diluicao = 0.0;
  bool aleatorio() const { return estrutura != Estrutura::Fixo; }
  std::size_t n_coef() const { return base.empty() ? 1 : base.size(); }
};

// Um grupo de covariancia: os termos que dividem uma matriz C_g.
struct Grupo {
  std::string nome;
  std::vector<std::size_t> termos;   // indices em Modelo::termos
  std::size_t dim = 0;               // soma dos coeficientes dos termos
  std::size_t offset = 0;            // onde os parametros comecam em theta
  std::size_t nparam = 0;            // dim (dim + 1) / 2
  Estrutura estrutura = Estrutura::Diagonal;

  // Indice em theta do elemento (i,j) da covariancia, guardada por colunas no triangulo
  // inferior. Exige i >= j; a chamada com i < j e trocada.
  std::size_t theta_idx(std::size_t i, std::size_t j) const {
    if (i < j) std::swap(i, j);
    return offset + (j * (2 * dim - j + 1)) / 2 + (i - j);
  }
};

struct Modelo {
  std::string alvo;
  std::vector<Termo> termos;
  std::vector<Grupo> grupos;
  // Componentes PRESOS pelo usuario (kernel(..., fixed = v)): 1 = nao anda. O caso que
  // motiva e a covariancia de erro CONHECIDA, Var(y) = s2a A + s2env I + V_e com V_e de
  // coeficiente fixo em 1: com a escala livre, V_e e s2env nao sao simultaneamente
  // identificaveis quando os v_i variam pouco (heterogeneidade aditiva contra
  // multiplicativa, Thompson & Sharp 1999), e o ajuste devolve h2 = 1 contra verdade 0,6.
  // Vazio quando nada esta preso, que e o caminho de sempre.
  // guarda o VALOR em que prender; NaN = livre. Guardar so um sinalizador
  // congelaria o componente no valor de PARTIDA e nao no pedido.
  std::vector<double> theta_fixo;
  bool preso(std::size_t k) const {
    return !theta_fixo.empty() && theta_fixo[k] == theta_fixo[k];
  }
  std::size_t offset_residual = 0;
  std::size_t ntheta = 0;
  bool tem_ausente = false;
  double codigo_ausente = 0.0;

  // Nomes dos parametros, na ordem de theta. Sem isto o usuario recebe um vetor de numeros
  // e tem de adivinhar qual e qual.
  std::vector<std::string> nomes_theta() const;
};

// Constroi o layout: valida, resolve os grupos, e numera theta.
//
// Um termo aleatorio sem grupo declarado ganha um grupo proprio — e o comportamento que
// mantem compativel o modelo comum, em que cada efeito tem a sua variancia.
Modelo monta_modelo(const std::string& alvo, std::vector<Termo> termos,
                    const std::vector<std::pair<std::string, std::vector<std::string>>>& grupos,
                    bool tem_ausente, double codigo_ausente);

// Uma coluna de dados: numerica ou textual, decidido na leitura.
struct Coluna {
  bool texto = false;
  std::vector<double> num;
  std::vector<std::string> txt;
  std::size_t tamanho() const { return texto ? txt.size() : num.size(); }
  // Rotulo de uma linha. Um id lido como 1758.0 e um escrito 1758 tem de dar o MESMO rotulo,
  // senao o termo com parentesco rejeita o conjunto inteiro por nao achar nivel em comum.
  std::string rotulo(std::size_t i) const;
};

struct Tabela {
  std::vector<std::string> nomes;
  std::vector<Coluna> colunas;
  std::size_t nlin = 0;
  const Coluna* acha(const std::string& n) const;
  std::vector<std::string> rotulos(const std::string& n) const;
  std::vector<double> numerico(const std::string& n) const;
};

// A matriz de incidencia de um termo, com os niveis que a geraram.
struct DesenhoTermo {
  std::size_t termo = 0;
  std::string nome;
  std::vector<std::string> niveis;
  std::size_t n_niveis = 0;
  std::size_t n_coef = 1;
  Csc z;                        // nlin x (n_niveis * n_coef)
  // false na linha cujo nivel nao existe no conjunto de niveis do termo
  std::vector<char> casou;
};

// `niveis_fixos` vem do pedigree para termos com parentesco; vazio deixa os niveis virem dos
// proprios dados, na ordem em que aparecem.
DesenhoTermo monta_termo(const Modelo& m, std::size_t k, const Tabela& t,
                         const std::vector<std::string>* niveis_fixos);

// Colunas de posto completo, por Gram-Schmidt modificado. Devolve as que ficam e as que saem.
void posto_completo(const Densa& x, double tol, std::vector<std::size_t>& fica,
                    std::vector<std::size_t>& sai);

}  // namespace br

#endif
