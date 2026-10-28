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

#ifndef UAITREM_MODELO_H
#define UAITREM_MODELO_H

#include "estruturas.h"

namespace br {

enum class Efeito { Classe, Covariavel };
// Declarada: o termo traz a PROPRIA matriz de covariancia (kernel(id, K=)), D de
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
  // Diluicao do efeito indireto com o tamanho da baia (Bijma 2010, Genetics 186:1029-1031):
  // a entrada de Z_S de cada companheiro vale (n_i - 1)^(-diluicao), com n_i o numero de
  // animais DISTINTOS na baia do registro i. diluicao = 0 reproduz a soma do livro
  // (coeficiente 1 por companheiro, Mrode e Pocrnic 2023 cap. 9, n fixo); diluicao = 1 e a
  // media dos companheiros. So tem sentido num termo social; validado em monta_modelo.
  double diluicao = 0.0;
  // Modelo pai / avo materno (Quaas e Pollak; Mrode e Pocrnic cap. 3): a linha i leva 1 no
  // PAI e 1/2 no AVO MATERNO, os dois niveis do MESMO efeito de touro. A coluna do avo vem em
  // `aninhado` (o R a manda com o prefixo "mgs:", que modelo_do_R tira); avo desconhecido
  // ("0", "", NA) deixa so a entrada do pai.
  bool mgs = false;
  // So no multicaracter: os caracteres em que o termo existe (pe(id, traits = "y1")). Vazio
  // = todos. Um registro de repeticao junto com um caracter de registro unico pede o pe so no
  // primeiro: no segundo ele se confunde com o residuo e a matriz de informacao fica
  // singular. Num grupo, cada termo pode ter o seu (expande_layout_mt).
  std::vector<std::string> tracos;
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
  // So no multicaracter: os indices (em ordem crescente) dos caracteres de cada termo,
  // paralelo a `termos`; uma lista vazia (ou nenhuma) = todos os t, o caso de sempre. O
  // bloco do termo q tem nt(q, t) caracteres, e a posicao pt e o caracter traco(q, pt).
  std::vector<std::vector<std::size_t>> tracos;
  std::size_t nt(std::size_t q, std::size_t t) const {
    return q < tracos.size() && !tracos[q].empty() ? tracos[q].size() : t;
  }
  std::size_t traco(std::size_t q, std::size_t pt) const {
    return q < tracos.size() && !tracos[q].empty() ? tracos[q][pt] : pt;
  }

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
  // O RESIDUO POR CLASSE (lista de correcoes, item 35): os rotulos das classes, na ordem dos
  // componentes var(residual|classe) a partir de offset_residual. Vazio = um residuo so, o
  // caminho de sempre
  std::vector<std::string> classes_residuo;
  std::size_t n_residuos() const { return classes_residuo.empty() ? 1 : classes_residuo.size(); }
  bool tem_ausente = false;
  double codigo_ausente = 0.0;

  // Nomes dos parametros, na ordem de theta. Sem isto o usuario recebe um vetor de numeros
  // e tem de adivinhar qual e qual.
  std::vector<std::string> nomes_theta() const;
};

// Constroi o layout: valida, resolve os grupos, e numera theta.
//
// Um termo aleatorio sem grupo declarado ganha um grupo proprio, e o comportamento que
// mantem compativel o modelo comum, em que cada efeito tem a sua variancia.
Modelo monta_modelo(const std::string& alvo, std::vector<Termo> termos,
                    const std::vector<std::pair<std::string, std::vector<std::string>>>& grupos,
                    bool tem_ausente, double codigo_ausente);

// O ROTULO DE UM NUMERO, o unico formato de id numerico do pacote. Inteiro exato abaixo de
// 2^53 sai com todos os digitos e sem expoente ("100000", nunca "1e+05"); qualquer outro
// numero finito sai na menor escrita de 15 a 17 digitos significativos que volta ao mesmo
// double ("1.5", "1234567.5"), entao dois numeros distintos nunca tem o mesmo rotulo. E o
// rotulo da coluna numerica dos dados (Coluna::rotulo) e o que o lado R usa, por R_rotulos,
// para escrever ids numericos do pedigree, dos genotipos e das chaves (rotulo_motor() em
// R/pedigree.R): um formato so, entao os dois lados nao divergem.
std::string rotulo_numero(double x);
// rotulo_numero(x) e o inteiro com todos os digitos: inteiro finito abaixo de 2^53
bool rotulo_exato(double x);
// "1e+05", "1.5e+07": a forma que o as.character() e o factor() do R dao a um double inteiro
// redondo. Devolve o rotulo inteiro do mesmo numero ("100000"), ou "" quando s nao e um
// inteiro exato escrito como o as.character() o escreveria ("1.23457e+06" nao e: o R escreve
// 1234570 como "1234570").
std::string inteiro_de_cientifico(const std::string& s);
// O complemento de uma mensagem de id sem par: quando o par existe em `outros` escrito do
// outro jeito ("1e+05" e "100000"), a frase que diz isso; "" quando nao existe. Varre
// `outros` inteiro, entao e so para a hora do erro.
std::string dica_cientifica(const std::string& id, const std::vector<std::string>& outros);

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

// ROTULO AUSENTE numa coluna que da nivel (classe fixa, id de termo aleatorio, classe de
// aninhamento, avo materno, baia de indirect(), sujeito do AR(1)). Uma regra so, a mesma de
// tipos_rotulo_ausente()/sem_rotulo() (R/indirect.R), que os ajustadores em R, pegs() e
// associative_matrix() usam, com os mesmos nomes de tipo nas mensagens:
//   coluna numerica: NA ou NaN ("NA"), Inf ("Inf"), -Inf ("-Inf");
//   coluna textual (e fator, que cruza como texto): texto vazio depois de tirar espacos,
//   tabulacoes e quebras de linha ("blank text"), e o texto "NaN" que factor() e
//   as.character() fazem de um NaN numerico. O NA textual ja para em tabela_do_R.
// Os textos "NA" e "Inf" ficam de fora: podem ser o nome de uma classe de verdade. A
// comparacao e por bytes, entao vale para texto em qualquer codificacao.
// `branco_e_codigo` tira o texto vazio da conta: e a coluna do avo materno, onde o vazio e o
// avo DESCONHECIDO, como o pai vazio no pedigree. `onde` troca "the column '<c>'" na
// mensagem (a baia de indirect() diz "the pen column '<c>' of indirect()"). Devolve "" quando
// toda linha tem rotulo (ou a coluna nao existe, e quem chama para no erro de sempre); senao
// o comeco da mensagem, "<o que> in the column '<c>' in n row(s), the first at row k".
std::string rotulos_ausentes(const Tabela& t, const std::string& col,
                             bool branco_e_codigo = false, const std::string& onde = "");
// VALOR AUSENTE numa coluna numerica que entra como valor (covariavel, coluna da base): NA,
// NaN, Inf e -Inf, no mesmo formato e com os mesmos nomes, com `papel` no lugar de "column"
// ("the covariate column 'x'"). Coluna textual devolve "": Tabela::numerico a recusa.
std::string valores_ausentes(const Tabela& t, const std::string& col,
                             const std::string& papel = "column");
// o texto em branco: vazio depois de tirar espacos, tabulacoes e quebras de linha
bool texto_em_branco(const std::string& s);
// o codigo de animal ou pai DESCONHECIDO numa coluna de termo com parentesco: "0" (o nivel 0
// do BLUPF90, efeito ausente no registro) e o texto "NA", os mesmos do pedigree. O registro
// fica na analise com incidencia zero naquele termo (monta_termo).
bool nivel_desconhecido(const std::string& s);
// A coluna que um termo com parentesco le (o animal, a mae, o pai, o avo materno de
// sire(mgs =)). Ali o "0" e o animal DESCONHECIDO para todo termo aleatorio que a le: o pe(dam)
// de maternal(dam) + pe(dam), o pe(id) de animal(id) + pe(id). Sem isso, as linhas de maes
// desconhecidas distintas formavam UM nivel "0" do pe, como se 60 animais de base tivessem a
// mesma mae. Medido (simulate_breeding semente 5, 660 registros, 60 de mae 0, REML de
// animal(id) + maternal(dam) + pe(dam)): var(pe) 0,0029 com a incidencia zero que fica
// aqui, 0,0031 com o nivel "0" comum (o pe numa copia da coluna) e 0,0075 com um nivel
// proprio para cada registro de mae desconhecida. A diferenca medida e pequena; o motivo da
// regra e o modelo, e nao o tamanho dela neste dado.
bool coluna_de_parentesco(const Modelo& m, const std::string& col);
// o termo ALEATORIO cujo "0" e desconhecido: o de parentesco e o que le a coluna de um deles
bool zera_desconhecido(const Modelo& m, const Termo& tm);

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
  // So no multicaracter: copia de Grupo::tracos do grupo do termo (monta_desenho_mt). O
  // termo ocupa z.ncol * nt(t) colunas, e a coluna da posicao pt e do caracter traco(pt).
  std::vector<std::size_t> tracos;
  std::size_t nt(std::size_t t) const { return tracos.empty() ? t : tracos.size(); }
  std::size_t traco(std::size_t pt) const { return tracos.empty() ? pt : tracos[pt]; }
};

// `niveis_fixos` e o conjunto de niveis que o termo indexa quando ele vem de fora da coluna:
// os ids do pedigree num termo com parentesco, os ids da K (ja reduzida) num kernel(), e a
// uniao ordenada dos rotulos num grupo sem estrutura de varios termos (monta_aleatorios, em
// mme.cpp). Nulo deixa os niveis virem dos proprios dados, na ordem em que aparecem: o termo
// fixo e o termo iid de grupo proprio. Recusa o rotulo ausente (rotulos_ausentes) na coluna
// que da o nivel e o valor nao finito (valores_ausentes) na covariavel e na base, em TODA
// linha.
DesenhoTermo monta_termo(const Modelo& m, std::size_t k, const Tabela& t,
                         const std::vector<std::string>* niveis_fixos);

// Colunas de posto completo, por Gram-Schmidt modificado DENSO. Devolve as que ficam e as que
// saem. E a referencia do posto esparso (o portao de R_posto_fixo): os desenhos ja nao a
// usam, porque pede o X denso e guarda cada coluna ortogonalizada como um vetor de n, O(n p^2).
// A rota densa de reserva de posto_esparso faz a mesma conta sem o X denso.
void posto_completo(const Densa& x, double tol, std::vector<std::size_t>& fica,
                    std::vector<std::size_t>& sai);

// ---- fixos.cpp: o desenho dos efeitos fixos, ESPARSO
// As colunas candidatas de X: o intercepto e cada termo fixo expandido, na ordem do modelo,
// com os nomes "termo=nivel". Por coluna (CSC nlin x p), com os MESMOS nao-zeros do desenho
// denso de antes: uma classe poe um 1 por linha, a covariavel aninhada o valor no nivel da
// linha, a covariavel simples os seus valores nao nulos. Um valor nao finito numa covariavel
// nao chega aqui: monta_termo o recusa em toda linha (valores_ausentes); o posto ainda
// derruba a coluna com NaN ou Inf que lhe chegue por outra via (R_posto_fixo).
struct FixosCandidatos {
  Csc x;
  std::vector<std::string> nomes;
  // por coluna: 0 o intercepto ou o nivel de uma classe, 1 a covariavel simples, 2 uma coluna
  // de covariavel aninhada; termo e o indice do termo de origem (o intercepto leva o maior
  // size_t); par_classe, numa coluna aninhada, e a coluna candidata do MESMO nivel da classe
  // de aninhamento quando essa coluna tambem e um termo fixo de classe, senao SEM_PAR. centro
  // e escala sao os de padroniza_covariaveis(): a coluna guardada e (x - centro) / escala, e a
  // intocada tem centro 0 e escala 1.
  std::vector<int> covar;
  std::vector<std::size_t> termo, par_classe;
  std::vector<double> centro, escala;
};
constexpr std::size_t SEM_PAR = static_cast<std::size_t>(-1);
FixosCandidatos colunas_fixas(const Modelo& m, const Tabela& t);
// A PADRONIZACAO DAS COVARIAVEIS, antes do posto, como o WOMBAT faz com as covariaveis
// (Meyer, manual 2017, secao 4.8.1.1, desvia registros e covariaveis das medias antes da
// analise). A covariavel simples vira (x - m) / s, com a media m e o dp s das linhas usadas
// (usa), e o posto e as MME passam a viver nessa escala. Um deslocamento grande (datas,
// carimbos de tempo) ou uma escala minuscula estragavam as duas coisas: com idade + 1e9 o theta
// de partida saia INADMISSIVEL, com + 3e8 a regressao saia 5,90 contra 0,0102, com + 1e11 ou
// x 1e-12 a coluna saia em dropped_x como dependente. Centrar e reparametrizacao porque o
// intercepto esta sempre no modelo. A aninhada em c vira (x - m) I_c / s quando c tambem e
// termo fixo de classe (o deslocamento cai nos niveis de c), e x I_c / s quando nao e, ou
// quando o termo esta em sem_centro (a volta de b nao fecha, ver aninhadas_sem_volta()). A
// coluna sem variacao real (dp ate 100 eps do maior |x|) fica como esta.
void padroniza_covariaveis(FixosCandidatos& fx, const std::vector<char>& usa,
                           const std::vector<std::size_t>& sem_centro);
// A VOLTA DE b a unidade do dado, de uma caracteristica: por coluna de x (a numeracao do
// posto), a escala e as colunas que absorvem o centro, com o coeficiente: b_j = b'_j / s_j e,
// para cada (alvo, coef), b_alvo += coef b_j. A covariavel simples manda -m ao intercepto; a
// aninhada no nivel k manda -m a coluna do nivel k da classe quando ela ficou, e quando o nivel
// k e o que caiu (r, a referencia), I_r = 1 - soma dos outros niveis manda -m ao intercepto e
// +m a cada outro nivel que ficou. Os alvos nao sao covariaveis, entao a ordem nao importa.
struct VoltaB {
  std::vector<double> escala;
  std::vector<std::vector<std::pair<std::size_t, double> > > absorve;
  bool vazia() const { return escala.empty(); }
};
// mantida: por coluna CANDIDATA, se ela ficou e tem equacao nesta caracteristica; linhas: as
// linhas em que a caracteristica entra. Devolve os termos aninhados cuja centragem nao volta
// pela classe: mais de um nivel da classe com registro fora de mantida (a classe confundida com
// outro efeito fixo). Esses sao refeitos sem centro.
std::vector<std::size_t> aninhadas_sem_volta(const FixosCandidatos& fx,
                                             const std::vector<char>& mantida,
                                             const std::vector<char>& linhas);
// a volta de b na numeracao de fica (as colunas que ficaram), para uma caracteristica
VoltaB monta_volta(const FixosCandidatos& fx, const std::vector<std::size_t>& fica,
                   const std::vector<char>& mantida, const std::vector<char>& linhas);
// b de uma caracteristica entre a escala padronizada e a do dado. pos(j) da a posicao em b da
// coluna j de x, ou -1 quando o par nao tem equacao. para_dado = true leva a solucao do motor a
// unidade do dado; false faz o caminho de volta, para um b do ajuste que vai reentrar num
// desenho remontado (b'_alvo = b_alvo - coef b_j, b'_j = s_j b_j).
template <class F>
void desfaz_escala_b(const VoltaB& v, F pos, double* b, bool para_dado) {
  for (std::size_t j = 0; j < v.escala.size(); j++) {
    if (v.escala[j] == 1.0 && v.absorve[j].empty()) continue;
    const int ij = pos(j);
    if (ij < 0) continue;
    if (para_dado) b[ij] /= v.escala[j];
    for (const auto& ac : v.absorve[j]) {
      const int ia = pos(ac.first);
      if (ia < 0) continue;
      if (para_dado) b[ia] += ac.second * b[ij];
      else b[ia] -= ac.second * b[ij];
    }
    if (!para_dado) b[ij] *= v.escala[j];
  }
}
// 2 soma log(escala) das colunas padronizadas com equacao: o -2logL do REML na escala do dado
// e o da padronizada mais isto (log|X'V^-1 X| muda por 2 log|det T| com X' = X T, e o centro
// e um cisalhamento, de determinante 1)
double logdet_escala(const std::vector<double>& escala);
// As colunas de posto completo sobre as linhas ativas, na ordem dada: a coluna j fica quando
// a distancia dela ao espaco das que ficaram antes passa de tol * max(||x_j||, 1), o criterio
// do Gram-Schmidt de antes (posto_completo), entao a coluna que sai de um conjunto dependente
// e a ULTIMA dele, e o posto nunca passa do numero de linhas ativas. A distancia sai da
// Cholesky de X'X na ordem das colunas quando ela e confiavel, e do residuo EXPLICITO
// x_j - X_K c pelas equacoes semi-normais corrigidas (Bjorck 1987) quando nao e, com cada
// decisao certificada; se uma nao se certifica (desenho numericamente quase singular), o posto
// inteiro sai do Gram-Schmidt denso de antes, montado coluna a coluna (ver fixos.cpp).
// relativo = true troca o limiar por tol * ||x_j||, o criterio do QR com pivo limitado do R
// (dqrdc2, o qr() padrao) que os ajustadores do limiar e da sobrevivencia usavam. info, se
// dado, recebe por qual rota as colunas foram decididas; nomes da o nome da coluna a mensagem.
struct PostoInfo {
  std::size_t rapidas = 0;      // decididas pela Cholesky
  std::size_t explicitas = 0;   // decididas pelo residuo explicito
  std::size_t correcoes = 0;    // correcoes das equacoes semi-normais, no todo
  bool densa = false;           // o posto saiu da rota densa de reserva
  std::size_t coluna_densa = 0; // a coluna que a rota esparsa nao certificou
};
void posto_esparso(const Csc& x, const std::vector<char>& ativa, double tol,
                   std::vector<std::size_t>& fica, std::vector<std::size_t>& sai,
                   bool relativo = false, PostoInfo* info = nullptr,
                   const std::vector<std::string>* nomes = nullptr);
// O posto dos desenhos multicaracter (model_mt(), model_ar1()): o de posto_esparso numa ordem
// em que as colunas com registro em TODAS as caracteristicas vem depois das incompletas (o
// intercepto na frente), para a referencia que cai ser um nivel completo sempre que existe um.
// obs e a mascara linha * t + caracteristica das celulas observadas. fica e sai na numeracao
// de fx.x, crescentes; motivo, por coluna de sai, "no record" ou "aliased" (ver fixos.cpp).
void posto_referencia_completa(const FixosCandidatos& fx, const std::vector<char>& usa,
                               const std::vector<char>& obs, std::size_t t, double tol,
                               std::vector<std::size_t>& fica, std::vector<std::size_t>& sai,
                               std::vector<std::string>& motivo);
// X so com as colunas que ficam, por linhas (a forma que os desenhos guardam)
Csr fixas_por_linha(const Csc& x, const std::vector<std::size_t>& fica);
// a mesma matriz por colunas, para o posto por caracteristica dos desenhos multicaracter
Csc fixas_por_coluna(const Csr& x);
// uma linha de X em denso, para as rotas de referencia pela forma V
void linha_densa(const Csr& x, std::size_t i, double* out);

}  // namespace br

#endif
