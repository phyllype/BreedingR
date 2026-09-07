// Estruturas numericas do pacote. Nada aqui depende do R nem de biblioteca externa.
//
// Sem Rcpp e sem RcppEigen de proposito: os dois sao GPL e amarrariam a licenca do pacote.
// A API C do proprio R basta para a travessia, e a numerica e escrita aqui.
//
// As escolhas de layout tem motivo medido:
//
//  - denso e PLANO, `linha * ncol + coluna`. Um vector<vector<double>> custa um salto de
//    ponteiro por linha e impede qualquer laco em bloco. O bloco denso final do fator tem
//    milhares de colunas num modelo de passo unico; ali isso nao e estilo.
//  - indice de linha em uint32_t. Metade da memoria de indice. O que decide a velocidade
//    numa esparsa e cache, e um indice de 8 bytes gasta o dobro de linhas de cache por
//    nao-zero. O limite de 4,29 bilhoes de linhas esta longe de qualquer pedigree.
//  - triangulo inferior nas simetricas. A^-1 e a matriz de coeficientes sao simetricas;
//    guardar metade poupa memoria e trafego.

#ifndef BREEDINGR_ESTRUTURAS_H
#define BREEDINGR_ESTRUTURAS_H

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace br {

// Erro do dominio. Nunca escapa para o R como excecao: a casca converte em Rf_error, porque
// uma excecao de C++ atravessando a fronteira do R deixa o interpretador em estado indefinido.
struct Erro : std::runtime_error {
  explicit Erro(const std::string& m) : std::runtime_error(m) {}
};

// ------------------------------------------------------------------------------ densa

struct Densa {
  std::size_t nlin = 0, ncol = 0;
  std::vector<double> dados;

  Densa() = default;
  Densa(std::size_t l, std::size_t c) : nlin(l), ncol(c), dados(l * c, 0.0) {}

  inline double  at(std::size_t i, std::size_t j) const { return dados[i * ncol + j]; }
  inline double& at(std::size_t i, std::size_t j)       { return dados[i * ncol + j]; }
  inline const double* linha(std::size_t i) const { return dados.data() + i * ncol; }
  inline double*       linha(std::size_t i)       { return dados.data() + i * ncol; }

  double maior_assimetria() const {
    double w = 0.0;
    for (std::size_t i = 0; i < std::min(nlin, ncol); i++)
      for (std::size_t j = 0; j < i; j++)
        w = std::max(w, std::fabs(at(i, j) - at(j, i)));
    return w;
  }
};

// ------------------------------------------------------------------------------ esparsa

// Compressed Sparse Column, com as linhas CRESCENTES dentro de cada coluna.
struct Csc {
  std::size_t nlin = 0, ncol = 0;
  std::vector<std::size_t> colptr;   // ncol + 1
  std::vector<std::uint32_t> linha;
  std::vector<double> valor;

  Csc() = default;
  Csc(std::size_t l, std::size_t c) : nlin(l), ncol(c), colptr(c + 1, 0) {}

  std::size_t nnz() const { return valor.size(); }

  // Elemento (i, j), zero se ausente. Busca binaria dentro da coluna.
  double get(std::size_t i, std::size_t j) const {
    if (i >= nlin || j >= ncol) return 0.0;
    std::size_t lo = colptr[j], hi = colptr[j + 1];
    const std::uint32_t alvo = static_cast<std::uint32_t>(i);
    auto it = std::lower_bound(linha.begin() + lo, linha.begin() + hi, alvo);
    if (it != linha.begin() + hi && *it == alvo)
      return valor[static_cast<std::size_t>(it - linha.begin())];
    return 0.0;
  }

  Densa densa() const {
    Densa d(nlin, ncol);
    for (std::size_t j = 0; j < ncol; j++)
      for (std::size_t k = colptr[j]; k < colptr[j + 1]; k++)
        d.at(linha[k], j) = valor[k];
    return d;
  }

  // Espelha o triangulo guardado, para quem precisa da matriz cheia.
  Densa densa_simetrica() const {
    Densa d = densa();
    for (std::size_t j = 0; j < ncol; j++)
      for (std::size_t k = colptr[j]; k < colptr[j + 1]; k++) {
        std::size_t i = linha[k];
        if (i != j) d.at(j, i) = valor[k];
      }
    return d;
  }
};

// Constroi a CSC a partir de triplos, SOMANDO duplicados.
//
// Somar e obrigatorio e nao conveniencia: o A^-1 de Henderson escreve varias vezes na mesma
// posicao — o d/4 entre um casal chega por CADA filho — e a matriz so fica certa se essas
// contribuicoes acumularem. Uma montagem em que o ultimo vence apaga parte do G^-1 em
// silencio e ainda converge, para o lugar errado.
//
// Duas passagens de contagem, O(nnz), sem ordenacao por comparacao: primeiro estabiliza por
// linha, depois por coluna; como a entrada ja sai ordenada por linha e a contagem e estavel,
// dentro de cada coluna as linhas saem crescentes.
Csc de_triplos(std::size_t nlin, std::size_t ncol,
               const std::vector<std::uint32_t>& li,
               const std::vector<std::uint32_t>& cj,
               const std::vector<double>& v);

// ------------------------------------------------------------------------ densa: fatoracao

// Cholesky densa em bloco, no lugar. O triangulo inferior de `s` vira L. Devolve false se a
// matriz nao for positiva-definida — o que, neste engine, e informacao sobre theta e nao
// falha a reportar como tal.
bool chol_densa(Densa& s);

// inv(L L') a partir de L triangular inferior, sem nunca formar L L'.
Densa inv_do_fator(const Densa& l);

// Inversa de uma simetrica positiva-definida, pela Cholesky.
Densa inv_pd(const Densa& s);

// Inversa geral por Gauss-Jordan com pivotamento. Existe para os testes de exatidao: e o
// caminho independente contra o qual a rota da Cholesky e conferida.
Densa inv_geral(const Densa& a);

// log|A| de uma positiva-definida, ou NaN se nao for.
double logdet_pd(const Densa& a);
// Decomposicao espectral de uma simetrica pequena (Jacobi): a rota para quando a matriz
// pode ser singular e ainda assim precisa ser usada. Ver o corpo em estruturas.cpp.
void jacobi_sim(const Densa& a, std::vector<double>& ev, Densa& u);

}  // namespace br

#endif
