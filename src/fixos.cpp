// O desenho dos efeitos fixos, esparso.
//
// O desenho de antes era DENSO: a expansao densa de cada termo, n x p na tabela inteira, outra
// copia nas linhas usadas, o Gram-Schmidt guardando cada coluna ortogonalizada como um vetor de
// n, e o X final, perto de 5 n p doubles e O(n p^2) numa thread so. Medido em dado real:
// 209 646 linhas x 675 colunas fixas (um grupo contemporaneo de 673 niveis e duas covariaveis)
// = 4,95 GB de conjunto de trabalho num model(maxiter = 0).
//
// Aqui X vive por colunas (CSC) enquanto se decide o posto e por linhas (CSR) depois. Uma linha
// tem um nao-zero por termo fixo, entao X custa O(n x termos) e nao O(n p). O posto sai de um
// fator R de p x p (denso: p^2 / 2 doubles) e de vetores de n temporarios; uma base de n x
// posto so existe na rota densa de reserva, quando o desenho e numericamente quase singular
// (ver posto_esparso).

#include "mme.h"

#include <cstdio>
#include <limits>

namespace br {

FixosCandidatos colunas_fixas(const Modelo& m, const Tabela& t) {
  FixosCandidatos f;
  const std::size_t n = t.nlin;
  Csc& x = f.x;
  x.nlin = n;
  x.colptr.assign(1, 0);
  // intercepto: 1 em toda linha da tabela
  x.linha.reserve(n);
  x.valor.reserve(n);
  for (std::size_t i = 0; i < n; i++) {
    x.linha.push_back(static_cast<std::uint32_t>(i));
    x.valor.push_back(1.0);
  }
  x.colptr.push_back(x.linha.size());
  f.nomes.push_back("intercept");
  for (std::size_t k = 0; k < m.termos.size(); k++) {
    if (m.termos[k].aleatorio()) continue;
    const DesenhoTermo dt = monta_termo(m, k, t, nullptr);
    x.linha.reserve(x.linha.size() + dt.z.nnz());
    x.valor.reserve(x.valor.size() + dt.z.nnz());
    for (std::size_t j = 0; j < dt.z.ncol; j++) {
      for (std::size_t p = dt.z.colptr[j]; p < dt.z.colptr[j + 1]; p++) {
        x.linha.push_back(dt.z.linha[p]);
        x.valor.push_back(dt.z.valor[p]);
      }
      x.colptr.push_back(x.linha.size());
      f.nomes.push_back(dt.nome + "=" + dt.niveis[j % std::max<std::size_t>(dt.n_niveis, 1)]);
    }
  }
  x.ncol = f.nomes.size();
  return f;
}

namespace {

// R triangular superior m x m guardada por colunas compactadas: a coluna q ocupa as q + 1
// posicoes a partir de q (q + 1) / 2, linhas 0..q.
// c <- R^-T c (R' y = c, para frente)
void frente(const std::vector<double>& R, std::size_t m, std::vector<double>& c) {
  for (std::size_t q = 0; q < m; q++) {
    const double* col = R.data() + q * (q + 1) / 2;
    double s = c[q];
    for (std::size_t i = 0; i < q; i++) s -= col[i] * c[i];
    c[q] = s / col[q];
  }
}
// c <- R^-1 c (R x = y, para tras, por colunas)
void tras(const std::vector<double>& R, std::size_t m, std::vector<double>& c) {
  for (std::size_t q = m; q-- > 0;) {
    const double* col = R.data() + q * (q + 1) / 2;
    c[q] /= col[q];
    const double cq = c[q];
    for (std::size_t i = 0; i < q; i++) c[i] -= col[i] * cq;
  }
}

std::string gb(double bytes) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.1f GB", bytes / 1e9);
  return buf;
}

std::string nome_coluna(const std::vector<std::string>* nomes, std::size_t j) {
  if (nomes && j < nomes->size()) return "'" + (*nomes)[j] + "'";
  return "number " + std::to_string(j + 1);
}

// A ROTA DENSA DE RESERVA: o Gram-Schmidt modificado de antes (posto_completo, em modelo.cpp),
// conta por conta e na mesma ordem de soma, sobre as linhas ativas e da coluna 0 em diante.
// So a forma muda: cada coluna e montada da CSC quando chega (o X denso n x p nunca existe) e
// so as ortogonalizadas que ficam sao guardadas. Com o invariante do posto (nao passa do
// numero de linhas ativas), que o de antes nao tinha: num desenho com menos linhas do que
// colunas o ruido da ortogonalizacao podia manter uma coluna a mais do que o posto possivel.
void posto_denso_reserva(const Csc& x, const std::vector<std::uint32_t>& nova, std::size_t na,
                         double tol, bool relativo, std::vector<std::size_t>& fica,
                         std::vector<std::size_t>& sai) {
  fica.clear();
  sai.clear();
  const std::uint32_t fora = std::numeric_limits<std::uint32_t>::max();
  std::vector<std::vector<double>> ortog;
  for (std::size_t j = 0; j < x.ncol; j++) {
    std::vector<double> v(na, 0.0);
    for (std::size_t k = x.colptr[j]; k < x.colptr[j + 1]; k++) {
      const std::uint32_t i = nova[x.linha[k]];
      if (i != fora) v[i] = x.valor[k];
    }
    double norma0 = 0.0;
    for (double a : v) norma0 += a * a;
    norma0 = std::sqrt(norma0);
    if (ortog.size() >= na) {
      sai.push_back(j);
      continue;
    }
    for (const auto& b : ortog) {
      double proj = 0.0;
      for (std::size_t i = 0; i < na; i++) proj += v[i] * b[i];
      for (std::size_t i = 0; i < na; i++) v[i] -= proj * b[i];
    }
    double norma = 0.0;
    for (double a : v) norma += a * a;
    norma = std::sqrt(norma);
    if (norma > tol * (relativo ? norma0 : std::max(norma0, 1.0))) {
      for (double& a : v) a /= norma;
      ortog.push_back(std::move(v));
      fica.push_back(j);
    } else {
      sai.push_back(j);
    }
  }
}

}  // namespace

// O CRITERIO e o do Gram-Schmidt modificado de antes (posto_completo, em modelo.cpp), coluna a
// coluna na ordem dada: a coluna j fica quando a distancia dela ao espaco das que ficaram
// antes passa de tol * max(||x_j||, 1). O conjunto que fica e o primeiro conjunto
// independente na ordem das colunas (num conjunto dependente sai a ULTIMA). Toda decisao que
// a rota esparsa toma e CERTIFICADA (abaixo); quando ela nao consegue certificar uma, o posto
// inteiro sai da rota densa de reserva, que e o Gram-Schmidt de antes. Nas duas o posto nunca
// passa do numero de linhas ativas.
//
// A CONTA RAPIDA. Com Z_K as colunas que ficaram, escaladas a norma 1 (o espaco gerado nao muda
// e o fator fica melhor condicionado), R o fator triangular de Z_K (R'R = Z_K'Z_K) e z_j a
// coluna j escalada, g = Z_K' z_j, l = R^-T g e c = R^-1 l = (R'R)^-1 g (os coeficientes de
// minimos quadrados), a distancia relativa ao quadrado e
//
//   d = ||z_j||^2 - ||l||^2,
//
// o passo da Cholesky de Z'Z na ordem das colunas: g sai das linhas do suporte de z_j (uma
// coluna de classe toca so os registros do nivel) e o resto e O(|K|^2). A Cholesky e estavel
// para tras: o d calculado e o d exato de Z'Z + E, |E| da ordem de |K| eps, e em primeira ordem
// o d de Z'Z fica a |K| eps (1 + ||c||)^2 dele. Mas d e a diferenca de dois numeros perto de 1,
// entao ele so DECIDE acima de 1e-6 ||z_j||^2 e de 1e6 vezes a estimativa do erro,
// |K| eps (||z_j||^2 + ||c||^2), com a distancia sqrt(d) ||x_j|| longe do limiar (fora de
// 1 +- 1e-3). Num fator mal condicionado o c calculado cresce, a estimativa com ele, e a coluna
// vai a conta explicita. A coluna que fica entra com [l; sqrt(d)] no fator.
//
// A CONTA EXPLICITA, para o resto (as dependentes, as quase dependentes e as de c grande, que e
// o mal condicionado): o residuo r = z_j - Z_K c, partindo do c acima, com as correcoes das
// equacoes semi-normais (Bjorck 1987), c += R^-1 h, h = R^-T Z_K' r. O que cada passo
// certifica:
//   - ||r|| e a distancia de z_j a UM ponto do espaco, entao e um teto da distancia: com ||r||
//     no limiar ou abaixo dele a coluna SAI;
//   - ||h|| = ||P r|| (com R'R = Z_K'Z_K) e a parte de r que ainda esta no espaco: com
//     ||h||^2 <= 1e-6 ||r||^2 a distancia e sqrt(||r||^2 - ||h||^2) a 1e-6 relativo, a coluna
//     fica (ou sai) por ela e entra no fator como [R c + h; sqrt(||r||^2 - ||h||^2)], a coluna
//     que o QR de [Z_K z_j] teria.
// Ate quatro correcoes. Sem certificado, ou com uma correcao que nao reduz ||h|| (o fator mal
// condicionado demais para corrigir: as correcoes divergem ou param, e com o menor residuo
// como distancia uma coluna dependente ficava e arrastava as seguintes), a rota esparsa nao
// decide: o posto inteiro sai da rota densa de reserva, que custa o tempo do Gram-Schmidt de
// antes, O(n p^2) numa thread, e menos memoria (8 n posto bytes, nenhuma copia de X); acima de
// 4 GB (contados com o teto kmax, nao com o posto) o ajuste para com a mensagem. O que leva a
// reserva: colunas anteriores quase colineares (o fator mal condicionado) que geram a coluna so
// com coeficientes enormes, como x1, x2 = x1 + 1e-8 z e x3 = x1 + 1e8 (x2 - x1). Polinomios sem
// escala (ate o grau 8 em 1e3 + U(0, 1)) e quase copias sozinhas sao certificados pela rota
// esparsa (portao em test-desenho-fixo-esparso.R).
//
// O CUSTO da rota esparsa. A rapida custa O(|supp(z_j)| x termos + |K|^2) por coluna; a
// explicita, um produto Z_K c e um Z_K' r de O(nnz(X)) cada por passo (uma coluna dependente
// de um desenho bem condicionado sai no primeiro residuo). No todo O(nnz(X) x termos + p^3 / 3
// + e nnz(X)), com e os passos explicitos, o fator de p^2 / 2 doubles e dois vetores de n
// temporarios. O fator e DENSO em p e as substituicoes sao sequenciais: com 10 mil colunas
// fixas sao 0,4 GB e minutos numa thread; acima de 4 GB o ajuste para com a mensagem antes do
// trabalho. Os produtos da rota explicita rodam em paralelo com dono unico (um bloco fixo de
// linhas de r com a sua soma de quadrados, ou uma entrada de Z_K' r, por thread) e as somas em
// ordem fixa, e a rota rapida e sequencial: a decisao e o fator sao os mesmos com qualquer
// numero de threads.
void posto_esparso(const Csc& x, const std::vector<char>& ativa, double tol,
                   std::vector<std::size_t>& fica, std::vector<std::size_t>& sai,
                   bool relativo, PostoInfo* info, const std::vector<std::string>* nomes) {
  fica.clear();
  sai.clear();
  PostoInfo conta;
  if (info) *info = conta;
  const std::size_t p = x.ncol;
  if (p == 0) return;
  const std::uint32_t fora = std::numeric_limits<std::uint32_t>::max();
  std::vector<std::uint32_t> nova(x.nlin, fora);
  std::size_t na = 0;
  for (std::size_t i = 0; i < x.nlin; i++)
    if (ativa[i]) nova[i] = static_cast<std::uint32_t>(na++);

  // as colunas nas linhas ativas, renumeradas, e a norma de cada uma: a soma dos quadrados na
  // ordem das linhas, a mesma conta (e o mesmo numero) do Gram-Schmidt denso, cujos zeros nao
  // mudam a soma
  std::vector<std::size_t> cp(p + 1, 0);
  std::vector<std::uint32_t> cl;
  std::vector<double> cv;
  std::vector<double> norma(p, 0.0);
  std::size_t p_cheias = 0;
  for (std::size_t j = 0; j < p; j++) {
    double s = 0.0;
    for (std::size_t k = x.colptr[j]; k < x.colptr[j + 1]; k++) {
      const std::uint32_t i = nova[x.linha[k]];
      if (i == fora) continue;
      cl.push_back(i);
      cv.push_back(x.valor[k]);
      s += x.valor[k] * x.valor[k];
    }
    cp[j + 1] = cl.size();
    norma[j] = std::sqrt(s);
    if (norma[j] > 0.0 && std::isfinite(norma[j])) p_cheias++;
  }
  // o fator R guarda |K| (|K| + 1) / 2 doubles, |K| <= kmax = min(linhas ativas, colunas nao
  // nulas nas linhas ativas). A guarda e a reserva de R contam esse TETO, nao o posto (que so
  // se sabe no fim): um desenho com muitas candidatas dependentes (uma classe aninhada grande)
  // para aqui, e reserva memoria, pelo numero de candidatas. A guarda vem ANTES do trabalho
  // O(p^3), nao depois dele
  const double kmax = static_cast<double>(std::min(na, p_cheias));
  const double bytes_r = 4.0 * kmax * (kmax + 1.0);
  if (bytes_r > 4.0e9)
    throw Erro("the fixed-effect design has " + std::to_string(p_cheias) + " candidate "
               "columns with an entry on the " + std::to_string(na) + " used records. Finding "
               "its rank sizes a dense triangular factor for up to " +
               std::to_string(static_cast<std::size_t>(kmax)) + " kept columns (the smaller "
               "of the two counts), " + gb(bytes_r) + " (more than 4 GB), before O(p^3) work, "
               "even when many of these columns would be dropped as dependent. Fit the class "
               "with the most levels as a random effect, or merge its levels");
  // escaladas a norma 1 (so as de norma finita e positiva sao lidas), por colunas e por linhas
  for (std::size_t j = 0; j < p; j++) {
    const bool usavel = norma[j] > 0.0 && std::isfinite(norma[j]);
    for (std::size_t k = cp[j]; k < cp[j + 1]; k++) cv[k] = usavel ? cv[k] / norma[j] : 0.0;
  }
  std::vector<std::size_t> rp(na + 1, 0);
  for (std::uint32_t i : cl) rp[i + 1]++;
  for (std::size_t i = 0; i < na; i++) rp[i + 1] += rp[i];
  std::vector<std::uint32_t> rc(cl.size());
  std::vector<double> rv(cl.size());
  {
    std::vector<std::size_t> prox(rp.begin(), rp.end() - 1);
    for (std::size_t j = 0; j < p; j++)
      for (std::size_t k = cp[j]; k < cp[j + 1]; k++) {
        rc[prox[cl[k]]] = static_cast<std::uint32_t>(j);
        rv[prox[cl[k]]++] = cv[k];
      }
  }

  std::vector<char> mantida(p, 0);
  std::vector<std::size_t> K, pos(p, 0);
  std::vector<double> R;
  R.reserve(static_cast<std::size_t>(kmax * (kmax + 1.0) / 2.0));
  std::vector<double> v, r, cheio(p, 0.0), l, c, g;
  const int nth = threads();
  const bool paralelo = nth > 1 && cl.size() > 50000;
  (void) paralelo;   // sem OpenMP so os pragmas o leem
  const double eps = std::numeric_limits<double>::epsilon();
  // blocos FIXOS de linhas para a soma de ||r||^2: a mesma soma com qualquer numero de threads
  const std::size_t bloco = 4096;
  const std::size_t nblocos = (na + bloco - 1) / bloco;
  std::vector<double> parcial(nblocos, 0.0);

  // r = v - Z_K c, com c em cheio (indexado pela coluna), e ||r||^2; um bloco de linhas por dono
  auto residuo = [&]() -> double {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth) if (paralelo)
#endif
    for (long bb = 0; bb < static_cast<long>(nblocos); bb++) {
      const std::size_t i0 = static_cast<std::size_t>(bb) * bloco;
      const std::size_t i1 = std::min(na, i0 + bloco);
      double acc = 0.0;
      for (std::size_t i = i0; i < i1; i++) {
        double s = v[i];
        for (std::size_t e = rp[i]; e < rp[i + 1]; e++)
          if (mantida[rc[e]]) s -= rv[e] * cheio[rc[e]];
        r[i] = s;
        acc += s * s;
      }
      parcial[static_cast<std::size_t>(bb)] = acc;
    }
    double rr = 0.0;
    for (double a : parcial) rr += a;
    return rr;
  };
  // out = Z_K' r; uma entrada por dono
  auto transposto = [&](std::vector<double>& out) {
    const std::size_t m = K.size();
    out.assign(m, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nth) if (paralelo)
#endif
    for (long qq = 0; qq < static_cast<long>(m); qq++) {
      const std::size_t q = static_cast<std::size_t>(qq), j = K[q];
      double s = 0.0;
      for (std::size_t k = cp[j]; k < cp[j + 1]; k++) s += cv[k] * r[cl[k]];
      out[q] = s;
    }
  };

  for (std::size_t j = 0; j < p; j++) {
    const double s0 = norma[j];
    const double lim = tol * (relativo ? s0 : std::max(s0, 1.0));
    // nula, abaixo do limiar, NaN ou Inf: o Gram-Schmidt tambem a derrubava
    if (!(s0 > lim) || !std::isfinite(s0)) {
      sai.push_back(j);
      continue;
    }
    const std::size_t m = K.size();
    // o posto nao passa do numero de linhas ativas: com |K| = na toda coluna esta no espaco
    // (e nenhuma se certificaria: o residuo e so arredondamento, e iria a reserva densa)
    if (m >= na) {
      sai.push_back(j);
      continue;
    }
    // g = Z_K' z_j pelas linhas do suporte de z_j, e ||z_j||^2
    l.assign(m, 0.0);
    double gjj = 0.0;
    for (std::size_t k = cp[j]; k < cp[j + 1]; k++) {
      gjj += cv[k] * cv[k];
      const std::size_t i = cl[k];
      for (std::size_t e = rp[i]; e < rp[i + 1]; e++)
        if (mantida[rc[e]]) l[pos[rc[e]]] += rv[e] * cv[k];
    }
    frente(R, m, l);
    double ll = 0.0;
    for (double a : l) ll += a * a;
    c = l;
    tras(R, m, c);
    double cc = 0.0;
    for (double a : c) cc += a * a;
    const double dd = gjj - ll;
    const double erro = static_cast<double>(m + 1) * eps * (gjj + cc);
    double rho = 0.0;
    int decide = 0;   // 1 fica pela rota rapida, 2 pela explicita, -1 sai
    if (dd > 1e-6 * gjj && dd > 1e6 * erro) {
      rho = std::sqrt(dd);
      if (rho * s0 > lim * (1.0 + 1e-3)) decide = 1;
      else if (rho * s0 < lim * (1.0 - 1e-3)) decide = -1;
    }
    if (decide != 0) {
      conta.rapidas++;
    } else {
      if (v.empty()) {
        v.assign(na, 0.0);
        r.assign(na, 0.0);
      }
      for (std::size_t k = cp[j]; k < cp[j + 1]; k++) v[cl[k]] = cv[k];
      // o limiar ao quadrado na escala da coluna normalizada
      const double lim2 = (lim / s0) * (lim / s0);
      double hh_antes = std::numeric_limits<double>::infinity();
      for (int it = 0;; it++) {
        for (std::size_t q = 0; q < m; q++) cheio[K[q]] = c[q];
        const double rr = residuo();
        if (!std::isfinite(rr)) break;
        // ||r|| e um teto da distancia: no limiar ou abaixo dele a coluna sai
        if (!(rr > lim2)) {
          decide = -1;
          break;
        }
        double hh = 0.0;
        if (m > 0) {
          transposto(g);
          frente(R, m, g);
          for (double a : g) hh += a * a;
          if (!std::isfinite(hh)) break;
        }
        // ||P r|| <= 1e-3 ||r||: a distancia esta certificada
        if (hh <= 1e-6 * rr) {
          rho = std::sqrt(rr - hh);
          decide = rho * s0 > lim ? 2 : -1;
          break;
        }
        if (it == 4 || !(hh < hh_antes)) break;
        hh_antes = hh;
        tras(R, m, g);
        for (std::size_t q = 0; q < m; q++) c[q] += g[q];
        conta.correcoes++;
      }
      for (std::size_t k = cp[j]; k < cp[j + 1]; k++) v[cl[k]] = 0.0;
      if (decide == 0) {
        // sem certificado a rota esparsa nao decide: o posto inteiro sai da densa
        const double bytes_d = 8.0 * static_cast<double>(na) * (kmax + 1.0);
        if (bytes_d > 4.0e9)
          throw Erro("the fixed-effect design is numerically near-singular: the sparse rank "
                     "cannot settle fixed column " + nome_coluna(nomes, j) + " against the "
                     "columns before it, and the dense check that settles such a design would "
                     "need " + gb(bytes_d) + " (more than 4 GB) for " + std::to_string(na) +
                     " records. It happens when columns before it are nearly collinear (such as "
                     "a covariate and a near copy of it) and span this column only with very "
                     "large coefficients. Drop the near copy, or the column built from it");
        conta.densa = true;
        conta.coluna_densa = j;
        posto_denso_reserva(x, nova, na, tol, relativo, fica, sai);
        if (info) *info = conta;
        return;
      }
      conta.explicitas++;
    }
    if (decide < 0) {
      sai.push_back(j);
      continue;
    }
    // a coluna nova do fator: [l; sqrt(d)] pela rota rapida, [R c + h; sqrt(rr - hh)] pela
    // explicita (h, que esta em g, e R (c* - c): o que falta a R c para os minimos quadrados)
    const std::size_t ini = R.size();
    R.resize(ini + m + 1, 0.0);
    if (decide == 1) {
      for (std::size_t q = 0; q < m; q++) R[ini + q] = l[q];
    } else {
      for (std::size_t q = 0; q < m; q++) {
        const double* col = R.data() + q * (q + 1) / 2;
        for (std::size_t i = 0; i <= q; i++) R[ini + i] += col[i] * c[q];
      }
      for (std::size_t q = 0; q < m; q++) R[ini + q] += g[q];
    }
    R[ini + m] = rho;
    mantida[j] = 1;
    pos[j] = m;
    K.push_back(j);
    fica.push_back(j);
  }
  if (info) *info = conta;
}

Csr fixas_por_linha(const Csc& x, const std::vector<std::size_t>& fica) {
  Csr out;
  out.nlin = x.nlin;
  out.ncol = fica.size();
  out.linptr.assign(x.nlin + 1, 0);
  for (std::size_t jj = 0; jj < fica.size(); jj++)
    for (std::size_t k = x.colptr[fica[jj]]; k < x.colptr[fica[jj] + 1]; k++)
      out.linptr[x.linha[k] + 1]++;
  for (std::size_t i = 0; i < x.nlin; i++) out.linptr[i + 1] += out.linptr[i];
  out.coluna.assign(out.linptr[x.nlin], 0);
  out.valor.assign(out.linptr[x.nlin], 0.0);
  std::vector<std::size_t> prox(out.linptr.begin(), out.linptr.end() - 1);
  // colunas em ordem crescente: dentro de cada linha elas saem crescentes
  for (std::size_t jj = 0; jj < fica.size(); jj++)
    for (std::size_t k = x.colptr[fica[jj]]; k < x.colptr[fica[jj] + 1]; k++) {
      const std::size_t i = x.linha[k];
      out.coluna[prox[i]] = static_cast<std::uint32_t>(jj);
      out.valor[prox[i]++] = x.valor[k];
    }
  return out;
}

Csc fixas_por_coluna(const Csr& x) {
  Csc out(x.nlin, x.ncol);
  for (std::uint32_t j : x.coluna) out.colptr[j + 1]++;
  for (std::size_t j = 0; j < x.ncol; j++) out.colptr[j + 1] += out.colptr[j];
  out.linha.assign(x.nnz(), 0);
  out.valor.assign(x.nnz(), 0.0);
  std::vector<std::size_t> prox(out.colptr.begin(), out.colptr.end() - 1);
  for (std::size_t i = 0; i < x.nlin; i++)
    for (std::size_t e = x.linptr[i]; e < x.linptr[i + 1]; e++) {
      const std::size_t j = x.coluna[e];
      out.linha[prox[j]] = static_cast<std::uint32_t>(i);
      out.valor[prox[j]++] = x.valor[e];
    }
  return out;
}

void linha_densa(const Csr& x, std::size_t i, double* out) {
  for (std::size_t j = 0; j < x.ncol; j++) out[j] = 0.0;
  for (std::size_t e = x.linptr[i]; e < x.linptr[i + 1]; e++) out[x.coluna[e]] = x.valor[e];
}

}  // namespace br
