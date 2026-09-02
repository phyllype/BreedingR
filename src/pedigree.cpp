// Pedigree: ordenacao topologica, endogamia de Meuwissen e Luo (1992), e A^-1 de Henderson (1976).
//
// Este e o modulo onde um erro passa despercebido com mais facilidade, porque A^-1 sai
// simetrica, positiva-definida e plausivel mesmo quando esta errada. Dois pontos concentram
// isso:
//
//  1. a VARIANCIA MENDELIANA depende de quantos pais sao conhecidos — 1/2 com os dois,
//     3/4 com um, 1 com nenhum — e a endogamia dos pais entra nela. Usar 1/2 sempre e o erro
//     classico, e ele so aparece num pedigree com endogamia de verdade.
//  2. as contribuicoes ACUMULAM na mesma posicao: o d/4 entre um casal chega por cada filho.
//     Uma montagem em que o ultimo valor vence produz uma matriz que ainda converge, para o
//     lugar errado.

#include "mme.h"
#include <unordered_map>
#include <queue>

namespace br {

// Constroi o pedigree em ordem topologica.
//
// A ordenacao e por profundidade ITERATIVA e nao recursiva: um pedigree real tem dezenas de
// milhares de animais e a recursao estoura a pilha muito antes disso.
Pedigree constroi_pedigree(const std::vector<std::string>& id0,
                           const std::vector<std::string>& pa0,
                           const std::vector<std::string>& ma0,
                           const std::vector<std::string>& mf,
                           const std::vector<double>& gama) {
  if (mf.size() != gama.size())
    throw Erro("metafounders and gamma with different lengths");
  for (double g : gama)
    if (!(g > 0.0) || !(g < 2.0))
      throw Erro("gamma outside (0, 2): a metafounder needs positive base variance and a positive Mendelian variance for its offspring");

  // metafundadores viram linhas-base VIRTUAIS, prefixadas ao pedigree
  std::vector<std::string> id(mf);
  id.insert(id.end(), id0.begin(), id0.end());
  std::vector<std::string> pa(mf.size(), "0"), ma(mf.size(), "0");
  pa.insert(pa.end(), pa0.begin(), pa0.end());
  ma.insert(ma.end(), ma0.begin(), ma0.end());

  const std::size_t n = id.size();
  if (pa.size() != n || ma.size() != n) throw Erro("pedigree with columns of different lengths");
  if (n == 0) throw Erro("empty pedigree");

  std::unordered_map<std::string, std::size_t> pos;
  pos.reserve(n * 2);
  for (std::size_t i = 0; i < n; i++) {
    if (id[i].empty() || id[i] == "0") throw Erro("there is an animal with an empty or '0' identifier");
    if (!pos.emplace(id[i], i).second)
      throw Erro(i < mf.size()
                     ? "metafounder '" + id[i] + "' repeats or collides with an animal id"
                     : "repeated animal in the pedigree: '" + id[i] + "'");
  }

  auto liga = [&](const std::string& s) -> std::int64_t {
    if (s.empty() || s == "0" || s == "NA") return -1;
    auto it = pos.find(s);
    // Um pai citado e ausente NAO pode virar desconhecido em silencio: isso muda a variancia
    // mendeliana do filho e o parentesco de toda a descendencia.
    if (it == pos.end())
      throw Erro("parent '" + s + "' is cited and has no line of its own in the pedigree");
    return static_cast<std::int64_t>(it->second);
  };
  std::unordered_map<std::string, double> gmap;
  for (std::size_t k = 0; k < mf.size(); k++) gmap[mf[k]] = gama[k];

  std::vector<std::int64_t> p(n), m(n);
  for (std::size_t i = 0; i < n; i++) { p[i] = liga(pa[i]); m[i] = liga(ma[i]); }

  // profundidade por varredura iterativa, com deteccao de ciclo
  std::vector<int> estado(n, 0);         // 0 nao visitado, 1 na pilha, 2 pronto
  std::vector<std::size_t> ordem;
  ordem.reserve(n);
  std::vector<std::pair<std::size_t, int>> pilha;
  for (std::size_t s = 0; s < n; s++) {
    if (estado[s]) continue;
    pilha.push_back({s, 0});
    while (!pilha.empty()) {
      auto& [v, fase] = pilha.back();
      if (fase == 0) {
        if (estado[v] == 1) throw Erro("the pedigree has a cycle at '" + id[v] + "'");
        if (estado[v] == 2) { pilha.pop_back(); continue; }
        estado[v] = 1;
        fase = 1;
        if (p[v] >= 0) { pilha.push_back({static_cast<std::size_t>(p[v]), 0}); continue; }
      }
      if (fase == 1) {
        fase = 2;
        if (m[v] >= 0) { pilha.push_back({static_cast<std::size_t>(m[v]), 0}); continue; }
      }
      estado[v] = 2;
      ordem.push_back(v);
      pilha.pop_back();
    }
  }

  std::vector<std::size_t> novo(n);
  for (std::size_t k = 0; k < ordem.size(); k++) novo[ordem[k]] = k;

  Pedigree out;
  out.ids.resize(n);
  out.pai.resize(n);
  out.mae.resize(n);
  out.eh_mf.assign(n, 0);
  out.gama.assign(n, 0.0);
  for (std::size_t k = 0; k < ordem.size(); k++) {
    const std::size_t v = ordem[k];
    out.ids[k] = id[v];
    out.pai[k] = p[v] >= 0 ? static_cast<std::int64_t>(novo[p[v]]) : -1;
    out.mae[k] = m[v] >= 0 ? static_cast<std::int64_t>(novo[m[v]]) : -1;
    auto it = gmap.find(id[v]);
    if (it != gmap.end()) { out.eh_mf[k] = 1; out.gama[k] = it->second; }
  }
  return out;
}

// Endogamia por Meuwissen e Luo (1992).
//
// Para cada animal com os dois pais conhecidos, sobe o ramo dos ancestrais comuns em vez de
// formar A. O custo e proporcional ao numero de ancestrais, nao a n^2.
std::vector<double> endogamia(const Pedigree& p) {
  const std::size_t n = p.ids.size();
  std::vector<double> f(n, 0.0);
  std::vector<double> l(n, 0.0);
  std::vector<double> d(n, 0.0);

  // variancia mendeliana de cada animal, ja com a endogamia dos pais
  auto mendel = [&](std::size_t i) {
    const std::int64_t s = p.pai[i], t = p.mae[i];
    if (s >= 0 && t >= 0) return 0.5 - 0.25 * (f[s] + f[t]);
    if (s >= 0) return 0.75 - 0.25 * f[s];
    if (t >= 0) return 0.75 - 0.25 * f[t];
    return 1.0;
  };

  std::priority_queue<std::size_t> fila;   // maior indice primeiro: o pedigree esta ordenado
  std::vector<char> na_fila(n, 0);

  for (std::size_t i = 0; i < n; i++) {
    if (!p.eh_mf.empty() && p.eh_mf[i]) {
      // metafundador: linha-base com a_ff = gamma, logo F = gamma - 1 e D = gamma.
      // Todas as formulas abaixo leem so f dos pais, entao ESTE pre-carregamento e a
      // unica mudanca que a generalizacao inteira pede.
      f[i] = p.gama[i] - 1.0;
      d[i] = p.gama[i];
      continue;
    }
    d[i] = mendel(i);
    const std::int64_t s = p.pai[i], t = p.mae[i];
    if (s < 0 || t < 0) { f[i] = 0.0; continue; }

    std::fill(l.begin(), l.end(), 0.0);
    l[i] = 1.0;
    fila.push(i);
    na_fila[i] = 1;
    double soma = 0.0;

    while (!fila.empty()) {
      const std::size_t j = fila.top();
      fila.pop();
      na_fila[j] = 0;
      const double lj = l[j];
      l[j] = 0.0;
      soma += lj * lj * d[j];
      const std::int64_t sj = p.pai[j], tj = p.mae[j];
      if (sj >= 0) {
        l[sj] += 0.5 * lj;
        if (!na_fila[sj]) { fila.push(static_cast<std::size_t>(sj)); na_fila[sj] = 1; }
      }
      if (tj >= 0) {
        l[tj] += 0.5 * lj;
        if (!na_fila[tj]) { fila.push(static_cast<std::size_t>(tj)); na_fila[tj] = 1; }
      }
    }
    f[i] = soma - 1.0;
    if (f[i] < 0.0) f[i] = 0.0;   // so arredondamento pode levar abaixo de zero
    d[i] = mendel(i);
  }
  return f;
}

// A^-1 de Henderson, no triangulo INFERIOR.
//
//   A^-1 = sum_i d_i^-1 v_i v_i',  v = (1, -1/2, -1/2) sobre (animal, pai, mae)
//
// O produto externo e montado a partir dos coeficientes JA FUNDIDOS por indice, e nao
// termo a termo. A razao e um caso que quebrou a primeira versao deste arquivo: quando pai
// e mae sao o MESMO animal, o vetor colapsa de (1, -1/2, -1/2) para (1, -1), e a
// contribuicao na diagonal do progenitor passa a ser 1/d e nao 3/4 de 1/d.
//
// Guardando as duas metades da matriz, somar (s,d) e (d,s) separadamente acerta esse caso
// por acidente. Guardando so o triangulo inferior, as duas caem na mesma posicao e uma
// delas se perde. Fundir os coeficientes antes acerta os dois casos sem ramo nenhum.
//
// Autofecundacao nao acontece em suino, mas acontece em planta e, sobretudo, acontece como
// ERRO DE DADO. Uma A^-1 errada por causa disso continua simetrica e positiva-definida.
Csc a_inversa(const Pedigree& p, const std::vector<double>& f) {
  const std::size_t n = p.ids.size();
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  li.reserve(n * 6); cj.reserve(n * 6); v.reserve(n * 6);

  // (indice, coeficiente) do animal e dos progenitores, com repetidos somados
  std::vector<std::pair<std::size_t, double>> vi;
  vi.reserve(3);

  for (std::size_t i = 0; i < n; i++) {
    const std::int64_t s = p.pai[i], t = p.mae[i];
    double di;
    if (!p.eh_mf.empty() && p.eh_mf[i])
      di = p.gama[i];
    else if (s >= 0 && t >= 0) di = 0.5  - 0.25 * (f[s] + f[t]);
    else if (s >= 0)           di = 0.75 - 0.25 * f[s];
    else if (t >= 0)           di = 0.75 - 0.25 * f[t];
    else                       di = 1.0;
    if (!(di > 0.0))
      throw Erro("non-positive Mendelian variance at '" + p.ids[i] +
                 "'; the pedigree has impossible inbreeding");
    const double a = 1.0 / di;

    vi.clear();
    auto junta = [&](std::size_t idx, double c) {
      for (auto& e : vi) if (e.first == idx) { e.second += c; return; }
      vi.emplace_back(idx, c);
    };
    junta(i, 1.0);
    if (s >= 0) junta(static_cast<std::size_t>(s), -0.5);
    if (t >= 0) junta(static_cast<std::size_t>(t), -0.5);

    for (std::size_t x = 0; x < vi.size(); x++)
      for (std::size_t y = 0; y <= x; y++) {
        std::size_t r = vi[x].first, c = vi[y].first;
        double val = a * vi[x].second * vi[y].second;
        if (r < c) std::swap(r, c);
        li.push_back(static_cast<std::uint32_t>(r));
        cj.push_back(static_cast<std::uint32_t>(c));
        v.push_back(val);
      }
  }
  return de_triplos(n, n, li, cj, v);
}

}  // namespace br
