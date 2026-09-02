#include "modelo.h"
#include <unordered_map>
#include <unordered_set>
#include <cstdio>

namespace br {

std::string Coluna::rotulo(std::size_t i) const {
  if (texto) return txt[i];
  const double x = num[i];
  // Um float que representa um inteiro exato perde o ".0". O limite 2^53 e onde o double
  // deixa de representar todo inteiro: acima dele o arredondamento ja nao e fiel e o rotulo
  // sai em decimal, menos util porem honesto.
  if (std::isfinite(x) && x == std::floor(x) && std::fabs(x) < 9007199254740992.0) {
    char b[32];
    std::snprintf(b, sizeof b, "%lld", static_cast<long long>(x));
    return b;
  }
  char b[64];
  std::snprintf(b, sizeof b, "%g", x);
  return b;
}

const Coluna* Tabela::acha(const std::string& n) const {
  for (std::size_t k = 0; k < nomes.size(); k++)
    if (nomes[k] == n) return &colunas[k];
  return nullptr;
}

std::vector<std::string> Tabela::rotulos(const std::string& n) const {
  const Coluna* c = acha(n);
  if (!c) throw Erro("no column '" + n + "' in the data");
  std::vector<std::string> out(nlin);
  for (std::size_t i = 0; i < nlin; i++) out[i] = c->rotulo(i);
  return out;
}

std::vector<double> Tabela::numerico(const std::string& n) const {
  const Coluna* c = acha(n);
  if (!c) throw Erro("no column '" + n + "' in the data");
  // Uma coluna textual pedida como numerica e ERRO declarado, nao conversao silenciosa que
  // inventaria niveis.
  if (c->texto) throw Erro("column '" + n + "' is text and was requested as numeric");
  return c->num;
}

std::vector<std::string> Modelo::nomes_theta() const {
  std::vector<std::string> out(ntheta);
  for (const Grupo& g : grupos) {
    // nomes dos coeficientes do grupo, na ordem dos slots
    std::vector<std::string> rot;
    for (std::size_t t : g.termos) {
      const Termo& tm = termos[t];
      for (std::size_t c = 0; c < tm.n_coef(); c++) {
        if (tm.n_coef() == 1) rot.push_back(tm.nome);
        else rot.push_back(tm.nome + "[" + std::to_string(c) + "]");
      }
    }
    for (std::size_t j = 0; j < g.dim; j++)
      for (std::size_t i = j; i < g.dim; i++) {
        const std::size_t k = g.theta_idx(i, j);
        out[k] = (i == j) ? "var(" + rot[i] + ")" : "cov(" + rot[i] + "," + rot[j] + ")";
      }
  }
  out[offset_residual] = "var(residual)";
  return out;
}

Modelo monta_modelo(const std::string& alvo, std::vector<Termo> termos,
                    const std::vector<std::pair<std::string, std::vector<std::string>>>& grupos,
                    bool tem_ausente, double codigo_ausente) {
  if (alvo.empty()) throw Erro("the model needs a trait");
  if (termos.empty()) throw Erro("the model needs at least one effect");

  std::unordered_map<std::string, std::size_t> idx;
  for (std::size_t k = 0; k < termos.size(); k++) {
    const Termo& t = termos[k];
    if (t.nome.empty()) throw Erro("there is a term without a name");
    if (t.coluna.empty()) throw Erro("term '" + t.nome + "' without a column");
    if (!idx.emplace(t.nome, k).second)
      throw Erro("term declared twice: '" + t.nome + "'");
    // Uma base num termo fixo nao tem semantica definida aqui, e um aninhamento num termo
    // aleatorio e ambiguo: os niveis viriam da classe e a covariancia do termo, e os dois
    // parariam de casar com kron(C, K). Recusados na declaracao.
    if (!t.aleatorio() && !t.base.empty())
      throw Erro("fixed term '" + t.nome + "' with base: base is for random terms");
    if (t.aleatorio() && !t.aninhado.empty() && !t.social)
      throw Erro("random term '" + t.nome + "' nested: the levels would come from the class and the covariance from the term");
    if (t.social && t.aninhado.empty())
      throw Erro("indirect term '" + t.nome + "' without the pen column: without knowing who lives with whom there is no indirect effect");
    if (t.social && !t.aleatorio())
      throw Erro("indirect term '" + t.nome + "' fixed: the indirect effect is genetic, and therefore random");
    // A diluicao reescreve a incidencia SOCIAL e nada mais: fora de um termo social ela
    // nao descreve coisa alguma, e um valor negativo amplificaria com o tamanho da baia.
    if (t.diluicao < 0.0 || !std::isfinite(t.diluicao))
      throw Erro("term '" + t.nome + "': dilution must be finite and >= 0");
    if (t.diluicao != 0.0 && !t.social)
      throw Erro("term '" + t.nome + "': dilution only applies to an indirect() term");
  }

  Modelo m;
  m.alvo = alvo;
  m.termos = std::move(termos);
  m.tem_ausente = tem_ausente;
  m.codigo_ausente = codigo_ausente;

  std::unordered_set<std::size_t> em_grupo;
  for (const auto& [nome, membros] : grupos) {
    Grupo g;
    g.nome = nome;
    if (membros.empty()) throw Erro("group '" + nome + "' without terms");
    for (const std::string& mn : membros) {
      auto it = idx.find(mn);
      if (it == idx.end()) throw Erro("group '" + nome + "' cites a nonexistent term '" + mn + "'");
      const std::size_t k = it->second;
      const Termo& t = m.termos[k];
      if (!t.aleatorio()) throw Erro("fixed term '" + mn + "' in a covariance group");
      if (!em_grupo.insert(k).second) throw Erro("term '" + mn + "' in two groups");
      // Todos os termos de um grupo dividem a MESMA estrutura, porque a penalidade e um
      // unico kron(C_g^-1, K^-1): misturar parentesco com diagonal nao tem K comum.
      if (g.termos.empty()) g.estrutura = t.estrutura;
      else if (g.estrutura != t.estrutura)
        throw Erro("group '" + nome + "' mixes structures: the penalty is a single kron(C,K)");
      g.termos.push_back(k);
      g.dim += t.n_coef();
    }
    m.grupos.push_back(std::move(g));
  }

  // termo aleatorio sem grupo declarado ganha grupo proprio
  for (std::size_t k = 0; k < m.termos.size(); k++) {
    if (!m.termos[k].aleatorio() || em_grupo.count(k)) continue;
    Grupo g;
    g.nome = m.termos[k].nome;
    g.termos = {k};
    g.dim = m.termos[k].n_coef();
    g.estrutura = m.termos[k].estrutura;
    m.grupos.push_back(std::move(g));
  }
  if (m.grupos.empty()) throw Erro("the model has no random effect at all");

  // numeracao de theta: os grupos, depois o residual
  std::size_t off = 0;
  for (Grupo& g : m.grupos) {
    g.offset = off;
    g.nparam = g.dim * (g.dim + 1) / 2;
    off += g.nparam;
  }
  m.offset_residual = off;
  m.ntheta = off + 1;
  return m;
}

DesenhoTermo monta_termo(const Modelo& m, std::size_t k, const Tabela& t,
                         const std::vector<std::string>* niveis_fixos) {
  const Termo& tm = m.termos[k];
  DesenhoTermo d;
  d.termo = k;
  d.nome = tm.nome;
  d.n_coef = tm.n_coef();

  const std::size_t nlin = t.nlin;
  d.casou.assign(nlin, 1);

  // niveis: do pedigree quando ha parentesco, senao dos dados na ordem de aparicao.
  // Ordem de aparicao e nao alfabetica: e o que mantem a saida comparavel com dado ja
  // renumerado, onde o nivel "10" vem depois do "9" e nao entre "1" e "2".
  std::vector<std::string> rot = t.rotulos(tm.coluna);
  if (niveis_fixos) {
    d.niveis = *niveis_fixos;
  } else {
    std::unordered_set<std::string> visto;
    for (const std::string& r : rot)
      if (visto.insert(r).second) d.niveis.push_back(r);
  }
  d.n_niveis = d.niveis.size();
  std::unordered_map<std::string, std::size_t> pos;
  for (std::size_t i = 0; i < d.niveis.size(); i++) pos.emplace(d.niveis[i], i);

  // termo SOCIAL: a incidencia da linha i marca os companheiros de baia, nao o animal.
  //
  // A baia vem de `aninhado`; os companheiros sao os animais DISTINTOS que aparecem naquela
  // baia nos dados. Um companheiro fora do conjunto de niveis (fora do pedigree, num termo
  // com parentesco) e ERRO declarado: soma-lo como zero afirmaria que o efeito social dele e
  // nulo, e descarta-lo mudaria o grupo de convivencia em silencio.
  if (tm.social) {
    const std::vector<std::string> baia = t.rotulos(tm.aninhado);
    // baia -> animais distintos nela
    std::unordered_map<std::string, std::vector<std::size_t>> membros;
    {
      std::unordered_map<std::string, std::unordered_set<std::string>> visto;
      for (std::size_t i = 0; i < nlin; i++) {
        auto it = pos.find(rot[i]);
        if (it == pos.end()) {
          throw Erro("indirect term '" + tm.nome + "': animal '" + rot[i] +
                     "' is not in the level set (pedigree)");
        }
        if (visto[baia[i]].insert(rot[i]).second)
          membros[baia[i]].push_back(it->second);
      }
    }
    std::vector<std::uint32_t> li2, cj2;
    std::vector<double> v2;
    for (std::size_t i = 0; i < nlin; i++) {
      const std::size_t proprio = pos.at(rot[i]);
      const std::vector<std::size_t>& mem = membros[baia[i]];
      // Diluicao (Bijma 2010, Genetics 186:1013-1028): cada companheiro entra com
      // (n_i - 1)^(-d), onde n_i - 1 e o numero de companheiros do registro i (o proprio
      // animal esta em `mem`, dai o -1). d = 0 e a soma do livro, coeficiente 1. A baia
      // de tamanho 1 nunca chega ao expoente: sem companheiro o laco abaixo nao executa
      // e a linha fica zero, entao o 0^-d jamais e avaliado.
      const std::size_t ncomp = mem.size() - 1;
      const double peso = (tm.diluicao != 0.0 && ncomp > 1)
                              ? std::pow(static_cast<double>(ncomp), -tm.diluicao)
                              : 1.0;
      for (std::size_t m2 : mem) {
        if (m2 == proprio) continue;
        li2.push_back(static_cast<std::uint32_t>(i));
        cj2.push_back(static_cast<std::uint32_t>(m2));
        v2.push_back(peso);
      }
    }
    d.z = de_triplos(nlin, d.n_niveis, li2, cj2, v2);
    return d;
  }

  // valores da base: 1 para escalar, ou as colunas declaradas
  std::vector<std::vector<double>> base;
  if (tm.base.empty()) {
    if (tm.efeito == Efeito::Covariavel && tm.aninhado.empty() && !tm.aleatorio()) {
      // covariavel simples: uma coluna, um coeficiente, nivel unico
      base.push_back(t.numerico(tm.coluna));
      d.niveis = {tm.coluna};
      d.n_niveis = 1;
      pos.clear();
      pos.emplace(tm.coluna, 0);
      // toda linha casa: o "nivel" e a propria covariavel
      Csc z(nlin, 1);
      std::vector<std::uint32_t> li, cj;
      std::vector<double> v;
      for (std::size_t i = 0; i < nlin; i++)
        if (base[0][i] != 0.0) {
          li.push_back(static_cast<std::uint32_t>(i));
          cj.push_back(0);
          v.push_back(base[0][i]);
        }
      d.z = de_triplos(nlin, 1, li, cj, v);
      return d;
    }
    base.push_back(std::vector<double>(nlin, 1.0));
  } else {
    for (const std::string& b : tm.base) base.push_back(t.numerico(b));
  }

  // covariavel aninhada: o nivel vem da classe de aninhamento, o valor da propria coluna
  const std::vector<std::string>* rot_nivel = &rot;
  std::vector<std::string> rot_aninhado;
  if (!tm.aninhado.empty()) {
    rot_aninhado = t.rotulos(tm.aninhado);
    rot_nivel = &rot_aninhado;
    if (!niveis_fixos) {
      d.niveis.clear();
      std::unordered_set<std::string> visto;
      for (const std::string& r : rot_aninhado)
        if (visto.insert(r).second) d.niveis.push_back(r);
      d.n_niveis = d.niveis.size();
      pos.clear();
      for (std::size_t i = 0; i < d.niveis.size(); i++) pos.emplace(d.niveis[i], i);
    }
    base.clear();
    base.push_back(t.numerico(tm.coluna));
  }

  const std::size_t ncol = d.n_niveis * d.n_coef;
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  li.reserve(nlin * d.n_coef);
  cj.reserve(nlin * d.n_coef);
  v.reserve(nlin * d.n_coef);

  for (std::size_t i = 0; i < nlin; i++) {
    auto it = pos.find((*rot_nivel)[i]);
    if (it == pos.end()) {
      // Linha sem nivel no conjunto NAO entra com incidencia zero: isso afirmaria que o
      // efeito dela e exatamente zero e puxaria a estimativa para baixo em silencio. Ela e
      // marcada e quem monta o desenho a exclui, reportando a contagem.
      d.casou[i] = 0;
      continue;
    }
    const std::size_t nivel = it->second;
    for (std::size_t c = 0; c < d.n_coef; c++) {
      const double x = base[c][i];
      if (x == 0.0) continue;
      // coluna = coeficiente * n_niveis + nivel: o NIVEL varia mais rapido, e e essa
      // convencao que faz a penalidade sair exatamente kron(C^-1, K^-1)
      li.push_back(static_cast<std::uint32_t>(i));
      cj.push_back(static_cast<std::uint32_t>(c * d.n_niveis + nivel));
      v.push_back(x);
    }
  }
  d.z = de_triplos(nlin, ncol, li, cj, v);
  return d;
}

void posto_completo(const Densa& x, double tol, std::vector<std::size_t>& fica,
                    std::vector<std::size_t>& sai) {
  fica.clear();
  sai.clear();
  const std::size_t n = x.nlin, p = x.ncol;
  std::vector<std::vector<double>> ortog;
  for (std::size_t j = 0; j < p; j++) {
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; i++) v[i] = x.at(i, j);
    double norma0 = 0.0;
    for (double a : v) norma0 += a * a;
    norma0 = std::sqrt(norma0);
    for (const auto& b : ortog) {
      double proj = 0.0;
      for (std::size_t i = 0; i < n; i++) proj += v[i] * b[i];
      for (std::size_t i = 0; i < n; i++) v[i] -= proj * b[i];
    }
    double norma = 0.0;
    for (double a : v) norma += a * a;
    norma = std::sqrt(norma);
    if (norma > tol * std::max(norma0, 1.0)) {
      for (double& a : v) a /= norma;
      ortog.push_back(std::move(v));
      fica.push_back(j);
    } else {
      sai.push_back(j);
    }
  }
}

}  // namespace br
