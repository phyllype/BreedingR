#include "modelo.h"
#include <unordered_map>
#include <unordered_set>
#include <cstdio>
#include <cstdlib>

namespace br {

bool rotulo_exato(double x) {
  return std::isfinite(x) && x == std::floor(x) && std::fabs(x) < 9007199254740992.0;
}

std::string rotulo_numero(double x) {
  // Um float que representa um inteiro exato perde o ".0". O limite 2^53 e onde o double
  // deixa de representar todo inteiro: acima dele o arredondamento ja nao e fiel e o rotulo
  // sai em decimal, menos util porem honesto.
  if (rotulo_exato(x)) {
    char b[32];
    std::snprintf(b, sizeof b, "%lld", static_cast<long long>(x));
    return b;
  }
  // Fora disso, a MENOR escrita com 15, 16 ou 17 digitos significativos que volta ao mesmo
  // double: "1.5" continua "1.5", e dois doubles distintos nunca dividem o rotulo. O %g de
  // antes guardava 6 digitos, e 123456.7 saia "123457", o rotulo do inteiro 123457: na
  // coluna de dados os dois viravam UM nivel sem aviso, e 1234567.5 saia "1.23457e+06", que
  // parecia o inteiro 1234570 escrito pelo R. Com 17 digitos todo double finito volta a si
  // mesmo. snprintf e strtod leem o mesmo separador decimal, entao a volta nao depende do
  // locale. NaN e Inf seguem o %g ("nan", "inf"); como nivel de classe sao recusados em
  // monta_termo.
  char b[48];
  if (!std::isfinite(x)) {
    std::snprintf(b, sizeof b, "%g", x);
    return b;
  }
  for (int p = 15; p <= 17; p++) {
    std::snprintf(b, sizeof b, "%.*g", p, x);
    if (std::strtod(b, nullptr) == x) break;
  }
  return b;
}

std::string Coluna::rotulo(std::size_t i) const {
  if (texto) return txt[i];
  return rotulo_numero(num[i]);
}

std::string inteiro_de_cientifico(const std::string& s) {
  // a forma do R: sinal opcional, um digito de 1 a 9, talvez um ponto com mais digitos, "e+"
  // e o expoente. A conta e feita nos DIGITOS, sem strtod, que depende do separador decimal
  // do locale.
  const std::size_t e = s.find("e+");
  if (e == std::string::npos || e + 2 >= s.size()) return "";
  const std::size_t sinal = (s[0] == '-') ? 1 : 0;
  std::size_t i = sinal;
  if (i >= e || s[i] < '1' || s[i] > '9') return "";
  std::string dig(1, s[i]);
  std::size_t frac = 0;
  i++;
  if (i < e) {
    if (s[i] != '.' || i + 1 >= e) return "";
    for (i++; i < e; i++) {
      if (s[i] < '0' || s[i] > '9') return "";
      dig.push_back(s[i]);
      frac++;
    }
    // o R nao escreve zero final na mantissa
    if (dig.back() == '0') return "";
  }
  // o expoente do R tem dois digitos no minimo, e zero a esquerda so para chegar a dois
  const std::size_t ne = s.size() - (e + 2);
  if (ne < 2 || (ne > 2 && s[e + 2] == '0')) return "";
  std::size_t expo = 0;
  for (std::size_t k = e + 2; k < s.size(); k++) {
    if (s[k] < '0' || s[k] > '9') return "";
    expo = expo * 10 + static_cast<std::size_t>(s[k] - '0');
    if (expo > 32) return "";
  }
  // inteiro so quando o expoente cobre as casas decimais
  if (expo < frac) return "";
  dig.append(expo - frac, '0');
  if (dig.size() > 16) return "";
  // E SO O QUE O as.character() ESCREVERIA. O R escolhe a notacao cientifica quando ela e
  // MAIS CURTA que a fixa (empate fica na fixa): 1e+05 (5 contra 6 de "100000"), 1.2e+07 (7
  // contra 8), mas 1200000 e nao 1.2e+06. Um texto como "1.23457e+06" le como 1234570, e o
  // as.character() desse numero e "1234570": e um nao inteiro arredondado a 6 digitos (o %g
  // de 1234567.5), e dizer que ele e o mesmo numero que "1234570" seria fundir dois animais
  // na mensagem.
  if (dig.size() <= s.size() - sinal) return "";
  long long v = 0;
  for (char c : dig) v = v * 10 + (c - '0');
  const double x = static_cast<double>(s[0] == '-' ? -v : v);
  if (!rotulo_exato(x)) return "";
  return rotulo_numero(x);
}

std::string dica_cientifica(const std::string& id, const std::vector<std::string>& outros) {
  if (id.empty()) return "";
  const std::string c = inteiro_de_cientifico(id);
  for (const std::string& o : outros) {
    const bool par = (!c.empty() && o == c) ||
                     (o.find("e+") != std::string::npos && inteiro_de_cientifico(o) == id);
    if (par)
      return " ('" + id + "' and '" + o + "' are the same number written two ways, the one "
             "with 'e+' as as.character() and factor() write a round number: write the ids "
             "the same way on both sides, as numbers or with "
             "format(x, scientific = FALSE, trim = TRUE))";
  }
  return "";
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
    if (t.aleatorio() && !t.aninhado.empty() && !t.social && !t.mgs)
      throw Erro("random term '" + t.nome + "' nested: the levels would come from the class and the covariance from the term");
    if (t.mgs && (!t.aleatorio() || !t.base.empty() || t.efeito != Efeito::Classe))
      throw Erro("sire term '" + t.nome + "' with mgs =: the maternal grandsire enters a random class term without base");
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

  // NIVEL AUSENTE e ERRO, no termo fixo de classe como no aleatorio. Coluna::rotulo faz de
  // qualquer numero um rotulo: o NA de uma coluna numerica (NaN no motor) virava o nivel
  // "nan", e as linhas sem nivel formavam UM nivel so, uma classe "cg=nan" no fixo ou um
  // animal "nan" no aleatorio, com efeito estimado e sem aviso (num grupo iid, "nan" ainda
  // pareava com o "nan" do outro termo). Num termo com parentesco o "nan" nao casava com o
  // pedigree e a linha saia calada. O NA de coluna textual ja para em tabela_do_R
  // (entrada.cpp), em qualquer linha, e esta e a mesma regra para a coluna numerica; Inf e
  // -Inf tambem nao sao nivel. A coluna conferida e a que da o nivel: a classe, o animal do
  // termo social, o pai do modelo pai / avo materno E a coluna do avo, ou a classe de
  // aninhamento de uma covariavel aninhada. O avo AUSENTE (NA numerico) virava o avo "nan",
  // que parava adiante com "maternal grandsire 'nan' is not in the level set", sem a linha;
  // o avo DESCONHECIDO e o "0" (ou texto vazio), que deixa so o pai. A covariavel simples
  // nao tem nivel, e a observacao ausente continua tirando so a linha (monta_desenho).
  {
    const bool so_valor = tm.base.empty() && tm.efeito == Efeito::Covariavel &&
                          tm.aninhado.empty() && !tm.aleatorio();
    std::vector<std::string> cols_nivel;
    if (!so_valor)
      cols_nivel.push_back((!tm.aninhado.empty() && !tm.social && !tm.mgs) ? tm.aninhado
                                                                          : tm.coluna);
    if (tm.mgs) cols_nivel.push_back(tm.aninhado);
    for (const std::string& col_nivel : cols_nivel) {
      const Coluna* cn = t.acha(col_nivel);
      if (!cn || cn->texto) continue;
      std::size_t n_aus = 0, primeira = 0;
      bool viu_na = false, viu_inf = false;
      for (std::size_t i = 0; i < nlin; i++) {
        const double x = cn->num[i];
        if (std::isfinite(x)) continue;
        if (n_aus++ == 0) primeira = i;
        (std::isnan(x) ? viu_na : viu_inf) = true;
      }
      if (n_aus > 0)
        throw Erro("term '" + tm.nome + "': " +
                   (viu_na ? std::string(viu_inf ? "NA or Inf" : "NA") : std::string("Inf")) +
                   " in the column '" + col_nivel + "' in " + std::to_string(n_aus) +
                   " row(s), the first at row " + std::to_string(primeira + 1) +
                   ": a record without a level has no known effect, and those rows would "
                   "share one level; drop those rows or fill them" +
                   (tm.mgs && col_nivel == tm.aninhado
                        ? " (an unknown maternal grandsire is 0)" : ""));
    }
  }

  // niveis: o conjunto que vem de fora da coluna (pedigree, ids da K, ou a uniao ordenada
  // num grupo iid de varios termos, ver monta_aleatorios em mme.cpp), senao os dados na
  // ordem de aparicao. Ordem de aparicao e nao alfabetica: e o que mantem a saida comparavel
  // com dado ja renumerado, onde o nivel "10" vem depois do "9" e nao entre "1" e "2".
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
    // Baia AUSENTE e ERRO. Coluna::rotulo faz de qualquer valor um rotulo: o NaN numerico
    // vira "nan", o Inf vira "inf", e o texto em branco (o que read.csv() e fread() dao a uma
    // celula vazia de coluna textual) fica "". Em todos os casos as linhas sem baia viravam
    // UMA baia, e animais sem relacao entravam como companheiros uns dos outros sem aviso.
    // Tratar cada uma como baia de um tambem seria inventar: o registro TEM companheiros, so
    // nao se sabe quais. O que conta como ausente:
    //   coluna numerica: NA ou NaN, Inf, -Inf;
    //   coluna textual (e fator, que cruza como texto): texto vazio depois de tirar espacos,
    //   tabulacoes e quebras de linha, e o texto "NaN", que factor() e as.character() fazem
    //   de um NaN numerico. O NA textual ja para antes, em tabela_do_R.
    // Os textos "NA" e "Inf" NAO entram: podem ser o nome de uma baia de verdade. A mesma
    // regra esta em sem_rotulo() (R/indirect.R), que indirect_residual(),
    // associative_matrix() e competition_strength() usam; mudar uma e mudar a outra. A
    // recusa fica aqui, onde a incidencia social e montada, para valer em todo ajustador.
    {
      const Coluna* cb = t.acha(tm.aninhado);
      if (!cb) throw Erro("no column '" + tm.aninhado + "' in the data");
      static const char* const tipos_num[] = {"NA", "Inf", "-Inf"};
      static const char* const tipos_txt[] = {"blank text", "the text 'NaN'", ""};
      const char* const* nome_tipo = cb->texto ? tipos_txt : tipos_num;
      // 0 quando a linha tem baia; senao 1 + o indice do tipo em nome_tipo
      auto tipo = [&](std::size_t i) -> int {
        if (!cb->texto) {
          const double x = cb->num[i];
          if (std::isnan(x)) return 1;
          if (std::isinf(x)) return x > 0 ? 2 : 3;
          return 0;
        }
        const std::string& s = cb->txt[i];
        const std::size_t a = s.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) return 1;
        const std::size_t b = s.find_last_not_of(" \t\r\n");
        return s.compare(a, b - a + 1, "NaN") == 0 ? 2 : 0;
      };
      std::size_t n_aus = 0, primeira = 0;
      bool viu[3] = {false, false, false};
      for (std::size_t i = 0; i < nlin; i++) {
        const int c = tipo(i);
        if (c == 0) continue;
        if (n_aus++ == 0) primeira = i;
        viu[c - 1] = true;
      }
      if (n_aus > 0) {
        std::string oque;
        for (int c = 0; c < 3; c++)
          if (viu[c]) oque += (oque.empty() ? "" : " or ") + std::string(nome_tipo[c]);
        throw Erro("indirect term '" + tm.nome + "': " + oque + " in the pen column '" +
                   tm.aninhado + "' of indirect() in " + std::to_string(n_aus) +
                   " row(s), the first at row " + std::to_string(primeira + 1) +
                   ": a record without a pen has no known pen mates, and grouping those "
                   "rows would make them mates of each other; drop those rows or assign "
                   "them a pen");
      }
    }
    const std::vector<std::string> baia = t.rotulos(tm.aninhado);
    // baia -> animais distintos nela
    std::unordered_map<std::string, std::vector<std::size_t>> membros;
    {
      std::unordered_map<std::string, std::unordered_set<std::string>> visto;
      for (std::size_t i = 0; i < nlin; i++) {
        auto it = pos.find(rot[i]);
        if (it == pos.end()) {
          throw Erro("indirect term '" + tm.nome + "': animal '" + rot[i] +
                     "' is not in the level set (pedigree)" +
                     dica_cientifica(rot[i], d.niveis));
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
      // Diluicao (Bijma 2010, Genetics 186:1029-1031): cada companheiro entra com
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

  // modelo pai / avo materno: 1 no pai e 1/2 no avo, no mesmo efeito. Pai fora do conjunto
  // de niveis exclui o registro (como qualquer classe); avo citado e fora dos niveis e ERRO,
  // como o companheiro de baia do termo social; avo desconhecido deixa so o pai.
  if (tm.mgs) {
    const std::vector<std::string> avo = t.rotulos(tm.aninhado);
    std::vector<std::uint32_t> li, cj;
    std::vector<double> v;
    for (std::size_t i = 0; i < nlin; i++) {
      auto it = pos.find(rot[i]);
      if (it == pos.end()) { d.casou[i] = 0; continue; }
      li.push_back(static_cast<std::uint32_t>(i));
      cj.push_back(static_cast<std::uint32_t>(it->second));
      v.push_back(1.0);
      const std::string& a = avo[i];
      if (a.empty() || a == "0" || a == "NA") continue;
      auto ja = pos.find(a);
      if (ja == pos.end())
        throw Erro("sire term '" + tm.nome + "': maternal grandsire '" + a +
                   "' is not in the level set (pedigree)" + dica_cientifica(a, d.niveis));
      li.push_back(static_cast<std::uint32_t>(i));
      cj.push_back(static_cast<std::uint32_t>(ja->second));
      v.push_back(0.5);
    }
    d.z = de_triplos(nlin, d.n_niveis, li, cj, v);
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
