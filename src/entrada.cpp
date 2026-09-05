// A travessia para o R. Converte objetos R em estruturas do pacote e de volta.
//
// Nenhuma conta acontece aqui. Se um dia aparecer aritmetica neste arquivo, ela estara fora
// dos testes que validam a numerica, e sera por ai que os numeros vao divergir.
//
// Uma excecao de C++ atravessando a fronteira do R deixa o interpretador em estado
// indefinido, entao TODO ponto de entrada embrulha o corpo e converte a excecao em Rf_error.

#include "mme.h"
#include <cstdio>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


// Le triplos do R (base 1, triangulo qualquer) e devolve a Csc do triangulo INFERIOR.
static br::Csc csc_do_R(SEXP i, SEXP j, SEXP x, SEXP n) {
  const int nn = Rf_asInteger(n);
  if (nn <= 0) Rf_error("invalid dimension");
  const R_xlen_t nz = XLENGTH(x);
  if (XLENGTH(i) != nz || XLENGTH(j) != nz) Rf_error("i, j and x with different lengths");
  std::vector<std::uint32_t> li, cj;
  std::vector<double> v;
  li.reserve(static_cast<std::size_t>(nz));
  cj.reserve(static_cast<std::size_t>(nz));
  v.reserve(static_cast<std::size_t>(nz));
  for (R_xlen_t k = 0; k < nz; k++) {
    int a = INTEGER(i)[k] - 1, b = INTEGER(j)[k] - 1;
    if (a < 0 || b < 0 || a >= nn || b >= nn) Rf_error("triplet outside the matrix");
    if (a < b) std::swap(a, b);
    li.push_back(static_cast<std::uint32_t>(a));
    cj.push_back(static_cast<std::uint32_t>(b));
    v.push_back(REAL(x)[k]);
  }
  return br::de_triplos(static_cast<std::size_t>(nn), static_cast<std::size_t>(nn), li, cj, v);
}

static SEXP csc_para_R(const br::Csc& m) {
  const R_xlen_t nz = static_cast<R_xlen_t>(m.nnz());
  SEXP li = PROTECT(Rf_allocVector(INTSXP, nz));
  SEXP cj = PROTECT(Rf_allocVector(INTSXP, nz));
  SEXP vv = PROTECT(Rf_allocVector(REALSXP, nz));
  R_xlen_t k = 0;
  for (std::size_t c = 0; c < m.ncol; c++)
    for (std::size_t t = m.colptr[c]; t < m.colptr[c + 1]; t++, k++) {
      INTEGER(li)[k] = static_cast<int>(m.linha[t]) + 1;
      INTEGER(cj)[k] = static_cast<int>(c) + 1;
      REAL(vv)[k] = m.valor[t];
    }
  const char* campos[] = {"i", "j", "x", "n"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
  for (int q = 0; q < 4; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
  SET_VECTOR_ELT(out, 0, li);
  SET_VECTOR_ELT(out, 1, cj);
  SET_VECTOR_ELT(out, 2, vv);
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(static_cast<int>(m.ncol)));
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(5);
  return out;
}

// Converte a excecao do dominio em condicao do R.
//
// Variadica de proposito: chaves NAO protegem virgulas no preprocessador, so parenteses. Com
// um parametro so, qualquer virgula dentro do corpo — e ha uma em cada Rf_allocVector — vira
// separador de argumento e a macro deixa de existir. Foi exatamente esse o erro na primeira
// versao deste arquivo.
//
// Rf_error faz longjmp, entao nada que precise de destrutor pode estar vivo quando ele
// dispara: a mensagem e copiada para um buffer ANTES de sair do bloco try.
#define GUARDA(...)                                                              char msg__[512];                                                               msg__[0] = 0;                                                                  try { __VA_ARGS__ }                                                            catch (const std::exception& e) {                                                std::snprintf(msg__, sizeof msg__, "%s", e.what());                          }                                                                              catch (...) { std::snprintf(msg__, sizeof msg__, "unknown error"); }       Rf_error("%s", msg__);

static std::vector<std::string> textos(SEXP v, const char* quem) {
  if (TYPEOF(v) != STRSXP) Rf_error("%s: expected a text vector", quem);
  const R_xlen_t n = XLENGTH(v);
  std::vector<std::string> out;
  out.reserve(static_cast<std::size_t>(n));
  for (R_xlen_t i = 0; i < n; i++) {
    SEXP e = STRING_ELT(v, i);
    if (e == NA_STRING)
      Rf_error("%s: there is NA. A missing identifier would become an animal named 'NA'", quem);
    out.emplace_back(CHAR(e));
  }
  return out;
}

extern "C" {

// Devolve o pedigree em ORDEM TOPOLOGICA, com a endogamia de cada animal.
//
// A ordem importa e por isso ela sai: o A^-1 e os efeitos genéticos vêm nessa ordem, e
// juntar de volta pela posicao original seria trocar animal em silencio.
SEXP R_pedigree(SEXP id, SEXP pai, SEXP mae, SEXP mfx, SEXP gmx) {
  GUARDA(
    auto i = textos(id, "id");
    auto p = textos(pai, "sire");
    auto m = textos(mae, "dam");
    br::Pedigree ped = br::constroi_pedigree(i, p, m, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)));
    std::vector<double> f = br::endogamia(ped);

    const R_xlen_t n = static_cast<R_xlen_t>(ped.ids.size());
    SEXP ids = PROTECT(Rf_allocVector(STRSXP, n));
    SEXP vp  = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP vm  = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP vf  = PROTECT(Rf_allocVector(REALSXP, n));
    for (R_xlen_t k = 0; k < n; k++) {
      SET_STRING_ELT(ids, k, Rf_mkChar(ped.ids[k].c_str()));
      // base 1 para o R; NA para progenitor desconhecido, que e o que o R espera ver
      INTEGER(vp)[k] = ped.pai[k] >= 0 ? static_cast<int>(ped.pai[k] + 1) : NA_INTEGER;
      INTEGER(vm)[k] = ped.mae[k] >= 0 ? static_cast<int>(ped.mae[k] + 1) : NA_INTEGER;
      REAL(vf)[k] = f[k];
    }
    const char* campos[] = {"id", "sire", "dam", "F"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
    for (int j = 0; j < 4; j++) SET_STRING_ELT(nms, j, Rf_mkChar(campos[j]));
    SET_VECTOR_ELT(out, 0, ids);
    SET_VECTOR_ELT(out, 1, vp);
    SET_VECTOR_ELT(out, 2, vm);
    SET_VECTOR_ELT(out, 3, vf);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(6);
    return out;
  )
}

// A^-1 no formato de triplos do triangulo inferior, pronta para virar dgCMatrix no R.
SEXP R_a_inversa(SEXP id, SEXP pai, SEXP mae, SEXP mfx, SEXP gmx) {
  GUARDA(
    auto i = textos(id, "id");
    auto p = textos(pai, "sire");
    auto m = textos(mae, "dam");
    br::Pedigree ped = br::constroi_pedigree(i, p, m, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)));
    std::vector<double> f = br::endogamia(ped);
    br::Csc a = br::a_inversa(ped, f);

    const R_xlen_t nz = static_cast<R_xlen_t>(a.nnz());
    SEXP li = PROTECT(Rf_allocVector(INTSXP, nz));
    SEXP cj = PROTECT(Rf_allocVector(INTSXP, nz));
    SEXP vv = PROTECT(Rf_allocVector(REALSXP, nz));
    R_xlen_t k = 0;
    for (std::size_t c = 0; c < a.ncol; c++)
      for (std::size_t t = a.colptr[c]; t < a.colptr[c + 1]; t++, k++) {
        INTEGER(li)[k] = static_cast<int>(a.linha[t]) + 1;
        INTEGER(cj)[k] = static_cast<int>(c) + 1;
        REAL(vv)[k] = a.valor[t];
      }
    SEXP ids = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(ped.ids.size())));
    for (std::size_t t = 0; t < ped.ids.size(); t++)
      SET_STRING_ELT(ids, static_cast<R_xlen_t>(t), Rf_mkChar(ped.ids[t].c_str()));

    const char* campos[] = {"i", "j", "x", "n", "id"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 5));
    for (int q = 0; q < 5; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, li);
    SET_VECTOR_ELT(out, 1, cj);
    SET_VECTOR_ELT(out, 2, vv);
    SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(static_cast<int>(a.ncol)));
    SET_VECTOR_ELT(out, 4, ids);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(6);
    return out;
  )
}

// Inversa de uma simetrica positiva-definida, pela Cholesky em bloco.
SEXP R_inv_pd(SEXP m) {
  GUARDA(
    SEXP dim = Rf_getAttrib(m, R_DimSymbol);
    if (TYPEOF(m) != REALSXP || dim == R_NilValue || XLENGTH(dim) != 2)
      Rf_error("expected a numeric matrix");
    const int nl = INTEGER(dim)[0], nc = INTEGER(dim)[1];
    if (nl != nc) Rf_error("the matrix must be square");
    br::Densa s(static_cast<std::size_t>(nl), static_cast<std::size_t>(nc));
    // R guarda por COLUNA e a Densa por LINHA
    for (int j = 0; j < nc; j++)
      for (int i = 0; i < nl; i++)
        s.at(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) =
            REAL(m)[static_cast<R_xlen_t>(j) * nl + i];
    br::Densa r = br::inv_pd(s);
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, nl, nc));
    for (int j = 0; j < nc; j++)
      for (int i = 0; i < nl; i++)
        REAL(out)[static_cast<R_xlen_t>(j) * nl + i] =
            r.at(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
    UNPROTECT(1);
    return out;
  )
}

// Fatoracao esparsa completa: ordena por grau minimo, permuta, fatora, e devolve L com a
// permutacao e o log-determinante.
SEXP R_chol_esparsa(SEXP i, SEXP j, SEXP x, SEXP n, SEXP ordenar) {
  GUARDA(
    br::Csc a = csc_do_R(i, j, x, n);
    std::vector<std::size_t> perm(a.ncol);
    if (Rf_asLogical(ordenar) == TRUE) perm = br::grau_minimo(a);
    else for (std::size_t k = 0; k < a.ncol; k++) perm[k] = k;
    br::Csc pa = br::permuta_sim(a, perm);
    br::Simbolica sb = br::simbolica(pa);
    br::Csc L;
    if (!br::cholesky(pa, sb, L))
      throw br::Erro("the matrix is not positive-definite");

    SEXP vperm = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(perm.size())));
    for (std::size_t k = 0; k < perm.size(); k++)
      INTEGER(vperm)[k] = static_cast<int>(perm[k]) + 1;
    SEXP fac = PROTECT(csc_para_R(L));
    const char* campos[] = {"L", "perm", "logdet", "dense_block"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
    for (int q = 0; q < 4; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, fac);
    SET_VECTOR_ELT(out, 1, vperm);
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(br::logdet(L)));
    SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(static_cast<int>(br::bloco_denso_final(L))));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(4);
    return out;
  )
}

// Inversa seletiva. `bloco` de 0 detecta a cauda densa; 1 forca a recorrencia pura.
SEXP R_inv_seletiva(SEXP i, SEXP j, SEXP x, SEXP n, SEXP bloco) {
  GUARDA(
    br::Csc a = csc_do_R(i, j, x, n);
    std::vector<std::size_t> perm = br::grau_minimo(a);
    br::Csc pa = br::permuta_sim(a, perm);
    br::Simbolica sb = br::simbolica(pa);
    br::Csc L;
    if (!br::cholesky(pa, sb, L)) throw br::Erro("the matrix is not positive-definite");
    br::SelInv z = br::inversa_seletiva(L, static_cast<std::size_t>(Rf_asInteger(bloco)));

    // devolve na numeracao ORIGINAL, senao quem le PEV le a diagonal trocada
    const R_xlen_t nz = static_cast<R_xlen_t>(z.valor.size());
    SEXP li = PROTECT(Rf_allocVector(INTSXP, nz));
    SEXP cj = PROTECT(Rf_allocVector(INTSXP, nz));
    SEXP vv = PROTECT(Rf_allocVector(REALSXP, nz));
    R_xlen_t k = 0;
    for (std::size_t c = 0; c < z.n; c++)
      for (std::size_t t = z.colptr[c]; t < z.colptr[c + 1]; t++, k++) {
        INTEGER(li)[k] = static_cast<int>(perm[z.linha[t]]) + 1;
        INTEGER(cj)[k] = static_cast<int>(perm[c]) + 1;
        REAL(vv)[k] = z.valor[t];
      }
    const char* campos[] = {"i", "j", "x", "n"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
    for (int q = 0; q < 4; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, li);
    SET_VECTOR_ELT(out, 1, cj);
    SET_VECTOR_ELT(out, 2, vv);
    SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(static_cast<int>(z.n)));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(5);
    return out;
  )
}

SEXP R_resolve(SEXP i, SEXP j, SEXP x, SEXP n, SEXP b) {
  GUARDA(
    br::Csc a = csc_do_R(i, j, x, n);
    std::vector<std::size_t> perm = br::grau_minimo(a);
    br::Csc pa = br::permuta_sim(a, perm);
    br::Simbolica sb = br::simbolica(pa);
    br::Csc L;
    if (!br::cholesky(pa, sb, L)) throw br::Erro("the matrix is not positive-definite");
    const std::size_t nn = a.ncol;
    if (static_cast<std::size_t>(XLENGTH(b)) != nn) Rf_error("b with the wrong length");
    std::vector<double> pb(nn);
    for (std::size_t k = 0; k < nn; k++) pb[k] = REAL(b)[perm[k]];
    std::vector<double> px = br::resolve(L, pb);
    SEXP out = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(nn)));
    for (std::size_t k = 0; k < nn; k++) REAL(out)[perm[k]] = px[k];
    UNPROTECT(1);
    return out;
  )
}


// ---------------------------------------------------------------- o ajuste completo

// Monta Modelo + Tabela + Pedigree a partir dos objetos do R.
//
// `tdil` e a diluicao por termo, paralela a tsoc; o padrao R_NilValue existe porque so os
// pontos de entrada que ja transportam dilution= (R_ajustar, R_avaliar) a enviam. Os
// demais ajustadores recusam dilution > 0 no lado R (recusa_dilution em R/model.R) e
// continuam chamando com a aridade antiga.
static br::Modelo modelo_do_R(SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov, SEXP test,
                              SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP ausente,
                              SEXP usa_ausente, SEXP tdil = R_NilValue,
                              SEXP tfix = R_NilValue) {
  const R_xlen_t nt = XLENGTH(tnome);
  std::vector<br::Termo> termos;
  std::vector<std::pair<std::string, std::vector<std::string>>> grupos;
  for (R_xlen_t k = 0; k < nt; k++) {
    br::Termo t;
    t.nome = CHAR(STRING_ELT(tnome, k));
    t.coluna = CHAR(STRING_ELT(tcol, k));
    t.efeito = LOGICAL(tcov)[k] ? br::Efeito::Covariavel : br::Efeito::Classe;
    const int e = INTEGER(test)[k];
    t.estrutura = e == 3 ? br::Estrutura::Declarada
                : e == 2 ? br::Estrutura::Parentesco
                : e == 1 ? br::Estrutura::Diagonal : br::Estrutura::Fixo;
    const char* nest = CHAR(STRING_ELT(tnest, k));
    if (*nest) t.aninhado = nest;
    t.social = LOGICAL(tsoc)[k] != 0;
    if (tdil != R_NilValue) {
      if (TYPEOF(tdil) != REALSXP || XLENGTH(tdil) != nt)
        Rf_error("dilution: expected one numeric value per term");
      t.diluicao = REAL(tdil)[k];
    }
    // base: colunas separadas por virgula. E o que faz um termo virar regressao aleatoria —
    // um termo com m coeficientes e uma covariancia m x m, nao m termos independentes.
    const char* base = CHAR(STRING_ELT(tbase, k));
    if (*base) {
      std::string b = base;
      std::size_t ini = 0;
      while (ini <= b.size()) {
        std::size_t fim = b.find(',', ini);
        if (fim == std::string::npos) fim = b.size();
        std::string col = b.substr(ini, fim - ini);
        if (!col.empty()) t.base.push_back(col);
        ini = fim + 1;
      }
    }
    const char* g = CHAR(STRING_ELT(tgrp, k));
    if (*g && t.estrutura != br::Estrutura::Fixo) {
      bool achou = false;
      for (auto& gn : grupos)
        if (gn.first == g) { gn.second.push_back(t.nome); achou = true; break; }
      if (!achou) grupos.push_back({g, {t.nome}});
    }
    termos.push_back(std::move(t));
  }
  br::Modelo m = br::monta_modelo(CHAR(STRING_ELT(alvo, 0)), std::move(termos), grupos,
                                  Rf_asLogical(usa_ausente) == TRUE, Rf_asReal(ausente));
  // kernel(..., fixed = v): PRENDE a variancia daquele termo em v. tfix chega por TERMO e
  // theta_fixo vive por POSICAO EM THETA, entao o mapeamento e feito aqui, onde os grupos
  // ja existem. Um kernel forma grupo escalar proprio, entao e a diagonal (0,0) dele.
  if (tfix != R_NilValue) {
    if (TYPEOF(tfix) != REALSXP || XLENGTH(tfix) != nt)
      Rf_error("fixed: expected one numeric value per term");
    for (R_xlen_t k = 0; k < nt; k++) {
      if (!R_finite(REAL(tfix)[k])) continue;
      const std::string nome = CHAR(STRING_ELT(tnome, k));
      for (const br::Grupo& gr : m.grupos) {
        bool meu = false;
        for (std::size_t t : gr.termos) if (m.termos[t].nome == nome) meu = true;
        if (!meu) continue;
        if (gr.dim != 1)
          Rf_error("kernel(fixed=) holds a SCALAR component: term '%s' shares a "
                   "covariance group, where holding one entry of the matrix and "
                   "estimating the rest is a different problem", nome.c_str());
        if (m.theta_fixo.empty())
          m.theta_fixo.assign(m.ntheta, std::numeric_limits<double>::quiet_NaN());
        m.theta_fixo[gr.offset] = REAL(tfix)[k];
      }
    }
  }
  return m;
}

static br::Tabela tabela_do_R(SEXP dados, SEXP nomes) {
  br::Tabela t;
  const R_xlen_t nc = XLENGTH(dados);
  for (R_xlen_t k = 0; k < nc; k++) {
    SEXP col = VECTOR_ELT(dados, k);
    br::Coluna c;
    if (TYPEOF(col) == STRSXP) {
      c.texto = true;
      const R_xlen_t n = XLENGTH(col);
      c.txt.reserve(static_cast<std::size_t>(n));
      for (R_xlen_t i = 0; i < n; i++) {
        SEXP e = STRING_ELT(col, i);
        if (e == NA_STRING) Rf_error("column with text NA at position %d", (int) k + 1);
        c.txt.emplace_back(CHAR(e));
      }
    } else if (TYPEOF(col) == REALSXP) {
      c.num.assign(REAL(col), REAL(col) + XLENGTH(col));
    } else {
      Rf_error("column of unsupported type at position %d", (int) k + 1);
    }
    t.nomes.emplace_back(CHAR(STRING_ELT(nomes, k)));
    if (t.nlin == 0) t.nlin = c.tamanho();
    else if (t.nlin != c.tamanho())
      Rf_error("column with different length at position %d", (int) k + 1);
    t.colunas.push_back(std::move(c));
  }
  return t;
}

// Kernels declarados (kernel(id, K=)): lista PARALELA aos termos, NULL para termo comum e
// list(ids, K) para termo kernel. A matriz cruza como Densa (o R guarda por COLUNA e a
// Densa por LINHA), e a validacao de forma acontece aqui, antes de qualquer conta.
static std::vector<br::KernelDecl> kernels_do_R(SEXP kern) {
  std::vector<br::KernelDecl> out;
  if (Rf_isNull(kern)) return out;
  if (TYPEOF(kern) != VECSXP) Rf_error("kernels: expected a list parallel to the terms");
  const R_xlen_t nt = XLENGTH(kern);
  out.resize(static_cast<std::size_t>(nt));
  for (R_xlen_t k = 0; k < nt; k++) {
    SEXP e = VECTOR_ELT(kern, k);
    if (Rf_isNull(e)) continue;
    if (TYPEOF(e) != VECSXP || XLENGTH(e) != 2)
      Rf_error("kernel of term %d: expected list(ids, K)", (int) k + 1);
    auto ids = textos(VECTOR_ELT(e, 0), "kernel ids");
    SEXP km = VECTOR_ELT(e, 1);
    SEXP dim = Rf_getAttrib(km, R_DimSymbol);
    if (TYPEOF(km) != REALSXP || dim == R_NilValue || XLENGTH(dim) != 2)
      Rf_error("kernel of term %d: K must be a numeric matrix", (int) k + 1);
    const int nl = INTEGER(dim)[0], nc = INTEGER(dim)[1];
    if (nl != nc || (std::size_t) nl != ids.size())
      Rf_error("kernel of term %d: K of %d x %d for %d ids", (int) k + 1, nl, nc,
               (int) ids.size());
    br::Densa m((std::size_t) nl, (std::size_t) nc);
    for (int j2 = 0; j2 < nc; j2++)
      for (int i2 = 0; i2 < nl; i2++)
        m.at((std::size_t) i2, (std::size_t) j2) = REAL(km)[(R_xlen_t) j2 * nl + i2];
    out[(std::size_t) k].ids = std::move(ids);
    out[(std::size_t) k].k = std::move(m);
  }
  return out;
}

extern "C++" {
template <class DES>
std::string genomica_no_desenho(DES& d, const br::Pedigree* pp, br::Pedigree& ped,
                                SEXP gid, SEXP gm, SEXP mistura, SEXP anucleo,
                                SEXP vk, std::vector<double>* dh = nullptr,
                                std::vector<std::size_t>* dr = nullptr) {
  if (XLENGTH(gid) == 0) return "";
  if (!pp) Rf_error("genotypes without a pedigree: H^-1 needs A^-1");
  auto gids = textos(gid, "genotypes");
  SEXP dim = Rf_getAttrib(gm, R_DimSymbol);
  if (TYPEOF(gm) != REALSXP || dim == R_NilValue || XLENGTH(dim) != 2)
    Rf_error("the genotype matrix must be numeric");
  const int nl = INTEGER(dim)[0], nm2 = INTEGER(dim)[1];
  if ((std::size_t) nl != gids.size())
    Rf_error("%d identifiers for %d genotype rows", (int) gids.size(), nl);
  // R guarda por COLUNA e a Densa por LINHA: transposicao aqui, uma vez
  br::Densa mg((std::size_t) nl, (std::size_t) nm2);
  for (int j2 = 0; j2 < nm2; j2++)
    for (int i2 = 0; i2 < nl; i2++)
      mg.at((std::size_t) i2, (std::size_t) j2) = REAL(gm)[(R_xlen_t) j2 * nl + i2];
  std::vector<std::string> nuc;
  if (XLENGTH(anucleo) > 0) nuc = textos(anucleo, "APY core");
  const std::size_t kv = (std::size_t) std::max(0, Rf_asInteger(vk));
  br::RelatorioG rel = br::aplica_genomica(d, ped, gids, mg, Rf_asReal(mistura), nuc, kv);
  // a priori dos genotipados, para accuracy(): so quem pede recebe
  if (dh) *dh = rel.diag_gstar;
  if (dr) *dr = rel.linha_ped;
  std::string nota = "single-step: " + std::to_string(gids.size()) + " genotyped x " +
                     std::to_string(nm2) + " markers, " +
                     std::to_string(rel.n_imputados) + " missing value(s) imputed by the " +
                     "mean, " + std::to_string(rel.n_monomorficos) +
                     " monomorphic marker(s) left out of G";
  if (!nuc.empty())
    nota += "; G* inverted by APY with a core of " + std::to_string(nuc.size()) +
            " (an approximation; the core size is part of the result)";
  if (kv > 0)
    nota += "; G* inverted by the Vecchia recursion with k = " + std::to_string(kv) +
            " (an approximation; k is part of the result)";
  return nota;
}
}  // extern "C++"

// Avalia -2logL pelas DUAS vias, mais score, AI e EM num theta dado. Existe para os gates:
// a identidade MME <-> forma V, e o score contra diferencas finitas centrais.
SEXP R_avaliar(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
               SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid, SEXP ppai,
               SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP theta, SEXP com_densa, SEXP mfx, SEXP gmx, SEXP pesos, SEXP gid, SEXP gm, SEXP mistura, SEXP anucleo, SEXP vk, SEXP kern, SEXP tdil) {
  GUARDA(
    br::Modelo m = modelo_do_R(alvo, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc, ausente, usa_ausente, tdil);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)));
      pp = &ped;
    }
    std::vector<double> pw;
    if (XLENGTH(pesos) > 0) pw.assign(REAL(pesos), REAL(pesos) + XLENGTH(pesos));
    std::vector<br::KernelDecl> ks = kernels_do_R(kern);
    br::Desenho d = br::monta_desenho(m, t, pp, pw.empty() ? nullptr : &pw,
                                      ks.empty() ? nullptr : &ks);
    std::vector<double> th(REAL(theta), REAL(theta) + XLENGTH(theta));
    if (th.size() != m.ntheta) Rf_error("theta with %d entries; the layout asks for %d",
                                        (int) th.size(), (int) m.ntheta);

    genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk);
    br::Avaliacao a = br::avalia(d, th);
    if (!a.ok) Rf_error("theta INADMISSIBLE: some covariance is not positive-definite");
    const double dv = Rf_asLogical(com_densa) == TRUE ? br::neg2logl_densa_V(d, th)
                                                      : std::nan("");

    const std::size_t nt2 = m.ntheta;
    SEXP score = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP em    = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP ai    = PROTECT(Rf_allocMatrix(REALSXP, (int) nt2, (int) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(score)[k] = a.score[k];
      REAL(em)[k] = a.em_theta[k];
      for (std::size_t j = 0; j < nt2; j++)
        REAL(ai)[j * nt2 + k] = a.ai.at(k, j);
    }
    const char* campos[] = {"neg2logl", "neg2logl_V", "score", "em", "ai", "off_pattern"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 6));
    for (int q = 0; q < 6; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(a.neg2logl));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(dv));
    SET_VECTOR_ELT(out, 2, score);
    SET_VECTOR_ELT(out, 3, em);
    SET_VECTOR_ELT(out, 4, ai);
    SET_VECTOR_ELT(out, 5, Rf_ScalarInteger((int) a.fora_do_padrao));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(5);
    return out;
  )
}

// Conversao + aplicacao da genomica, comum aos tres ajustadores (o desenho e template
// porque uni, multi e AR carregam os mesmos campos que o passo unico toca). Devolve a
// nota do relatorio, vazia sem genotipos. O extern "C++" existe porque este arquivo vive
// num bloco extern "C" — e template nao tem linkage de C.

SEXP R_ajustar(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
               SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid, SEXP ppai,
               SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP maxiter, SEXP tol, SEXP n_em,
               SEXP gid, SEXP gm, SEXP mistura, SEXP anucleo, SEXP vk, SEXP verb, SEXP mfx, SEXP gmx, SEXP pesos, SEXP inicio, SEXP kern, SEXP tdil, SEXP tfix) {
  GUARDA(
    br::Modelo m = modelo_do_R(alvo, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc, ausente, usa_ausente, tdil, tfix);
    br::Tabela t = tabela_do_R(dados, nomes);
    // a priori dos genotipados sob passo unico, para accuracy(). Vazio sem genomica, e
    // nesse caso accuracy() segue com o 1 + F do pedigree.
    std::vector<double> diag_h_geno;
    std::vector<std::size_t> linha_h_geno;
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)));
      pp = &ped;
    }
    std::vector<double> pw;
    if (XLENGTH(pesos) > 0) pw.assign(REAL(pesos), REAL(pesos) + XLENGTH(pesos));
    std::vector<br::KernelDecl> ks = kernels_do_R(kern);
    br::Desenho d = br::monta_desenho(m, t, pp, pw.empty() ? nullptr : &pw,
                                      ks.empty() ? nullptr : &ks);

    std::string nota = genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk,                                           &diag_h_geno, &linha_h_geno);

    std::vector<double> th0;
    if (XLENGTH(inicio) > 0) {
      th0.assign(REAL(inicio), REAL(inicio) + XLENGTH(inicio));
      if (th0.size() != m.ntheta)
        Rf_error("start with %d entries; the layout asks for %d", (int) th0.size(), (int) m.ntheta);
    }
    br::Ajuste r = br::ajusta(d, th0.empty() ? nullptr : &th0, (std::size_t) Rf_asInteger(n_em),
                              (std::size_t) Rf_asInteger(maxiter), Rf_asReal(tol),
                              Rf_asLogical(verb) == TRUE);
    if (!nota.empty()) r.mensagem = r.mensagem.empty() ? nota : nota + "; " + r.mensagem;

    std::vector<std::string> nomes_th = m.nomes_theta();
    const std::size_t nt2 = m.ntheta;
    SEXP theta = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP se    = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP nms_t = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(theta)[k] = r.theta.empty() ? NA_REAL : r.theta[k];
      REAL(se)[k] = r.se.empty() ? NA_REAL : r.se[k];
      SET_STRING_ELT(nms_t, (R_xlen_t) k, Rf_mkChar(nomes_th[k].c_str()));
    }
    Rf_setAttrib(theta, R_NamesSymbol, nms_t);
    Rf_setAttrib(se, R_NamesSymbol, nms_t);

    // EBV e PEV por grupo: fatias da solucao e da inversa seletiva, com os niveis como nomes
    //
    // Um ajuste que parou ANTES da primeira avaliacao volta sem theta: montar as MME com
    // esse vetor vazio lia fora dele. Sem theta nao ha solucao para fatiar, e o objeto ja
    // carrega a mensagem — as fatias saem NA, como no resto do caminho de falha.
    br::Montado M;
    if (!r.theta.empty()) M = br::monta_mme(d, r.theta);
    if (M.offset_grupo.empty()) M.offset_grupo.assign(m.grupos.size(), 0);
    SEXP ebv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) m.grupos.size()));
    SEXP pev = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) m.grupos.size()));
    SEXP ebv_nomes = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) m.grupos.size()));
    for (std::size_t g = 0; g < m.grupos.size(); g++) {
      const std::size_t larg = d.largura(g);
      SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
      SEXP pv = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
      for (std::size_t c = 0; c < larg; c++) {
        REAL(v)[c] = r.solucao.empty() ? NA_REAL : r.solucao[M.offset_grupo[g] + c];
        REAL(pv)[c] = r.pev.empty() ? NA_REAL : r.pev[M.offset_grupo[g] + c];
      }
      for (const br::DesenhoTermo& a : d.aleatorios)
        if (a.termo == m.grupos[g].termos[0]) {
          SEXP nv = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) larg));
          for (std::size_t c = 0; c < larg; c++)
            SET_STRING_ELT(nv, (R_xlen_t) c, Rf_mkChar(a.niveis[c % a.n_niveis].c_str()));
          Rf_setAttrib(v, R_NamesSymbol, nv);
          Rf_setAttrib(pv, R_NamesSymbol, nv);
          UNPROTECT(1);
          break;
        }
      SET_VECTOR_ELT(ebv, (R_xlen_t) g, v);
      SET_VECTOR_ELT(pev, (R_xlen_t) g, pv);
      SET_STRING_ELT(ebv_nomes, (R_xlen_t) g, Rf_mkChar(m.grupos[g].nome.c_str()));
      UNPROTECT(2);
    }
    Rf_setAttrib(ebv, R_NamesSymbol, ebv_nomes);
    Rf_setAttrib(pev, R_NamesSymbol, Rf_duplicate(ebv_nomes));

    SEXP sc = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP vc = PROTECT(Rf_allocMatrix(REALSXP, (int) nt2, (int) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(sc)[k] = r.score.empty() ? NA_REAL : r.score[k];
      for (std::size_t j = 0; j < nt2; j++)
        REAL(vc)[j * nt2 + k] = r.vcov.empty() ? NA_REAL : r.vcov[k * nt2 + j];
    }
    Rf_setAttrib(sc, R_NamesSymbol, Rf_duplicate(nms_t));

    // Solucoes dos efeitos fixos: as primeiras x.ncol posicoes da solucao, nomeadas pelas
    // colunas de X — os MESMOS nomes termo=nivel que dropped_x usa. Nada e recalculado: o
    // solver ja resolveu o sistema inteiro, isto e uma fatia. Sem theta nao ha solucao e a
    // fatia sai NA, como no resto do caminho de falha.
    SEXP bfix = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) d.x.ncol));
    SEXP bnms = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) d.x.ncol));
    for (std::size_t k = 0; k < d.x.ncol; k++) {
      REAL(bfix)[k] = r.solucao.empty() ? NA_REAL : r.solucao[k];
      SET_STRING_ELT(bnms, (R_xlen_t) k, Rf_mkChar(d.nomes_x[k].c_str()));
    }
    Rf_setAttrib(bfix, R_NamesSymbol, bnms);

    const char* campos[] = {"theta", "se", "neg2logl", "converged", "iters", "reldelta",
                            "message", "n_used", "n_columns", "ebv", "dropped_x", "pev",
                            "score", "vcov", "b", "newton_dec",
                            "h_prior", "h_prior_row"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 18));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 18));
    for (int q = 0; q < 18; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SEXP saiu = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_x.size()));
    for (std::size_t k = 0; k < d.saiu_x.size(); k++)
      SET_STRING_ELT(saiu, (R_xlen_t) k, Rf_mkChar(d.saiu_x[k].c_str()));
    SET_VECTOR_ELT(out, 0, theta);
    SET_VECTOR_ELT(out, 1, se);
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(r.neg2logl));
    SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(r.convergiu));
    SET_VECTOR_ELT(out, 4, Rf_ScalarInteger((int) r.iters));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(r.reldelta));
    SET_VECTOR_ELT(out, 6, Rf_mkString(r.mensagem.c_str()));
    SET_VECTOR_ELT(out, 7, Rf_ScalarInteger((int) d.n_usadas()));
    SET_VECTOR_ELT(out, 8, Rf_ScalarInteger((int) d.total_colunas()));
    SET_VECTOR_ELT(out, 9, ebv);
    SET_VECTOR_ELT(out, 10, saiu);
    SET_VECTOR_ELT(out, 11, pev);
    SET_VECTOR_ELT(out, 12, sc);
    SET_VECTOR_ELT(out, 13, vc);
    SET_VECTOR_ELT(out, 14, bfix);
    SET_VECTOR_ELT(out, 15, Rf_ScalarReal(r.decremento));
    // Sob passo unico a variancia a priori de um genotipado e a diagonal de H, que
    // naquele bloco e a de G*, e nao 1 + F do pedigree. Vao as duas coisas: o valor e a
    // LINHA do pedigree a que ele pertence, porque accuracy() indexa por linha. Vazio
    // quando nao houve genomica, e ai accuracy() segue com o F do pedigree.
    SEXP dh = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) diag_h_geno.size()));
    SEXP dr = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) linha_h_geno.size()));
    for (std::size_t q = 0; q < diag_h_geno.size(); q++) REAL(dh)[q] = diag_h_geno[q];
    for (std::size_t q = 0; q < linha_h_geno.size(); q++)
      INTEGER(dr)[q] = (int) linha_h_geno[q] + 1;
    SET_VECTOR_ELT(out, 16, dh);
    SET_VECTOR_ELT(out, 17, dr);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(15);
    return out;
  )
}

// A22^-1 pelo Schur esparso, exposta para os gates. `geno` em base 1 na ordem do pedigree
// TOPOLOGICO (a mesma de pedigree()).
SEXP R_a22_inversa(SEXP pid, SEXP ppai, SEXP pmae, SEXP geno) {
  GUARDA(
    auto i = textos(pid, "id");
    auto p = textos(ppai, "sire");
    auto ma = textos(pmae, "dam");
    br::Pedigree ped = br::constroi_pedigree(i, p, ma);
    std::vector<double> f = br::endogamia(ped);
    br::Csc ainv = br::a_inversa(ped, f);
    std::vector<std::size_t> idx;
    for (R_xlen_t k = 0; k < XLENGTH(geno); k++) {
      const int v = INTEGER(geno)[k];
      if (v < 1 || (std::size_t) v > ped.ids.size()) Rf_error("genotyped index outside the pedigree");
      idx.push_back((std::size_t) v - 1);
    }
    br::Densa a22i = br::a22_inversa(ainv, idx);
    const std::size_t n2 = idx.size();
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, (int) n2, (int) n2));
    for (std::size_t j2 = 0; j2 < n2; j2++)
      for (std::size_t i2 = 0; i2 < n2; i2++)
        REAL(out)[j2 * n2 + i2] = a22i.at(i2, j2);
    UNPROTECT(1);
    return out;
  )
}


// ---------------------------------------------------------------- multi-caracteristica

SEXP R_avaliar_mt(SEXP dados, SEXP nomes, SEXP alvos, SEXP tnome, SEXP tcol, SEXP tcov,
                  SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                  SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP theta,
                  SEXP com_densa, SEXP mfx, SEXP gmx, SEXP kern) {
  GUARDA(
    SEXP alvo1 = PROTECT(Rf_mkString(CHAR(STRING_ELT(alvos, 0))));
    br::Modelo m = modelo_do_R(alvo1, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc,
                               ausente, usa_ausente);
    UNPROTECT(1);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p2 = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p2, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)));
      pp = &ped;
    }
    std::vector<std::string> alv = textos(alvos, "traits");
    std::vector<br::KernelDecl> kd = kernels_do_R(kern);
    br::DesenhoMT d = br::monta_desenho_mt(m, alv, t, pp, &kd);
    std::vector<double> th(REAL(theta), REAL(theta) + XLENGTH(theta));
    if (th.size() != d.modelo.ntheta)
      Rf_error("theta with %d entries; the layout asks for %d", (int) th.size(),
               (int) d.modelo.ntheta);

    br::AvaliacaoMT a = br::avalia_mt(d, th);
    if (!a.ok) Rf_error("theta INADMISSIBLE");
    const double dv = Rf_asLogical(com_densa) == TRUE ? br::neg2logl_densa_V_mt(d, th)
                                                      : std::nan("");
    const std::size_t nt2 = d.modelo.ntheta;
    SEXP score = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP ai    = PROTECT(Rf_allocMatrix(REALSXP, (int) nt2, (int) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(score)[k] = a.score[k];
      for (std::size_t j = 0; j < nt2; j++) REAL(ai)[j * nt2 + k] = a.ai.at(k, j);
    }
    const char* campos[] = {"neg2logl", "neg2logl_V", "score", "ai", "off_pattern"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 5));
    for (int q = 0; q < 5; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(a.neg2logl));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(dv));
    SET_VECTOR_ELT(out, 2, score);
    SET_VECTOR_ELT(out, 3, ai);
    SET_VECTOR_ELT(out, 4, Rf_ScalarInteger((int) a.fora_do_padrao));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(4);
    return out;
  )
}

SEXP R_ajustar_mt(SEXP dados, SEXP nomes, SEXP alvos, SEXP tnome, SEXP tcol, SEXP tcov,
                  SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                  SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP maxiter,
                  SEXP tol, SEXP gid, SEXP gm, SEXP mistura, SEXP anucleo, SEXP vk, SEXP verb, SEXP mfx, SEXP gmx, SEXP kern) {
  GUARDA(
    SEXP alvo1 = PROTECT(Rf_mkString(CHAR(STRING_ELT(alvos, 0))));
    br::Modelo m = modelo_do_R(alvo1, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc,
                               ausente, usa_ausente);
    UNPROTECT(1);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p2 = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p2, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)));
      pp = &ped;
    }
    std::vector<std::string> alv = textos(alvos, "traits");
    std::vector<br::KernelDecl> kd = kernels_do_R(kern);
    br::DesenhoMT d = br::monta_desenho_mt(m, alv, t, pp, &kd);

    std::string nota = genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk);

    br::AjusteMT r = br::ajusta_mt(d, (std::size_t) Rf_asInteger(maxiter), Rf_asReal(tol),
                                  Rf_asLogical(verb) == TRUE);
    if (!nota.empty()) r.mensagem = r.mensagem.empty() ? nota : nota + "; " + r.mensagem;

    std::vector<std::string> nomes_th = br::nomes_theta_do_mt(d);
    const std::size_t nt2 = d.modelo.ntheta;
    SEXP theta = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP se    = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP nms_t = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(theta)[k] = r.theta.empty() ? NA_REAL : r.theta[k];
      REAL(se)[k] = r.se.empty() ? NA_REAL : r.se[k];
      SET_STRING_ELT(nms_t, (R_xlen_t) k, Rf_mkChar(nomes_th[k].c_str()));
    }
    Rf_setAttrib(theta, R_NamesSymbol, nms_t);
    Rf_setAttrib(se, R_NamesSymbol, nms_t);

    // EBV e PEV por grupo: a fatia do grupo, com nomes "nivel|caracteristica" — a coluna
    // (tau * n_coef + ct) * n_niveis + nv e caracteristica-major no coeficiente
    SEXP ebv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) d.modelo.grupos.size()));
    SEXP pev = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) d.modelo.grupos.size()));
    SEXP ebv_nomes = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.modelo.grupos.size()));
    {
      std::size_t off = d.x.ncol * d.t;
      for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
        const std::size_t larg = d.largura(g);
        SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
        SEXP pv = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
        SEXP nv = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) larg));
        for (std::size_t c = 0; c < larg; c++) {
          REAL(v)[c] = r.solucao.empty() ? NA_REAL : r.solucao[off + c];
          REAL(pv)[c] = r.pev.empty() ? NA_REAL : r.pev[off + c];
        }
        std::size_t c0 = 0;
        for (std::size_t tm : d.modelo.grupos[g].termos)
          for (const br::DesenhoTermo& a : d.aleatorios)
            if (a.termo == tm) {
              for (std::size_t tau = 0; tau < d.t; tau++)
                for (std::size_t ct = 0; ct < a.n_coef; ct++)
                  for (std::size_t l = 0; l < a.n_niveis; l++) {
                    std::string nm = a.niveis[l] + "|" + d.alvos[tau];
                    if (a.n_coef > 1) nm += "[" + std::to_string(ct) + "]";
                    SET_STRING_ELT(nv, (R_xlen_t) (c0 + (tau * a.n_coef + ct) * a.n_niveis + l),
                                   Rf_mkChar(nm.c_str()));
                  }
              c0 += a.z.ncol * d.t;
            }
        Rf_setAttrib(v, R_NamesSymbol, nv);
        Rf_setAttrib(pv, R_NamesSymbol, Rf_duplicate(nv));
        SET_VECTOR_ELT(ebv, (R_xlen_t) g, v);
        SET_VECTOR_ELT(pev, (R_xlen_t) g, pv);
        SET_STRING_ELT(ebv_nomes, (R_xlen_t) g, Rf_mkChar(d.modelo.grupos[g].nome.c_str()));
        UNPROTECT(3);
        off += larg;
      }
    }
    Rf_setAttrib(ebv, R_NamesSymbol, ebv_nomes);
    Rf_setAttrib(pev, R_NamesSymbol, Rf_duplicate(ebv_nomes));

    // Solucoes dos efeitos fixos: a coluna j da caracteristica tau vive em j*t + tau (a
    // convencao do bloco fixo da montagem), nomeada "coluna|caracteristica" como o ebv.
    const std::size_t nb = d.x.ncol * d.t;
    SEXP bfix = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nb));
    SEXP bnms = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) nb));
    for (std::size_t j = 0; j < d.x.ncol; j++)
      for (std::size_t tau = 0; tau < d.t; tau++) {
        const std::size_t k = j * d.t + tau;
        REAL(bfix)[k] = r.solucao.empty() ? NA_REAL : r.solucao[k];
        const std::string nm = d.nomes_x[j] + "|" + d.alvos[tau];
        SET_STRING_ELT(bnms, (R_xlen_t) k, Rf_mkChar(nm.c_str()));
      }
    Rf_setAttrib(bfix, R_NamesSymbol, bnms);
    SEXP saiu = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_x.size()));
    for (std::size_t k = 0; k < d.saiu_x.size(); k++)
      SET_STRING_ELT(saiu, (R_xlen_t) k, Rf_mkChar(d.saiu_x[k].c_str()));

    const char* campos[] = {"theta", "se", "neg2logl", "converged", "iters", "reldelta",
                            "message", "n_used", "n_columns", "ebv", "pev", "b",
                            "dropped_x", "newton_dec"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 14));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 14));
    for (int q = 0; q < 14; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, theta);
    SET_VECTOR_ELT(out, 1, se);
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(r.neg2logl));
    SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(r.convergiu));
    SET_VECTOR_ELT(out, 4, Rf_ScalarInteger((int) r.iters));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(r.reldelta));
    SET_VECTOR_ELT(out, 6, Rf_mkString(r.mensagem.c_str()));
    SET_VECTOR_ELT(out, 7, Rf_ScalarInteger((int) d.n_usadas()));
    SET_VECTOR_ELT(out, 8, Rf_ScalarInteger((int) d.total_colunas()));
    SET_VECTOR_ELT(out, 9, ebv);
    SET_VECTOR_ELT(out, 10, pev);
    SET_VECTOR_ELT(out, 11, bfix);
    SET_VECTOR_ELT(out, 12, saiu);
    SET_VECTOR_ELT(out, 13, Rf_ScalarReal(r.decremento));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(11);
    return out;
  )
}


// ---------------------------------------------------------------- residuo AR(1)/CAR(1)

SEXP R_avaliar_ar1(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
                   SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                   SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP sujeito,
                   SEXP tempo, SEXP theta, SEXP com_densa, SEXP mfx, SEXP gmx, SEXP kern) {
  GUARDA(
    SEXP alvo1 = PROTECT(Rf_mkString(CHAR(STRING_ELT(alvo, 0))));
    br::Modelo m = modelo_do_R(alvo1, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc,
                               ausente, usa_ausente);
    UNPROTECT(1);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p2 = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p2, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)));
      pp = &ped;
    }
    std::vector<std::string> alv = textos(alvo, "traits");
    std::vector<br::KernelDecl> kd = kernels_do_R(kern);
    br::DesenhoAR d = br::monta_desenho_ar1(m, alv, t, pp,
                                            CHAR(STRING_ELT(sujeito, 0)),
                                            CHAR(STRING_ELT(tempo, 0)), &kd);
    std::vector<double> th(REAL(theta), REAL(theta) + XLENGTH(theta));
    if (th.size() != d.modelo.ntheta)
      Rf_error("theta with %d entries; the layout asks for %d", (int) th.size(),
               (int) d.modelo.ntheta);
    br::AvaliacaoAR a = br::avalia_ar1(d, th);
    if (!a.ok) Rf_error("theta INADMISSIBLE (s2e <= 0 or |rho| >= 1)");
    const double dv = Rf_asLogical(com_densa) == TRUE ? br::neg2logl_densa_V_ar1(d, th)
                                                      : std::nan("");
    const std::size_t nt2 = d.modelo.ntheta;
    SEXP score = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    for (std::size_t k = 0; k < nt2; k++) REAL(score)[k] = a.score[k];
    const char* campos[] = {"neg2logl", "neg2logl_V", "score", "off_pattern"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
    for (int q = 0; q < 4; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(a.neg2logl));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(dv));
    SET_VECTOR_ELT(out, 2, score);
    SET_VECTOR_ELT(out, 3, Rf_ScalarInteger((int) a.fora_do_padrao));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(3);
    return out;
  )
}

SEXP R_ajustar_ar1(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
                   SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                   SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP sujeito,
                   SEXP tempo, SEXP maxiter, SEXP tol, SEXP gid, SEXP gm, SEXP mistura,
                   SEXP anucleo, SEXP vk, SEXP verb, SEXP mfx, SEXP gmx, SEXP kern) {
  GUARDA(
    SEXP alvo1 = PROTECT(Rf_mkString(CHAR(STRING_ELT(alvo, 0))));
    br::Modelo m = modelo_do_R(alvo1, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc,
                               ausente, usa_ausente);
    UNPROTECT(1);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p2 = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p2, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)));
      pp = &ped;
    }
    std::vector<std::string> alv = textos(alvo, "traits");
    std::vector<br::KernelDecl> kd = kernels_do_R(kern);
    br::DesenhoAR d = br::monta_desenho_ar1(m, alv, t, pp,
                                            CHAR(STRING_ELT(sujeito, 0)),
                                            CHAR(STRING_ELT(tempo, 0)), &kd);

    std::string nota = genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk);

    br::AjusteMT r = br::ajusta_ar1(d, (std::size_t) Rf_asInteger(maxiter), Rf_asReal(tol),
                                   Rf_asLogical(verb) == TRUE);
    if (!nota.empty()) r.mensagem = r.mensagem.empty() ? nota : nota + "; " + r.mensagem;
    std::vector<std::string> nomes_th = br::nomes_theta_ar1(d);
    const std::size_t nt2 = d.modelo.ntheta;
    SEXP theta = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP se    = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nt2));
    SEXP nms_t = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) nt2));
    for (std::size_t k = 0; k < nt2; k++) {
      REAL(theta)[k] = r.theta.empty() ? NA_REAL : r.theta[k];
      REAL(se)[k] = r.se.empty() ? NA_REAL : r.se[k];
      SET_STRING_ELT(nms_t, (R_xlen_t) k, Rf_mkChar(nomes_th[k].c_str()));
    }
    Rf_setAttrib(theta, R_NamesSymbol, nms_t);
    Rf_setAttrib(se, R_NamesSymbol, Rf_duplicate(nms_t));
    // EBV e PEV por grupo, como no uni: fatias com os niveis como nomes
    SEXP ebv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) d.modelo.grupos.size()));
    SEXP pev = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) d.modelo.grupos.size()));
    SEXP ebv_nomes = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.modelo.grupos.size()));
    {
      std::size_t off = d.x.ncol;
      for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
        const std::size_t larg = d.largura(g);
        SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
        SEXP pv = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
        for (std::size_t c = 0; c < larg; c++) {
          REAL(v)[c] = r.solucao.empty() ? NA_REAL : r.solucao[off + c];
          REAL(pv)[c] = r.pev.empty() ? NA_REAL : r.pev[off + c];
        }
        {
          SEXP nv = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) larg));
          std::size_t c0 = 0;
          for (std::size_t tm : d.modelo.grupos[g].termos)
            for (const br::DesenhoTermo& a : d.aleatorios)
              if (a.termo == tm) {
                for (std::size_t tau = 0; tau < d.t; tau++)
                  for (std::size_t ct = 0; ct < a.n_coef; ct++)
                    for (std::size_t l = 0; l < a.n_niveis; l++) {
                      std::string nm = a.niveis[l];
                      if (d.t > 1) nm += "|" + d.alvos[tau];
                      if (a.n_coef > 1) nm += "[" + std::to_string(ct) + "]";
                      SET_STRING_ELT(nv,
                          (R_xlen_t) (c0 + (tau * a.n_coef + ct) * a.n_niveis + l),
                          Rf_mkChar(nm.c_str()));
                    }
                c0 += a.z.ncol * d.t;
              }
          Rf_setAttrib(v, R_NamesSymbol, nv);
          Rf_setAttrib(pv, R_NamesSymbol, Rf_duplicate(nv));
          UNPROTECT(1);
        }
        SET_VECTOR_ELT(ebv, (R_xlen_t) g, v);
        SET_VECTOR_ELT(pev, (R_xlen_t) g, pv);
        SET_STRING_ELT(ebv_nomes, (R_xlen_t) g, Rf_mkChar(d.modelo.grupos[g].nome.c_str()));
        UNPROTECT(2);
        off += larg;
      }
    }
    Rf_setAttrib(ebv, R_NamesSymbol, ebv_nomes);
    Rf_setAttrib(pev, R_NamesSymbol, Rf_duplicate(ebv_nomes));

    // Solucoes dos efeitos fixos, na convencao do bloco fixo (identica a da multi): a
    // coluna j da caracteristica tau em j*t + tau. Com t = 1 o nome e so termo=nivel.
    const std::size_t nb = d.x.ncol * d.t;
    SEXP bfix = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) nb));
    SEXP bnms = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) nb));
    for (std::size_t j = 0; j < d.x.ncol; j++)
      for (std::size_t tau = 0; tau < d.t; tau++) {
        const std::size_t k = j * d.t + tau;
        REAL(bfix)[k] = r.solucao.empty() ? NA_REAL : r.solucao[k];
        std::string nm = d.nomes_x[j];
        if (d.t > 1) nm += "|" + d.alvos[tau];
        SET_STRING_ELT(bnms, (R_xlen_t) k, Rf_mkChar(nm.c_str()));
      }
    Rf_setAttrib(bfix, R_NamesSymbol, bnms);
    SEXP saiu = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_x.size()));
    for (std::size_t k = 0; k < d.saiu_x.size(); k++)
      SET_STRING_ELT(saiu, (R_xlen_t) k, Rf_mkChar(d.saiu_x[k].c_str()));

    const char* campos[] = {"theta", "se", "neg2logl", "converged", "iters", "reldelta",
                            "message", "n_used", "n_columns", "n_subjects", "ebv", "pev",
                            "b", "dropped_x", "newton_dec"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 15));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 15));
    for (int q = 0; q < 15; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, theta);
    SET_VECTOR_ELT(out, 1, se);
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(r.neg2logl));
    SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(r.convergiu));
    SET_VECTOR_ELT(out, 4, Rf_ScalarInteger((int) r.iters));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(r.reldelta));
    SET_VECTOR_ELT(out, 6, Rf_mkString(r.mensagem.c_str()));
    SET_VECTOR_ELT(out, 7, Rf_ScalarInteger((int) d.n_usadas()));
    SET_VECTOR_ELT(out, 8, Rf_ScalarInteger((int) d.total_colunas()));
    SET_VECTOR_ELT(out, 9, Rf_ScalarInteger((int) d.sujeitos.size()));
    SET_VECTOR_ELT(out, 10, ebv);
    SET_VECTOR_ELT(out, 11, pev);
    SET_VECTOR_ELT(out, 12, bfix);
    SET_VECTOR_ELT(out, 13, saiu);
    SET_VECTOR_ELT(out, 14, Rf_ScalarReal(r.decremento));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(11);
    return out;
  )
}

SEXP R_gibbs(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
             SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid, SEXP ppai,
             SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP gid, SEXP gm, SEXP mistura,
             SEXP anucleo, SEXP n_iter, SEXP burnin, SEXP thin, SEXP loc_fixa,
             SEXP theta_fixo, SEXP vk, SEXP verb, SEXP mfx, SEXP gmx) {
  GUARDA(
    br::Modelo m = modelo_do_R(alvo, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc, ausente, usa_ausente);
    br::Tabela t = tabela_do_R(dados, nomes);
    br::Pedigree ped;
    const br::Pedigree* pp = nullptr;
    if (XLENGTH(pid) > 0) {
      auto i = textos(pid, "id");
      auto p = textos(ppai, "sire");
      auto ma = textos(pmae, "dam");
      ped = br::constroi_pedigree(i, p, ma, textos(mfx, "metafounders"), std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)));
      pp = &ped;
    }
    br::Desenho d = br::monta_desenho(m, t, pp);
    std::string nota = genomica_no_desenho(d, pp, ped, gid, gm, mistura, anucleo, vk);

    std::vector<double> thf;
    const std::vector<double>* pthf = nullptr;
    if (XLENGTH(theta_fixo) > 0) {
      thf.assign(REAL(theta_fixo), REAL(theta_fixo) + XLENGTH(theta_fixo));
      if (thf.size() != m.ntheta)
        Rf_error("theta with %d entries; the layout asks for %d", (int) thf.size(),
                 (int) m.ntheta);
      pthf = &thf;
    }
    br::GibbsSaida S = br::gibbs(d, (std::size_t) Rf_asInteger(n_iter),
                                 (std::size_t) Rf_asInteger(burnin),
                                 (std::size_t) Rf_asInteger(thin),
                                 Rf_asLogical(loc_fixa) == TRUE, pthf,
                                 Rf_asLogical(verb) == TRUE);
    if (!nota.empty()) S.mensagem = S.mensagem.empty() ? nota : nota + "; " + S.mensagem;

    std::vector<std::string> nomes_th = m.nomes_theta();
    SEXP amostras = PROTECT(Rf_allocMatrix(REALSXP, (int) S.n_amostras, (int) S.ntheta));
    for (std::size_t i = 0; i < S.n_amostras; i++)
      for (std::size_t k = 0; k < S.ntheta; k++)
        REAL(amostras)[k * S.n_amostras + i] = S.amostras[i * S.ntheta + k];
    SEXP nms_t = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) S.ntheta));
    for (std::size_t k = 0; k < S.ntheta; k++)
      SET_STRING_ELT(nms_t, (R_xlen_t) k, Rf_mkChar(nomes_th[k].c_str()));

    // medias e desvios das localizacoes por grupo, com os niveis como nomes
    br::Montado M = br::monta_mme(d, br::partida(d));
    SEXP ebv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) m.grupos.size()));
    SEXP esd = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) m.grupos.size()));
    SEXP ebv_nomes = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) m.grupos.size()));
    for (std::size_t g = 0; g < m.grupos.size(); g++) {
      const std::size_t larg = d.largura(g);
      SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
      SEXP w = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
      for (std::size_t c = 0; c < larg; c++) {
        const std::size_t k = M.offset_grupo[g] + c;
        REAL(v)[c] = S.media_loc.empty() ? NA_REAL : S.media_loc[k];
        REAL(w)[c] = S.var_loc.empty() ? NA_REAL : std::sqrt(S.var_loc[k]);
      }
      for (const br::DesenhoTermo& a : d.aleatorios)
        if (a.termo == m.grupos[g].termos[0]) {
          SEXP nv = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) larg));
          for (std::size_t c = 0; c < larg; c++)
            SET_STRING_ELT(nv, (R_xlen_t) c, Rf_mkChar(a.niveis[c % a.n_niveis].c_str()));
          Rf_setAttrib(v, R_NamesSymbol, nv);
          Rf_setAttrib(w, R_NamesSymbol, Rf_duplicate(nv));
          UNPROTECT(1);
          break;
        }
      SET_VECTOR_ELT(ebv, (R_xlen_t) g, v);
      SET_VECTOR_ELT(esd, (R_xlen_t) g, w);
      SET_STRING_ELT(ebv_nomes, (R_xlen_t) g, Rf_mkChar(m.grupos[g].nome.c_str()));
      UNPROTECT(2);
    }
    Rf_setAttrib(ebv, R_NamesSymbol, ebv_nomes);
    Rf_setAttrib(esd, R_NamesSymbol, Rf_duplicate(ebv_nomes));

    // Media e dp a posteriori dos efeitos fixos: as primeiras x.ncol posicoes do vetor de
    // localizacao, com os mesmos nomes termo=nivel que dropped_x usa.
    SEXP bfix = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) d.x.ncol));
    SEXP bsd  = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) d.x.ncol));
    SEXP bnms = PROTECT(Rf_allocVector(STRSXP,  (R_xlen_t) d.x.ncol));
    for (std::size_t k = 0; k < d.x.ncol; k++) {
      REAL(bfix)[k] = S.media_loc.empty() ? NA_REAL : S.media_loc[k];
      REAL(bsd)[k]  = S.var_loc.empty() ? NA_REAL : std::sqrt(S.var_loc[k]);
      SET_STRING_ELT(bnms, (R_xlen_t) k, Rf_mkChar(d.nomes_x[k].c_str()));
    }
    Rf_setAttrib(bfix, R_NamesSymbol, bnms);
    Rf_setAttrib(bsd, R_NamesSymbol, Rf_duplicate(bnms));
    SEXP saiu = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_x.size()));
    for (std::size_t k = 0; k < d.saiu_x.size(); k++)
      SET_STRING_ELT(saiu, (R_xlen_t) k, Rf_mkChar(d.saiu_x[k].c_str()));

    const char* campos[] = {"samples", "names", "ebv", "ebv_sd", "message", "n_used",
                            "b", "b_sd", "dropped_x"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 9));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 9));
    for (int q = 0; q < 9; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SET_VECTOR_ELT(out, 0, amostras);
    SET_VECTOR_ELT(out, 1, nms_t);
    SET_VECTOR_ELT(out, 2, ebv);
    SET_VECTOR_ELT(out, 3, esd);
    SET_VECTOR_ELT(out, 4, Rf_mkString(S.mensagem.c_str()));
    SET_VECTOR_ELT(out, 5, Rf_ScalarInteger((int) d.n_usadas()));
    SET_VECTOR_ELT(out, 6, bfix);
    SET_VECTOR_ELT(out, 7, bsd);
    SET_VECTOR_ELT(out, 8, saiu);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(11);
    return out;
  )
}

SEXP R_snp_blup(SEXP dados, SEXP nomes, SEXP alvo, SEXP tnome, SEXP tcol, SEXP tcov,
                SEXP test, SEXP tgrp, SEXP tnest, SEXP tbase, SEXP tsoc, SEXP pid,
                SEXP ppai, SEXP pmae, SEXP ausente, SEXP usa_ausente, SEXP gid, SEXP gm,
                SEXP rpg, SEXP theta, SEXP tol, SEXP maxiter, SEXP verb, SEXP mfx, SEXP gmx) {
  GUARDA(
    br::Modelo m = modelo_do_R(alvo, tnome, tcol, tcov, test, tgrp, tnest, tbase, tsoc, ausente, usa_ausente);
    br::Tabela t = tabela_do_R(dados, nomes);
    if (XLENGTH(pid) == 0) Rf_error("snp_blup needs a pedigree: the model is single step");
    auto i = textos(pid, "id");
    auto p = textos(ppai, "sire");
    auto ma = textos(pmae, "dam");
    br::Pedigree ped = br::constroi_pedigree(i, p, ma, textos(mfx, "metafounders"),
        std::vector<double>(REAL(gmx), REAL(gmx) + XLENGTH(gmx)));
    br::Desenho d = br::monta_desenho(m, t, &ped);

    std::vector<double> th(REAL(theta), REAL(theta) + XLENGTH(theta));
    if (th.size() != m.ntheta)
      Rf_error("theta with %d entries; the layout asks for %d", (int) th.size(),
               (int) m.ntheta);
    auto gids = textos(gid, "genotypes");
    if (gids.empty()) Rf_error("snp_blup without genotypes has nothing to solve");
    SEXP dim = Rf_getAttrib(gm, R_DimSymbol);
    if (TYPEOF(gm) != REALSXP || dim == R_NilValue || XLENGTH(dim) != 2)
      Rf_error("the genotype matrix must be numeric");
    const int nl2 = INTEGER(dim)[0], nm2 = INTEGER(dim)[1];
    if ((std::size_t) nl2 != gids.size())
      Rf_error("%d identifiers for %d genotype rows", (int) gids.size(), nl2);
    br::Densa mg((std::size_t) nl2, (std::size_t) nm2);
    for (int j2 = 0; j2 < nm2; j2++)
      for (int i2 = 0; i2 < nl2; i2++)
        mg.at((std::size_t) i2, (std::size_t) j2) = REAL(gm)[(R_xlen_t) j2 * nl2 + i2];

    br::SnpBlup S = br::snp_blup(d, mg, gids, th, Rf_asReal(rpg),
                                 (std::size_t) Rf_asInteger(maxiter), Rf_asReal(tol),
                                 Rf_asLogical(verb) == TRUE);

    // fatias como no ajuste exato: b nomeado pelas colunas de X, ebv por grupo
    SEXP bfix = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) S.n_fixo));
    SEXP bnms = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) S.n_fixo));
    for (std::size_t k = 0; k < S.n_fixo; k++) {
      REAL(bfix)[k] = S.solucao[k];
      SET_STRING_ELT(bnms, (R_xlen_t) k, Rf_mkChar(d.nomes_x[k].c_str()));
    }
    Rf_setAttrib(bfix, R_NamesSymbol, bnms);

    SEXP ebv = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) m.grupos.size()));
    SEXP ebv_nomes = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) m.grupos.size()));
    for (std::size_t g = 0; g < m.grupos.size(); g++) {
      const std::size_t larg = d.largura(g);
      SEXP v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) larg));
      for (std::size_t c = 0; c < larg; c++)
        REAL(v)[c] = S.solucao[S.offset_grupo[g] + c];
      for (const br::DesenhoTermo& a : d.aleatorios)
        if (a.termo == m.grupos[g].termos[0]) {
          SEXP nv = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) larg));
          for (std::size_t c = 0; c < larg; c++)
            SET_STRING_ELT(nv, (R_xlen_t) c, Rf_mkChar(a.niveis[c % a.n_niveis].c_str()));
          Rf_setAttrib(v, R_NamesSymbol, nv);
          UNPROTECT(1);
          break;
        }
      SET_VECTOR_ELT(ebv, (R_xlen_t) g, v);
      SET_STRING_ELT(ebv_nomes, (R_xlen_t) g, Rf_mkChar(m.grupos[g].nome.c_str()));
      UNPROTECT(1);
    }
    Rf_setAttrib(ebv, R_NamesSymbol, ebv_nomes);

    SEXP ef = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) S.efeitos.size()));
    for (std::size_t k = 0; k < S.efeitos.size(); k++)
      REAL(ef)[k] = std::isnan(S.efeitos[k]) ? NA_REAL : S.efeitos[k];

    std::string nota = "ssSNPBLUP: " + std::to_string(gids.size()) + " genotyped x " +
        std::to_string(nm2) + " markers, " + std::to_string(S.n_imputados) +
        " missing value(s) imputed by the mean, " + std::to_string(S.n_monomorficos) +
        " monomorphic marker(s) left out";

    const char* campos[] = {"b", "ebv", "g", "converged", "iters", "resnorm", "message",
                            "n_used", "n_columns", "dropped_x"};
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 10));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 10));
    for (int q = 0; q < 10; q++) SET_STRING_ELT(nms, q, Rf_mkChar(campos[q]));
    SEXP saiu = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) d.saiu_x.size()));
    for (std::size_t k = 0; k < d.saiu_x.size(); k++)
      SET_STRING_ELT(saiu, (R_xlen_t) k, Rf_mkChar(d.saiu_x[k].c_str()));
    SET_VECTOR_ELT(out, 0, bfix);
    SET_VECTOR_ELT(out, 1, ebv);
    SET_VECTOR_ELT(out, 2, ef);
    SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(S.convergiu));
    SET_VECTOR_ELT(out, 4, Rf_ScalarInteger((int) S.iters));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(S.residuo));
    SET_VECTOR_ELT(out, 6, Rf_mkString(nota.c_str()));
    SET_VECTOR_ELT(out, 7, Rf_ScalarInteger((int) d.n_usadas()));
    SET_VECTOR_ELT(out, 8, Rf_ScalarInteger((int) d.total_colunas()));
    SET_VECTOR_ELT(out, 9, saiu);
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(8);
    return out;
  )
}

SEXP R_versao(void) { return Rf_mkString("0.1.0"); }

static const R_CallMethodDef metodos[] = {
  {"R_pedigree",   (DL_FUNC) &R_pedigree,   5},
  {"R_a_inversa",  (DL_FUNC) &R_a_inversa,  5},
  {"R_inv_pd",       (DL_FUNC) &R_inv_pd,       1},
  {"R_chol_esparsa", (DL_FUNC) &R_chol_esparsa, 5},
  {"R_inv_seletiva", (DL_FUNC) &R_inv_seletiva, 5},
  {"R_resolve",      (DL_FUNC) &R_resolve,      5},
  {"R_avaliar",      (DL_FUNC) &R_avaliar,     28},
  {"R_ajustar",      (DL_FUNC) &R_ajustar,     32},
  {"R_a22_inversa",  (DL_FUNC) &R_a22_inversa,  4},
  {"R_avaliar_mt",   (DL_FUNC) &R_avaliar_mt,  21},
  {"R_ajustar_mt",   (DL_FUNC) &R_ajustar_mt,  27},
  {"R_avaliar_ar1",  (DL_FUNC) &R_avaliar_ar1, 23},
  {"R_ajustar_ar1",  (DL_FUNC) &R_ajustar_ar1, 29},
    {"R_gibbs", (DL_FUNC) &R_gibbs, 29},
  {"R_snp_blup",   (DL_FUNC) &R_snp_blup,  25},
{"R_versao",     (DL_FUNC) &R_versao,     0},
  {NULL, NULL, 0}
};

void R_init_BreedingR(DllInfo* dll) {
  R_registerRoutines(dll, NULL, metodos, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

}  // extern "C"
