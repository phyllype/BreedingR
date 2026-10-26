// Residuos do modelo misto: o esqueleto das linhas de W e o bootstrap parametrico sob o
// modelo ajustado.
//
// ## O residuo e a alavanca
//
// As MME do model() estao na escala de s2e: C_u = W~'W~ + kron(C_g^-1 s2e, K^-1), com W~ a
// linha ja escalada por sqrt(w). No otimo, o residuo condicional e^ = y~ - W~ s^ = s2e P y~ e
//
//   Var(e^_i) = s2e (1 - h_ii) / w_i,   h_ii = w~_i' C_u^-1 w~_i
//
// (Hodges 1998, os minimos quadrados aumentados). Os pares (j, k) de nao-zeros de UMA linha
// estao no padrao de W'W, logo no padrao do fator, e a inversa seletiva que a PEV ja monta os
// entrega exatamente: a alavanca nao custa fatoracao nenhuma. A forma Var(e^) =
// s2e (I - s2e W C^-1 W') que circulou no desenho nao fecha em dimensao; as corretas sao a de
// cima e R - W C*^-1 W' com as MME absolutas (ar1.cpp, multitrait.cpp).
//
// ## O bootstrap parametrico
//
// Com theta fixo no estimado: y* = X b^ + Z u* + e*, u*_g ~ N(0, C_g (x) K_g) e e* ~ N(0, R),
// resolvido com UM fator das MME em theta^ (o mesmo C_u para todas as replicas, porque C_u nao
// depende de y). O que sai por replica e o residuo condicional e^* e o ajustado fixo X b^*,
// exatamente o que o R calcula nos dados; uma estatistica qualquer calculada pelo MESMO codigo
// nos dados e nas replicas tem p de Monte Carlo valido sob o modelo ajustado (condicional em
// theta^). O sorteio de u*_g usa a Cholesky esparsa de K_g^-1 = L L' (na ordem de grau minimo):
// x = L^-T z tem covariancia K_g, e u*_{a} = soma_b chol(C_g)[a, b] x_b nos slots do grupo.
//
// O RNG e o do R (GetRNGstate / norm_rand), em ordem fixa (replica, grupo, slot, nivel, depois
// as linhas), entao set.seed() governa e o resultado nao depende do numero de threads.
//
// O preparo (prepara_nulos) monta o desenho, fatora as MME e prepara o sorteio UMA vez, e o
// objeto vive num ponteiro externo do R; cada lote so sorteia e resolve, e lotes seguidos
// continuam a mesma sequencia do RNG. O preparo resolve tambem os DADOS no mesmo fator: o
// residuo que sai e o do ajuste se, e so se, desenho, pedigree e genotipos forem os dele, e e
// assim que residual_diagnostics() confere o que recebeu. E le a diagonal de cada K_g pela
// inversa seletiva do fator de K_g^-1 que o sorteio ja tem, a variancia a priori exata de cada
// nivel para os valores geneticos padronizados.

#include "mme.h"
#include <R_ext/Random.h>
#include <Rmath.h>

namespace br {

// O erro padronizado (variancia 1) de uma replica: normal (familia 0), Laplace (familia < 0,
// curtose em excesso 3, todos os momentos finitos) ou t de Student com `familia` graus de
// liberdade (> 4), escalada. O normal e um norm_rand() so, como antes: a sequencia do RNG do
// caso normal nao muda.
double erro_padrao(double familia) {
  if (familia == 0.0) return norm_rand();
  // os sorteios em sequencia explicita: a ordem de avaliacao dos operandos nao e fixada
  const double a = familia < 0.0 ? exp_rand() : norm_rand();
  const double b = familia < 0.0 ? exp_rand() : Rf_rchisq(familia);
  if (familia < 0.0) return (a - b) / std::sqrt(2.0);
  return a * std::sqrt((familia - 2.0) / b);
}

LinhasW linhas_w(const Desenho& d, const std::vector<std::size_t>& offset_grupo) {
  LinhasW W;
  W.ptr.assign(d.nlin + 1, 0);
  // passada 1: quantas entradas por linha
  for (std::size_t i = 0; i < d.nlin; i++) {
    if (!d.usa[i]) continue;
    for (std::size_t e = d.x.linptr[i]; e < d.x.linptr[i + 1]; e++)
      if (d.x.valor[e] != 0.0) W.ptr[i + 1]++;
  }
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++)
    for (std::size_t t : d.modelo.grupos[g].termos)
      for (const DesenhoTermo& a : d.aleatorios)
        if (a.termo == t)
          for (std::size_t p = 0; p < a.z.valor.size(); p++) {
            const std::size_t i = a.z.linha[p];
            if (d.usa[i] && a.z.valor[p] != 0.0) W.ptr[i + 1]++;
          }
  for (std::size_t i = 0; i < d.nlin; i++) W.ptr[i + 1] += W.ptr[i];
  W.col.assign(W.ptr[d.nlin], 0);
  W.val.assign(W.ptr[d.nlin], 0.0);
  // passada 2: X primeiro, depois os grupos na ordem das MME (termos na ordem do grupo)
  std::vector<std::size_t> pos(W.ptr.begin(), W.ptr.end() - 1);
  for (std::size_t i = 0; i < d.nlin; i++) {
    if (!d.usa[i]) continue;
    for (std::size_t e = d.x.linptr[i]; e < d.x.linptr[i + 1]; e++) {
      const double v = d.x.valor[e];
      if (v == 0.0) continue;
      W.col[pos[i]] = d.x.coluna[e];
      W.val[pos[i]++] = v;
    }
  }
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++) {
    std::size_t col0 = offset_grupo[g];
    for (std::size_t t : d.modelo.grupos[g].termos)
      for (const DesenhoTermo& a : d.aleatorios)
        if (a.termo == t) {
          for (std::size_t c = 0; c < a.z.ncol; c++)
            for (std::size_t p = a.z.colptr[c]; p < a.z.colptr[c + 1]; p++) {
              const std::size_t i = a.z.linha[p];
              if (!d.usa[i] || a.z.valor[p] == 0.0) continue;
              W.col[pos[i]] = static_cast<std::uint32_t>(col0 + c);
              W.val[pos[i]++] = a.z.valor[p];
            }
          col0 += a.z.ncol;
        }
  }
  return W;
}

// x tal que L' x = z, com L triangular inferior em CSC (a coluna j de L e a linha j de L')
static void resolve_lt_no_lugar(const Csc& L, std::vector<double>& x) {
  const std::size_t n = L.ncol;
  for (std::size_t jj = n; jj-- > 0;) {
    double s = x[jj];
    double diag = 0.0;
    for (std::size_t p = L.colptr[jj]; p < L.colptr[jj + 1]; p++) {
      const std::size_t i = L.linha[p];
      if (i == jj) diag = L.valor[p];
      else s -= L.valor[p] * x[i];
    }
    x[jj] = s / diag;
  }
}

SorteioAleatorio prepara_sorteio(const Modelo& mo, const std::vector<Csc>& kinv,
                                 const std::vector<DesenhoTermo>& aleatorios,
                                 const std::vector<double>& theta,
                                 const std::vector<std::size_t>& offset_grupo) {
  SorteioAleatorio S;
  const std::size_t ng = mo.grupos.size();
  S.offset = offset_grupo;
  S.dim.assign(ng, 0);
  S.nl.assign(ng, 0);
  S.chol_c.resize(ng);
  S.perm_k.resize(ng);
  S.lk.resize(ng);
  for (std::size_t g = 0; g < ng; g++) {
    const Grupo& gr = mo.grupos[g];
    S.dim[g] = gr.dim;
    for (const DesenhoTermo& a : aleatorios)
      if (a.termo == gr.termos[0]) { S.nl[g] = a.n_niveis; break; }
    Densa c(gr.dim, gr.dim);
    for (std::size_t j = 0; j < gr.dim; j++)
      for (std::size_t i = j; i < gr.dim; i++) {
        const double v = theta[gr.theta_idx(i, j)];
        c.at(i, j) = v;
        c.at(j, i) = v;
      }
    // C_g semidefinida (uma correlacao em +-1, uma variancia em zero) e um ponto que o
    // ajuste aceita: o fator e o de Cholesky com pivo nulo zerado, que reproduz C_g
    Densa l = c;
    if (!chol_densa(l)) {
      // fallback espectral: C = U diag(ev) U', fator U diag(sqrt(max(ev, 0)))
      std::vector<double> ev;
      Densa u;
      jacobi_sim(c, ev, u);
      l = Densa(gr.dim, gr.dim);
      for (std::size_t i = 0; i < gr.dim; i++)
        for (std::size_t k = 0; k < gr.dim; k++)
          l.at(i, k) = u.at(i, k) * std::sqrt(std::max(ev[k], 0.0));
      S.chol_c[g] = l;
    } else {
      for (std::size_t i = 0; i < gr.dim; i++)
        for (std::size_t j = i + 1; j < gr.dim; j++) l.at(i, j) = 0.0;
      S.chol_c[g] = l;
    }
    if (kinv[g].ncol > 0) {
      if (kinv[g].ncol != S.nl[g])
        throw Erro("group '" + gr.nome + "': K^-1 of " + std::to_string(kinv[g].ncol) +
                   " rows for " + std::to_string(S.nl[g]) + " levels");
      S.perm_k[g] = grau_minimo(kinv[g]);
      const Csc pk = permuta_sim(kinv[g], S.perm_k[g]);
      const Simbolica sb = simbolica(pk);
      if (!cholesky(pk, sb, S.lk[g]))
        throw Erro("group '" + gr.nome + "': the precision K^-1 is not positive-definite, "
                   "so its levels cannot be drawn for the parametric bootstrap");
    }
  }
  return S;
}

void sorteia_aleatorios(const SorteioAleatorio& S, std::vector<double>& u) {
  for (std::size_t g = 0; g < S.dim.size(); g++) {
    const std::size_t dim = S.dim[g], nl = S.nl[g], off = S.offset[g];
    for (std::size_t a = 0; a < dim; a++)
      for (std::size_t l = 0; l < nl; l++) u[off + a * nl + l] = 0.0;
    std::vector<double> x(nl), z(nl);
    for (std::size_t b = 0; b < dim; b++) {
      for (std::size_t l = 0; l < nl; l++) z[l] = norm_rand();
      if (S.lk[g].ncol > 0) {
        resolve_lt_no_lugar(S.lk[g], z);
        for (std::size_t k = 0; k < nl; k++) x[S.perm_k[g][k]] = z[k];
      } else {
        x = z;
      }
      for (std::size_t a = 0; a < dim; a++) {
        const double f = S.chol_c[g].at(a, b);
        if (f == 0.0) continue;
        double* ua = u.data() + off + a * nl;
        for (std::size_t l = 0; l < nl; l++) ua[l] += f * x[l];
      }
    }
  }
}

std::vector<double> diagonal_de_k(const Csc& lk, const std::vector<std::size_t>& perm) {
  const SelInv z = inversa_seletiva(lk, 0);
  const std::vector<double> dz = z.diagonal();
  std::vector<double> out(dz.size());
  for (std::size_t k = 0; k < dz.size(); k++) out[perm[k]] = dz[k];
  return out;
}

std::vector<std::vector<double>> diagonais_k(const std::vector<Csc>& kinv) {
  std::vector<std::vector<double>> out(kinv.size());
  for (std::size_t g = 0; g < kinv.size(); g++) {
    if (kinv[g].ncol == 0) continue;
    const std::vector<std::size_t> perm = grau_minimo(kinv[g]);
    const Csc pk = permuta_sim(kinv[g], perm);
    const Simbolica sb = simbolica(pk);
    Csc L;
    if (!cholesky(pk, sb, L))
      throw Erro("a precision K^-1 of the fit is not positive-definite, so the prior "
                 "variance of its levels cannot be read");
    out[g] = diagonal_de_k(L, perm);
  }
  return out;
}

namespace {

// O preparo do model(): o desenho (dono dele), o fator das MME em theta na ordem perm, as
// linhas de W~, o sorteio de u* e X~ b^ por linha. Tudo o que nao muda entre replicas.
class PreparoUni : public PreparoNulos {
 public:
  Desenho d;
  std::vector<std::size_t> perm, usadas;
  Csc L;
  LinhasW W;
  SorteioAleatorio SA;
  std::vector<double> xb;
  double sd_e = 0.0;
  std::size_t q_tot = 0;

  NulosSaida lote(std::size_t nsim, bool com_u, double familia) override {
    NulosSaida S;
    const std::size_t nu = usadas.size(), nx = d.x.ncol;
    S.n_usadas = nu;
    S.nsim = nsim;
    S.residuo.assign(nu * nsim, 0.0);
    S.fixo.assign(nu * nsim, 0.0);
    if (com_u) S.u.assign((q_tot - nx) * nsim, 0.0);
    std::vector<double> u(q_tot, 0.0), ys(d.nlin, 0.0), rhs(q_tot), pb(q_tot), s(q_tot);
    GetRNGstate();
    try {
      for (std::size_t r = 0; r < nsim; r++) {
        checa_interrupcao();
        sorteia_aleatorios(SA, u);
        for (std::size_t i : usadas) {
          double v = xb[i];
          for (std::size_t p = W.ptr[i]; p < W.ptr[i + 1]; p++)
            if (W.col[p] >= nx) v += W.val[p] * u[W.col[p]];
          ys[i] = v + sd_e * erro_padrao(familia);
        }
        std::fill(rhs.begin(), rhs.end(), 0.0);
        for (std::size_t i : usadas)
          for (std::size_t p = W.ptr[i]; p < W.ptr[i + 1]; p++) rhs[W.col[p]] += W.val[p] * ys[i];
        for (std::size_t k = 0; k < q_tot; k++) pb[k] = rhs[perm[k]];
        const std::vector<double> px = resolve(L, pb);
        for (std::size_t k = 0; k < q_tot; k++) s[perm[k]] = px[k];
        double* er = S.residuo.data() + r * nu;
        double* fr = S.fixo.data() + r * nu;
        for (std::size_t k = 0; k < nu; k++) {
          const std::size_t i = usadas[k];
          double fa = 0.0, fx = 0.0;
          for (std::size_t p = W.ptr[i]; p < W.ptr[i + 1]; p++) {
            const double v = W.val[p] * s[W.col[p]];
            fa += v;
            if (W.col[p] < nx) fx += v;
          }
          er[k] = ys[i] - fa;
          fr[k] = fx;
        }
        if (com_u)
          for (std::size_t k = nx; k < q_tot; k++) S.u[r * (q_tot - nx) + (k - nx)] = u[k];
      }
    } catch (...) {
      PutRNGstate();
      throw;
    }
    PutRNGstate();
    return S;
  }
};

}  // namespace

std::unique_ptr<PreparoNulos> prepara_nulos(Desenho&& d0, const std::vector<double>& theta,
                                            const std::vector<double>& b) {
  std::unique_ptr<PreparoUni> P(new PreparoUni());
  P->d = std::move(d0);
  const Desenho& d = P->d;
  if (b.size() != d.x.ncol)
    throw Erro("b with " + std::to_string(b.size()) + " entries for " +
               std::to_string(d.x.ncol) + " fixed columns");
  // C, o rhs dos dados e o fator: C e a copia permutada saem de cena antes do sorteio (que
  // fatora as K^-1), e so o fator fica; num passo unico cada um deles tem o bloco denso
  std::vector<std::size_t> offset_grupo;
  {
    Montado M = monta_mme(d, theta);
    if (!M.ok) throw Erro("theta INADMISSIBLE: some covariance is not positive-definite");
    P->perm = grau_minimo(M.c);
    {
      const Csc pc = permuta_sim(M.c, P->perm);
      const Simbolica sb = simbolica(pc);
      if (!cholesky(pc, sb, P->L))
        throw Erro("the mixed model equations are not positive-definite");
    }
    P->q_tot = M.total;
    P->sd_e = std::sqrt(M.s2e);
    P->W = linhas_w(d, M.offset_grupo);
    offset_grupo = M.offset_grupo;
    for (std::size_t i = 0; i < d.nlin; i++)
      if (d.usa[i]) P->usadas.push_back(i);
    // o residuo dos dados no fator de agora: o mesmo numero do ajuste se tudo for o dele
    std::vector<double> pb(M.total), s(M.total);
    for (std::size_t k = 0; k < M.total; k++) pb[k] = M.rhs[P->perm[k]];
    const std::vector<double> px = resolve(P->L, pb);
    for (std::size_t k = 0; k < M.total; k++) s[P->perm[k]] = px[k];
    const LinhasW& W = P->W;
    P->residuo_dados.assign(P->usadas.size(), 0.0);
    for (std::size_t k = 0; k < P->usadas.size(); k++) {
      const std::size_t i = P->usadas[k];
      double fa = 0.0;
      for (std::size_t p = W.ptr[i]; p < W.ptr[i + 1]; p++) fa += W.val[p] * s[W.col[p]];
      P->residuo_dados[k] = d.y[i] - fa;
    }
  }
  P->SA = prepara_sorteio(d.modelo, d.kinv, d.aleatorios, theta, offset_grupo);
  // X~ b^ do ajuste, a media de y* em todas as replicas
  P->xb.assign(d.nlin, 0.0);
  for (std::size_t i : P->usadas)
    for (std::size_t p = P->W.ptr[i]; p < P->W.ptr[i + 1]; p++)
      if (P->W.col[p] < d.x.ncol) P->xb[i] += P->W.val[p] * b[P->W.col[p]];
  P->diag_k.resize(d.modelo.grupos.size());
  for (std::size_t g = 0; g < d.modelo.grupos.size(); g++)
    if (P->SA.lk[g].ncol > 0) P->diag_k[g] = diagonal_de_k(P->SA.lk[g], P->SA.perm_k[g]);
  P->niveis = niveis_grupos(d);
  P->usa = d.usa;
  return std::unique_ptr<PreparoNulos>(P.release());
}

}  // namespace br
