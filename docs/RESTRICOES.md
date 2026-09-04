# As 19 restricoes declaradas, e o estado de cada uma

Levantadas por varredura do fonte em 2026-09-03 (`grep` por "declared limit", "nesta
versao", "not implemented", "por enquanto" em `R/` e `src/`), mais os itens ABERTO do
CHECKLIST. Cada linha tem arquivo e linha de origem, para que ninguem precise confiar na
memoria de quem escreveu.

Estado: ABERTO / EM CURSO / FEITO / PROJETO (a restricao e matematicamente correta e o
entregavel e mensagem e documentacao, nao codigo).

## Ajustadores espelho: model_mt() e model_ar1()

| # | origem | restricao | estado |
|---|---|---|---|
| 1 | `src/multitrait2.cpp:118` | `kernel()` nao existe no multicaracter | ABERTO |
| 2 | `src/ar1b.cpp:73` | `kernel()` nao existe no AR(1) | ABERTO |
| 3 | `src/multitrait2.cpp:220` | theta cru, sem pisos nem log-Cholesky | ABERTO |
| 4 | `src/ar1b.cpp:274` | idem | ABERTO |
| 5 | `CHECKLIST:27` | posto de X sobre a tabela inteira, nao sobre as linhas usadas | ABERTO |

## Dados incompletos

| # | origem | restricao | estado |
|---|---|---|---|
| 6 | `src/ar1b.cpp:26`, `R/ar1.R:13` | AR(1) multicaracter exige registro COMPLETO entre tracos | ABERTO |

## ssSNPBLUP

| # | origem | restricao | estado |
|---|---|---|---|
| 7 | `src/sssnp.cpp:73` | um so grupo de parentesco | ABERTO |
| 8 | `src/sssnp.cpp:79` | termo genomico tem de ser grupo escalar | ABERTO |

## Covariancia declarada

| # | origem | restricao | estado |
|---|---|---|---|
| 9 | `src/mme.cpp:433` | `kernel()` so na rota `model()` (mesma familia de 1 e 2) | ABERTO |
| 10 | `CHECKLIST:743` | nao ha como FIXAR o componente de um `kernel` em 1 | ABERTO |

## Limiar e sobrevivencia

| # | origem | restricao | estado |
|---|---|---|---|
| 11 | `R/threshold.R:40` | componentes DADAS, nao estimadas: `start=` obrigatorio | ABERTO |
| 12 | `R/threshold.R:96`, `R/model.R:584` | `predict()` e PEV indisponiveis no modo conjunto | ABERTO |
| 13 | `R/survival.R:45` | sem covariaveis dependentes do tempo | ABERTO |

## Metafundadores

| # | origem | restricao | estado |
|---|---|---|---|
| 14 | `src/pedigree.cpp` | estimar Gamma dos genotipos; inversa generalizada para Gamma singular | ABERTO |
| 15 | `src/genomica.cpp`, `src/sssnp.cpp` | nao chegam ao lado genomico (DEFEITO, o unico nao declarado) | ABERTO |

## Restantes

| # | origem | restricao | estado |
|---|---|---|---|
| 16 | `R/model.R:575` | `accuracy()` le F do PEDIGREE mesmo em passo unico | ABERTO |
| 17 | `R/indirect.R:39` | `indirect_residual()` trata so heterogeneidade de VARIANCIA | ABERTO |
| 18 | `CHECKLIST:90` | pedigree pai/avo-materno aceito em silencio | ABERTO |
| 19 | `R/nonadditive.R` | D densa, sem inversa de Hoeschele & VanRaden, so epistasia A x A | ABERTO |

## Fora da lista, pedido a parte

- **Dados incompletos e desbalanceados em TODA funcao.** Varredura propria: o que cada
  ajustador faz com animal medido num traco e nao no outro, serie longitudinal com buraco,
  baia de tamanho desigual, nivel fixo com um registro so, e pedigree em que a maioria nao
  tem fenotipo. A restricao 6 e um caso particular disto.
- **Verbose com estado das estimativas.** Toda funcao que itera tem de imprimir, a cada
  iteracao, onde as estimativas estao e nao so o passo. Hoje algumas imprimem -2logL e
  passo relativo; falta o vetor de componentes.
