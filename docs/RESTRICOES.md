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
| 1 | `src/multitrait2.cpp:118` | `kernel()` nao existe no multicaracter | **FEITO** 2026-09-04 |
| 2 | `src/ar1b.cpp:73` | `kernel()` nao existe no AR(1) | **FEITO** 2026-09-04 |
| 3 | `src/multitrait2.cpp:220` | theta cru, sem pisos nem log-Cholesky | ABERTO, com custo MEDIDO |
| 4 | `src/ar1b.cpp:274` | idem | ABERTO |
| 5 | `CHECKLIST:27` | posto de X sobre a tabela inteira, nao sobre as linhas usadas | **FEITO** 2026-09-04 |

Os itens 1 e 2 foram fechados montando a K declarada uma unica vez, em
`kinv_declarada()` (`src/mme.cpp`), que os tres ajustadores passaram a chamar. Nao ha
segunda implementacao para divergir depois. O portao e `test-kernel-espelhos.R`, pela
identidade de escala: trocar K por cK e C por C/c tem de deixar -2logL intacto, porque o
c cancela entre `nl*log|C|` e `dim*log|K|`. Ele reprova se a K for ignorada, se a
identidade for posta no lugar dela, ou se for invertida errado.

Fechar 1 e 2 nao foi so ligar a inversao: o portao encontrou dois defeitos que ja
estavam la e ninguem tinha exercitado.

O primeiro DERRUBAVA A SESSAO. Os dois espelhos passavam `d.modelo` para a montagem da K
enquanto ainda iteravam `m.grupos`, e `d.modelo` so e preenchido bem depois
(`src/multitrait2.cpp:175`, `src/ar1b.cpp:158`): naquele ponto o modelo esta vazio. A unica
linha que o lia era a que monta o NOME do termo para a mensagem de erro, entao com uma K
boa nada acontecia e com uma K nao positiva-definida o indice caia fora do vetor e o R
morria com segmentation fault, sem erro, sem stack, no meio da suite. O conserto foi passar
o modelo vivo; a funcao compartilhada tambem passou a conferir o indice antes de usa-lo,
porque foi exatamente esse o tipo de leitura que escapou.

O segundo trocava o resultado por uma falha de montagem. Num termo com K declarada os
niveis tem de vir DA K, e nao da tabela: e assim que todo nivel da estrutura ganha equacao,
com registro ou sem. O univariado fazia isso; os espelhos passavam `nullptr` e tiravam os
niveis dos dados, de modo que uma K de 12 ids contra menos niveis observados terminava em
"triplet outside the matrix". Junto vinha a reducao por LINHA NULA da K (a convencao das
matrizes parciais multirraca, Mrode & Pocrnic 2023, p.243-244), que os espelhos tambem nao
tinham. As duas coisas viraram `reduz_kernels()` e `casa_niveis_nulos()` em `src/mme.cpp`,
chamadas pelos tres.

O custo do item 3, medido nesse mesmo portao. A identidade vale na VEROSSIMILHANCA nos
tres ajustadores, ao zero absoluto:

| ajustador | -2logL em (K, theta) | -2logL em (4K, theta/4) | diferenca |
|---|---|---|---|
| univariado | 139.952428 | 139.952428 | 0 |
| AR(1) | -20.951439 | -20.951439 | 0 |
| multicaracter | 288.19906066 | 288.19906066 | 0 |

No multicaracter a conta foi refeita tambem pela forma V densa, que nao compartilha
nenhuma linha de montagem com a esparsa: 288.19906058 contra 288.19906058. A
verossimilhanca esta certa.

O OTIMIZADOR nao esta. Partindo dos valores iniciais automaticos, o `ajusta_mt` para em
295.396991 no caso 4K contra 288.199061 no caso K, ou seja 7.2 unidades de -2logL pior, e
reporta `converged = TRUE`. Univariado e AR(1) chegam ao mesmo ponto nos dois casos. E o
preco exato de andar em theta cru: o passo nao tem piso, a convergencia nao e certificada
pelo decremento de Newton e a partida (`partida_mt()`, `src/multitrait2.cpp:185`) sai de
`var(y)` por traco, sem olhar a escala da K. Por isso o portao trava a identidade na
verossimilhanca e nao no otimo ajustado: o otimo ajustado mediria o item 3, nao o 1 e o 2.

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
| 10 | `CHECKLIST:743` | nao ha como FIXAR o componente de um `kernel` em 1 | **FEITO** 2026-09-04 |

## Limiar e sobrevivencia

| # | origem | restricao | estado |
|---|---|---|---|
| 11 | `R/threshold.R:40` | componentes DADAS, nao estimadas: `start=` obrigatorio | ABERTO |
| 12 | `R/threshold.R:96`, `R/model.R:584` | `predict()` e PEV indisponiveis no modo conjunto | ABERTO |
| 13 | `R/survival.R:45` | sem covariaveis dependentes do tempo | ABERTO |

## Metafundadores

| # | origem | restricao | estado |
|---|---|---|---|
| 14 | `src/pedigree.cpp` | estimar Gamma dos genotipos; inversa generalizada para Gamma singular | **PARCIAL** 2026-09-04: pseudo-inversa FEITA; estimar Gamma segue aberto |
| 15 | `src/genomica.cpp`, `src/sssnp.cpp` | nao chegam ao lado genomico (DEFEITO, o unico nao declarado) | ABERTO |

## Restantes

| # | origem | restricao | estado |
|---|---|---|---|
| 16 | `R/model.R:575` | `accuracy()` le F do PEDIGREE mesmo em passo unico | **FEITO** 2026-09-04 |
| 17 | `R/indirect.R:39` | `indirect_residual()` trata so heterogeneidade de VARIANCIA | ABERTO |
| 18 | `CHECKLIST:90` | pedigree pai/avo-materno aceito em silencio | ABERTO |
| 19 | `R/nonadditive.R` | D densa, sem inversa de Hoeschele & VanRaden, so epistasia A x A | **PARCIAL** 2026-09-04: A x D, D x D e ordens superiores feitos; a D densa e a inversa de Hoeschele & VanRaden seguem abertas |

## Fora da lista, pedido a parte

- **Dados incompletos e desbalanceados em TODA funcao.** Varredura propria: o que cada
  ajustador faz com animal medido num traco e nao no outro, serie longitudinal com buraco,
  baia de tamanho desigual, nivel fixo com um registro so, e pedigree em que a maioria nao
  tem fenotipo. A restricao 6 e um caso particular disto.
- **Verbose com estado das estimativas: FEITO em 2026-09-04.** Os quatro ajustadores
  iterativos (, , , ) imprimem o vetor de
  componentes a cada iteracao, com os nomes que o ajuste ja usa, quebrado em linhas de ate
  80 colunas. O impressor e um so,  em , para que os quatro
  nao divirjam. Portao: , 11 asserts, incluindo
  o silencio com  e a quebra de linha.

## Quantos tracos, e a que custo (MEDIDO em 2026-09-04)

Nao ha limite declarado em `model_mt()`. O unico limite duro do pacote e `t > 32` no
AR(1) (`src/ar1b.cpp:24`), e ele e CAP DEFENSIVO e nao estrutural: nao ha bitmask por
tras, entao o numero podia ser outro.

Segundos por iteracao, nesta maquina, modelo `cbind(...) ~ cg + animal(id)`:

| registros | 3 tracos | 6 tracos | 10 tracos |
|---|---|---|---|
| 350 | 0,01 | 0,09 | 0,51 |
| 1.050 | 0,11 | 0,83 | 4,27 |
| 2.800 | 1,27 | 10,07 | 50,73 |

Com poucos registros, empurrando so o numero de tracos (dois passos): 15 tracos 5,2 s;
20 tracos 20,4 s; 25 tracos 67,6 s; 30 tracos 310 s. Memoria nunca passou de 30 MB: o
gargalo e TEMPO, nao RAM.

**Teto pratico.** Um ajuste real pede 20 a 50 iteracoes. Ate 6 tracos e trivial; 10 tracos
em alguns milhares de registros e meia hora; 15 a 20 sao horas; 25 a 30 roda e nao e
usavel. Nos tracos o crescimento e pior que quadratico porque a AI e `ntheta x ntheta` com
`ntheta = t(t+1)/2` por grupo e e invertida a cada iteracao.

**ABERTO, e o que mais incomoda:** nos REGISTROS o custo medido cresce ~n^2,2, quando um
modelo misto esparso devia crescer bem melhor. Investigar se o caminho multicaracter faz
algo denso em n. Isto NAO estava registrado em lugar nenhum.
