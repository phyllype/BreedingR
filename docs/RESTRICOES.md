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
| 3 | `src/multitrait2.cpp:220` | theta cru, sem pisos nem log-Cholesky | **FEITO** 2026-09-04 |
| 4 | `src/ar1b.cpp:274` | idem | **FEITO** 2026-09-04 |
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

### O item 3 e o 4, e uma correcao do que eu mesmo tinha escrito aqui

A primeira medicao que registrei nesta secao dizia que o multicaracter parava 7.2 unidades
de -2logL acima do otimo com a K escalada por 4, reportando `converged = TRUE`. O numero
estava certo e a LEITURA estava errada, e a correcao importa porque muda o que se conclui.
Aquele ajuste usava `y2 = 0.6*y + ruido`, que da correlacao genetica exatamente 1: o otimo
REML fica ON a fronteira `det(C) = 0`, e num otimo de fronteira o escore NAO zera, por
construcao. Perturbar o componente em 1e-06 ja devolvia theta inadmissivel. Nada do que se
media ali separava passo bom de passo ruim.

Refeita a medida num bivariado com correlacao genetica INTERIOR (r_g = 0.56) e com um
efeito de kernel simulado de verdade, `u ~ N(0, K (x) C)`, a mesma pergunta tem resposta
limpa. A verossimilhanca e invariante sob `(K, C) -> (cK, C/c)`, entao o otimo AJUSTADO tem
de ser tambem; isso e propriedade do caminhante, e nao da verossimilhanca.

| escala de K | antes | so com o passo em z | com passo, pisos, conjunto ativo e certificado em z |
|---|---|---|---|
| c = 1 | 307.195 (27 it) | 355.146 (1000 it) | **354.571004** (13 it) |
| c = 4 | 310.876 (37 it) | 355.975 (1000 it) | **354.571101** (15 it) |
| c = 25 | 453.552 (30 it) | 356.041 (1000 it) | **354.571126** (35 it) |
| dispersao | 146 unidades | 0.895 | **1.2e-04** |
| converged | TRUE nos tres | FALSE nos tres | TRUE nos tres |

Foram tres coisas, e nenhuma sozinha resolvia.

A primeira: o passo anda em log-Cholesky por bloco (Pinheiro e Bates, 1996), com o residuo
como bloco de t x t no multicaracter, e o rho do AR(1) em atanh. A segunda: a partida
divide pela media geometrica dos autovalores da K, que sai do `log|K^-1|` que o desenho ja
guarda, sem o que c = 1 e c = 25 partem de pontos diferentes da mesma superficie.

A terceira foi a que faltou por mais tempo, e era o defeito de verdade. O CERTIFICADO ficou
em theta quando o passo foi para z, com conjunto ativo proprio. O teste que ele usava,
`lmin < 1e-3 lmax`, e a exata desigualdade que o grampo do passo torna FALSA por
construcao: no ponto em que o laco para, o grampo deixa `lmin = 1e-3 lmax` e o `<` estrito
nunca dispara. A direcao grampeada ficava congelada no passo e cobrada integralmente no
certificado. Medido: decremento parado em 3.25e+03 com -2logL identico por 980 iteracoes,
num ponto que uma descida por coordenada melhorava em 0.200 unidade.

Duas coisas a mais sairam da mesma investigacao, as duas medidas. Uma direcao cuja
`dtheta/dz` colapsou (a diagonal presa tem `dtheta/dz = 2L^2`, e L = 4.0e-04 dava 3.2e-07)
nao carrega curvatura e nao pode ficar nem no passo nem no certificado; deixa-la devolvia
decremento 14.87 onde a folga real era 0.0063. E a parede prende o BLOCO, nao a coordenada:
o piso relativo fixa a RAZAO entre as diagonais de Cholesky, entao um bloco encostado nele
esta confinado a uma face de dimensao menor e a coordenada fora da diagonal, que nao tem
piso proprio, tambem nao e livre. Exigir score positivo para excluir reprovava o proprio
otimo restrito quando o multiplicador era numericamente zero: com `sz = -1e-04` o
decremento ficava em 0.034 num ponto onde 4000 sorteios multivariados e uma busca em linha
nao acharam melhora nenhuma, e onde perturbar o bloco em 1e-06 ja da theta inadmissivel.

O que sobrou fora dos espelhos e o RESGATE EM. O univariado tem um, os espelhos nao, e a
consequencia e visivel: numa celula de fronteira o AR(1) precisa de 726 iteracoes para
certificar um ponto 0.48 unidade melhor que o antigo, onde o univariado andaria. Por isso o
`maxiter` padrao dos espelhos e 1000 e nao 300. Os itens 14 a 16 do levantamento dizem que
o porte e viavel: `u`, `tr(K^-1 [C^-1]_ab)` e `q` ja existem nos dois espelhos com os mesmos
nomes, e sem o fator `s2e` do univariado, porque as MME dos espelhos ja sao absolutas.

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
