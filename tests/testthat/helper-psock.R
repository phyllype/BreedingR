# Os portoes de serie e paralelo conferem a cadeia, nao o socket: com a maquina carregada o
# cluster PSOCK de gibbs(cores = ) pode nao subir ("Cluster setup failed. 2 of 2 workers failed
# to connect"), o que nao e defeito do pacote. com_psock() avalia `expr`; se ela para com esse
# erro, avalia de novo uma vez e, se parar outra vez, pula o teste com a mensagem. Qualquer
# outro erro, e a falha de uma expectativa que nao venha do cluster, segue.
com_psock <- function(expr) {
  ex <- substitute(expr)
  env <- parent.frame()
  for (tentativa in 1:2) {
    v <- tryCatch(eval(ex, env), error = function(e) e)
    if (!inherits(v, "error")) return(v)
    if (!grepl("Cluster setup failed|failed to connect", conditionMessage(v))) stop(v)
  }
  testthat::skip(paste("o cluster PSOCK nao subiu duas vezes:", conditionMessage(v)))
}
