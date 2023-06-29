#' Converting periodic hsmm to periodic hmm
#'
#' @param omega Matrix of conditional transition probabilities
#' @param dm a list (in time) of lists (of states) of pmf vectors
#' @param eps accuracy parameter
#'
#' @return List of Sparse Gamma matrices
#' @export
#'
#' @examples
phsmm2phmm2 = function(omega,dm,eps=1e-25){
  # dm is now a list (in time) of lists (over states)
  L = length(dm) # length of one cycle
  mv = sapply(dm[[1]],length) # lengths of the pmf vectors
  m = length(mv) # number of state aggregates
  G_all = vector("list")
  F = vector("list") # computing all cdfs
  for(k in 1:L){
    Fk = vector("list")
    for (i in 1:m){ Fk[[i]] = cumsum(c(0,dm[[k]][[i]][-mv[i]])) }
    F[[k]] = Fk
  }
  for (t in 1:L){ # loop over time
    G = matrix(0,0,sum(mv)) # constructing an empty matrix with number of columns equal to final dimension (sum_i m_i)
    for (i in 1:m){ # loop over states
      mi = mv[i] # length of pmf for state aggregate i
      ci = numeric(mi)
      for (k in 0:(mi-1)){
        l = ifelse(t==k, L, (t-k)%%L)
        l = ifelse(l==0, L, l)
        ci[k+1] = ifelse(abs(1-F[[l]][[i]][k+1])>eps,dm[[l]][[i]][k+1]/(1-F[[l]][[i]][k+1]),1)
      }
      cim = ifelse(1-ci>0,1-ci,0)
      Gi = matrix(0,mi,0)
      for (j in 1:m){
        if(i==j) {
          if(mi==1){
            Gi = cbind(Gi,c(rep(0,mv[[j]]-1),cim))
          } else{
            Gi = cbind(Gi,rbind(cbind(rep(0,mi-1),diag(cim[-mi],mi-1,mi-1)),
                                c(rep(0,mi-1),cim[[mi]])))}
        } else   { if(mi==1)
        { Gi = cbind(Gi,matrix(c(omega[[i,j]]*ci,rep(0,mv[[j]]-1)),1))} else
        { Gi = cbind(Gi,cbind(omega[[i,j]]*ci,matrix(0,mv[[i]],mv[[j]]-1)))}
        }
      }
      G = rbind(G,Gi)
    }
    G_all[[t]] = SparseM::as.matrix.csr(G)
  }
  G_all
}
