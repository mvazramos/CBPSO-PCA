# Fit function
fitF <- function(A,Xs){
  #F<-sum(diag(t(Xs-Xs%*%A%*%t(A))%*%(Xs-Xs%*%A%*%t(A))))/(sum(diag(t(Xs)%*%Xs)) ) 
  
  F <-(norm(Xs-Xs%*%A%*%t(A), "F")^2/norm(Xs,"F")^2)
  
  return(F)
}


RandMat <- function(J, Q) {
  U0 <- matrix(0, J, Q)
  ###########################################################Changes made to implement the situation dim1=dim2 
  cc=1
  for ( i in sample(1:Q,Q) ){
    U0[cc,] <- diag(Q)[i,]
    cc=cc+1
  }  # end for    
  if (Q<J){
    for (j in (Q+1):J){
      p <- sample(1:Q, 1)
      U0[j, p] <- 1
    } # end for
  }  # end if
  U0
}  # end RandMat function


PMlargestEigen <- function(M){
  # Power Method
  # M covariance matrix restricted to some attributes
  #   (symmetric)
  tol=10^(-5)
  maxiteig <- 1000
  #x0 <- cbind(rep(1,nrow(M))) 
  x0 <- cbind(runif(nrow(M))) #; x0
  it <- 0 
  repeat  
  { 
    x1 = M %*% x0 ; 
    x1 = x1 / norm(x1, "F") #; x1; x1[1] <- -1 ; x1 = x1 / norm(x1, "F"); x1
    
    if (norm(abs(x1 - x0), "F") <= tol ) break
    x0 <- x1 
    it <- it + 1 
    if (it == maxiteig) break 
  }  
  lambda = sum((M %*% x1) * x1) 
  list(vector = x1, value = lambda, iterations = it) 
} 


compute_A <- function(I,J,Q,V,Xs){
  
  A <- matrix(0, J, Q)
  
  for (j in 1: Q) {
    
    Jxx <- which(V[, j] == 1)
    len.Jxx <- length(Jxx)
    
    if (sum(V[, j])>1) {
      if (I >= len.Jxx) {
        # CASE 1: I >= J (more observations than variables)
        S.group <- t(Xs[, Jxx])%*%Xs[, Jxx] 
        #
        #NEW
        #          
        #         A[Jxx, j] <- eigen(S.group)$vectors[, 1]
        A[Jxx, j] <- PMlargestEigen(S.group)$vector
        
        
      }else{
        # CASE 2: I < J ( more variables than observations)
        SS.group <- Xs[, Jxx]%*%t(Xs[, Jxx])         # O que faz isto??  May 7 2021 Question
        ##          A[Jxx, j] <- t(X.group[, Jxx])%*%eigen(SS.group)$vectors[, 1]/sqrt(eigen(SS.group)$values[1])  #Esta linha do CRAN ? substitu?da pela seguinte:
        PMinitial <- PMlargestEigen(SS.group)
        #PMinitial <- partial_eigen(SS.group,1)
        
        A[Jxx, j] <- t(Xs[, Jxx])%*%PMinitial$vector/sqrt(PMinitial$value)
      }  # end if	
    }else{
      A[Jxx, j] <- 1
    }  # end if
  }  
  
  return(A)
}


################ General PSO 
"
Parameters
  P: number of particles
  nIter: number of iterations
  minIner: minimum of intertia value
  maxIner: maximum of inertia value
  wC: cognitive weight
  wS: social weight
"



"
Each Particle is composed by a matrix :
  - V:  with de variable configurations   #Position
  - vel: with the velocities              #Velocity
"


################################################# Operator to get Valid Solutions 
randVel <- function(J,Q){
  "
  Initiates random velocity matrix with entries between [-1,1]
  "
  Vr <- matrix(0, J,Q)
  for(i in 1:J){
    for(j in 1:Q){
      rand <- runif(1,-1,1)
      Vr[i,j]<- rand
      
      
    }
  }
  return(Vr)
}


operatorO <- function(Mat,Vel,J,Q){
  
  "
                This is the O() operator that converts V(t+1) = V(t) + velocity 
                to a feasable solution!
  "
  tempMat= matrix(0,J,Q)
  
  Mat=abs(Mat)
  for (j in 1:J) {
    k <- which(Mat[j,]==max(Mat[j,]))
    tempMat[j,k]=1
  }
  tempMat <- checkColZeros(tempMat,Vel)
  return(tempMat)
}


operatorL <- function(Mat){
  M=Mat
  
  "
  To guarantee that the entries of the newvelocity is in the interval [-1,1]
  "
  Lz <- function(x) {
    return(2*(1+exp(1)^(-x))^(-1)-1)
  } 
  newMat=apply(M, 2, Lz)
  return(newMat)
}

#A=operatorO(Mat=randVel(J,Q),J,Q)
#Vel=A[[1]]

#TestMat = A[[2]]
#print(TestMat)

checkColZeros <- function(test,Vel){
  "
                  Function that makes sure that V is a valid matrix: binary-row stochastic 
                  Inputs:
                  - test: V matrix to check and correct
  "  
  Mat=test
  
  colunas <- colSums(Mat)   # Da a soma de cada coluna
  
  
  
  if(is.element(0,colunas)){
    
    lista <- which(colunas==0)   # Diz-me a coluna de zeros
    
    
    if(length(lista)==1){
      
      # Caso com apenas uma coluna cheia de zeros  
      
      col_zeros=lista  # indice da coluna com mais zeros
      max <- which(colunas==max(colunas))[1] # Diz-me a coluna com mais 1's
      lista_uns <- which(Mat[,max]==1) # Diz as linhas que tem 1's 
      i <- which(Vel[lista_uns,max]==min(abs(Vel[lista_uns,max])))[1]  # Diz-me a linha associada a menor velocidade 
      
      # Correção dos zeros
      Mat[i,col_zeros]=1
      Mat[i,max]=0
    } 
    else{
      
      # Caso com mais do que uma coluna com zeros
      temp_list=lista
      count=0
      
      while(length(temp_list)!=0){
        count=count+1
        col_zeros=temp_list[1]   # indice da coluna com mais zeros 
        max <- which(colunas==max(colunas))[1] # Diz-me a coluna com mais 1's
        lista_uns <- which(Mat[,max]==1) # Diz as linhas que tem 1's 
        
        i <- which(Vel[lista_uns,max]==min(abs(Vel[lista_uns,max])))  # Diz-me a linha associada a menor velocidade 
        
        # Correção dos zeros
        Mat[i,col_zeros]=1
        Mat[i,max]=0
        
        temp_list=temp_list[-1]
      }
      
      
    }
    
    
    
  }
  
  return(Mat)
}
