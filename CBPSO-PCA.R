source("CBPSOPCA_functions.R")
"""
Parameters:
  - X: Data matrix.
  - Xs: 0 if the Standardization is the one used by Macedo and Freitas, 1 if used like in Ramirez-Figueroa et al. and Vigneau and Qannari, 
        2 if data is only centered.
  - Q: Number of clusters of variables (disjoint principal components)
  - nPar: Number of particles in the swarm,
  - nIter: Maximum number of iterations ( this is the stopping criterion)
  - minIner: Minimum value of the Inertia function
  - maxIner: Maximum value of the Inertia function 
  - wCognition: Cognitive Weight Value (individual behaviour)
  - wSocial: Social Weight Value (collective behaviour)

Returns a list where:
  - the first element corresponds to the best V matrix found
  - the second element corresponds to the corresponding loadings matrix A
  - the third element corresponds to the value of the fit error
  - the fourth element corresponds to the value of the variances for each column of A ( in decreasing order)
  - the fifth element corresponds to the scores matrix Y=XA.
"""

PSO <- function(X,sX=0,Q,n_Par=300,nIter=30,minIner=0.5,maxIner=3.5,wCognition=1.5,wSocial=1.5){
  
  I <- dim(X)[1]
  J <- dim(X)[2]
  
  if(sX==0){
    "
    This performs the standardization as it is defined in Macedo and Freitas (2015).
    "
    X <- scale(X,center=T,scale=T)
    X <- as.matrix(X)
    X <- as.matrix(X*sqrt((I-1)/(I))) 
  }
  if(sX==1){
    "
    Standardization as assumed in Vigneau and Qannari, and Ramirez-Figueroa et al.
    "
    X <- scale(X,center=T,scale=T)
    X <- as.matrix(X)
  }
  
  if(sX==2){
    "
    Only centers the data, as it can also be considered in the works of Ramirez-Figueroa et al. 
    "
    X <- scale(X,center=T,scale=F)
    X <- as.matrix(X)
  }
  
  #################
  Vel=list()      # all the velocities for each particle
  
  Vp=list()   # all the V matrix for each particle 
  Ap=list()   # all the A matrix for each particle 
  
  Vp_best=list()  # the best V matrix for each particle
  Ap_best=list()  # the best A matrix for each particle
 
  
  Ap_Fit=c() # All the fit for each particle
  Ap_bestFit=c() # All the best fit for each particle
  
  # Step 1 - Initilization
  
  for (i in 1:n_Par) {
    # Compute loading matrix A_p
    # Compute fit of Vp ---> Fit(A_p)
    #Local + Local bests
    Vp[[i]]=RandMat(J,Q) # Matrix
    Vp_best[[i]]=Vp[[i]]
    Ap[[i]]=compute_A(I,J,Q,Vp_best[[i]],X) # Matrix
    Ap_best[[i]]=Ap[[i]]
    
    # Compute each fit 
    Ap_Fit[i]= fitF(Ap[[i]],X)  # double
    Ap_bestFit[i]=Ap_Fit[i]
    
    Vel[[i]]=randVel(J,Q) # Matrix
  }
  
  
  #Global 
  
  # Get index with best fit:
  p <- which(Ap_bestFit==min(Ap_bestFit))[1]    # the index that identifies the best particle 

  # Compute global best
   
  
  Vstar <- Vp_best[[p]]    # the V for the best particle
  Astar <-Ap_best[[p]]     # the A for the best particle
  FVstar <- fitF(Ap_best[[p]],X)  # the Fit for the best Particle 
  
  
# Iteration and Update 
  inertia=maxIner
  for (k in 1:nIter) {
    #print(FVstar)
    for (p in 1:n_Par) {
      #Iteration
      # newVelocity;
      rho1=runif(1,0,1)
      rho2=runif(1,0,1)
      newVelp= inertia*Vel[[p]] + rho1*wCognition*(Ap_best[[p]]-Ap[[p]])+rho2*wSocial*(Astar-Ap[[p]])
      
      Vel[[p]]= operatorL(M=newVelp) # Make velocity between [-1,1]
      Vtemp= Ap[[p]]+ Vel[[p]]       # New Position Vtemp
      Vtemp=operatorL(M=Vtemp)
      
      # Compute valid feasible solution
      Vp[[p]]=operatorO(Vtemp,Vel=Vel[[p]],J,Q) # Make Vtemp a valid solution with operator O
      
      # Compute A
      Ap[[p]]=compute_A(I,J,Q,Vp[[p]],X)
      # Compute fit o A
      Ap_Fit[p]= fitF(Ap[[p]],X)
      
      #Updating
      
      if (Ap_Fit[p]<Ap_bestFit[p]) {  # Faz update do Best de cada particula
        
        Vp_best[[p]]=Vp[[p]]       # Matrix V
        Ap_best[[p]]=Ap[[p]]       # Matrix A
        Ap_bestFit[p]=Ap_Fit[p]
      }
      
      if (Ap_Fit[p]<FVstar ) {       # Faz update do Best Global 
        Vstar <- Vp[[p]]    # the V for the best particle
        Astar <-Ap[[p]]     # the A for the best particle
        FVstar <- Ap_Fit[p]  # the best Fit value for the best particle
        
      }
      
    }
    
    inertia= maxIner- ((maxIner-minIner)/nIter)*k
    
  }
  
  
  if(sX==0){ 
    "
      In the case of Macedo & Freitas(2015) standardization.
      
    "
  Y <- X%*%Astar
  varPSO<- var(Y)
  d <- round(diag(varPSO)*100/J,2)
  dorder <- d[order(d, decreasing = TRUE)]
  Resultorder <- t(t(rbind(d,Astar ))[order(d, decreasing = TRUE), ])[-1, ]
  }
  else{
    "
    Accordingly to the Ramirez-Figueroa et al. considerations.
    "
    
    Y <- X%*%Astar
    varPSO<- var(Y)
    d <-round(diag(varPSO)*100/(sum(diag(var(X)))),2)
    dorder <- d[order(d, decreasing = TRUE)]
    Resultorder <- t(t(rbind(d,Astar ))[order(d, decreasing = TRUE), ])[-1, ]
  }
  
  return(list(as.matrix(Vstar),as.matrix(Resultorder),FVstar,dorder,Y))  
}


"""
In this section of the code, the results from Ramirez-Figueora et al. paper were reproduced, for the environmental data
that is referred in the paper and freely available: http://weppi.gtk.fi/publ/foregsatlas/ForegsData.php ( file C_XRF_data_2v6_8Feb06.xls in the subsoil data)
The results were successfully reproduced.
"""

# Q=3
# library(readxl)
# df <- read_excel("C_XRF_data_2v6_8Feb06.xls", sheet = "Sheet2")
# Result_XRF<- PSO(X=df,sX=1,Q=Q,n_Par=100,nIter=30,minIner=0.3,maxIner=2.5,wCognition=1.5,wSocial=2.5)
# Result_XRF[[2]] # Loadings matrix
# Result_XRF[[4]] # Explained variances by disjoint component (cluster)
# sum(Result_XRF[[4]]) # Total explained variance
