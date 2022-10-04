#' @title  Non-Homogeneous Poisson Process simulation with different algorithms
#' @description This function simulates Non-Homogeneous Poisson Process with 
#' @description 1- Ogata's modified thinning algorithm (Ogata, 1981, p.24, Algorithm 1)
#' @description 2- A proposed algorithm that combines Inversion, Order Statistics and Transformation method
#' @description Two inversion methods are proposed
#' @description a- The closed form of the inverse cumulative intensity function is available ("cf Inversion")
#' @description b- The closed form of the inverse cumulative intensity function is not available ("Numerical Inversion")
#' @description In the later case, the numerical inversion is used
#' 
#' @param n The number of independent replications for the Hawkes Process simulation
#' @param T The simulation duration
#' @param lam0  The base intensity function of Hawkes Process
#' @param Lam0  The cumulative base intensity function of Hawkes 
#' @param Lam0Inv The inverse cumulative base intensity function of Hawkes Process
#' @param method The simulation method, "Thinning", "Numerical Inversion","cf Inversion"
#'
#' @return The event times of n replications and the replication number of each event
#' @export  
#' @examples rNHPP(n=10000,T=1000,lam0=function(x) {2000*(x+250)^(-2.33)},Lam0=NULL,Lam0Inv=NULL,method="Thinning")
#' @examples rNHPP(n=10000,T=1000,lam0=function(x) {2000*(x+250)^(-2.33)},Lam0=NULL,Lam0Inv=NULL,method="Numerical Inversion")
#' @examples rNHPP(n=10000,T=1000,lam0=function(x) {2000*(x+250)^(-2.33)},Lam0=function(x) {(-2000/1.33)*(x+250)^(-1.33)+(2000/1.33*250^(-1.33))},Lam0Inv=function(x) {(-(x-(2000/1.33*250^(-1.33)))*(1.33/2000))^(-(1/1.33))-250},method="cf Inversion")
#'
rNHPP <- function (n,T,lam0,Lam0,Lam0Inv,method=c("Thinning","cf Inversion","Numerical Inversion"))
  {
  if (method == "Thinning") return (NHPP_Modified_Thinning(n,T,lam0))
  if (method == "cf Inversion") return (NHPP_Closed_form_INV_Order_VECTOR(n,T,lam0,Lam0,Lam0Inv))
  if (method == "Numerical Inversion") return (NHPP_INV_Order_VECTOR(n,T,lam0))
  if (method != "Thinning" & method != "Numerical Inversion" & method != "cf Inversion") 
    return (print("Please choose one of the method, thinning, closed form inversion or numerical inversion  for the simulation algorithm"))
}
#--------------------------------------------------------------------------------
NHPP_INV_Order_VECTOR <- function(n,T,lam0)
{
  # computes inverse of cumulative base intensity function, Lam0Inv(t)
  Lam0Inv                <- Runuran::pinv.new(pdf=lam0, lb=0, ub=T,uresolution=1e-010)
  # computes the total area under the base intensity function
  area.base <- vector()
  UNU.RAN.object.Lam0Inv    <- Runuran::unuran.details(Lam0Inv,show=FALSE, return.list=TRUE)
  area.base                 <- UNU.RAN.object.Lam0Inv$area.pdf
  # computes the total area under the excitation function
  
  # stores the event times of all n replications
  NHPP <- vector()
  # stores the replication numbers of the event times of all n replications
  replication_no <- vector()
  
  # for n independent NHPP 
  # NHPP simulation Start-----------------------------------------------
  immigrant <- vector()
  # immigrant stores the number off springs for each immigrant
  immigrant <- rpois(n,area.base)
  # area base = total area under the base intensity function from 0 to T
  # n = the total number of independent replications of the NHPP simulation  
  
  # sum(immigrant)           indicates the total number off the off springs in the all replications
  # runif(sum(immigrant) ,min = 0, max = area.base) 
  #         generates standard uniform random numbers for each offspring in the all replications
  # Runuran::uq(Lam0Inv,runif(sum(immigrant) ,min = 0, max = area.base)/area.base)
  # generates the events of NHPP
  NHPP   <- Runuran::uq(Lam0Inv,runif(sum(immigrant) ,min = 0, max = area.base)/area.base)
  # stores the replication number of the events 
  replication_no  <- rep(which (immigrant > 0), times=immigrant[immigrant>0])
  
  # NHPP simulation End-----------------------------------------------------
  
  # sorts the event times according to the cascade number(y) and then event time(x)
  dataframe_replication <- data.frame(  x = unlist(NHPP),y =unlist(replication_no))
  sorted_dataframe_dataframe_replication <- dataframe_replication[with(dataframe_replication, order(y,x)),]
  
  output_event_replication <- list()
  # stores the sorted event times of all n replications
  output_event_replication[[1]]  <- sorted_dataframe_dataframe_replication[[1]] 
  # stores the sorted replication number of the events of all n replications
  output_event_replication[[2]]  <- sorted_dataframe_dataframe_replication[[2]]
  
  return(output_event_replication)
  
}
#--------------------------------------------------------------------------------
NHPP_Modified_Thinning <- function(n,T,lam0)
{
  # STEP 1.	Compute maximum value of lam0 (t) between time 0 and T, maxlam0
  #maxlam0     <- optimize(lam0, interval=c(0,T), maximum=TRUE)[[2]]
  maxlam0     <- lam0(0)
  # STEP 2.	lstar = maxlam0		# ... upper bound for the thinning algorithm
  lstar       <- maxlam0
  
  # stores the event times of all n replications
  NHPP <- list()
  # stores the replication numbers of the event times of all n replications
  replication_no <- list()
  
  for (i in 1:n) {
    lstar       <- maxlam0
    # STEP 3. 	NHPPEv =  empty set 	# ... Non-Homogeneous Poisson Process(NHHP) event set
    NHPPEv    <- vector()
    # STEP 4.	 n = 0 						      # . number of events in the set NHPPEv
    j <- 0
    # STEP 5.	curT = 0					      # . current time
    curT <- 0
    # STEP 6.	Repeat
    repeat
    {
      # STEP 6.a.	Generate standard Exponential, E
      # STEP 6.b 	t = E\lstar				# . inter event time
      t <- rexp(1)/lstar
      # STEP 6.c.	curT = curT + t				# update current time
      curT <- curT + t
      # STEP 6.c.	If (curT > T), {Go to Step 7}
      if (curT > T) {break}
      # STEP 6.d.	curInt = lam0(curT)		# .  current intensity
      curInt <- lam0(curT)
      # STEP 6.e.	Generate standard Uniform, U
      # STEP 6.f 	If U L.E. curInt\lstar
      if (runif(1, min = 0, max = 1) <= curInt / lstar)
        # STEP 6.f.i.	then
      {
        # STEP 6.f.ii.	j = j + 1
        j <- j + 1
        # STEP 6.f.iii.	NHPP[n] = curT	# Add the new event into the NHHP event set
        NHPPEv[j] <- curT
      }
      lstar       <- curInt
    }  
    # stores the events of the ith  replications
    NHPP[[i]]   <- NHPPEv
    # stores the replication number of the events of NNHP
    replication_no[[i]]  <- rep(i, times=length(NHPPEv))
  }
  
  # NHPP  End-----------------------------------------------------
  
  output_event_replication <- list()
  output_event_replication[[1]]  <- unlist(NHPP)
  output_event_replication[[2]]  <- unlist(replication_no)
  return(output_event_replication)
}  
#--------------------------------------------------------------------------------
NHPP_Closed_form_INV_Order_VECTOR <- function(n,T,lam0,Lam0,Lam0Inv)
{
  
  # computes the total area under the base intensity function
  area.base <- vector()
  
  area.base                 <- Lam0(T)
  # computes the total area under the excitation function
  
  # stores the event times of all n replications
  NHPP <- vector()
  # stores the replication numbers of the event times of all n replications
  replication_no <- vector()
  
  # for n independent NHPP 
  # NHPP simulation Start-----------------------------------------------
  immigrant <- vector()
  # immigrant stores the number off springs for each immigrant
  immigrant <- rpois(n,area.base)
  # area base = total area under the base intensity function from 0 to T
  # n = the total number of independent replications of the NHPP simulation  
  
  # sum(immigrant)           indicates the total number off the off springs in the all replications
  # runif(sum(immigrant) ,min = 0, max = area.base) 
  #         generates standard uniform random numbers for each offspring in the all replications
  # Lam0Inv(runif(sum(immigrant) ,min = 0, max = area.base))
  # generates the events of NHPP
  NHPP   <- Lam0Inv(runif(sum(immigrant) ,min = 0, max = area.base))
  # stores the replication number of the events 
  replication_no  <- rep(which (immigrant > 0), times=immigrant[immigrant>0])
  
  # NHPP simulation End-----------------------------------------------------
  
  # sorts the event times according to the cascade number(y) and then event time(x)
  dataframe_replication <- data.frame(  x = unlist(NHPP),y =unlist(replication_no))
  sorted_dataframe_dataframe_replication <- dataframe_replication[with(dataframe_replication, order(y,x)),]
  
  output_event_replication <- list()
  # stores the sorted event times of all n replications
  output_event_replication[[1]]  <- sorted_dataframe_dataframe_replication[[1]] 
  # stores the sorted replication number of the events of all n replications
  output_event_replication[[2]]  <- sorted_dataframe_dataframe_replication[[2]]
  
  return(output_event_replication)
  
}  
#--------------------------------------------------------------------------------
