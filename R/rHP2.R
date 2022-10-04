#' @title  Hawkes Process simulation with different algorithms
#' @description This function simulates Hawkes Process with 
#' @description 1- Ogata's modified thinning algorithm (Ogata, 1981, p.25, Algorithm 2)
#' @description 2- Generating Hawkes Process by clusters (Hawkes and Oakes,1974)
#' @description The clusters of Hawkes Process, the Non-Homogeneous Poisson Processes are simulated with
#' @description an proposed algorithm that combines, inversion, order statistics and transformation method.
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
#' @param phi The excitation function of Hawkes Process
#' @param Phi The cumulative excitation function of Hawkes Process
#' @param PhiInv The inverse cumulative excitation function of Hawkes Process
#' @param method The simulation method, "Thinning", "Numerical Inversion","cf Inversion"
#'
#' @return The event times of n replications and the replication number of each event
#' @export  
#' @examples rHP(n=10000,T=1000,lam0=function(x) {2000*(x+250)^(-2.33)},phi=function(x) {2000*(x+250)^(-2.33)},Lam0=NULL,Phi=NULL,Lam0Inv=NULL,PhiInv=NULL,method="Thinning")
#' @examples rHP(n=10000,T=1000,lam0=function(x) {2000*(x+250)^(-2.33)},phi=function(x) {2000*(x+250)^(-2.33)},Lam0=NULL,Phi=NULL,Lam0Inv=NULL,PhiInv=NULL,method="Numerical Inversion")
#' @examples rHP(n=10000,T=1000,lam0=function(x) {2000*(x+250)^(-2.33)},phi=function(x) {2000*(x+250)^(-2.33)},Lam0=function(x) {(-2000/1.33)*(x+250)^(-1.33)+(2000/1.33*250^(-1.33))},Phi=function(x) {(-2000/1.33)*(x+250)^(-1.33)+(2000/1.33*250^(-1.33))},Lam0Inv=function(x) {(-(x-(2000/1.33*250^(-1.33)))*(1.33/2000))^(-(1/1.33))-250},PhiInv=function(x) {(-(x-(2000/1.33*250^(-1.33)))*(1.33/2000))^(-(1/1.33))-250},method="cf Inversion")
#'
rHP <- function (n,T,lam0,phi,Lam0,Phi,Lam0Inv,PhiInv,method=c("Thinning","Numerical Inversion","cf Inversion"))
{
  if (method == "Thinning") return (HP.thinnig (n,T,lam0,phi))

  if (method == "Numerical Inversion") return (HP_Composition_VECTOR(n,T,lam0,phi))
  # numerical inversion with rejected events

  if (method == "cf Inversion") return (HP_Composition_VECTOR_Closed_form (n,T,lam0,phi,Lam0,Phi,Lam0Inv,PhiInv))
  
  if (method != "Thinning" & method != "Numerical Inversion" & method != "cf Inversion")
    return (print("Please choose one of the method, thinning, numerical inversion or closed form inversion for the simulation algorithm"))
}
#--------------------------------------------------------------------------------
HP_Composition_VECTOR <- function(n,T,lam0,phi)
# with rejected events
{
  # computes the inverse of cumulative base intensity function, Lam0Inv(t) Using numerical inversion
  Lam0Inv                <- Runuran::pinv.new(pdf=lam0, lb=0, ub=T,uresolution=1e-010)
  # computes the inverse of cumulative excitation function, PhiInv(t) Using numerical inversion
  PhiInv                 <- Runuran::pinv.new(pdf=phi, lb=0, ub=T,uresolution=1e-010)
  # computes the total area under the base intensity function
  area.base <- vector()
  UNU.RAN.object.Lam0Inv    <- Runuran::unuran.details(Lam0Inv,show=FALSE, return.list=TRUE)
  area.base                 <- UNU.RAN.object.Lam0Inv$area.pdf
  # computes the total area under the excitation function
  area.exci <- vector()
  UNU.RAN.object.PhiInv    <-Runuran::unuran.details(PhiInv,show=FALSE, return.list=TRUE)
  area.exci                <-UNU.RAN.object.PhiInv$area.pdf
  
  # stores the event times of all n replications
  HawkesGen <- list()
  # stores the replication numbers of the event times of all n replications
  replication_no <- list()
  
  # NHPP simulation for Generation 1 (using base intensity - immigrant process) 
  # for n independent replications 
  # NHPP simulation for generation 1, Start-----------------------------------------------
  immigrant <- vector()
  # immigrant stores the number off springs for each immigrant
  immigrant <- rpois(n,area.base)
  # area base = total area under the base intensity function from 0 to T
  # n = the total number of independent replications of the Hawkes Process simulation
  # each of these n numbers represents the number of immigrant events in the replication (first generation events)
  # immigrant[immigrant>0]   indicates that the number off springs of the immigrants that have off springs
  # which (immigrant > 0)    indicates the index(replication number) of parent who have off springs
  # sum(immigrant)           indicates the total number off the off springs in the generation
  # runif(sum(immigrant) ,min = 0, max = area.base) 
  #         generates standard uniform random numbers for each offspring in the generation
  # Runuran::uq(Lam0Inv,(runif(sum(immigrant) ,min = 0, max = area.base)/area.base))
  #         generates the inter event times between immigrants and their off springs 
  start <- 0
  # i indicates the generation number
  i <- 1
  # generates the events of the generation 1 for n replications
  HawkesGen[[i]]   <- start + Runuran::uq(Lam0Inv, runif(sum(immigrant) ,min = 0, max = area.base)/area.base)
  # stores the replication number of the events of the generation 1
  replication_no[[i]]  <- rep(which (immigrant > 0), times=immigrant[immigrant>0])
  
  reject    <- list()
  # stores the index (replication number) of the event times that are outside of the simulation period (0,T) for generation i
  # There is no immigrant events created that is outside the simulation period (0,T) for generation 1.
  
  # NHPP simulation for generation 1, End-----------------------------------------------------
  
  # NHPP simulation for generations 2,3,... (using excitation function - offspring process) 
  # for n independent replications 
  # NHPP simulation for generations 2,3,..., Start-----------------------------------------------
  
  while (length(HawkesGen[[i]]) > 0 ) {
    offspring <- vector()
    offspring <- rpois(length(HawkesGen[[i]]),area.exci)
    # area.exci = total area under the excitation function from 0 to T
    # length(HawkesGen[[i]]) = the total number of living replications that can generate further off springs
    if(sum(offspring) == 0) {break}
    # checks for living replications
    i <- i+1
    # the generation number, i is updated
    # rep(HawkesGen[[i-1]], times= offspring) indicates that the parent event time of each offspring generated in this generation 
    # runif(sum(offspring),min = 0, max = area.exci)
    #            generates standard uniform random numbers for each offspring in the generation
    # Runuran::uq(PhiInv,runif(sum(offspring),min = 0, max = area.exci)/area.exci)
    #           generates the inter event times between parents and their off springs 
    # generates the events of the generation i for n replications
    HawkesGen[[i]]   <- rep(HawkesGen[[i-1]], times= offspring) + Runuran::uq(PhiInv,runif(sum(offspring),min = 0, max = area.exci)/area.exci)
    # stores the replication number of the events of the generation i
    replication_no[[i]]  <- rep(replication_no[[(i-1)]], times=offspring)
    # stores the index of the event time that is greater than simulation duration, T
    reject[[i]] <- which (HawkesGen[[i]] > T)
    if (length(reject[[i]])>0){
      HawkesGen[[i]] <- HawkesGen[[i]][-(reject[[i]])]
      replication_no[[i]] <- replication_no[[i]][-(reject[[i]])]      
    }
  }
  # NHPP simulation for generations 2,3,..., End-----------------------------------------------------
  
  # sorts the event times according to the replication number(y) and then event time(x)
  dataframe_replication <- data.frame(  x = unlist(HawkesGen),y =unlist(replication_no))
  sorted_dataframe_dataframe_replication <- dataframe_replication[with(dataframe_replication, order(y,x)),]
  output_event_replication <- list()
  # stores the sorted event times of all n replications
  output_event_replication[[1]]  <- sorted_dataframe_dataframe_replication[[1]] 
  # stores the replication number of the sorted events of all n replications
  output_event_replication[[2]]  <- sorted_dataframe_dataframe_replication[[2]]
  return(output_event_replication)
}  

# How to sort a vector according to a second vector
# https://chartio.com/resources/tutorials/how-to-sort-a-data-frame-by-multiple-columns-in-r/

#--------------------------------------------------------------------------------
HP.thinnig <- function(n,T,lam0,phi)
{
  # STEP 1.	Compute maximum value of lam0 (t) between time 0 and T, maxlam0
  maxlam0_initial     <- optimize(lam0, interval=c(0,T), maximum=TRUE)[[2]]
  # STEP 2.	Compute maximum value of phi(t) between time 0 and infinity, maxphi
  maxphi_initial      <- optimize(phi, interval=c(0,100), maximum=TRUE)[[2]]
  
  # stores the event times of all n replications
  HawkesReplicationEvent <- list()
  # stores the replication numbers of the event times of all n replications
  replication_no <- list()
  
  for (i in 1:n) {
    maxlam0 <- maxlam0_initial 
    maxphi  <- maxphi_initial
    # STEP 3.	lstar = maxlam0		# ... upper bound for the thinning algorithm
    lstar   <- maxlam0 
    
    # STEP 4. 	HawkesEv =  empty set 					# ... Hawkes Process(HP) event set
    HawkesEv    <- vector()
    # STEP 5.	 n = 0 						# . number of events in the set HawkesEv
    n <- 0
    # STEP 6.	curT = 0					# . current time
    curT <- 0
    # STEP 7.	Repeat
    repeat
    {
      # STEP 7.a.	Generate standard Exponential, E
      # STEP 7.b 	t = E\lstar				# . inter event time
      t <- rexp(1)/lstar
      # STEP 7.c.	curT = curT + t				# update current time
      curT <- curT + t
      # STEP 7.d.	If (curT > T), {Go to Step 8}
      if (curT > T) {break}
      # STEP 7.e.	sumphi = 0	# . sum of excitations triggered by the events in the process history
      sumphi <- 0
      # STEP 7.f. If (n > 0)
      if (n > 0)
      {
        # STEP 7.f.i. For i = 1 to n {sumphi = sumphi + phi (curT - HawkesEv[n])}
        for (n in 1:length(HawkesEv))
        {
          sumphi <- sumphi + phi (curT - HawkesEv[n])
        } 
      }
      # STEP 7.g.	curInt = lam0(curT) + sumphi		# .  current intensity
      curInt <- lam0(curT) + sumphi
      # STEP 7.h.	Generate standard Uniform, U
      # STEP 7.i 	If U L.E. curInt\lstar
      if (runif(1, min = 0, max = 1) <= curInt / lstar)
        # STEP 7.i.i.i.	then 
      {
        # STEP 7.i.i.1	n = n + 1
        n <- n + 1
        # STEP 7.i.i.2	HawkesEv[n] = curT	# Add the new event into the HP event set
        HawkesEv[n] <- curT
        # STEP 7.i.i.3.	Update maximum value of lam0(t) between time curT and T, maxlam0
        maxlam0     <- optimize(lam0, interval=c(curT,T), maximum=TRUE)[[2]]
        # STEP 7.i.i.4.	lambda star = maxlam0 + maxphi + sumphi
        lstar <- maxlam0 + maxphi + sumphi
      }
      # STEP 7.i.i.ii.	else 
      else
      {
        # STEP 7.i.i.ii.1 Update maximum value of lam0(t) between time curT and T, maxlam0
        maxlam0     <- optimize(lam0, interval=c(curT,T), maximum=TRUE)[[2]]
        # STEP 7.i.i.ii.2.	lambda star = maxlam0 + sumphi
        lstar <- maxlam0 + sumphi
      }
    }
    
    # stores the event times of all n replications
    HawkesReplicationEvent[[i]] <- HawkesEv
    # stores the replication numbers of the event times of all n replications
    replication_no[[i]]  <- rep(i, times=length(HawkesReplicationEvent[[i]]))
  }
  # STEP 8.	Return Hawkes Process event set and replication number
  output_event_replication <- list()
  # stores the event times of all n replications
  output_event_replication[[1]]  <- unlist(HawkesReplicationEvent) 
  # stores the replication number of the sorted events of all n replications
  output_event_replication[[2]]  <- unlist(replication_no)
  return(output_event_replication)
}
#--------------------------------------------------------------------------------
HP_Composition_VECTOR_Closed_form <- function(n,T,lam0,phi,Lam0,Phi,Lam0Inv,PhiInv)
{
  # computes the total area under the base intensity function
  area.base              <- Lam0(T)
  # computes the total area under the excitation function
  area.exci              <- Phi(T)
  
  # stores the event times of all n replications
  HawkesGen <- list()
  # stores the replication numbers of the event times of all n replications
  replication_no <- list()
  
  # NHPP for Generation 1 (using base intensity - immigrant event) 
  # for n independent cascades 
  # NHPP for generation for 1 Start-----------------------------------------------
  immigrant <- vector()
  # immigrant stores the number off springs for each immigrant
  immigrant <- rpois(n,area.base)
  # immigrant[immigrant>0]      indicates that the number off springs of the immigrants who have off springs
  # which (immigrant > 0)    indicates the index of parent who have off springs
  # sum(immigrant)           indicates the total number off the off springs in the generation
  # runif(sum(immigrant) ,min = 0, max = area.base) 
  #         generates standard uniform random numbers for each offspring in the generation
  
  start <- 0
  # i indicates the generation number
  i <- 1
  # generates the events of the generation 1 for n replications
  HawkesGen[[i]]   <- start + Lam0Inv(runif(sum(immigrant) ,min = 0, max = area.base))
  # stores the replication number of the events of the generation 1
  replication_no[[i]]  <- rep(which (immigrant > 0), times=immigrant[immigrant>0])
  
  reject    <- list()
  # stores the index of the event times that are outside of the simulation period (0,T) for generation i
  reject[[i]] <- which (HawkesGen[[i]] > T)
  if (length(reject[[i]])>0){
    HawkesGen[[i]] <- HawkesGen[[i]][-(reject[[i]])]
    replication_no[[i]] <- replication_no[[i]][-(reject[[i]])]      
  }
  
  # NHPP for generation 1 End-----------------------------------------------------
  
  while (length(HawkesGen[[i]]) > 0 ) {
    offspring <- vector()
    offspring <- rpois(length(HawkesGen[[i]]),area.exci)
    if(sum(offspring) == 0) {break}
    
    # NHPP for Generation 2,3,4,... (using excitation function - offspring event) 
    # for n independent cascades 
    # NHPP for generations 2,3,... Start--------------------------------------------     
    
    # the generation number, i is updated
    i <- i+1
    
    # rep(HawkesGen[[i-1]], times= offspring) indicates that the parent event time of each offspring generated in this generation 
    # runif(sum(offspring),min = 0, max = area.exci)
    #            generates standard uniform random numbers for each offspring in the generation
    
    # generates the events of the generation i for n replications
    HawkesGen[[i]]   <- rep(HawkesGen[[i-1]], times= offspring) + PhiInv(runif(sum(offspring),min = 0, max = area.exci))
    # stores the replication number of the events of the generation i
    replication_no[[i]]  <- rep(replication_no[[(i-1)]], times=offspring)
    # stores the index of the event time that is greater than simulation duration, T
    reject[[i]] <- which (HawkesGen[[i]] > T)
    if (length(reject[[i]])>0){
      HawkesGen[[i]] <- HawkesGen[[i]][-(reject[[i]])]
      replication_no[[i]] <- replication_no[[i]][-(reject[[i]])]      
    }
    # NHPP for generations 2,3,... END----------------------------------------------    
  }
  # sorts the event times according to the cascade number(y) and then event time(x)
  dataframe_replication <- data.frame(  x = unlist(HawkesGen),y =unlist(replication_no))
  sorted_dataframe_dataframe_replication <- dataframe_replication[with(dataframe_replication, order(y,x)),]
  
  output_event_replication <- list()
  # stores the sorted event times of all n replications
  output_event_replication[[1]]  <- sorted_dataframe_dataframe_replication[[1]] 
  # stores the sorted replication number of the events of all n replications
  output_event_replication[[2]]  <- sorted_dataframe_dataframe_replication[[2]]
  
  return(output_event_replication)
  #return(HawkesGen)
}  
#--------------------------------------------------------------------------------


