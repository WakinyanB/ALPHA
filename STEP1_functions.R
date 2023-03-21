#########################################################################
#######           Phenotypic evolution of SARS-CoV-2:             #######
#######             a statistical inference approach              #######
#########################################################################
#######            STEP 1 - R Script with functions               #######
#########################################################################

  
# Packages
# (may first require installations: install.packages())

library(tidyverse) # for data manipulation
library(plyr) # for data manipulation
library(ggplot2) # for graphic visualizations
library(gridExtra) # to show multiple plots
library(cowplot) # to show multiple plots
library(deSolve) # to solve a system of ordinary differential equations (ODE integration, function 'ode')
library(boot) # for the inverse logit function 'inv.logit'
library(scales) # to manipulate the internal scaling infrastructure used by ggplot2
library(knitr) # for function 'kable' (a table generator)
library(crayon) # to colour some messages in the terminal
library(R.utils)
library(optimx)
library(investr)
library(MASS) # for matrix inversions (function 'ginv')
# For parallelized computations:
library(parallel)
library(bigparallelr)
library(foreach)

# Mathematical model (system of ODEs to describe the epidemiological dynamics of phase 1)

ODE_phase1 <- function(t, y, parms){
  
  # Parameters
  kappa = parms[["kappa"]] # 1/kappa is the average latency period
  gamma = parms[["gamma"]] # recovery rate
  alpha = parms[["alpha"]] # virulence (additional mortality due to disease)
  omega = parms[["omega"]] # once infected, probability of developing symptoms
  p = parms[["p"]] # probability to die of the disease when infected symptomatically
  if(any(names(parms) == "R0")){ # R0 is the basic reproduction number
                                 # beta is the transmission rate
    beta = parms[["R0"]]*gamma # beta = R0*gamma*alpha/(alpha+omega*p*(gamma-alpha))
                           # <=> beta \approx R0*gamma
  }else{
    beta = parms[["beta"]] 
  }
  if(any(names(parms) == "stringency_index")){
    # The Stringency Index reflects the severity of non-pharmaceutical interventions
    psi = parms[["stringency_index"]]
    k = parms[["k"]] # maximum achievable efficiency
    a = parms[["a"]] # shape parameter
    c = k*(psi[floor(t)+1]/100)^a
  }else{
    c = 0 # no control
  }
  # State variables
  S <- y[1] # Susceptible
  E <- y[2] # Exposed
  IA <- y[3] # Infectious (Asymptomatic)
  ISr <- y[4] # Infectious (Symptomatic who will recover)
  ISd <- y[5] # Infectious (Symptomatic who will die)
  R <- y[6] # Recovered
  D <- y[7] # Deceased
  
  IS_cumul <- y[8] # cumulated IS (ISr + ISd)
  
  I <- IA+ISr+ISd # Total number of I
  Nt <- S+E+I+R # total number of living individuals at time t
  Ni <- (1-c)*beta*S*I/Nt # New infections
  Nc <- kappa*E # New contagious individuals
  Ns <- omega*Nc # New symptomatic individuals
  Nd <- alpha*ISd # New deaths
  
  # Temporal derivatives
  dS <- -Ni
  dE <- Ni - Nc
  dIA <- Nc-Ns - gamma*IA
  dISr <- (1-p)*Ns - gamma*ISr
  dISd <- p*Ns - Nd
  dR <- gamma*(IA+ISr)
  dD <- Nd

  dIS_cumul <- Ns
  
  # Return results
  list(c(dS, dE, dIA, dISr, dISd, dR, dD, dIS_cumul))
}

# Alternative model where the ISd individuals does not participate in the force of infection

ODE_phase1.v2 <- function(t, y, parms){
  
  # Parameters
  kappa = parms[["kappa"]] # 1/kappa is the average latency period
  gamma = parms[["gamma"]] # recovery rate
  alpha = parms[["alpha"]] # virulence (additional mortality due to disease)
  omega = parms[["omega"]] # once infected, probability of developing symptoms
  p = parms[["p"]] # probability to die of the disease when infected symptomatically
  if(any(names(parms) == "R0")){ # R0 is the basic reproduction number
    # beta is the transmission rate
    beta = parms[["R0"]]*gamma # beta = R0*gamma*alpha/(alpha+omega*p*(gamma-alpha))
    # <=> beta \approx R0*gamma
  }else{
    beta = parms[["beta"]] 
  }
  if(any(names(parms) == "stringency_index")){
    # The Stringency Index reflects the severity of non-pharmaceutical interventions
    psi = parms[["stringency_index"]]
    k = parms[["k"]] # maximum achievable efficiency
    a = parms[["a"]] # shape parameter
    c = k*(psi[floor(t)+1]/100)^a
  }else{
    c = 0 # no control
  }
  # State variables
  S <- y[1] # Susceptible
  E <- y[2] # Exposed
  IA <- y[3] # Infectious (Asymptomatic)
  ISr <- y[4] # Infectious (Symptomatic who will recover)
  ISd <- y[5] # Infectious (Symptomatic who will die)
  R <- y[6] # Recovered
  D <- y[7] # Deceased
  
  IS_cumul <- y[8] # cumulated IS (ISr + ISd)
  
  I <- IA+ISr+ISd # Total number of I
  Nt <- S+E+I+R # total number of living individuals at time t
  Ni <- (1-c)*beta*S*(I-ISd)/Nt # New infections (ISd individuals are removed from the force of infection)
  Nc <- kappa*E # New contagious individuals
  Ns <- omega*Nc # New symptomatic individuals
  Nd <- alpha*ISd # New deaths
  
  # Temporal derivatives
  dS <- -Ni
  dE <- Ni - Nc
  dIA <- Nc-Ns - gamma*IA
  dISr <- (1-p)*Ns - gamma*ISr
  dISd <- p*Ns - Nd
  dR <- gamma*(IA+ISr)
  dD <- Nd
  
  dIS_cumul <- Ns
  
  # Return results
  list(c(dS, dE, dIA, dISr, dISd, dR, dD, dIS_cumul))
}

# Numerical simulations

simul_and_plot <- function(init, times, t = NULL, N, parms, func = ODE_phase1, solver = "ode45"){
  
  # init: initial values for state variables 
  # times: vector of time points at which the model will be numerically integrated
  # t: time points with change in the Stringency Index (optional)
  # N: size of the population
  # parms: parameters of the model
  # func: the function giving the ODE system of the model
  # solver: argument 'method' for function 'ode' ('ode45' by default)
  
  if(!any(names(parms) == "stringency_index")){
    cat(crayon::cyan("MESSAGE: no stringency Index detected \n"))
    # warns that there are no values for the Stringency Index (~ no control (default))
  }
  # Simulation
  SEIR_tab <- ode(y = init, times = times, parms, func = func, method = solver)
  
  # Data formatting for ggplot2
  data_plot <- data.frame("Time" = rep(SEIR_tab[,1],4),
                          "Numbers" = c(SEIR_tab[,2],
                                        SEIR_tab[,3],
                                        apply(SEIR_tab[,4:6], 1, sum),
                                        SEIR_tab[,7]),
                          "Compartment" = c(rep("S", length(times)), rep("E", length(times)),
                                            rep("I", length(times)), rep("R", length(times))))
  data_plot$Compartment <- factor(data_plot$Compartment, levels = c("S", "E", "I", "R"))
  
  Fig_epidemic <- ggplot(data_plot, aes(x = Time, y = Numbers, color = Compartment)) +
    geom_line(cex = 0.8) +
    theme_bw() +
    scale_color_manual(values = c("#157be2", "#b116bb", "#c21919", "#1fad13")) +
    theme(axis.text.x = element_text(size=11),
          axis.text.y = element_text(size=11),
          axis.title.x = element_blank(),
          legend.position = "top") +
    scale_y_continuous(trans='log10', limits = c(1e-5,N))

  if(!is.null(t)){
    Fig_epidemic <- Fig_epidemic + geom_vline(xintercept=t, linetype="dashed", col = "#5b5b5b")
  }
  return(list("tab" = SEIR_tab, "Fig" = Fig_epidemic))
}

# Statistical inference

## Functions for transformations of parameters

log_ab <- function(x, a, b){ # ]a,b[ -> ]-Inf, +Inf[
  # x is a number (or a vector of numbers) between a (min) and b (max)
  # for the specific case (a,b) = (0,1), use logit function instead
  
  if(length(x) == 1){
    if(is.na(x)){
      return(NA)
    }else if(x < a | x > b){
      stop(paste0("x (= ", x, ") is not between a (= ", a, ") and b (= ", b,")..."))
    }else{
      return(log((x-a)/(b-x)))
    }
  }else if(length(x) > 1){
    return(sapply(x, FUN=log_ab, a=a, b=b))
  }else{
    return(numeric())
  }
}

inv.log_ab <- function(y, a, b){ # R (]-Inf, +Inf[) -> ]a,b[
  # for the specific case (a,b) = (0,1), use inverse logit function instead
  
  if(length(y) == 1){
    
    num <- a+b*exp(y)
    denom <- 1+exp(y)
    
    if(is.infinite(num) | is.infinite(denom)){
      if(y %>% names %>% is.null){
        return(b)
      }else{
        return(setNames(object = b, nm = names(y)))
      }
    }else{
      return(num/denom)
    }
  }else if(length(y) > 1){
    return(sapply(y, FUN = inv.log_ab, a=a, b=b))
  }else{
    return(numeric())
  }
}

## Function defining Weighted Least Squares (WLS) value

WLS <- function(par, constants, model, approx_init_I,
                deaths_approx_init_IS, stringency_index, times, N, data_obs, solver){
  
  # par: initial values of parameters to be optimized (with potential log- or logit-transformation)
  # constants: fixed parameters, constants, etc... (NULL by default)
  # model: function specifying the ODEs of the model - here, ODE_phase1() -
  # approx_init_I: Do we approximate the initial states of compartment I ? (Boolean)
  # deaths_approx_init_IS: whether approx_init_I is TRUE, number of deaths for approximating
  #                        initialization of compartment I_S (by default: first number of deaths in data_obs)
  #                        NB : this is only useful if the simulation starts more than a day before the first data
  # times: vector of time points (observations start at time t=1)
  # N: size of the population
  # data_obs: matrix of observed data. Columns (/!\ order) : (1) daily negative tests,
  #                                                          (2) daily positive tests,
  #                                                          (3) daily fatality cases
  # solver: argument 'method' for function 'ode' ('ode45' by default)
  
  parms <- c(par[names(par) == "beta"] %>% exp,
             par[names(par) == "gamma"] %>% exp,
             par[names(par) == "alpha"] %>% inv.logit,
             par[names(par) == "kappa"] %>% exp,
             par[names(par) == "k"] %>% inv.logit,
             par[names(par) == "a"] %>% exp,
             par[names(par) == "omega"] %>% inv.logit,
             par[names(par) == "p"] %>% inv.log_ab(a=0,b=0.2),
             constants)
  
  parms$stringency_index <- stringency_index
  
  # Initial states
  if(any(names(constants) == "pS")){ # pS: initial proportion of S -> S(t=t0)
    init_simul <- constants[["pS"]]*N
  }else if(any(names(par) == "pS")){
    init_simul <- inv.logit(par[["pS"]])*N
  }else{
    return(cat(crayon::red("ERROR: initialisation - a value is missing for S(t=t0)...")))
  }
  if(any(names(constants) == "pE")){ # pE: initial proportion of E -> E(t=t0)
    init_simul <- c(init_simul, constants[["pE"]]*N)
  }else if(any(names(par) == "pE")){
    init_simul <- c(init_simul, inv.log_ab(par[["pE"]],a=0,b=0.01)*N)
  }else{
    return(cat(crayon::red("ERROR: initialisation - a value is missing for E(t=t0)...")))
  }
  if(approx_init_I){ # Boolean: approx_init_I = TRUE or FALSE
    if(is.null(deaths_approx_init_IS)){
      deaths_approx_init_IS <- data_obs[1,3]
    }
    ISd_0 <- deaths_approx_init_IS/parms[["alpha"]]
    ISr_0 <- ISd_0*(1-parms[["p"]])/parms[["p"]]
    init_simul <- c(init_simul, (ISr_0+ISd_0)*(1-parms[["omega"]])/parms[["omega"]], ISr_0, ISd_0)
  }else{
    if(any(names(par) == "pIA")){
      init_simul <- c(init_simul, inv.log_ab(par[["pIA"]],a=0,b=0.01)*N)
    }else if(any(names(constants) == "pIA")){
      init_simul <- c(init_simul, constants[["pIA"]]*N)
    }else{
      return(cat(crayon::red("ERROR: initialisation - a value is missing for I_A(t=t0)...")))
    }
    if(any(names(par) == "pISr")){
      init_simul <- c(init_simul, inv.log_ab(par[["pISr"]],a=0,b=0.01)*N)
    }else if(any(names(constants) == "pISr")){
      init_simul <- c(init_simul, constants[["pISr"]]*N)
    }else{
      return(cat(crayon::red("ERROR: initialisation - a value is missing for I_Sr(t=t0)...")))
    }
    if(any(names(par) == "pISd")){
      init_simul <- c(init_simul, inv.log_ab(par[["pISd"]],a=0,b=0.01)*N)
    }else if(any(names(constants) == "pISd")){
      init_simul <- c(init_simul, constants[["pISd"]]*N)
    }else{
      return(cat(crayon::red("ERROR: initialisation - a value is missing for I_Sd(t=t0)...")))
    }
  }
  init_simul <- c(init_simul, N-sum(init_simul), 0, 0)
  
  # Simulation
  simul <- ode(y = init_simul, times = times, parms = parms, func = model, method = solver)
  simul_obs <- simul[times >= 0,]
  d <- dim(simul_obs)[1]
  
  # rho, reporting rate of individuals S and IA
  if(any(names(par) == "rho")){ # if rho is constant through time
    rho <- par[["rho"]] %>% inv.logit
  }else if(any(names(constants) == "rho")){
    rho <- constants[["rho"]]
  }else{ # if rho is linear through time: rho(t) = mu + eta.t
    if(any(names(par) == "eta")){
      eta <- par[["eta"]] %>% inv.logit
    }else if(any(names(constants) == "eta")){
      eta <- constants[["eta"]]
    }else{
      return(cat(crayon::red("ERROR: Missing value for reporting rate...")))
    }
    if(any(names(par) == "mu")){
      mu <- par[["mu"]] %>% inv.logit
    }else if(any(names(constants) == "mu")){
      mu <- constants[["mu"]]
    }else{
      return(cat(crayon::red("ERROR: Missing value for reporting rate...")))
    }
    rho <- mu + eta*times[times > 0]
  }
  tests_neg_simul <- rho*simul_obs[-1,2]
  tests_pos_simul <- (simul_obs[-1,9]-simul_obs[-d,9]) + rho*simul_obs[-1,4]
  deaths_simul <- simul_obs[-1,8]-simul_obs[-d,8]
  
  data_exp <- cbind(tests_neg_simul, tests_pos_simul, deaths_simul) # expected
  
  # Preventing division by too small numbers
  min_threshold <- max(1e-3, min(data_obs))
  data_exp_cor <- data_exp
  data_exp_cor[data_exp_cor < min_threshold] <- min_threshold
  
  return(((data_obs - data_exp)^2/data_exp_cor) %>% sum) # (observed - expected)^2/expected
}

## Function for estimate parameters by optimizing (minimize) outputs of the previous function 'WLS' (using optim)

Estimate <- function(transf_parm_optim, constants = NULL, model = ODE_phase1, fn = WLS,
                     solver = "ode45", approx_init_I = FALSE, deaths_approx_init_IS = NULL,
                     times, N, stringency_index = NULL, data_obs, method = "Nelder-Mead",
                     maxit = 2000, abstol = 1e-6, reltol = 1e-6, save_fits = TRUE){

  # transf_parm_optim: initial values of parameters to be optimized (with potential log- or logit-transformation)
  # fn: function to minimized (here, the previous function 'WLS')
  # method, maxit, abstol and reltol are arguments for function 'optim' (package 'stats')
  # save_fits: logical indicating whether results of function 'optim' are to be saved in a RDS file
  # ...: for other arguments see previous function 'WLS'
  
  n_starts <- dim(transf_parm_optim)[1] # Number of set(s) of starting values
  if(is.null(n_starts)){ # special case where transf_parm_optim is a vector
    n_starts <- 1
    transf_parm_optim <- transf_parm_optim[1,] 
  }
  estim_tab <- data.frame() # data frame for saving results of non-linear optimization
  if(save_fits){
    RDS <- list() # list for saving results of 'optim' in a RDS file
  }
  
  if(is.null(stringency_index)){
    cat(crayon::cyan("MESSAGE: no stringency Index detected \n"))
  }
  for(i in 1:n_starts){

    par <- transf_parm_optim[i,]
    fit <- optim(par = par, fn = fn, constants = constants, model = model,
                 approx_init_I = approx_init_I, deaths_approx_init_IS = deaths_approx_init_IS,
                 stringency_index = stringency_index, times = times, N = N, data_obs = data_obs,
                 solver = solver, method = method, control = list(maxit = maxit,
                                                                  abstol = abstol, reltol = reltol))
  
    # Saving repeat N°, WLS value, convergence features and estimates of parameter
    estim_tab <- rbind(estim_tab, c(i, fit$value, fit$convergence,
                                    fit$par[names(fit$par) == "pS"] %>% inv.logit,
                                    fit$par[names(fit$par) == "pE"] %>% inv.log_ab(a=0,b=0.01),
                                    fit$par[names(fit$par) == "pIA"] %>% inv.log_ab(a=0,b=0.01),
                                    fit$par[names(fit$par) == "pISr"] %>% inv.log_ab(a=0,b=0.01),
                                    fit$par[names(fit$par) == "pISd"] %>% inv.log_ab(a=0,b=0.01),
                                    fit$par[names(fit$par) == "beta"] %>% exp,
                                    fit$par[names(fit$par) == "gamma"] %>% exp,
                                    fit$par[names(fit$par) == "alpha"] %>% inv.logit,
                                    fit$par[names(fit$par) == "kappa"] %>% exp,
                                    fit$par[names(fit$par) == "k"] %>% inv.logit,
                                    fit$par[names(fit$par) == "a"] %>% exp,
                                    fit$par[names(fit$par) == "omega"] %>% inv.logit,
                                    fit$par[names(fit$par) == "p"] %>% inv.log_ab(a=0,b=0.2),
                                    fit$par[names(fit$par) == "rho"] %>% inv.logit,
                                    fit$par[names(fit$par) == "eta"] %>% inv.logit,
                                    fit$par[names(fit$par) == "mu"] %>% inv.logit))
    if(save_fits){
      RDS[[i]] <- fit
    }
    print(paste("Progress:", i, "out of", n_starts, "completed"))
  }
  estim_tab[,3] <- gsub(10, "Degeneracy of the Nelder-Mead simplex", estim_tab[,3])
  estim_tab[,3] <- gsub(0, "Successful completion", estim_tab[,3])
  estim_tab[,3] <- gsub(1, "iteration limit had been reached", estim_tab[,3])
  
  names <- c("pS", "pE", "pIA", "pISr", "pISd", "beta", "gamma", "alpha", "kappa", "k", "a",
             "omega", "p", "rho", "eta", "mu")
  if(n_starts == 1){
    par_names <- names(transf_parm_optim)
  }else{
    par_names <- colnames(transf_parm_optim)
  }
  colnames(estim_tab) <- c("Repeat", "Value", "Convergence",
                           names[par_names %>% match(names) %>% sort])
  
  color_conv <- c("Successful completion" = "#00BFC4",
                  "iteration limit had been reached" = "#F8766D",
                  "Degeneracy of the Nelder-Mead simplex" = "#C77CFF")
  
  return(list("Results" = estim_tab,
              "fits" = RDS,
              "Fig" = ggplot(data = estim_tab, aes(x = Repeat, y = Value)) +
                geom_line(cex = 0.3) +
                geom_point(aes(col = Convergence), cex = 2, pch = 19) +
                theme_bw() +
                xlab("Optim() repeat") +
                ylab("Minimum WLS value") +
                ggtitle("Results (value and convergence features) returned by function optim()") +
                scale_color_manual(values = color_conv[estim_tab$Convergence %>% unique %>%
                                                         match(names(color_conv))]) +
                theme(axis.text.x = element_text(size=11),
                      axis.text.y = element_text(size=11),
                      legend.position = "bottom")))
}

# Same function but with parallelized computations

Estimate_parallel <- function(transf_parm_optim, constants = NULL, model = ODE_phase1, fn = WLS,
                              solver = "ode45", approx_init_I = FALSE, deaths_approx_init_IS = NULL,
                              times, N, stringency_index = NULL, data_obs, method = "Nelder-Mead",
                              maxit = 2000, abstol = 1e-6, reltol = 1e-6, Ncores = NA){
  
  # transf_parm_optim: initial values of parameters to be optimized (with potential log- or logit-transformation)
  # fn: function to minimized (here, the previous function 'WLS')
  # method, maxit, abstol and reltol are arguments for function 'optim' (package 'stats')
  # Ncores : number of cores for parallelized optimizations
  
  # see previous function 'WLS' for details for the other arguments
  
  n_parms <- ncol(transf_parm_optim) # Number of estimated parameter(s)
  n_starts <- nrow(transf_parm_optim) # Number of set(s) of starting values
  if(is.null(n_starts)){ # special case where transf_parm_optim is a vector
    n_starts <- 1
    transf_parm_optim <- transf_parm_optim[1,] 
  }
  estim_tab <- data.frame() # data frame for saving results of non-linear optimization
  
  if(is.null(stringency_index)){
    cat(crayon::cyan("MESSAGE: no stringency Index detected \n"))
  }
  
  # Parallelized optimizations
  
  if(is.na(Ncores)){Ncores <- bigparallelr::nb_cores()}
  print(paste("Number of cores:", Ncores))
  
  t <- format(Sys.time(), "%Y-%m-%d %Hh%M")
  folder_name <- paste0("ALPHA_Estim_optim_parallel_",
                        stringi::stri_replace_all_fixed(t, c("-", " "), "_", vectorize_all = F))
  dir.create(folder_name)
  lapply(list("PROJECT ALPHA\n", t, paste0("\nParallel - number of cores: ", Ncores),
              paste("\nODEs solver:", solver),
              paste("\nOptimization algorithm: optim() with method", method),
              "\nControl:", paste(c("Max iterations", "Absolute tolerance", "Relative tolerance"),
                                  c(maxit, abstol, reltol), sep = " = ", collapse = "\n"),
              paste("\nNumber of starts:", n_starts),
              paste0("\nParameters to estimate (", n_parms, "):"), colnames(transf_parm_optim) %>% paste(collapse = ", "),
              "\nConstants:", paste(names(constants), constants, sep = " = ", collapse = "\n")),
         cat, "\n", file=paste0(folder_name,"/README.txt"), append=TRUE)
  
  cl <- parallel::makeCluster(Ncores)
  doParallel::registerDoParallel(cl)
  
  estim_tab <- foreach::foreach(i = 1:n_starts, .combine = 'rbind', .inorder = FALSE,
                                .packages = c("tidyverse", "deSolve", "boot"), .export = ls(globalenv())) %dopar% {
                                  
                                  fit <- optim(par = transf_parm_optim[i,], constants = constants, fn = fn, model = model,
                                               approx_init_I = approx_init_I, deaths_approx_init_IS = deaths_approx_init_IS,
                                               stringency_index = stringency_index, times = times, N = N, data_obs = data_obs,
                                               solver = solver, method = method,
                                               control = list(maxit = maxit, abstol = abstol, reltol = reltol))
                                  
                                  # Saving repeat N°, WLS value, convergence features and estimates of parameter
                                  res <- c("Repeat" = i, "Value" = fit$value, "Convergence" = fit$convergence,
                                           fit$par[names(fit$par) == "pS"] %>% inv.logit,
                                           fit$par[names(fit$par) == "pE"] %>% inv.log_ab(a=0,b=0.01),
                                           fit$par[names(fit$par) == "pIA"] %>% inv.log_ab(a=0,b=0.01),
                                           fit$par[names(fit$par) == "pISr"] %>% inv.log_ab(a=0,b=0.01),
                                           fit$par[names(fit$par) == "pISd"] %>% inv.log_ab(a=0,b=0.01),
                                           fit$par[names(fit$par) == "beta"] %>% exp,
                                           fit$par[names(fit$par) == "gamma"] %>% exp,
                                           fit$par[names(fit$par) == "alpha"] %>% inv.logit,
                                           fit$par[names(fit$par) == "kappa"] %>% exp,
                                           fit$par[names(fit$par) == "k"] %>% inv.logit,
                                           fit$par[names(fit$par) == "a"] %>% exp,
                                           fit$par[names(fit$par) == "omega"] %>% inv.logit,
                                           fit$par[names(fit$par) == "p"] %>% inv.log_ab(a=0,b=0.2),
                                           fit$par[names(fit$par) == "rho"] %>% inv.logit,
                                           fit$par[names(fit$par) == "eta"] %>% inv.logit,
                                           fit$par[names(fit$par) == "mu"] %>% inv.logit)
                                  
                                  saveRDS(res, file = paste0(folder_name, "/result_task_", i, ".rds"))
                                  return(res)
                                }
  parallel::stopCluster(cl)
  
  estim_tab <- as.data.frame(estim_tab)
  estim_tab <- estim_tab[order(estim_tab$Repeat),]
  
  estim_tab$Convergence <- gsub(10, "Degeneracy of the Nelder-Mead simplex", estim_tab$Convergence)
  estim_tab$Convergence <- gsub(0, "Successful completion", estim_tab$Convergence)
  estim_tab$Convergence <- gsub(1, "iteration limit had been reached", estim_tab$Convergence)
  
  color_conv <- c("Successful completion" = "#00BFC4",
                  "iteration limit had been reached" = "#F8766D",
                  "Degeneracy of the Nelder-Mead simplex" = "#C77CFF")
  
  return(list("Results" = estim_tab,
              "Fig" = ggplot(data = estim_tab, aes(x = Repeat, y = Value)) +
                geom_line(cex = 0.3) +
                geom_point(aes(col = Convergence), cex = 2, pch = 19) +
                theme_bw() +
                labs(x="Optim() repeat", y="WLS value\n") +
                ggtitle("Results (value and convergence features) returned by function optim()") +
                scale_color_manual(values = color_conv[estim_tab$Convergence %>% unique %>%
                                                         match(names(color_conv))]) +
                theme(axis.text.x = element_text(size=11),
                      axis.text.y = element_text(size=11),
                      legend.position = "bottom")))
}

## Wild Bootstrap

### Generating random weights W_i for residual disturbance with E[W_i] = 0 and E[W_i²] = 1

random.disturbance <- function(size, distrib, seed){
  
  # distrib = distribution of the random disturbance {W_i}_(i=1...n) with E[W_i] = 0 and E[W_i²] = 1:
  #           - "Standard Normal" (default)
  #           - "Rademacher": random variable that puts probability 0.5 on the values 1 and -1
  #           - "Mammen2" (Mammen's two-point distribution): (1-sqrt(5))/2 with probability (sqrt(5)+1)/(2*sqrt(5))
  #                                                          (1+sqrt(5))/2 with probability (sqrt(5)-1)/(2*sqrt(5))
  # seed = an integer to set the seed for Random Number Generation (purpose of reproducibility).
  
  set.seed(seed)
  
  if(distrib %in% c("Standard Normal", "standard normal", "Standard normal", "standard Normal")){
    W <- rnorm(size)
    
  }else if(distrib %in% c("Rademacher", "rademacher")){
    W <- ifelse(runif(size) < 0.5, 1, -1)
    
  }else if(distrib %in% c("Mammen2", "mammen2")){
    W <- sample(c((1-sqrt(5))/2,(1+sqrt(5))/2), size, replace = TRUE,
                prob = c((sqrt(5)+1)/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5))))
  }else{
    return(cat(crayon::red("ERROR: argument 'distrib' must be 'Standard Normal' (default), 'Rademacher' or 'Mammen2'")))
  }
  return(W)
}

### Wild bootstrap procedure

Bootstrap.wild <- function(data_obs, simul_obs, nb, transf_parm_optim, constants = NULL, times, N,
                           approx_init_I = FALSE, deaths_approx_init_IS = NULL, model = ODE_phase1,
                           fn = WLS, stringency_index = NULL, solver = "ode45",
                           method = "Nelder-Mead", maxit = 2000, abstol = 1e-6, reltol = 1e-6,
                           random.distrib = "Standard Normal", seed = 123){
  
  # see: Patrick Kline & Andres Santos. A Score Based Approach to Wild Bootstrap Inference
  #      Journal of Econometric Methods 2012; 1(1): 23–4. (DOI 10.1515/2156-6674.10)
  # https://eml.berkeley.edu/~pkline/papers/ScoreFinal_web.pdf
  # 2. Wild Bootstrap Review (p.4)
  
  # data_obs = real observation data (/!\ order: daily negative tests, positive tests and fatality cases)
  # simul_obs = simulated observation data
  # nb = total number of bootstraps to perform
  # random.distrib = distribution of the random disturbance {W_i}_(i=1...n) with E[W_i] = 0 and E[W_i²] = 1:
  #                  - "Standard Normal" (default)
  #                  - "Rademacher": 1 or -1, both with probability 0.5
  #                  - "Mammen2" (Mammen's two-point distribution): (1-sqrt(5))/2 with probability (sqrt(5)+1)/(2*sqrt(5))
  #                                                                 (1+sqrt(5))/2 with probability (sqrt(5)-1)/(2*sqrt(5))
  # seed = an integer to set the seed for Random Number Generation (purpose of reproducibility) ; 123 by default.
  
  n_starts <- dim(transf_parm_optim)[1] # Number of set(s) of starting values
  n_parms <- dim(transf_parm_optim)[2] # Number of parameters
  n_time <- dim(data_obs)[1] # number of time points
  n_obs <- dim(data_obs)[2] # number of observable states
  n_points <- n_time*n_obs # total number of points = number of time points x number of observable sates
  
  if(is.null(stringency_index)){
    cat(crayon::cyan("MESSAGE: no Stringency Index detected \n"))
  }
  
  Residuals <- data_obs-simul_obs
  W <- random.disturbance(n_points*nb, random.distrib, seed) %>% matrix(ncol = n_points)
  
  estim_tab <- data.frame()
  
  for(b in 1:nb){
    
    bootstraped_values <- simul_obs + Residuals*W[b,]
    
    min_val <- +Inf
    
    for(i in 1:n_starts){
      par <- transf_parm_optim[i,]
      fit <- optim(par = par, fn = fn, constants = constants, model = model,
                   approx_init_I = approx_init_I, deaths_approx_init_IS = deaths_approx_init_IS,
                   stringency_index = stringency_index, times = times, N = N, data_obs = bootstraped_values,
                   solver = solver, method = method,
                   control = list(maxit = maxit, abstol = abstol, reltol = reltol))
      
      if(fit$value < min_val & fit$convergence == 0){
        best_fit <- fit
        min_val <- best_fit$value
      }
    }
    if(min_val != +Inf){
      
      estim_tab <- rbind(estim_tab, c(b, best_fit$value, "Successful completion", bootstraped_values[1,3],
                                      best_fit$par[names(best_fit$par) == "pS"] %>% inv.logit,
                                      best_fit$par[names(best_fit$par) == "pE"] %>% inv.log_ab(a=0,b=0.01),
                                      best_fit$par[names(best_fit$par) == "pIA"] %>% inv.log_ab(a=0,b=0.01),
                                      best_fit$par[names(best_fit$par) == "pISr"] %>% inv.log_ab(a=0,b=0.01),
                                      best_fit$par[names(best_fit$par) == "pISd"] %>% inv.log_ab(a=0,b=0.01),
                                      best_fit$par[names(best_fit$par) == "beta"] %>% exp,
                                      best_fit$par[names(best_fit$par) == "gamma"] %>% exp,
                                      best_fit$par[names(best_fit$par) == "alpha"] %>% inv.logit,
                                      best_fit$par[names(best_fit$par) == "kappa"] %>% exp,
                                      best_fit$par[names(best_fit$par) == "k"] %>% inv.logit,
                                      best_fit$par[names(best_fit$par) == "a"] %>% exp,
                                      best_fit$par[names(best_fit$par) == "omega"] %>% inv.logit,
                                      best_fit$par[names(best_fit$par) == "p"] %>% inv.log_ab(a=0,b=0.2),
                                      best_fit$par[names(best_fit$par) == "rho"] %>% inv.logit,
                                      best_fit$par[names(best_fit$par) == "eta"] %>% inv.logit,
                                      best_fit$par[names(best_fit$par) == "mu"] %>% inv.logit))
    }else{
      estim_tab <- rbind(estim_tab, c(b, NA, "No convergence", rep(NA,n_parms)))
    }
    print(paste("Progress:", b, "bootstrap(s) out of", nb, "completed"))
  }
  names <- c("pS", "pE", "pIA", "pISr", "pISd", "beta", "gamma", "alpha", "kappa", "k", "a",
             "omega", "p", "rho", "eta", "mu")
  par_names <- names(transf_parm_optim)
  colnames(estim_tab) <- c("Repeat", "Value", "Convergence", "Deaths_t1",
                           names[par_names %>% match(names) %>% sort])
  return(estim_tab)
}

# Graphical vizualisations

## Plotting Stringency Index and outbreak control

plot_control <- function(t0, tf, stringency_index, k, a){
  
  # t0 and tf represent respectively the first and the last time point.
  
  return(ggplot(data.frame("Time" = t0:tf, "stringency_Index" = stringency_index,
                           "control_efficiency" = 100*k*(stringency_index/100)^a),
                aes(x = Time)) +
           geom_point(aes(y = stringency_Index), cex = 1, pch = 19, col = "darkred") +
           geom_point(aes(y = control_efficiency), cex = 1, pch = 19, col = "#1c7ef6") +
           scale_y_continuous(name = "Stringency Index \n", limits = c(0,100),
                              sec.axis = sec_axis(~./100, name="Control efficiency \n")) +
           theme_bw() +
           theme(axis.text.x = element_text(size=11),
                 axis.text.y = element_text(size=11),
                 axis.title.y = element_text(color = "darkred"),
                 axis.title.y.right = element_text(color = "#1c7ef6")) +
           geom_vline(xintercept=t, linetype="dashed", col = "#5b5b5b"))
}

## Plotting estimates against WLS values for each fit

plot_estimates_vs_fit_value <- function(parm_name, tab, group = FALSE,
                                        log_scale_x = FALSE, log_scale_y = FALSE, cex = 1){
  
  # parm_name: name of the parameter
  # tab: a data frame with parameter estimates (1st column) and min. WLS value (2nd column)
  # group: TRUE if we want to color estimates differently according to the group to which they belong
  #                                         (optional third column 'group' in tab) ; FALSE by default.
  
  Fig <- ggplot(tab, aes(x = tab[,1], y = tab[,2])) +
    theme_bw() +
    labs(x = parm_name, y = "Minimum WLS value") +
    theme(axis.text.x = element_text(size=11),
          axis.text.y = element_text(size=11, angle = 45))
  
  if(group == TRUE){
    Fig <- Fig + geom_point(aes(col = group), cex = cex, pch = 19)
  }else{
    Fig <- Fig + geom_point(cex = cex, pch = 19)
  }
  if(log_scale_x){
    Fig <- Fig + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x)))
  }
  if(log_scale_y){
    Fig <- Fig + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x)))
  }
  return(Fig)
}

## Comparing simulation using optimized parameters with observed data

simul_and_compare <- function(estimates, real_values = NULL, constants = NULL, times = NULL, dates = NULL, N,
                              data_obs, stringency_index = NULL, approx_init_I = FALSE,
                              deaths_approx_init_IS = NULL, model = ODE_phase1, solver = "ode45"){
  if(is.null(dates)){
    dates <- times
  }
  if(is.null(times)){
    times <- 0:(dates[length(dates)]-dates[1])
  }
  parms <- c(estimates, constants)

  if(!is.null(stringency_index)){
    parms$stringency_index <- stringency_index
  }else{
    cat(crayon::cyan("MESSAGE: no stringency Index detected \n"))
  }
  
  init_estim <- c(parms[["pS"]], parms[["pE"]])*N
  if(approx_init_I){
    if(is.null(deaths_approx_init_IS)){
      deaths_approx_init_IS <- data_obs[1,3]
    }
    ISd_0 <- deaths_approx_init_IS/parms[["alpha"]]
    ISr_0 <- ISd_0*(1-parms[["p"]])/parms[["p"]]
    init_estim <- c(init_estim, (ISr_0+ISd_0)*(1-parms[["omega"]])/parms[["omega"]], ISr_0, ISd_0)
  }else{
    init_estim <- c(init_estim, parms[["pIA"]], parms[["pISr"]], parms[["pISd"]])*N
  }
  init_estim <- c(init_estim, N-sum(init_estim), 0, 0)
  
  simul_estim <- ode(y = init_estim, times = times, parms = parms, func = model, method = solver)
  simul_estim_obs <- simul_estim[times >= 0,]
  d <- dim(simul_estim_obs)[1]
  
  if(any(names(parms) == "rho")){
    rho <- parms[["rho"]]
  }else{
    rho <- parms[["mu"]] + parms[["eta"]]*times[times > 0]
  }
  T_neg_estim <- rho*simul_estim_obs[-1,2]
  T_pos_estim <- (simul_estim_obs[-1,9]-simul_estim_obs[-d,9]) + rho*simul_estim_obs[-1,4]
  Deaths_estim <- simul_estim_obs[-1,8]-simul_estim_obs[-d,8]
  
  Compare_data <- data.frame("Time" = dates[times > 0] %>% rep(2),
                             "T_neg" = c(data_obs[,1], T_neg_estim),
                             "T_pos" = c(data_obs[,2], T_pos_estim),
                             "Deaths" = c(data_obs[,3], Deaths_estim),
                             "Legend" = c(rep("Observed data", d-1),
                                          rep("With estimates", d-1)))
  
  Residuals <- data.frame("Time" = dates[times > 0],
                          "Residuals_T_neg" = c(data_obs[,1]-T_neg_estim),
                          "Residuals_T_pos" = c(data_obs[,2]-T_pos_estim),
                          "Residuals_Deaths" = c(data_obs[,3]-Deaths_estim))
  
  Fig_neg_tests <- ggplot(Compare_data, aes(x = Time, y = T_neg, col = Legend)) +
    geom_line(data = Compare_data[Compare_data$Legend == "With estimates",], cex = 1.2) +
    geom_point(data = Compare_data[Compare_data$Legend == "Observed data",], cex = 1.2, pch = 19) +
    ylab("Daily negative tests") +
    theme_bw() +
    scale_color_manual(values = c("#00BFC4", "#58ebef")) +
    theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11, angle = 45),
          legend.position = "right") +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE))
  
  Fig_resid_neg_tests <- ggplot(data = Residuals, aes(x = Time, y = Residuals_T_neg)) +
    geom_point(pch=1, cex = 1, col = "darkblue")+
    theme_bw() +
    ylab("Raw residuals") +
    theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11, angle = 45)) +
    geom_hline(yintercept = 0, lwd = 0.5) +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE))  

  Fig_pos_tests <- ggplot(Compare_data, aes(x = Time, y = T_pos, col = Legend)) +
    geom_line(data = Compare_data[Compare_data$Legend == "With estimates",], cex = 1.2) +
    geom_point(data = Compare_data[Compare_data$Legend == "Observed data",], cex = 1.2, pch = 19) +
    ylab("Daily positive tests") +
    theme_bw() +
    scale_color_manual(values = c("#F8766D", "#f584ea")) +
    theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11, angle = 45),
          legend.position = "right") +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE))
  
  Fig_resid_pos_tests <- ggplot(data = Residuals, aes(x = Time, y = Residuals_T_pos)) +
    geom_point(pch=1, cex = 1, col = "darkred")+
    theme_bw() +
    ylab("Raw residuals") +
    theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11, angle = 45)) +
    geom_hline(yintercept = 0, lwd = 0.5) +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE))
  
  Fig_deaths <- ggplot(Compare_data, aes(x = Time, y = Deaths, col = Legend)) +
    geom_line(data = Compare_data[Compare_data$Legend == "With estimates",], cex = 1.2) +
    geom_point(data = Compare_data[Compare_data$Legend == "Observed data",], cex = 1.2, pch = 19) +
    ylab("Daily fatality cases") +
    scale_color_manual(values = c("gray30", "gray60")) +
    theme_bw() +
    theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11),
          legend.position = "right")

  Fig_resid_deaths <- ggplot(data = Residuals, aes(x = Time, y = Residuals_Deaths)) +
    geom_point(pch=1, cex = 1)+
    theme_bw() +
    ylab("Raw residuals") +
    theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) +
    geom_hline(yintercept = 0, lwd = 0.5)

  if(!is.null(real_values)){
    Fig_control <- ggplot(data.frame("Time" = dates,
                                     "control" = c(real_values[["k"]]*(stringency_index/100)^real_values[["a"]],
                                                   parms[["k"]]*(stringency_index/100)^parms[["a"]]),
                                     "Legend" = c(rep("With real parameters", length(stringency_index)),
                                                  rep("With estimates", length(stringency_index)))),
                          aes(x = Time, y = control, color = Legend)) +
      geom_point(cex = 1, pch = 19) +
      scale_color_manual(values = c("#35d8f5", "#1c7ef6")) +
      ylab("Control efficiency \n") + ylim(c(0,1)) +
      theme_bw() + theme(axis.text.x = element_text(size=11), axis.text.y = element_text(size=11))
    
    return(list("simul" = simul_estim,
                "data_estim" = cbind("T_neg" = T_neg_estim, "T_pos" = T_pos_estim, "Deaths" = Deaths_estim),
                "Fig_positive_tests" = Fig_pos_tests, "Fig_resid_positive_tests" = Fig_resid_pos_tests,
                "Fig_negative_tests" = Fig_neg_tests, "Fig_resid_negative_tests" = Fig_resid_neg_tests,
                "Fig_deaths" = Fig_deaths, "Fig_resid_deaths" = Fig_resid_deaths,
                "Fig_control" = Fig_control))
  }else{
    return(list("simul" = simul_estim,
                "data_estim" = cbind("T_neg" = T_neg_estim, "T_pos" = T_pos_estim, "Deaths" = Deaths_estim),
                "Fig_positive_tests" = Fig_pos_tests, "Fig_resid_positive_tests" = Fig_resid_pos_tests,
                "Fig_negative_tests" = Fig_neg_tests, "Fig_resid_negative_tests" = Fig_resid_neg_tests,
                "Fig_deaths" = Fig_deaths, "Fig_resid_deaths" = Fig_resid_deaths))
  }
}

## Plot (1-alpha) % of the distribution (density) of a parameter

plot_density <- function(name, data, col = "red", alpha = 0.05){
  
  # name: name of the parameter (-> name of the x-axis)
  # data: vector containing all the values of the parameter
  # col: color of the curb (red by default)
  # alpha: order for defining quantiles (default: 0.05 / 5%)
  
  if(alpha < 0 | alpha >= 1){
    return(cat(crayon::red("ERROR: 'alpha' must belong to the interval [0 ; 1 [")))
  }
  data <- data[!is.na(data)]
  
  return(ggplot(as.data.frame(data) %>% setNames("parm"), aes(x = parm)) +
           geom_density(cex = 1, col = col, fill = col, alpha = 0.1) +
           xlim(quantile(data, probs = c(alpha/2,1-alpha/2))) +
           xlab(name) +
           theme_bw() +
           ylab(paste0("Density (", 100*(1-alpha), "%) \n")) +
           theme(axis.text.x = element_text(size=11),
                 axis.text.y = element_text(size=11)))
}

## Compute simulated observations

simul_obs_data <- function(estimates, constants = NULL, times, N,
                           stringency_index = NULL, approx_init_I = FALSE,
                           deaths_approx_init_IS = NULL, model = ODE_phase1, solver = "ode45"){
  
  parms <- c(estimates, constants)
  
  if(!is.null(stringency_index)){
    parms$stringency_index <- stringency_index
  }else{
    cat(crayon::cyan("MESSAGE: no stringency Index detected \n"))
  }
  
  init_simul <- c(parms[["pS"]], parms[["pE"]])*N
  if(approx_init_I){
    if(is.null(deaths_approx_init_IS)){
      return(cat(crayon::red("ERROR: the first mortality data is needed if you want to approximate the initial state of I_Sd \n")))
    }
    ISd_0 <- deaths_approx_init_IS/parms[["alpha"]]
    ISr_0 <- ISd_0*(1-parms[["p"]])/parms[["p"]]
    init_simul <- c(init_simul, (ISr_0+ISd_0)*(1-parms[["omega"]])/parms[["omega"]], ISr_0, ISd_0)
  }else{
    init_simul <- c(init_simul, parms[["pIA"]], parms[["pISr"]], parms[["pISd"]])*N
  }
  init_simul <- c(init_simul, N-sum(init_simul), 0, 0)
  
  simul <- ode(y = init_simul, times = times, parms = parms, func = model, method = solver)
  simul_obs <- simul[times >= 0,]
  d <- dim(simul_obs)[1]
  
  if(any(names(parms) == "rho")){
    rho <- parms[["rho"]]
  }else{
    rho <- parms[["mu"]] + parms[["eta"]]*times[times > 0]
  }
  T_neg_estim <- rho*simul_obs[-1,2]
  T_pos_estim <- (simul_obs[-1,9]-simul_obs[-d,9]) + rho*simul_obs[-1,4]
  Deaths_estim <- simul_obs[-1,8]-simul_obs[-d,8]
  
  return(cbind("T_neg" = T_neg_estim, "T_pos" =  T_pos_estim, "Deaths" = Deaths_estim))
}

# Identifiability analysis

## Compute the identifiability profile of a parameter of interest

Identifiability_profile <- function(parm_name, x, parm_optim, constants = NULL, model = ODE_phase1, fn = WLS,
                                    solver = "ode45", approx_init_I = FALSE, deaths_approx_init_IS = NULL,
                                    times, N, stringency_index = NULL, data_obs, method = "Nelder-Mead",
                                    maxit = 2000, abstol = 1e-6, reltol = 1e-6){

  # parm_name = name of the parameter for which the profile is built
  # x = x-axis values of the profile (sequential values to be set for the parameter for which the profile is built)
  # parm_optim = vector with initial values of parameters to optimized
  # fn = function to minimize
  # data = matrix of observed data. Columns: daily (1) negative tests, (2) positive tests, (3) fatality cases
  # ..., the other arguments are identical to those used in previous functions

  if(is.null(stringency_index)){
    cat(crayon::cyan("MESSAGE: no stringency Index detected \n"))
  }
  
  transf_parm_optim <- c(parm_optim[names(parm_optim) == "pS"] %>% logit,
                         parm_optim[names(parm_optim) == "pE"] %>% log_ab(a=0,b=0.01),
                         parm_optim[names(parm_optim) == "pIA"] %>% log_ab(a=0,b=0.01),
                         parm_optim[names(parm_optim) == "pISr"] %>% log_ab(a=0,b=0.01),
                         parm_optim[names(parm_optim) == "pISd"] %>% log_ab(a=0,b=0.01),
                         parm_optim[names(parm_optim) == "beta"] %>% log,
                         parm_optim[names(parm_optim) == "gamma"] %>% log,
                         parm_optim[names(parm_optim) == "alpha"] %>% logit,
                         parm_optim[names(parm_optim) == "kappa"] %>% log,
                         parm_optim[names(parm_optim) == "k"] %>% logit,
                         parm_optim[names(parm_optim) == "a"] %>% log,
                         parm_optim[names(parm_optim) == "omega"] %>% logit,
                         parm_optim[names(parm_optim) == "p"] %>% log_ab(a=0,b=0.2),
                         parm_optim[names(parm_optim) == "rho"] %>% logit,
                         parm_optim[names(parm_optim) == "eta"] %>% logit,
                         parm_optim[names(parm_optim) == "mu"] %>% logit)
  len <- length(x)
  values <- c()

  for(i in 1:len){

    csts <- c(constants, x[i])
    names(csts)[length(csts)] <- parm_name
    
    fit <- optim(par = transf_parm_optim, fn = fn, constants = csts, model = model,
                 approx_init_I = approx_init_I, deaths_approx_init_IS = deaths_approx_init_IS,
                 stringency_index = stringency_index, times = times, N = N, data_obs = data_obs,
                 solver = solver, method = method,
                 control = list(maxit = maxit, abstol = abstol, reltol = reltol))

    values <- values %>% c(ifelse(fit$convergence == 0, yes = fit$value, no = NA))

    print(paste("Progress:", round(100*i/len, digits = 2), "%"))
  }
  return(data.frame("parm" = x, "WLS" = values))
}


Identifiability_profile_parallel <- function(parm_name, x, parm_optim, constants = NULL, model = ODE_phase1, fn = WLS,
                                             solver = "ode45", approx_init_I = FALSE, deaths_approx_init_IS = NULL,
                                             times, N, stringency_index = NULL, data_obs, method = "Nelder-Mead",
                                             maxit = 2000, abstol = 1e-6, reltol = 1e-6, Ncores = NA){
  
  # parm_name = name of the parameter for which the profile is built
  # x = x-axis values of the profile (sequential values to be set for the parameter for which the profile is built)
  # parm_optim = vector with initial values of parameters to optimized
  # fn = function to minimize
  # data = matrix of observed data. Columns: daily (1) negative tests, (2) positive tests, (3) fatality cases
  # ..., the other arguments are identical to those used in previous functions
  
  if(is.null(stringency_index)){
    cat(crayon::cyan("MESSAGE: no stringency Index detected \n"))
  }
  
  transf_parm_optim <- c(parm_optim[names(parm_optim) == "pS"] %>% logit,
                         parm_optim[names(parm_optim) == "pE"] %>% log_ab(a=0,b=0.01),
                         parm_optim[names(parm_optim) == "pIA"] %>% log_ab(a=0,b=0.01),
                         parm_optim[names(parm_optim) == "pISr"] %>% log_ab(a=0,b=0.01),
                         parm_optim[names(parm_optim) == "pISd"] %>% log_ab(a=0,b=0.01),
                         parm_optim[names(parm_optim) == "beta"] %>% log,
                         parm_optim[names(parm_optim) == "gamma"] %>% log,
                         parm_optim[names(parm_optim) == "alpha"] %>% logit,
                         parm_optim[names(parm_optim) == "kappa"] %>% log,
                         parm_optim[names(parm_optim) == "k"] %>% logit,
                         parm_optim[names(parm_optim) == "a"] %>% log,
                         parm_optim[names(parm_optim) == "omega"] %>% logit,
                         parm_optim[names(parm_optim) == "p"] %>% log_ab(a=0,b=0.2),
                         parm_optim[names(parm_optim) == "rho"] %>% logit,
                         parm_optim[names(parm_optim) == "eta"] %>% logit,
                         parm_optim[names(parm_optim) == "mu"] %>% logit)
  len <- length(x)
  values <- c()
  
  # Parallelized optimizations
  
  if(is.na(Ncores)){Ncores <- bigparallelr::nb_cores()}
  print(paste("Number of cores:", Ncores))
  
  t <- format(Sys.time(), "%Y-%m-%d %Hh%M")
  folder_name <- paste0("ALPHA_Identifiability_optim_parallel_",
                        stringi::stri_replace_all_fixed(t, c("-", " "), "_", vectorize_all = F))
  dir.create(folder_name)
  lapply(list("PROJECT ALPHA\n", t, paste0("\nParallel - number of cores: ", Ncores),
              paste("\nODEs solver:", solver),
              paste("\nOptimization algorithm: optim() with method", method),
              "\nControl:", paste(c("Max iterations", "Absolute tolerance", "Relative tolerance"),
                                  c(maxit, abstol, reltol), sep = " = ", collapse = "\n"),
              paste("\nParameter of interest:", parm_name, "from", min(x), "to", max(x)),
              paste0("\nParameters to estimate (", length(parm_optim), "):"), names(parm_optim) %>% paste(collapse = ", "),
              "\nConstants:", paste(names(constants), constants, sep = " = ", collapse = "\n")),
         cat, "\n", file=paste0(folder_name,"/README.txt"), append=TRUE)
  
  cl <- parallel::makeCluster(Ncores)
  doParallel::registerDoParallel(cl)
  
  profile <- foreach::foreach(i = 1:n_starts, .combine = 'rbind', .inorder = FALSE,
                              .packages = c("tidyverse", "deSolve", "boot"), .export = ls(globalenv())) %dopar% {
                                
                                csts <- c(constants, x[i])
                                names(csts)[length(csts)] <- parm_name
                                
                                fit <- optim(par = transf_parm_optim, fn = fn, constants = csts, model = model,
                                             approx_init_I = approx_init_I, deaths_approx_init_IS = deaths_approx_init_IS,
                                             stringency_index = stringency_index, times = times, N = N, data_obs = data_obs,
                                             solver = solver, method = method,
                                             control = list(maxit = maxit, abstol = abstol, reltol = reltol))
                                
                                saveRDS(res, file = paste0(folder_name, "/result_task_", i, ".rds"))
                                
                                return(c("parm" = x[i], "WLS" = fit$value, "convergence" = fit$convergence))
                                } %>% as.data.frame
  parallel::stopCluster(cl)
  return(profile[order(profile$parm),])
}

## Plot the identifiability profile

plot_profile <- function(parm_name, profile, real_value = NULL, log = FALSE){
  
  # if real_value is not 'NULL', a vertical dashed line indicates this value on the x-axis
  # log: if TRUE, the log WLS value is plotted on the y-axis
  
  min <- profile[which(profile$WLS == min(profile$WLS, na.rm = TRUE)),]
  
  if(log == TRUE){
    Fig <- ggplot(data = profile, aes(x = parm, y = log(WLS))) +
      ylab("WLS value (log scale)\n")
  }else{
    Fig <- ggplot(data = profile, aes(x = parm, y = WLS)) +
      ylab("WLS value\n")
  }
  if(!is.null(real_value)){
    Fig <- Fig + geom_vline(xintercept = real_value, lty = 'dashed')
  }
  
  return(Fig + geom_line(cex = 0.7, color = "#29abd9") +
           geom_point(cex = 1.3, color = "#29abd9") +
           geom_point(data = min, pch = 1, cex = 7) +
           xlab(parm_name) +
           theme_classic() +
           theme(axis.title.x = element_text(size=13),
                 axis.title.y = element_text(size=11),
                 axis.text.x = element_text(size=10),
                 axis.text.y = element_text(size=10)))
}

# Miscellaneous

## For graphical purpose: Add a '+' before positive numbers ('-' before negative numbers, nothing before 0)

add_plus_sign <- function(x, digits = 2){
  if(length(x) > 1){
    return(sapply(x, add_plus_sign, digits = digits))
  }
  if(isTRUE(all.equal(x,0))){
    return(sprintf(paste0("% .",digits,"f"),x))
  }else{
    return(sprintf(paste0("%+.",digits,"f"),x))
  }
}
