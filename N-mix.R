#------ Function to simulate the data

# low_high
# alpha.lam = 1
# p = 0.8

# low_low
# alpha.lam = 1
# p = 0.2

# high_high
# alpha.lam = 10
# p = 0.8

# high_low
# alpha.lam = 10
# p = 0.2


data.fn <- function(# Sampling desing parameters
                    R = 300, # Number of sites
                    T = 8,   # Number of replicate surveys per site
                    J = 3,   # Number of qPCR runs per sample collected on a single host
                    K = 10,  # Number of seasons
                    
                    # First value = Uninfected; Second value= infected
                    alpha.lam = c(10, 10),   # Average number of hosts per site the first season
                    phi = c(0.8, 0.8),       # Apparent host survival probability of uninfected and infected hosts
                    gamma = c(2, 2),         # Average number of individuals arriving at each site for uninfected and infected hosts
                    # Average site-level infection intensity
                    infec = 2,
                    
                    # Host detection probability during each survey
                    p = c(0.8, 0.5),  

                    # Imperfect pathogen detection as a function of pathogen load
                    alpha_qPCR = -0.25,    # Intercept on logit scale
                    beta_qPCR = 0.5,       # Slope of logit scale
                    
                    # Transition paobabilities
                    recovery = 0.1,  # Recovery probability
                    infection = 0.5  # Infection probability
                    ){
  
# Empty matrices to hold the data    
  yN <- yI <- array(NA, dim = c(R, T, K))   # Pathogen-corrected abundance data
  NI <- NN <- array(NA, dim = c(R, K))      # True abundance data
  
#-------- First season 
  NN[,1] <- rpois(n = R, lambda = alpha.lam[1])   # Uninfected
  NI[,1] <- rpois(n = R, lambda = alpha.lam[2])   # Infected

#------ Second season  
# Empty matrices to hold latent abundance variables 
# i.e., number of hosts surviving, arriving, and transitioning
  
  SN <- SI <- GI <- GN <- TrN <- TrI <- array(0, dim = c(R, K-1))  
  
  for(k in 2:K){  
    for(i in 1:R){
      
      if(NN[i,k-1]>0){
        SN[i, k-1] <- rbinom(n=1, size=NN[i,k-1], prob=phi[1])  
             # Survival of uninfected
        TrN[i,k-1] <- rbinom(n=1, size=SN[i,k-1], prob=infection)    
             # Transition from uninfected to infected
      }
      if(NI[i,k-1]>0){
        SI[i, k-1] <-  rbinom(n=1, size=NI[i,k-1], prob=phi[2])   
             # Survival of infected 
        TrI[i, k-1] <- rbinom(n=1, size=SI[i,k-1], prob=recovery) 
             # Tranisition from infected to uninfected
      }
      # Recruitment
      GN[i, k-1] <- rpois(1, lambda = gamma[1])
        # Number of individuals gained for uninfected
      GI[i, k-1] <- rpois(1, lambda = gamma[2])
        # Number of individuals gained for infected
      
    }
    
# Total number of individuals the next season given the processes of survival, recruitment, and transition
    NI[,k] <-  SI[,k-1] + GI[,k-1] + TrN[,k-1] - TrI[,k-1]
    NN[,k] <-  SN[,k-1] + GN[,k-1] + TrI[,k-1] - TrN[,k-1]
    
  }

# Assign the true infection intensity to each individual

II_true_NN <- array(NA, dim = c(R, max(NN), T, K))
II_true_NI <- array(NA, dim = c(R, max(NI), T, K))

###------ True uninfected individuals
for(i in 1:R){
  for(l in 1:T){
    for(k in 1:K){
      for(j in 1:NN[i,k]){
        
        II_true_NN[i, j, l, k] <- 0
        
      }   
    }
  }
}

###----- True infection intensity
for(i in 1:R){
  for(l in 1:T){
    for(k in 1:K){
      for(j in 1:NI[i,k]){
        
        II_true_NI[i, j, l, k] <- rlnorm(1, meanlog = infec, sd = 1)
        
      }   
    }
  }
}
    
# Obervation process: 2 parts
  # 1. Imperfect host detection at the site level
  # 2. Imperfect pathogen detection on the host

# 1. Individual detection probability depends on number of individuals at the site
# Correcting for imperfect host detection probability

for(i in 1:R){
  for(j in 1:T){
    for(k in 1:K){
      
      yN[i, j, k] <- rbinom(n = 1, NN[i,k], p)
      yI[i, j, k] <- rbinom(n = 1, NI[i,k], p)
      
    }
  }
}

# 2. Correct for individuals that were misclassified as uninfected, when the pathogen was truly present 
# Correcting for imperfect pathogen detection

# Try making it a 4-dimensional matrix
  # ID, site, replicate site survey, season
  cN <- array(NA, dim = c(R, max(yN), J, T, K))
  cI <- array(NA, dim = c(R, max(yI), J, T, K))

p_qPCR <- plogis(alpha_qPCR + beta_qPCR * II_true_NI)
  
for(i in 1:R){
  for(j in 1:T){
    for(k in 1:K){
      for(l in 1:max(NN[i,k])){
        
      false_negative[i,k] <-  rbinom(n = 1, NN[i,k], p_qPCR[i,k])
        
      cN[i, l, j, k] <- rbinom(n = 1, II[i,k], p)
      cI[i, l, j, k] <- rbinom(n = 1, II[i,k], p)
      
    }
  }
}
      
  return(list(R = R, T = T, K = K,
              alpha.lam= alpha.lam,
              phi = phi,
              gamma = gamma,
              infection = infection,
              recovery = recovery,
              alpha_qPCR = alpha_qPCR,
              beta_qPCR = beta_qPCR,
              SN = SN,
              SI = SI,
              GN = GN,
              GI = GI,
              TrI = TrI,
              TrN = TrN,
              NN = NN, 
              NI = NI,
              II = II,
              p = p, 
              yN = yN,
              yI = yI,
              cI = cI,
              cN = cN))
}

# --------------- Code model in BUGS language

setwd("/Users/Cici/GitHub/2010to2014ElCope/Simulations_Detection_pop/")
sink("model.txt")
cat("
model{

# Priors
#------- NOT Infected
alpha.lamN  ~ dnorm(0,0.01)
pN     ~ dnorm(0, 0.368)
gammaN ~ dnorm(0,0.01) 
phiN   ~ dnorm(0, 0.368)
psi_NI ~ dnorm(0, 0.368)

#------- Infected
alpha.lamI  ~ dnorm(0,0.01)
pI     ~ dnorm(0, 0.368)
gammaI ~ dnorm(0,0.01)  
phiI   ~ dnorm(0, 0.368)
psi_IN ~ dnorm(0, 0.368)

#------------ Ecological model
#---- First season

for(i in 1:R){
  #------ Not infected
    NN[i, 1] ~ dpois(lambdaN[i])
      log(lambdaN[i]) <- alpha.lamN

  #----- Infected
    NI[i, 1] ~ dpois(lambdaI[i])  
      log(lambdaI[i]) <- alpha.lamI
}

#------ All other seasons
for(k in 2:K){
  for(i in 1:R){

    #------- Nost Infected
    SN[i,k] ~ dbin(phiN1[i,k], NN[i,k-1])  # Total survivors
        logit(phiN1[i,k]) <- phiN
    TN[i,k] ~ dbin(psi_NI1[i,k], SN[i,k])  # Survive, become infected
        logit(psi_NI1[i,k]) <- psi_NI
    GN[i,k] ~ dpois(GaN[i, k])   # Recruits
      log(GaN[i, k]) <- gammaN

    #------- Infected
    SI[i,k] ~ dbin(phiI1[i,k], NI[i,k-1] ) # Infecteds who survive
        logit(phiI1[i,k]) <- phiI
    TI[i,k] ~ dbin(psi_IN1[i,k], SI[i,k] ) # Get better, transition to uninfected
        logit(psi_IN1[i,k]) <- psi_IN

    GI[i,k] ~ dpois(GaI[i, k])  # Recruits
      log(GaI[i, k]) <- gammaI

# Totals
  NN[i, k] <- SN[i,k] - TN[i,k] + GN[i,k] + TI[i,k]  
  NI[i, k] <- SI[i,k] - TI[i,k] + GI[i,k] + TN[i,k]  

  }

}

#------------- Obervation model

for(i in 1:R){
  for(j in 1:T){
    for(k in 1:K){

      yN[i, j, k] ~ dbin(pN1[i,j,k], NN[i, k]) 
                logit(pN1[i,j,k]) <- pN

      yI[i, j, k] ~ dbin(pI1[i,j,k], NI[i, k]) 
                logit(pI1[i,j,k]) <- pI

    }
  }
}

}
", fill = TRUE)
sink()

# Simulate the data
sodata <- data.fn()

# Bundle data
win.data <- list(yN = sodata$yN, 
                 yI = sodata$yI,
                 R = dim(sodata$yN)[1], 
                 T = dim(sodata$yN)[2],
                 K = dim(sodata$yN)[3])

# Initial values
NIst <- apply(sodata$yI, c(1, 3), max, na.rm = TRUE)
NNst <- apply(sodata$yN, c(1, 3), max, na.rm = TRUE)

inits <- function() {list( 
  alpha.lamN = (sodata$alpha.lam[1]),
  pN = plogis(sodata$p),
  phiN =  plogis(sodata$phi[1]),
  gammaN = (sodata$gamma[1]),
  psi_NI =  plogis(sodata$infection),
  
  alpha.lamI = (sodata$alpha.lam[2]),
  pI =  plogis(sodata$p),
  phiI =  plogis(sodata$phi[2]),
  gammaI = (sodata$gamma[2]),                        
  psi_IN =  plogis(sodata$recovery)
  
)}

# Monitor Parameters
params <- c("alpha.lamN", "alpha.lamI",
            "pN", "pI",
            "phiN", "phiI",
            "gammaN", "gammaI",
            "psi_NI", "psi_IN"
)

# MCMC settings
ni <- 10000
nb <- 3000
nt <- 10
nc <- 3


library("jagsUI")

high_high <- jags(win.data, inits, params, "model.txt", 
            n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb, parallel = TRUE)


#setwd("/Users/Cici/GitHub/2010to2014ElCope/Simulations_Detection_pop/data/")
#save(low_high, file = "low_high.rda")

