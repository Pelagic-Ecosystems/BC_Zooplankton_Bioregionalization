# Create function that converts chl-a to phytoplankton size classes and functional types
# Enter a vector of chlorophyll-a values and return dataframes for the PSC and PFT for the global (Hirata et al. 2011)
# and the North Pacific (Zheng et al. 2018) versions.

# The Hirata et al. 2011 analysis showed that total chlorophyll-a can be used to derive the percentage of the
# phytoplankton community composed into size classes (Pico, Nano, Micro) and seven functional types (Picoeukaryotes,
# prokaryotes, Prochlorococcus, prymnesiophytes/haptophytes, chlorophytes, dinoflagellates, diatoms). The functional types
# though are interpretations of the diagnostic pigments and were not validated with new microscopy/cytometry samples. In 
# Zeng et al. 2018, the equations in Hirata et al. 2011 (also an update of the original Uitz et al. 2006 paper) were 
# validated for the northeast Pacific region with HPLC data and microscopy, confiriming the PFTs too. The Zeng et al. 2018
# paper also updated the coefficients of the Hirata et al. 2011 equations to better fit the NEP pigment dataset. This 
# improved the model for all RMSE and all r2 except Micro/Diatoms. Zeng et al. 2018 tested the improved model on satellite
# chl-a and calculated an error of 30% for estimated PSC/PFTs relative to the underway estimates.
# Using these equations with satellite chl-a is reasonable but interpret with caution.

# The calculates the fraction of the chla belonging to a group but not the chla of each group. At the end of this code,
# the fractions are also multiplied with chla so the values returned as PSC.chla & PFT.chla are the chla of each group.
# Mode is 1 = global, 2 = NEP
# Hirata et al. 2011 notes to set PFT fraction to 1 if >1 and 0 if <0.
# There is an error somewhere in PFT that the total of the row sums is not 1, maybe the confidence interval error?

# P.Pata - Sep 28, 2020

p_chla_SCFT <- function(chla, mode = 2) {
  x <- log10(chla)
  y <- chla
  z <- mode
  n <- length(x)
  
  # Phytoplankton size classes
  PSC <- data.frame( micro = numeric(n), nano = numeric(n) , pico = numeric(n)  )
  
  a0 <- c(0.9117, 1.172)
  a1 <- c(-2.7330,-2.423)
  a2 <- c(0.4003,-0.035)
  PSC$micro <- ( a0[z] + exp(a1[z]*x + a2[z]) ) ^ -1
  PSC$micro[PSC$micro < 0] <- 0
  PSC$micro[PSC$micro > 1] <- 1
  
  a0 <- c(0.1529, -0.061)
  a1 <- c(1.0306,-0.112)
  a2 <- c(-1.5576,-1.293)
  a3 <- c(-1.8597,0.503)
  a4 <- c(2.9954,4.947)
  PSC$pico <- -( a0[z] + exp(a1[z]*x + a2[z]) ) ^ -1 + a3[z]*x + a4[z]
  PSC$pico[PSC$pico < 0] <- 0
  PSC$pico[PSC$pico > 1] <- 1
  
  PSC$nano <- 1 - PSC$micro - PSC$pico
  PSC$nano[PSC$nano < 0] <- 0
  PSC$nano[PSC$nano > 1] <- 1
  
  # Phytoplankton functional types
  PFT <- data.frame( diatom = numeric(n), dinoflagellate  = numeric(n), chlorophyte = numeric(n), 
                     prymnesiophyte = numeric(n), prokaryote = numeric(n), eukaryote = numeric(n))
  a0 <- c(1.3272,1.317)
  a1 <- c(-3.9828,-2.388)
  a2 <- c(0.1953,0.129)
  PFT$diatom <- ( a0[z] + exp(a1[z]*x + a2[z]) ) ^ -1
  PFT$diatom[PFT$diatom < 0] <- 0
  PFT$diatom[PFT$diatom > 1] <- 1
  
  PFT$dinoflagellate <- PSC$micro - PFT$diatom
  PFT$dinoflagellate[PFT$dinoflagellate < 0] <- 0
  PFT$dinoflagellate[PFT$dinoflagellate > 1] <- 1
  
  a0 <- c(0.2490,1.112)
  a1 <- c(-1.2621,-0.418)
  a2 <- c(-0.5523,-2.290)
  PFT$chlorophyte <- (a0[z]/y) * exp(a1[z]*(x + a2[z]) ^ 2)
  PFT$chlorophyte[PFT$chlorophyte < 0] <- 0
  PFT$chlorophyte[PFT$chlorophyte > 1] <- 1
  
  PFT$prymnesiophyte <- PSC$nano - PFT$chlorophyte
  PFT$prymnesiophyte[PFT$prymnesiophyte < 0] <- 0
  PFT$prymnesiophyte[PFT$prymnesiophyte > 1] <- 1
  
  a0 <- c(0.0067,-0.253)
  a1 <- c(0.6154,-8.154)
  a2 <- c(-19.5190,-17.3)
  a3 <- c(0.9643,1.214)
  a4 <- c(0.1027,-0.02)
  a5 <- c(-0.1189,0.004)
  a6 <- c(0.0626,0.039)
  PFT$prokaryote <- (a0[z]/a1[z]/y) * exp ( (a2[z]*(x + a3[z])^2) / (a1[z]^2)) + a4[z]*x^2 + a5[z]*x + a6[z]
  PFT$prokaryote[PFT$prokaryote < 0] <- 0
  PFT$prokaryote[PFT$prokaryote > 1] <- 1
  
  PFT$eukaryote <- PSC$pico - PFT$prokaryote
  PFT$eukaryote[PFT$eukaryote < 0] <- 0
  PFT$eukaryote[PFT$eukaryote > 1] <- 1
  
  PSC.chla <- PSC * chla
  PFT.chla <- PFT * chla
  return(list(chla = chla, PSC = PSC,PFT = PFT,
              PSC.chla = PSC.chla, PFT.chla = PFT.chla))
  
}