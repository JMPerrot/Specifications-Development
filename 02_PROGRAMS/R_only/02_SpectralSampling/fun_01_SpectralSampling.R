
SpecSampl <- function(){
  ################################################################################
  # input output directories
  ################################################################################
  PathData <- '../../01_DATA'
  PathResults <- '../../03_RESULTS'
  SpectralSampling_Dir <- file.path(PathResults,'02_SpectralSampling')
  dir.create(path = SpectralSampling_Dir,
             showWarnings = F,recursive = T)
  
  ################################################################################
  # load leaf optics dataset
  ################################################################################
  dbName <- 'ANGERS'
  PathLOPdb <- file.path(PathData,dbName,'LeafOptics.RData')
  load(PathLOPdb)
  lambda <- Reflectance$V1
  Reflectance <- Reflectance[,-1]
  Transmittance <- Transmittance[,-1]
  
  ################################################################################
  # Define parameters for inversion
  ################################################################################
  Parms2Estimate <- c('EWT_LMA','CHL_CAR')
  Parms2Estimate_ind <- list('EWT_LMA'=c('EWT','LMA'),'CHL_CAR'=c('CHL','CAR'))
  # define spectral sampling
  SpectralSampling <- list()
  SpectralSampling$CHL_CAR <- as.list(c(seq(2,50,by=1),seq(60,100,by=10)))
  SpectralSampling$EWT_LMA <- as.list(c(seq(2,50,by=1),seq(60,100,by=10),seq(120,400,by=20)))
  
  # define spectral domain
  SpectralDomain <- list()
  SpectralDomain$CHL_CAR <- list('minWL' = 400, 'maxWL' = 900)
  SpectralDomain$EWT_LMA <- list('minWL' = 1300, 'maxWL' = 2400)
  
  # define parameters to estimate during inversion
  ParmsEstInv <- list()
  ParmsEstInv$CHL_CAR <- c('CHL', 'CAR', 'EWT', 'LMA', 'N')
  ParmsEstInv$EWT_LMA <- c('EWT', 'LMA', 'N')
  
  # define parameters to estimate during inversion
  InitValues <- list()
  InitValues$CHL_CAR <- data.frame(CHL=45, CAR=8, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
  InitValues$EWT_LMA <- data.frame(CHL=45, CAR=8, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
  
  # init NRMSE data frame
  NRMSE_df<- data.frame("NRMSE" = NA)
  ################################################################################
  # Perform inversion
  ################################################################################
  # for each parameter
  nbWorkers <- 4
  registerDoFuture()
  plan(multisession, workers = nbWorkers)
  Estimate_SpectralSampling <- list()
  SpectralSampling_subDir <- file.path(SpectralSampling_Dir,'_individualSteps')
  dir.create(path = SpectralSampling_subDir,
             showWarnings = F,recursive = T)
  
  for (parm in Parms2Estimate){
    print(parm)
    # apply multiprocessing for each spectral sampling
    Invert_SpectralSampling <- function() {
      foreach(subsample = SpectralSampling[[parm]]) %dopar% {
        if (subsample<10){subChar <- paste('00',as.character(subsample),sep = '')}
        else if (subsample<100){subChar <- paste('0',as.character(subsample),sep = '')}
        else if (subsample<1000){subChar <- as.character(subsample)}
        FileName <- file.path(SpectralSampling_subDir,
                              paste(parm,'_SpecSampling_',subChar,'.RData',sep = ''))
        if (!file.exists(FileName)){
          lambda_tmp <- seq(from = SpectralDomain[[parm]]$minWL,
                            to = SpectralDomain[[parm]]$maxWL,
                            by = subsample)
          
          DataFit <- FitSpectralData(SpecPROSPECT = SpecPROSPECT,
                                     lambda = lambda,
                                     Refl = Reflectance,
                                     Tran = Transmittance,
                                     UserDomain = lambda_tmp,
                                     UL_Bounds = FALSE)
          
          Invert_est <- prospect::Invert_PROSPECT(SpecPROSPECT = DataFit$SpecPROSPECT,
                                                  Refl = DataFit$Refl,
                                                  Tran = DataFit$Tran,
                                                  PROSPECT_version = 'D',
                                                  Parms2Estimate = ParmsEstInv[[parm]],
                                                  InitValues = InitValues[[parm]])
          NRMSE<- get_performances_inversion(target =  Factor[[parm]]*Biochemistry[[parm]],
                                             estimate = Factor[[parm]]*Estimate_SpectralSampling[[parm]][[col0]],
                                             categories = T)[3]
          NRMSE_df<-rbind(NRMSE_df,NRMSE)
          save(Invert_est,file = FileName)
        } else {
          load(file = FileName)
        }
        return(Invert_est)
      }
    }
    PROSPECT_EST <- Invert_SpectralSampling()
    for (subparm in Parms2Estimate_ind[[parm]]){
      Estimate_SpectralSampling[[subparm]] <- do.call(cbind,lapply(PROSPECT_EST,'[[',subparm))
      Estimate_SpectralSampling[[subparm]] <- data.frame(Estimate_SpectralSampling[[subparm]])
      colnames(Estimate_SpectralSampling[[subparm]]) <- unlist(SpectralSampling[[parm]])
    }
  }
  NRMSE_df<-NRMSE_df[-c(1),]
  plan(sequential)
  
  # save all spectral sampling results in a unique file per parameter
  for (parm in names(Estimate_SpectralSampling)){
    FileName <- file.path(SpectralSampling_Dir,paste(parm,'_SpecSampling.csv',sep = ''))
    write_delim(x = Estimate_SpectralSampling[[parm]],
                file = FileName,
                delim = '\t',
                col_names = T)
  }
  List_NRMSE<-data.frame("pas" = SpectralSampling[[parm]], "nrmse"<-NRMSE_df[[1]])
  return(List_NRMSE)
}
