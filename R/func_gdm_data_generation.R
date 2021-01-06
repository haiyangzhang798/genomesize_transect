## this R code is aimming for the GDM analysis
### writen by Haiyang Zhang (zhanghaiyang798@gmail.com)
##########################################################################################
##########################################################################################
# this funciton need to read community data, location of each sites, and environmental data. 
# as a return, it will give the dataset for making gdm figure.
# it can producing total and turnover component with corresponding modification.



# -------------------------------------------------------------------------


func_gdm_run <- function(comm, space, env, dis_select="turnover"){
  
## 1. read in community data and clean it.
  comm.data <- as.matrix(comm[, -c(1,2)])
  rownames(comm.data) <- comm$SiteID
  comm.data[comm.data > 0] <- 1 ## here we transfer abundance-based community to absent/present data.
  library(betapart)
  tdf_pa.core <- betapart.core(comm.data)
  #yield 3 dissimilarity matrices pairwise disimilaries for all sites
  tdf.dist <- beta.pair(tdf_pa.core, index.family="sor")
   
  if (dis_select=="turnover")
    Dist_Mat  <- as.matrix(tdf.dist$beta.sim)
  else if (dis_select=="total")
    Dist_Mat  <- as.matrix(tdf.dist$beta.sor)
  
## 2. read in location and environmental data
  Space <- space[, c("SiteID", "X", "Y")]
  Env <- env[, 2:7]
  Env <- merge(Env, Space, by = "SiteID")
  ## here we need to make sure comm, space and env data have all the same siteID.
  Env <- subset(Env, Env$SiteID %in% comm$SiteID)
  Space <- subset(Space, Space$SiteID %in% comm$SiteID)
  library(gdm)
  
    ##########################################################################################
    # 3.1 turnover from betapart

  ##########################################################################################
  
    comm$SiteID <- rownames(Dist_Mat)
    bio_Data = cbind(comm$SiteID, Dist_Mat)
    colnames(bio_Data)[1] <- "SiteID" ### this code is very important to make sure the colname right for the formatsitepair analysis 
    bio_Data <- data.frame(bio_Data)
    
    data_GDM_S_Env <- formatsitepair(bioData = bio_Data,bioFormat=3,siteColumn="SiteID",XColumn="X",YColumn="Y",predData=Env)
    data_GDM_S_Space <- formatsitepair(bio_Data,bioFormat=3,siteColumn="SiteID",XColumn="X",YColumn="Y",predData=Space)
    ##Fit GDM with only ENVIRONMENT (geo=FALSE), only Space and both environment and Space variables (geo=TRUE)
    tab_GDM_Env <- gdm(data_GDM_S_Env,geo=FALSE)
    tab_GDM_Space <- gdm(data_GDM_S_Space,geo=TRUE)
    tab_GDM_Both <- gdm(data_GDM_S_Env,geo=TRUE)
    
    ##...Here the partitioning of the variance (following Legendre 2008)
    var_env_interaction <- tab_GDM_Env$explained ##...Variance explained by Environment (A) + the interaction between Environment and Space (B).
    var_space_interaction <- tab_GDM_Space$explained ##...Variance explained by Space (C) + the interaction (B).
    var_env_space_interaction <- tab_GDM_Both$explained ##...Variance explained by Environment (A) + Space (C) + the interaction (B).
    
    var_Interaction <- round(var_env_interaction+var_space_interaction-var_env_space_interaction,4) ##...Variance explained by the interaction (B)
    var_env <- round(var_env_interaction-var_Interaction,4) ##...Variance explained by Environment (A)
    var_space <- round(var_space_interaction-var_Interaction,4) ##...Variance explained by Space (C)
    var_residual <- round(100-var_env_space_interaction,4) ##...Unexplained variance (residuals)
    var_explain <- c(var_Interaction,var_env,var_space,var_residual)
    var_class <- c('var_Interaction', 'var_env', 'var_space', 'var_residual')
    tab_GDM_Both$df_var_final <- data.frame(var_class,var_explain)
    
    IMP_SPECIES_ALL<-numeric()
    for (i in 1:length(tab_GDM_Both$predictor)) IMP_SPECIES_ALL[i]<-sum(tab_GDM_Both$coefficients[(1+(3*(i-1))):(i*3)])
    EXP<-tab_GDM_Both$explained/100
    IMP_SPECIES_ALL<-(IMP_SPECIES_ALL*EXP)/sum(IMP_SPECIES_ALL)
    names(IMP_SPECIES_ALL) <- tab_GDM_Both$predictor
    tab_GDM_Both$IMP_SPECIES_ALL <- data.frame(IMP_SPECIES_ALL)
    ##########################################################################################
  return(tab_GDM_Both)
  }
