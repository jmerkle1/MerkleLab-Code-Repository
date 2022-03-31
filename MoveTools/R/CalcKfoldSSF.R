#' Conduct a k-fold cross validation analysis for SSF
#'
#' Performs cross-validation and generates the cross-validation scores based on the number of used steps in each rank within each strata. Also calculates the proportion of used steps in a specified number of the top ranks. Based on the script from Roberts et al 2017 Ecography, updated by M. Laforge and J. Merkle to follow the methods of Fortin et al. 2009 (Ecology)
#'
#' @param data The dataset for which models should be run (should be a regular dataframe)
#' @param formula The formula (as a character string) of the model to be tested. Include response, explanatory variables, random effects, cluster term etc. Do NOT include a "data" argument,and do not include function name
#' @param modType The type of model you want to fit. One of: "clogit" (default), or "glmmTMB".
#' @param binVar The variable that you want to bin observations by. Typically animal ID. Must >= than K.
#' @param k  Number of training and testing bins to use. Default = 5.
#' @param resp  The name of the response variable in your data. Default = "case"
#' @param strat  The name of the column representing strata in your data. Default = "strata"
#' @param setSeed  Whether or not to set the randomization to the same each run (default = True)
#' @param TMBtheta   The value to pass along to "TMBStruc$parameters$theta[1]" before fitting a glmmTMB model. This argument is used for the Muff et al (2020) method to fit conditional mixed effects models. Defaults to log(1e3).
#' @param TMBmaparg  The value to pass along to "TMBStruc$mapArg" before fitting a glmmTMB model. This argument is used for the Muff et al (2020) method to fit conditional mixed effects models. Defaults to list(theta=factor(c(NA,1:1))), which is 1 random effect.
#' @param numb.top.ranks  How many top ranks to include in the calculate of proportion of used steps in top ranks Defaults to 1, which is the top bin.
#' @param update  Whether function should update on it's progress. Default is TRUE
#'
#'
#' @return Returns a dataframe of each fold, cross-validation scores (i.e., spearman rank correlations) for each fold, and the proportion of used steps in the specified number of the top ranks.
#'
#' @examples
#' #To come
#'
#' @export


CalcKfoldSSF <- function(
  data,
  formula,
  modType="clogit",
  binVar,
  k=5,
  resp="case",
  strat="strata",
  setSeed=TRUE,
  TMBtheta=log(1e3),
  TMBmaparg=list(theta=factor(c(NA,1:1))),
  update=TRUE,
  numb.top.ranks=1
)
{

  # some checks
  if("sf" %in% class(data))
    stop("Your data cannot be a sf dataframe")
  print("Warning: This method calculates predicted values based only on the fixed effects (i.e., setting all the random effects to 0)! It does not calculate based on the average of the random coefficients.")
  if(length(unique(tapply(data[,resp], data[,strat], sum)))>1)
    stop("You have some strata without a used step!")
  if(length(unique(table(data[,resp])))>1)
    print("Warning: You do not have the same number of available steps across all strata. Correlations might be biased!")

  # prep for breaking up into folds
  x<-1:length(unique(eval(parse(text=paste("data$",binVar,sep="")))))
  spl<-split(x,cut(x,k,labels=FALSE))

  # split randomly based on binVar
  newdata <- data.frame(binVar = unique(eval(parse(text=paste("data$",binVar,sep="")))))
  if(setSeed==T){
    set.seed(k)
  }

  random_sample <- data.frame(binVar = sample(newdata$binVar,
                                              length(unique(eval(parse(text=paste("data$",binVar,sep="")))))))
  random_sample$rand.vec <- 0

  for(i in 1:k){
    random_sample$rand.vec[min(spl[[i]]):max(spl[[i]])]<-i
  }

  data <- merge(data, random_sample, by.x = binVar, by.y = "binVar", all.x = T )

  if(update==T)
    print("Fitting the models. On fold: ")

  if(modType=="clogit"){
    require(survival)
    for(i in 1:k){  # HERE IS WHERE Parallel processing could be implemented for at least across the Ks!!!!!
      assign(paste("Mod",i,sep=""),clogit(eval(parse(text=paste(formula))),data[data$rand.vec != i,],method = "efron"))
      if(update==T)
        print(i)
    }
  }else{
    if(modType=="glmmTMB"){
      require(glmmTMB)
      for(i in 1:k){  # HERE IS WHERE Parallel processing could be implemented for at least across the Ks!!!!!
        TMBStruc = glmmTMB(eval(parse(text=paste(formula))), family=poisson,doFit=FALSE,data[data$rand.vec != i,])
        TMBStruc$mapArg = TMBmaparg
        TMBStruc$parameters$theta[1] = TMBtheta
        assign(paste("Mod",i,sep=""),glmmTMB:::fitTMB(TMBStruc))   # fit model!
        if(update==T)
          print(i)
      }
    }else{
      stop("error:model type not supported")
    }
  }

  if(update==T)
    print("Predicting values")

  if(modType=="clogit"){
    for(i in 1:k){
      data$SSFscores[data$rand.vec == i]<- predict(eval(parse(text=paste("Mod",i,sep=""))),
                                                   newdata=subset(data,rand.vec == i), type="risk",
                                                   reference="sample")
      if(update==TRUE)
        print(i)
    }
  }else{
    for(i in 1:k){

      # y <- predict(eval(parse(text=paste("Mod",i,sep=""))), newdata = subset(data,rand.vec == i),
      #              type = "link", re.form = NA, allow.new.levels = TRUE)

      # skirt around using predict() as it is finicky with glmmTMB for SSFs. Used some code from Bsmitty here
      # get the fixed coefficients anyway
      mm1 <- model.matrix(formula(eval(parse(text=paste("Mod",i,sep="")))$modelInfo$terms$cond$fixe),
                          data = subset(data,rand.vec == i)) # Manually type in fixed effects only
      betas <- fixef(eval(parse(text=paste("Mod",i,sep=""))))$cond
      data$SSFscores[data$rand.vec == i]<- mm1 %*% betas
      if(update==T)
        print(i)
    }
  }



  # Run the k-fold CV evaluation sensu Fortin et al. 2009 (Ecology)
  dataset <- data[complete.cases(data[,"SSFscores"]), c(strat, resp, "SSFscores","rand.vec")]

  u <- unique(dataset[,strat])
  #for each strata, rank SSFscores to determine the predicted rank of the target point
  dataset <- do.call(rbind, lapply(1:length(u), function(i){
    temp <- dataset[dataset[,strat] == u[i],]
    temp <- temp[order(temp$SSFscores),]
    temp$rank <- 1:nrow(temp)
    return(temp)
  }))

  # loop through each fold and figure out the spearman rank of the bin number and the number of used points in each bin
  toreturn <- do.call(rbind, lapply(1:k, function(w){
    fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[w])

    ranks <- data.frame(ranks=1:max(fold$rank))

    obs.rks <- as.data.frame(table(fold$rank[fold[,resp] == 1]))
    names(obs.rks) <- c("ranks","Freq")

    ranks <- merge(ranks, obs.rks, all.x = TRUE)
    ranks$Freq[is.na(ranks$Freq)==TRUE] <- 0

    return(data.frame(k=w, rho.spearman=cor(ranks$ranks,ranks$Freq,method="spearman"),
                      prop.in.top.ranks.obs=sum(ranks$Freq[(nrow(ranks)-numb.top.ranks+1):nrow(ranks)])/sum(ranks$Freq),
                      prop.in.top.ranks.expected=numb.top.ranks/max(fold$rank)))
  }))

  return(toreturn)
}  #   END  OF FUNCTION
