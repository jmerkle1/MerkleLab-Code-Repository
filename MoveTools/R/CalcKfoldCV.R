#' Conduct a k-fold cross validation analysis
#'
#' Performs cross-validation and generates the cross-validation score. Based on the script from Roberts et al 2017 Ecography, updated by M. Laforge and J. Merkle
#'
#' @param data The dataset for which models should be run
#' @param formula The formula of the model to be tested. Include response and explanatory variables, do NOT include a "data" argument,and do not include function name (e.g.: "Pres~x1+x2")
#' @param modType The type of model you want to fit. One of: "glmer" (default), "glm","lmer", "mclogit", and "glmmTMB".
#' @param binVar The variable that you want to bin observations by. Typically animal ID. Must >= than K.
#' @param k  Number of training and testing bins to use. Default = 5.
#' @param resp  The name of the response variable in your data. Default = Pres
#' @param modFam  The family of the model (default = binomial)
#' @param setSeed  The value to pass along to "TMBStruc$parameters$theta[1]" before running the model. This argument is used for the Muff et al (2020) method to fit conditional mixed effects models.
#' @param TMBtheta  Number of training and testing bins to use. Default = 5.
#' @param TMBmaparg  The value to pass along to "TMBStruc$mapArg" both when setting of the funtion and before fitting. This argument is used for the Muff et al (2020) method to fit conditional mixed effects models. Seems to be necessary to get "predict" to work when setting the variance structure manually for this method.
#' @param update  Whether function should update on it's progress. Default is T
#' @param numb.top.bins  How many top bins to include in the calculate of proportion of used points (in testing dataset) falling in top bins. Defaults to 1, which is the top bin out of 10.
#'
#'
#' @return Returns a vector of cross-validation scores (i.e., spearman rank correlations) for each fold. Uses the area-adjusted frequency of categories (bins). See Boyce et al. 2002 for details. Also provides the proportion of used points of the testing data that fall in the top bins (specified by numb.top.bins) identified from the quantiles of available data in the testing dataset, as well as the number of available points in the top bins.
#'
#' @examples
#' #To come
#'
#' @export

CalcKfoldCV<-function(
  data,
  formula,
  modType="glmer",
  binVar,
  k=5,
  resp="Pres",
  modFam="binomial",
  setSeed=FALSE,
  TMBtheta=NULL,
  TMBmaparg=NULL,
  update=TRUE,
  numb.top.bins=1
)
{
  x<-1:length(unique(eval(parse(text=paste("data$",binVar,sep="")))))

  spl<-split(x,cut(x,k,labels=FALSE))

  # split individuals randomly
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

  if(update==T){
    print("Fitting the models. On fold: ")
  }

  if(modType=="glmer"){
    require(lme4)
    for(i in 1:k){
      assign(paste("Mod",i,sep=""),glmer(eval(parse(text=paste(formula))),family=modFam,data[data$rand.vec != i,]))
      if(update==T){
        print(i)
      }else{
      }
    }
  }else{
    if(modType=="glm"){
      for(i in 1:k){
        assign(paste("Mod",i,sep=""),glm(eval(parse(text=paste(formula))),family=modFam,data[data$rand.vec != i,]))
        if(update==T){
          print(i)
        }else{
        }
      }
    }else{
      if(modType=="lmer"){
        require(lme4)
        for(i in 1:k){
          assign(paste("Mod",i,sep=""),lmer(eval(parse(text=paste(formula))),data[data$rand.vec != i,]))
          if(update==T){
            print(i)
          }else{
          }
        }
      }else{
        if(modType=="mclogit"){
          require(mclogit)
          for(i in 1:k){
            assign(paste("Mod",i,sep=""),mclogit(eval(parse(text=paste(formula))),data[data$rand.vec != i,]))
            if(update==T){
              print(i)
            }else{
            }
          }
        }else{
          if(modType=="glmmTMB"){
            require(glmmTMB)
            require(Matrix)
            for(i in 1:k){
              if(is.null(TMBmaparg)){
                TMBStruc = glmmTMB(eval(parse(text=paste(formula))), family=modFam,doFit=FALSE,data[data$rand.vec != i,])
              }else{
                TMBStruc = glmmTMB(eval(parse(text=paste(formula))), map=TMBmaparg,family=modFam,doFit=FALSE,data[data$rand.vec != i,])
                TMBStruc$mapArg = TMBmaparg
              }

              if(is.null(TMBtheta)){

              }else{
                TMBStruc$parameters$theta[1] = TMBtheta
              }

              assign(paste("Mod",i,sep=""),glmmTMB:::fitTMB(TMBStruc))
              if(update==T){
                print(i)
              }else{
              }
            }
          }else{
            print("error:model type not supported")
          }
        }
      }
    }
  }

  if(update==T){
    print("Predicting values")
  }

  for(i in 1:k){
    data$RSFscores[data$rand.vec == i]<-exp(predict(eval(parse(text=paste("Mod",i,sep=""))),
                                                    newdata=subset(data,rand.vec == i), type="link",
                                                    allow.new.levels=T))
    if(update==T){
      print(i)
    }
  }


  # Run the k-fold CV evaluation sensu Boyce et al. 2002
  dataset <- data[complete.cases(data[,"RSFscores"]),]

  toreturn <- do.call(rbind, lapply(1:k, function(w){
    fold <- subset(dataset,rand.vec==unique(dataset$rand.vec)[w])
    # grab the quantile of the RSF scores from the available data points
    q.pp <- quantile(fold$RSFscores[eval(parse(text=paste("fold$",resp,sep="")))==0],
                     probs=seq(0,1,.1)) ## computing quantiles of RSF scores
    bin <- rep(NA,length(fold$RSFscores))
    for (j in 1:10){
      bin[fold$RSFscores>=q.pp[j]& fold$RSFscores<q.pp[j+1]] = j  ## binning RSF scores (10 bins)
    }
    used<-eval(parse(text=paste("fold$",resp,sep="")))
    # --------------------------------------------------------
    a <- table(used,bin) ## area adjusted freq in used/available for each bin
    a <- t(a) #transpose the table
    a <- as.data.frame.matrix(a) ## the next few lines compute area-adjusted frequency of categories (bins) of RSF scores
    names(a) <- c("avail","used")
    a$areaadjusted <- rep(NA,length(10))
    sum0 <- sum(a[,1])
    sum1 <- sum(a[,2])
    a$areaadjusted <- (a[,2] / sum1 ) / (a[,1] / sum0)
    a$bins <- seq(1,10,by=1);a
    
    return(data.frame(k=w,
                      rho.spearman=with(a,cor.test(bins,areaadjusted,method="spearm"))$estimate,
                      prop.in.top.bins.obs=sum(a$used[(nrow(a)-numb.top.bins+1):nrow(a)])/sum(a$used),
                      prop.in.top.bins.expected=sum(a$avail[(nrow(a)-numb.top.bins+1):nrow(a)])/sum(a$avail)))
  }))
  return(toreturn)
}  #   END   ###################################
