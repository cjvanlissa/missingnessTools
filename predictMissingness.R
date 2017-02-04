#This function can be used to identify which (auxiliary) variables are related to missingness (in study variables, or dvs). It uses the Ranger algorithm to perform random forests. The output consists of a list of variable importance statistics per dv, and a simplified list, with the additive rank of each of the auxiliary variables as a predictor of missingness. 

#From Jochen Hardt, Max Herke and Rainer Leonhart (2012):
#"the imputation model should include all variables of the analysis, plus those highly correlated with responses or explanatory variables, and finally those variables that explain the mechanism leading to missing data e.g.[2, 14]."
#Rubin DB: Multiple imputations after 18 plus years. JASA. 1996, 91: 473-489.
#Schafer JL, Graham JW: Missing data: our view of the state of the art. Psychol Methods. 2002, 7: 147-177.

#In applied research contexts, it is often impractical to include all variables in a missingness model. Moreover, such an approach could lead to overfitting. This function is intended to help select a small subset of potentially important variables related to missingness in dvs, from a larger dataset.

predictMissingness<-function(data, dvs=NULL, auxiliary=NULL, ...){
  require(ranger)
  #require(missForest)
  
  #dvs=grep("soccom\\d+_1", names(data), value = TRUE)
  #auxiliary=c("lifesat2", "aandacht_t1", "aandacht_t2", "narcissism_t1", "narcissism_t2", "socmeddis_count_t1", "socmeddis_count_t2", "impulsiv_t1", "impulsiv_t2", "self_disclosure_t1", "self_disclosure_t2", "snsusers", "snsadopters", "snsdroppers", "sn_act3it_1", "sn_act3it_2", "snsdisc_tot_1", "snsdisc_tot_2", "snsdisc_norm_1", "snsdisc_norm_2", "sndisc_express_1", "sndisc_express_2", "snsdisc_sexy_1", "snsdisc_sexy_2", "attention_mn_1", "attention_mn_2", "regret_mn_1", "shame_mn_1", "negemo_mn_1", "impuls_mn_1", "impuls_mn_2", "fomo_mn_1", "fomo_mn_2", "publicsa_mn_1")

  if(is.null(dvs)&&!is.null(auxiliary)){
    dvs<-names(data)[!(names(data)%in% auxiliary)]
    }
  if(!is.null(dvs)&&is.null(auxiliary)){
    auxiliary<-names(data)[!(names(data)%in% dvs)]
    }
  if(length(which(c(is.null(dvs),is.null(auxiliary))))==2){ #if both are empty
    dvs<-names(data)
    auxiliary<-names(data)
  } 
  df<-data[,which(names(data) %in% c(dvs,auxiliary))]
  n<-nrow(df)
  xAttrib <- lapply(df, attributes)
  varType <- character(p)
  
  #Remove vars with all missing
  allmissing <- which(apply(is.na(df), 2, sum) == n)
  if(length(allmissing)){
    df <- df[, -allmissing]
    dvs<-dvs[!(dvs %in% names(allmissing))]
    auxiliary<-auxiliary[!(auxiliary %in% names(allmissing))]
    cat("  Removed variable(s) because all entries were missing:", names(df)[allmissing], "\n")
  }
  
  #Remove vars with constant values
  allsame <- which(apply(df, 2, function(x) length(rle(x)$lengths)) == 1)
  if(length(allmissing)){
    df <- df[, -allsame]
    dvs<-dvs[!(dvs %in% names(allsame))]
    auxiliary<-auxiliary[!(auxiliary %in% names(allsame))]
    cat("  Removed variable(s) because all entries were identical:", names(df)[allsame], "\n")
  }
  
  p <- ncol(df)
  df.imputed<-df
  
    for (t.co in 1:ncol(df.imputed)) {
    if (is.null(xAttrib[[t.co]])) {
      varType[t.co] <- "numeric"
      df.imputed[is.na(df[, t.co]), t.co] <- mean(df[, t.co], 
        na.rm = TRUE)
    }
    else {
      varType[t.co] <- "factor"
      max.level <- max(table(df.imputed[, t.co]))
      class.assign <- sample(names(which(max.level == table(df.imputed[, t.co]))), 1)
      if (class.assign != "NA's") {
        df.imputed[is.na(df[, t.co]), t.co] <- class.assign
      }
      else {
        while (class.assign == "NA's") {
          class.assign <- sample(names(which(max.level == table(df.imputed[, t.co]))), 1)
        }
        df.imputed[is.na(df[, t.co]), t.co] <- class.assign
      }
    }
  }
  

  df.ismissing<-as.data.frame(ifelse(is.na(df), 1,0))

  variableimportance<-lapply(dvs, function(x){
    data<-data.frame(df.ismissing[,x], df.imputed[,auxiliary[!auxiliary==x]])
    sort(ranger(as.formula(paste0(names(data)[1], "~.")), data=data, importance = "impurity", ...)$variable.importance, decreasing = TRUE)
  })
  names(variableimportance)<-dvs
  importance<-sort(sapply(auxiliary, function(x){
    sum(unlist(lapply(variableimportance, function(y){which(names(y)==x)})))
  }))
  plot(x=1:length(auxiliary), y=importance, xlab="Variable number", ylab="Additive rank (lower is better)", type="b")
  return(list(Importance=variableimportance, rank=importance))
}

