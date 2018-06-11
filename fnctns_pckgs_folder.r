# dm_hf_00.r
# Read in functions, packages and name folder

source("scripts/file_location_ignore_git.R")
## Functions ----

mypoisdiff <- function (x = c(10,20), y = c(100, 100), output = "arr", method = "rothman") {
    # Rothman formula for rate difference
    r_diff <- x[2]/y[2] - x[1]/y[1] 
    r_diff_var <- (x[1]/y[1]^2 + x[2]/y[2]^2 )
    r_diff_se <- r_diff_var^0.5
    r_diff_se_95 <- 1.96 *r_diff_se
    r_diff_lci <- r_diff - r_diff_se_95
    r_diff_uci <- r_diff + r_diff_se_95
    r_diff <- 1000*c(r_diff, r_diff_lci, r_diff_uci)
    # ROthman formula for rate ratio
    rr <- (x[2]/y[2]) /  (x[1]/y[1])
    rr_ln <- log(rr)
    rr_ln_var <- 1/x[1] + 1/x[2]
    rr_ln_se <- rr_ln_var^0.5
    
    rr_ln_lci <- rr_ln - 1.96*rr_ln_se
    rr_ln_uci <- rr_ln + 1.96*rr_ln_se
    rr_ln <- c(rr_ln, rr_ln_lci, rr_ln_uci)

    if(output == "arr") res <- r_diff else res <- exp(rr_ln)
    
    res <- format( round(res,2), digits = 2, nsmall = 2, trim = FALSE)
    paste0(res[1], " (", res[2], "-", res[3], ")")
}

CombineEst <- function (mymodel, myvars){
  covwant <- vcov(mymodel)[ myvars !=0, myvars !=0]
  coefwant <- coef(mymodel) [ myvars !=0]
  coefwant <- coefwant * myvars[ myvars!=0]
  std_er <- sqrt(t(coefwant) %*% covwant %*% coefwant)
  mydf<- data.frame (
    est = exp(sum(coefwant)),
    lower = exp(sum(coefwant) - 1.96*std_er),
    upper = exp(sum(coefwant) + 1.96*std_er)
  )
  mydf <- lapply(mydf, function (x) format(round(x,2), digit = 3, trim = TRUE))
  mydf$fin_line <- paste0(mydf$est, " (95% CI ", mydf$lower, " to ", mydf$upper, ")")
  mydf$fin_brac <- paste0(mydf$est, "; 95% CI ", mydf$lower, " to ", mydf$upper)
  mydf
}

CombineEstLgst <- function (mymodel, myvars){
  covwant <- vcov(mymodel)[ myvars !=0, myvars !=0]
  coefwant <- coef(mymodel) [ myvars !=0]
  coefwant <- coefwant * myvars[ myvars!=0]
  std_er <- sqrt(t(coefwant) %*% covwant %*% coefwant)
  mydf<- data.frame (
    est = exp(sum(coefwant)),
    lower = exp(sum(coefwant) - 1.96*std_er),
    upper = exp(sum(coefwant) + 1.96*std_er)
  )
  mydf <- lapply(mydf, function (x) format(round(x,2), digit = 2, trim = TRUE, nsmall = 2))
  mydf$fin_line <- paste0(mydf$est, " (95% CI ", mydf$lower, " to ", mydf$upper, ")")
  mydf$fin_brac <- paste0(mydf$est, "; 95% CI ", mydf$lower, " to ", mydf$upper)
  mydf
}

SimpleCI <- function (mymodel){
  coefwant <- coef(mymodel)
  ciwant   <- confint.default(mymodel)
  mydf <- data.frame (est = coefwant, lower = ciwant[,1], upper = ciwant[,2])
  rownames(mydf) <- rownames(ciwant)
  mydf[] <- lapply(mydf, function (x) formatC(round(exp(x),2), format = "f", digits = 2))
  mydf$fin_line <- paste0(mydf$est, " (95% CI ", mydf$lower, " to ", mydf$upper, ")")
  mydf$fin_brac <- paste0(mydf$est, "; 95% CI ", mydf$lower, " to ", mydf$upper)
  mydf
}

# Function to convert centred varaibles back to original scale and make prettier labels
UnTransform <- function (x, excludevars = "year") {
  within (x, {
    if (! "sex" %in% excludevars) sex <- factor(sex, levels = 0:1, labels = c("Men", "Women"))
    if (! "type" %in% excludevars) type <- factor(type, levels = c("t2dm", "t1dm", "pop"), labels = c("Type 2", "Type 1","No diabetes"))
    if (! "age_ten" %in% excludevars) age_ten <- (5 + age_ten) *10
    if (! "dep_two" %in% excludevars) dep_two <- 5 * dep_two
    if (! "year" %in% excludevars) year <- 2008 + year
  })
}


# Function to estimate modelled rates with 95% CI
MakePredict <- function(mymodel = mymodel, mynew = mydataframe, do_logit = FALSE){
  hf_p <- predict(mymodel, mynew, type = "link", se.fit = TRUE)
  hf_p <- as.data.frame(hf_p)
  hf_p$est <- hf_p$fit
  hf_p$lower <- hf_p$est - 1.96 * hf_p$se.fit
  hf_p$upper <- hf_p$est + 1.96 * hf_p$se.fit
  if (do_logit == TRUE)  hf_p [ c("est", "lower", "upper")] <- lapply ( hf_p [ c("est", "lower", "upper")], plogis)
  hf_p [ c("est", "lower", "upper")]
}

# function to make weighted sd
weighted.sd <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  a <- (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
                                            na.rm)
  a^0.5
}

# Summary results from final time model
GlmToTable <- function (mymodel = mod_time, transformation = exp) {
  mytable <- transformation(cbind (coef(mymodel), confint.default(mymodel)))
  mytable <- as.data.frame (format(round(mytable, 3), nsmall = 3, scientific = FALSE, trim = TRUE))
  names(mytable) <- c("est", "lower", "upper")
  mytable$fin_line <- paste0(mytable$est, " (95% CI ", mytable$lower, " to ", mytable$upper, ")")
  mytable$fin_brac <- paste0(mytable$est, "; 95% CI ", mytable$lower, " to ", mytable$upper)
  mytable$percent_decline <- 100*(1-as.numeric(as.character(mytable$est)))
  mytable$p_value <- summary(mymodel)$coef[,4]
  mytable$p_value <- ifelse(mytable$p_value <0.001, "<0.001", round(mytable$p_value, 3))
  mytable [ grep("year", row.names(mytable)), ]
}

# Simple count summarise function
Npercent <- function (x) {
  x <- na.omit(x) 
  n = sum(x)
  percent = round(100*mean(x), 1)
  paste0(n, " (", percent, "%)")
}

## packages ----

library(Hmisc)
library(mgcv)
library(tidyverse)
library(ggplot2)
# also require reshape2, but call with reshape2::dcast

