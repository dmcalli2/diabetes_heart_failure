## Estimate incidence rates by diabetes type

source("scripts/fnctns_pckgs_folder.r") 

## Read in files ----
# Reads data for main and sensitivity analyses
for (filename in c("HFbasic", "HFSensitivity", "HFExclIHD", "HfFirstPos")){
  hf_raw <- read.csv (file = paste0(datastored, "/", filename, ".csv"), as.is = TRUE)

## Data processing and checking plausibility of data ----
hf <- hf_raw
names(hf) <- tolower(names(hf))
hf <- hf [ hf$persontime != 0, ]
hf <- hf [ hf$deprivation !=11, ]

## Compare admissions in case fatality and admissions datasets
sum(hf$admission) 
sum(hf$death)

## As per JM, ignore heart failure deaths as PM diagnosis of heart failure considered unreliable
hf$count_death <- hf$death
hf$death <- 0

## Restrict to age 20 or older
hf <- dplyr::filter(hf, hf$age  >=20)
## Rescale and label variables for regression models
## Relabel types of diagnosis
hf$sex <- hf$sex - 1
hf$type <- factor(hf$type, levels = c(1,2,3,4), labels = c("other", "t1dm", "t2dm", "pop"))
hf <- dplyr::filter(hf, type != "other")
hf$type <- factor(hf$type, levels = c("pop", "t1dm", "t2dm"))
hf$age_ten <- hf$age / 10
## Change ordering of deprivation so that 1 means low and 10 means high to aid interpretability
hf$deprivation <- 11 - hf$deprivation
hf$dep_two <- hf$deprivation/5
hf$rate <- 1000 * hf$admission/ hf$persontime 
## Centre year around mid-point (2008)
hf$year_original <- hf$year
hf$year <- hf$year - 2008
## Centre age around 50
hf$age_ten <- hf$age_ten - 5
## Set prediction age for graphs
predict_age <- 50 

## Examine proportions with diabetes diagnoses ----
pt_age <- hf %>%
  group_by(type, age) %>%
  summarise(persontime_total = sum(persontime)) %>% 
  ungroup() %>%
  group_by(type) %>%
  mutate(persontime_prop = persontime_total/ sum(persontime_total))

age_plot <- ggplot(pt_age, aes(x = age, y = persontime_prop, colour = type)) +  geom_line() 
age_plot
tapply(hf$persontime, hf$type, function(x) format(sum(x), big.mark = ",", trim = TRUE))

## Calculate summary rates for table ---- 
# Using approximate method from Rothman, Epidemiology an Introduction, pages 137-8
sum_rate <- hf %>%
  filter (year_original == max(year_original))  %>%
  mutate (age_cut = cut2(age, cuts = c(30, 50, 70,90))) %>%
  group_by (age_cut, sex, type) %>%
  summarise (events = sum(admission + death),
             persontime = sum(persontime),
             rate = events/persontime)
rates <- reshape2::dcast (sum_rate, sex + age_cut ~ type, value.var = "rate")

diffs_t2 <- sum_rate %>%
  filter (type %in% c("pop", "t2dm")) %>%
  group_by (age_cut, sex) %>%
  summarise (rr = mypoisdiff(events, persontime, output = "rr"),
             arr = mypoisdiff(events, persontime, output = "arr"),
             type = "t2dm")

diffs_t1 <- sum_rate %>%
  filter (type %in% c("pop", "t1dm")) %>%
  group_by (age_cut, sex) %>%
  summarise (rr = mypoisdiff(events, persontime, output = "rr"),
             arr = mypoisdiff(events, persontime, output = "arr"),
             type = "t1dm")

## Pull into a table
diffs <- rbind (diffs_t1, diffs_t2)
rrs <- reshape2::dcast(diffs, sex + age_cut ~ type, value.var = "rr" )
arr <- reshape2::dcast(diffs, sex + age_cut ~ type, value.var = "arr" )
rrs$pop <- "-"
arr$pop <- "-"

sum_rate$rate <- round (1000* sum_rate$rate, 2)
sum_rate$persontime <- format (sum_rate$persontime, big.mark = ",", digits = 0, scientific = FALSE, trim = TRUE)
sum_rate$e_pt_rate <- with (sum_rate, paste0(rate, " (", events, "/", persontime, ")"))
e_pt_rate <- reshape2::dcast (sum_rate, sex + age_cut ~ type, value.var = "e_pt_rate")

table1 <- rbind (e_pt_rate, rrs, arr)
table1$myvar <- rep(c("e_pt_rate", "rrs", "arr"), each = length(table1$sex)/3)
table1$rate <- factor(table1$myvar, levels = c("e_pt_rate", "rrs", "arr"),
                      labels = c("Rate", "Rate ratio", "Rate diff."))
table1 <- table1 [ with(table1, order(sex, age_cut, rate )),]
table1$age_cut <- as.character(table1$age_cut)
table1$age_cut [table1$age_cut == "[20,30)"] <-"20-29"
table1$age_cut [table1$age_cut == "[30,50)"] <-"30-49"
table1$age_cut [table1$age_cut == "[50,70)"] <-"50-69"
table1$age_cut [table1$age_cut == "[70,90]"] <-"70-89"
table1$sex <- ifelse(table1$sex == 0, "Men", "Women")

# Rearrange table to improve readability
table1 <- table1 %>% 
  arrange(age_cut, sex)
table1$age_cut [duplicated(paste0(table1$sex, table1$age_cut))] <- ""
table1$sex [ duplicated(paste0(table1$sex, table1$age_cut))] <- ""
table1$sex[lag(table1$sex) == table1$sex] <- ""

table1 <- table1 [ c("sex", "age_cut", "rate", "pop", "t1dm", "t2dm")]
names(table1) <- c("Age", "Sex", "", "No diabetes", "Type 1", "Type 2") 
rm(sum_rate, arr, rrs, diffs, diffs_t1, diffs_t2)
table1[] <- map(table1, ~ ifelse(str_detect(.x, "NA|NaN"), "", .x))

## Cross sectional incidence rate models ----
# Incidence rate by age, sex and deprivation, without diabetes
mod_final <- glm( I(admission+ death) ~ age_ten*sex*dep_two - age_ten:sex:dep_two +  
                    offset(log(persontime)), family = quasipoisson, data = hf)
summary (mod_final)

## Age, sex and deprivation interactions for diabetes tested
mod_type <- update(mod_final, . ~ . + type + 
                      type:age_ten + type:sex + type:dep_two +
                      type:age_ten:sex + type:age_ten:dep_two)
summary (mod_type)
# This is final model for cross-sectional analysis

# Examine association for everyone aged over 60 to compare results with other paper
hf_cprd <- hf %>% 
  filter(age >= 60, type != "t1dm") %>% 
  mutate(typef = if_else(type == "t2dm" & sex ==1, 1, 0),
         typem = if_else(type == "t2dm" & sex ==0, 1, 0))

mod_final_over60 <- update(mod_final, . ~ . + typef + typem, data = hf_cprd)
cprd_coef <- exp(coef(mod_final_over60)) %>% round(2)
cprd_conf <- exp(confint(mod_final_over60)) %>% round(2)
cprd <- cbind(cprd_coef, cprd_conf) %>% 
  as_tibble() %>% 
  setNames(c("coef", "lci", "uci")) %>% 
  mutate(param = names(cprd_coef),
         res = paste0(coef, " 95% CI (", lci,"-", uci,")")) %>% 
  filter(param %in% c("typef", "typem")) %>% 
  select(param, res)

## Longitudinal incidence rate models ----
age_trend <- hf %>% 
  group_by (year_original) %>%
  summarise(age_mean = weighted.mean(age, persontime))

# Year only to examine trends
mod_time_unad <- update (mod_type, . ~  + year)
summary(mod_time_unad)

# covariates
mod_time_cov <- update(mod_final, . ~ age_ten*sex*dep_two - age_ten:sex:dep_two +
                     year +
                     year:age_ten +
                     offset(log(persontime)))
summary(mod_time_cov)
mod_time_cov_sex <- update (mod_time_cov, . ~ . + year:sex)
mod_time_cov_dep <- update (mod_time_cov, . ~ . + year:dep_two)

mod_time <- update (mod_time_cov, . ~ . +
                      type +
                      type:age_ten + type:sex + type:dep_two +
                      type:age_ten:sex + type:age_ten:dep_two + 
                      year:type)

# Smooth time model
mod_time_s <- gam (I(admission + death) ~ age_ten*sex*dep_two - age_ten:sex:dep_two + type + 
                     type:age_ten + type:sex + type:dep_two +
                     type:age_ten:sex + type:age_ten:dep_two + 
                     s(year, by = type) +
                     offset(log(persontime)),
                   family = quasipoisson, data = hf)
summary(mod_time_s)

# Replicate model from Shah et al to compare
hf_shah <- hf
hf_shah$predicted <- predict (mod_type, type = "response" )
hf_shah <- hf_shah %>% 
  group_by (age = age>=60, sex = sex == 0, type) %>%
  summarise(pred_events = sum(predicted),
            persontime = sum(persontime))
hf_shah$pred_rate <- hf_shah$pred_events / hf_shah$persontime
hf_shah <- reshape2::dcast (hf_shah, age + sex ~ type, value.var = "pred_rate" )
hf_shah$rr_t2 <- hf_shah$t2dm / hf_shah$pop
hf_shah [] <- lapply (hf_shah, round, 2)
hf_shah$rr_compare_shah <- hf_shah$rr_t2 / 2
paste(hf_shah$rr_t2, collapse = ", ")

## Produce model summaries ----
allres <- GlmToTable(mod_time)
time_unad <- GlmToTable(mod_time_unad)
time_cov <- GlmToTable(mod_time_cov)
time_cov_sex <- GlmToTable(mod_time_cov_sex)
time_cov_dep <- GlmToTable(mod_time_cov_dep)

## Plot estimates from models alongside aggregated data ----
## Rates by type, age and sex plot
# Create dataset for predictions
hf_smry <- expand.grid(age_ten = seq(-3,4,0.1/10), sex = 0:1, dep_two = (1:10)/5, 
                       type = levels(hf$type), persontime = 1000)
pred_res <- MakePredict(mymodel = mod_type, hf_smry)
hf_smry <- cbind(hf_smry, pred_res)
hf_smry <- UnTransform (hf_smry)

# Make datapoints for age plot
plot1_points <- hf %>% 
  group_by(age_ten, sex, type) %>%
  summarise (est = log(1000* sum(admission + death)/sum(persontime)),
             mysize = sum(admission + death))
plot1_points <- UnTransform(plot1_points, c("year", "dep_two"))
# Draw age, sex and type plot
plot1 <- ggplot(filter (hf_smry, dep_two == 5),
                aes(x = age_ten, y = est, colour = type,
                    group = interaction (type))) + 
  geom_line() +
  geom_ribbon(alpha = 0.1, mapping = aes(ymin = lower, ymax = upper, fill = type)) +
  geom_point(data = plot1_points, mapping = aes(size = mysize), alpha = 0.3) +
  facet_grid (~sex)  +
  scale_color_discrete("") + 
  scale_fill_discrete(("")) +
  scale_y_continuous("Rate per 1000 person-years", breaks = log(2^(-3:6)), labels = 2^(-3:6)) +
  scale_x_continuous("Age (years)") +
  scale_size(guide = FALSE)
plot1

figure1_diabetes_age_sex <- plot1
rm(plot1)

## Make table for appendix of data aggregated by age, sex and deprivation ----
hf_cross_sect <- filter(hf_smry, dep_two ==5  & (age_ten %in% seq(20, 80,5) ))
hf_cross_sect [ , c("est", "lower", "upper")] <- lapply(hf_cross_sect [ , c("est", "lower", "upper")],
                                                        function (x) round(exp(x),2))
hf_cross_sect <- select (hf_cross_sect, -persontime)
hf_cross_sect <-mutate (hf_cross_sect, 
                        fin_line = paste0(est, " (95% CI ", lower, " to ", upper, ")"),
                        fin_brac = paste0(est, "; 95% CI ", lower, " to ", upper))

hf_cross_sect <- reshape2::dcast(hf_cross_sect, sex  + age_ten ~ type, value.var = "fin_brac" )

## Plot age, sex, type on absolute scale ----
hf_smry_abs <- hf_smry
hf_smry_abs[ , c("est", "lower", "upper")] <- sapply (hf_smry_abs[ , c("est", "lower", "upper")] , exp)
plot1_points$est_abs <- exp(plot1_points$est)

plot_abs <- ggplot(filter (hf_smry_abs, dep_two == 5),
                aes(x = age_ten, y = est, colour = type,
                    group = interaction (type))) + 
  geom_line() +
  geom_ribbon(alpha = 0.1, mapping = aes(ymin = lower, ymax = upper, fill = type)) +
  geom_point(data = plot1_points, mapping = aes(size = mysize, y = est_abs), alpha = 0.3) +
  facet_grid (~sex)  +
  scale_color_discrete("") + 
  scale_fill_discrete(("")) +
  scale_y_continuous("Rate per 1000 person-years", limits = c(0,50)) +
  scale_x_continuous("Age (years)") +
  scale_size(guide = FALSE)
plot_abs
figure_sup_diabetes_age_sex_absolute <- plot_abs

## Plot rates by deprivation, sex and type ----
# Relative scale
plot1_dep <- ggplot(filter (hf_smry, age_ten == predict_age),
                aes(x = dep_two, y = est, colour = type,
                    group = interaction (type))) + 
  geom_line() +
  geom_ribbon(alpha = 0.1, mapping = aes(ymin = lower, ymax = upper, fill = type)) +
  facet_grid (~sex)  +
  scale_color_discrete("") + 
  scale_fill_discrete(("")) +
  scale_y_continuous("Rate per 1000 person-years", breaks = log(2^(-3:6)), labels = 2^(-3:6)) +
  scale_x_continuous("Deprivation deciles of SIMD") +
  scale_size(guide = FALSE)
plot1_dep
figure_sup_diabetes_dep_sex <- plot1_dep
rm(plot1_dep)
# Absolute scale
plot1_dep_abs <- ggplot(filter (hf_smry_abs, age_ten == predict_age),
                    aes(x = dep_two, y = est, colour = type,
                        group = interaction (type))) + 
  geom_line() +
  geom_ribbon(alpha = 0.1, mapping = aes(ymin = lower, ymax = upper, fill = type)) +
  facet_grid (~sex)  +
  scale_color_discrete("") + 
  scale_fill_discrete(("")) +
  scale_y_continuous("Rate per 1000 person-years") +
  scale_x_continuous("Deprivation deciles of SIMD") +
  scale_size(guide = FALSE)
plot1_dep_abs

## Plot of regression model by type and time ----
# Plot on absolute and relative scale
for (mymodel in c("mod_time", "mod_time_s")){
  for (myscale in c("relative", "absolute")){
    hf_smry_time <- expand.grid(age_ten = seq(-3,4,0.1), sex = 0:1, dep_two = (1:10)/5, 
                                type = levels(hf$type), 
                                year = -4:5,
                                persontime = 1000)
    hf_smry_time <- cbind (hf_smry_time, MakePredict(mymodel = get(mymodel), hf_smry_time))
    hf_smry_time <- UnTransform (hf_smry_time, excludevars = "")
    hf_smry_time$predicted <- ifelse(hf_smry_time$year <= 2014, "Predicted", "Estimated")
    hf_smry_time$predicted <- factor(hf_smry_time$predicted, levels = c("Predicted", "Estimated"))
    
    if (myscale == "absolute") hf_smry_time[ , c("est", "lower", "upper")] <- 
      lapply(hf_smry_time[ , c("est", "lower", "upper")], exp)
    
    plot_time <- ggplot(filter (hf_smry_time, dep_two == 5 & age_ten == predict_age),
                        aes(x = year, y = est, colour = type, group = interaction(type))) + 
      geom_line() +
      geom_ribbon(alpha = 0.1, mapping = aes(ymin = lower, ymax = upper, fill = type)) +
      facet_grid (~sex)  +
      scale_color_discrete("") + 
      scale_fill_discrete(("")) +
      scale_y_continuous("Rate per 1000 person-years", breaks = log(2^(-3:6)), labels = 2^(-3:6)) +
      scale_x_continuous("Calendar years") +
      scale_size(guide = FALSE) +
      scale_linetype(guide = FALSE)
    if (myscale == "absolute") plot_time <- plot_time + 
      scale_y_continuous("Rate per 1000 person-years")
    assign (paste("plot_time", mymodel, myscale, sep = "_"), plot_time)
    assign (paste("data_for_plot_time", mymodel, myscale, sep = "_"), hf_smry_time)

  }
}

## Summary stats for paper ----
strt_cohort <- min(hf$year) + 2000 + 8
end_cohort <- max(hf$year) + 2000 + 8

sum_data <- data.frame (
  events =  tapply( (hf$admission+ hf$death), hf$type, sum),
  deaths = tapply( hf$death, hf$type, sum),
  persontimes =  tapply( hf$persontime, hf$type, sum)
)
sum_data$rates <- 1000* sum_data$events / sum_data$persontimes
sum_data$prop_death <- 100*sum_data$deaths / sum_data$events
sum_data [, c("rates", "prop_death")] <- lapply (sum_data [ , c("rates", "prop_death")], round, digits = 1)
sum_data [] <- lapply (sum_data, function(x) format(x, scientific = FALSE, big.mark = ",", trim = FALSE))


sum_data_sa <- data.frame (
  events =  tapply( (hf$admission+ hf$count_death), hf$type, sum),
  admissions = tapply( hf$admission, hf$type, sum),
  persontimes =  tapply( hf$persontime, hf$type, sum)
)
sum_data_sa$rates <- 1000* sum_data_sa$events / sum_data_sa$persontimes
sum_data_sa$death_incrs <- sum_data_sa$events / sum_data_sa$admission
sum_data_sa$rates <- round(sum_data_sa$rates, digits = 1)
sum_data_sa$persontimes <- round(sum_data_sa$persontimes)
sum_data_sa [] <- lapply (sum_data_sa, function(x) format(x, scientific = FALSE, big.mark = ",", trim = FALSE))

incrs_if_incld_death <- round(sum(hf$admission+ hf$count_death)/ sum(hf$admission),1)

## Summarise model outputs for text ----
# Women in 2013
# Women and men, 2013, t2dm
selectvar <- rep(0, length(coef(mod_time)))
names(selectvar)<- rownames(vcov(mod_time))
selectvar [ c('typet2dm', 'sex:typet2dm', 'dep_two:typet2dm')] <- 1
selectvar [ c('year:typet2dm')] <- 2013-2008
# as set age to 50, can predict 50 with the age values being zero
selectvar [ c('age_ten:typet2dm', 'age_ten:sex:typet2dm', 'age_ten:dep_two:typet2dm')] <- 0
women <- CombineEst (mod_time, selectvar)
selectvar["sex:typet2dm"] <- 0
men <- CombineEst (mod_time, selectvar)
both_2013_t2dm <- rbind(women, men) 
rm(men, women)

# Women and men, 2013, t1dm
selectvar <- rep(0, length(coef(mod_time)))
names(selectvar)<- rownames(vcov(mod_time))
selectvar [ c('typet1dm', 'sex:typet1dm', 'dep_two:typet1dm')] <- 1
selectvar [ c('year:typet1dm')] <- 2013-2008
# as set age to 50, can predict 50 with the age values being zero
selectvar [ c('age_ten:typet1dm', 'age_ten:sex:typet1dm', 'age_ten:dep_two:typet1dm')] <- 0
women <- CombineEst (mod_time, selectvar)
selectvar["sex:typet1dm"] <- 0
men <- CombineEst (mod_time, selectvar)
both_2013_t1dm <- rbind(women, men) 
rm(men, women)

# 80-year old women t2dm
selectvar <- rep(0, length(coef(mod_type)))
names(selectvar)<- rownames(vcov(mod_type))
selectvar [ c("typet2dm",  "dep_two:typet2dm", "sex:typet2dm")] <- 1
selectvar [ c("age_ten:typet2dm","age_ten:sex:typet2dm", "age_ten:dep_two:typet2dm")] <- 8-5
women <- CombineEst (mod_type, selectvar)
selectvar [ c("sex:typet2dm", "age_ten:sex:typet2dm")] <- 0
men <- CombineEst (mod_type, selectvar)
both_age80 <- rbind (women, men)
rm(men, women)

# 30-year old women t2dm
selectvar <- rep(0, length(coef(mod_type)))
names(selectvar)<- rownames(vcov(mod_type))
selectvar [ c("typet2dm",  "dep_two:typet2dm", "sex:typet2dm")] <- 1
selectvar [ c("age_ten:typet2dm","age_ten:sex:typet2dm", "age_ten:dep_two:typet2dm")] <- 3-5
women <- CombineEst (mod_type, selectvar)
selectvar [ c("sex:typet2dm", "age_ten:sex:typet2dm")] <- 0
men <- CombineEst (mod_type, selectvar)
both_age30 <- rbind (women, men)
rm(men, women)

# 20-year old women t2dm
selectvar <- rep(0, length(coef(mod_type)))
names(selectvar)<- rownames(vcov(mod_type))
selectvar [ c("typet2dm",  "dep_two:typet2dm", "sex:typet2dm")] <- 1
selectvar [ c("age_ten:typet2dm","age_ten:sex:typet2dm", "age_ten:dep_two:typet2dm")] <- 2-5
women <- CombineEst (mod_type, selectvar)
selectvar [ c("sex:typet2dm", "age_ten:sex:typet2dm")] <- 0
men <- CombineEst (mod_type, selectvar)
both_age20 <- rbind (women, men)
rm(men, women)

## age comparison rate
nodep_men_pop <- filter(hf, sex == 0, dep_two == 1, type == "pop" )
mod_no_dep <- glm (I(admission + death) ~ age_ten + offset (log(persontime)), data = nodep_men_pop, family = quasipoisson)
out_no_dep <- cbind (coef(mod_no_dep), confint(mod_no_dep))
out_no_dep <- as.data.frame(exp(out_no_dep))
names(out_no_dep) <- c("est", "lower", "upper")
out_no_dep [] <- lapply(out_no_dep, function(x) format(round(x,2), digits = 3, nsmall = 2, trim = TRUE))
out_no_dep$fin_line <- paste0(out_no_dep$est, " (95% CI ", out_no_dep$lower, " to ", out_no_dep$upper, ")")
out_no_dep$fin_brac <- paste0(out_no_dep$est, "; 95% CI ", out_no_dep$lower, " to ", out_no_dep$upper)
out_no_dep <- out_no_dep[2,]

## deprivation association - not that intersting
dep_mod <- glm(I(admission + death) ~ age_ten + dep_two + sex + 
                 offset(log(persontime)), family = quasipoisson, data = hf)
out_dep <- cbind (coef(dep_mod), confint(dep_mod))
out_dep <- as.data.frame(exp(out_dep))
names(out_dep) <- c("est", "lower", "upper")
out_dep [] <- lapply(out_dep, function(x) format(round(x,2), digits = 3, nsmall = 2, trim = TRUE))
out_dep$fin_line <- paste0(out_dep$est, " (95% CI ", out_dep$lower, " to ", out_dep$upper, ")")
out_dep$fin_brac <- paste0(out_dep$est, "; 95% CI ", out_dep$lower, " to ", out_dep$upper)
out_dep <- out_dep["dep_two",]

## RR for men aged 50 for model with and without deprivation - type 2 DM
# Deprivation model
mod_type_no_dep <- update (mod_type, . ~ . - dep_two - age_ten:dep_two - sex:dep_two - dep_two:type -
                             age_ten:dep_two:type)
selectvar <- rep(0, length(coef(mod_type)))
names(selectvar)<- rownames(vcov(mod_type))
selectvar["typet2dm"] <- 1 
men_dep <- CombineEst (mod_type, selectvar)
# No deprivation model
selectvar <- rep(0, length(coef(mod_type_no_dep)))
names(selectvar)<- rownames(vcov(mod_type_no_dep))
selectvar["typet2dm"] <- 1 
men_no_dep <- CombineEst (mod_type_no_dep, selectvar)
dep_vs_no_dep_type2 <- rbind (men_dep, men_no_dep)
## RR for men aged 50 for model with and without deprivation - type 1 DM
# Deprivation model
selectvar <- rep(0, length(coef(mod_type)))
names(selectvar)<- rownames(vcov(mod_type))
selectvar["typet1dm"] <- 1 
men_dep <- CombineEst (mod_type, selectvar)
# No deprivation model
selectvar <- rep(0, length(coef(mod_type_no_dep)))
names(selectvar)<- rownames(vcov(mod_type_no_dep))
selectvar["typet1dm"] <- 1 
men_no_dep <- CombineEst (mod_type_no_dep, selectvar)
dep_vs_no_dep_type1 <- rbind (men_dep, men_no_dep)

dep_atten_t2 <- as.numeric(dep_vs_no_dep_type2[1,1])/ as.numeric(dep_vs_no_dep_type2[2,1])
dep_atten_t1 <- as.numeric(dep_vs_no_dep_type1[1,1])/ as.numeric(dep_vs_no_dep_type1[2,1])
dep_atten <- format(c(type1 = dep_atten_t1, type2 = dep_atten_t2), digits = 2, nsmall = 2)

## deprivation type interaction
mod_type_dep <- update (mod_type, . ~ . + dep_two:type)
summary(mod_type_dep)
out_dep_inter <- cbind (coef(mod_type_dep), confint(mod_type_dep))
out_dep_inter <- as.data.frame(exp(out_dep_inter))
names(out_dep_inter) <- c("est", "lower", "upper")
out_dep_inter [] <- lapply(out_dep_inter, function(x) format(round(x,2), digits = 3, nsmall = 2, trim = TRUE))
out_dep_inter$fin_line <- paste0(out_dep_inter$est, " (95% CI ", out_dep_inter$lower, " to ", out_dep_inter$upper, ")")
out_dep_inter$fin_brac <- paste0(out_dep_inter$est, "; 95% CI ", out_dep_inter$lower, " to ", out_dep_inter$upper)

## Supplementary table
supp_table <- hf
supp_table$age_bands <- cut2 (supp_table$age, cuts = c(50, 60, 70, 80))
lower_cut <- substr(levels(supp_table$age_bands), 2, 3)
upper_cut <- as.numeric(substr(levels(supp_table$age_bands), 5, 6)) -1
levels(supp_table$age_bands) <- paste0(lower_cut, " to ", upper_cut)

supp_table$dep_bands [supp_table$deprivation %in% 1:5] <- "1-5"
supp_table$dep_bands [supp_table$deprivation %in% 6:10] <- "6-10"

supp_table$year_cuts <- cut2 (supp_table$year_original, cuts = seq(2006,2014,2))
supp_table$sex <- factor(supp_table$sex, levels = 0:1, labels = c("men", "women"))

supp_table <- supp_table %>% 
  group_by (age_bands, sex, dep_bands, type) %>%
  summarise(admission = sum(admission),
            death = sum(death),
            persontime = round(sum(persontime),0))
supp_table [ , c("admission", "death")] <- lapply(supp_table [ , c("admission", "death")], 
                                                  function (x) ifelse(x <= 5, "<=5", x))
names(supp_table) <- c("Age", "Sex", "Deprivation", "Diagnosis", "Admissions", "Deaths", "Persontime")
## model coefficeints etc
summary_mod_type <- summary(mod_type)
vcov_mod_type <- vcov(mod_type)
summary_mod_time <- summary(mod_time)
vcov_mod_time <- vcov(mod_time)
summary_mod_time_s <- summary(mod_time_s)

## save objects for paper and supplement  ----
save(figure1_diabetes_age_sex, figure_sup_diabetes_dep_sex, 
     figure_sup_diabetes_age_sex_absolute, plot1_dep_abs,
     plot_time_mod_time_absolute, plot_time_mod_time_relative,
     plot_time_mod_time_s_absolute, plot_time_mod_time_s_relative,
     allres, time_unad, time_cov, time_cov_sex, time_cov_dep,
     strt_cohort, end_cohort, sum_data, sum_data_sa, incrs_if_incld_death,
     both_2013_t2dm, both_2013_t1dm, both_age80, both_age80,
     out_no_dep, out_dep, out_dep_inter, predict_age, 
     hf_cross_sect, supp_table,
     summary_mod_type, vcov_mod_type, summary_mod_time, vcov_mod_time,
     summary_mod_time_s, table1,
     dep_vs_no_dep_type1, dep_vs_no_dep_type2, dep_atten,
     file = paste0("data/", filename, "_plots and data.Rdata"))

# Save objects for interactive figure in Shiny app ----
save(data_for_plot_time_mod_time_absolute,
     data_for_plot_time_mod_time_relative,
     data_for_plot_time_mod_time_s_absolute,
     data_for_plot_time_mod_time_s_relative,
     file = paste0("data/", filename, "_for_shiny_app.Rdata"))
}



