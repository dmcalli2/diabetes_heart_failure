# Prevalence analysis
# Load functions and packages and foldername
source("scripts/fnctns_pckgs_folder.r")

hf_prev <- read.csv(file = paste0(datastored, "/HFPrevalence.csv"), as.is = TRUE)

# Data processsing ----
hf <- hf_prev
names(hf) <- tolower(names(hf))
hf <- filter (hf, hf != 0 | nohf != 0, !is.na(nohf), hf != - 1)
# Rescale and label variables for regression models
# Relabel types
hf <- filter(hf, age  >=20)
hf$sex <- hf$sex - 1
hf$type <- factor(hf$type, levels = c(1,2,3,4), labels = c("other", "t1dm", "t2dm", "pop"))
hf <- filter(hf, type != "other")
hf$type <- factor(hf$type, levels = c("pop", "t1dm", "t2dm"))
hf$age_ten <- hf$age / 10
# Change ordering of deprivation so that 1 means low and 10 means high
hf$deprivation <- 11 - hf$deprivation
hf$dep_two <- hf$deprivation/5
hf$n <- hf$hf + hf$nohf
hf$prev <- 100 * hf$hf/ hf$n 

## Summary statistics for 2004 ---- 
pt_age_2004 <- hf %>%
  filter(year == 2004) %>% 
  group_by(type) %>%
  summarise(age_m = weighted.mean (age, n), age_sd = weighted.sd (age, n),
            sex_t = sum (sex* n), sex_prop =  100 * sum (sex * n)/ sum(n),
            dep_m = weighted.mean (deprivation, n), dep_sd = weighted.sd (deprivation, n)) 
pt_age_2004 [,!names(pt_age_2004) %in% c("type", "sex_t")] <- 
  lapply(pt_age_2004[,!names(pt_age_2004) %in% c("type", "sex_t")] 
         , format, digits = 1, nsmall = 1, scientific = FALSE, trim = TRUE)
pt_age_2004 <- transmute (pt_age_2004, 
                          type = type,
                          age = paste0(age_m, " (", age_sd, ")"),
                          dep = paste0(dep_m, " (", dep_sd, ")"),
                          women = paste0(format(round(sex_t,0), big.mark = ","), " (", sex_prop, "%)"))
row.names(pt_age_2004) <- pt_age_2004$type

## Centre age around 50 for estimating associations 
hf$age_ten <- hf$age_ten - 5
## Set prediction age for graphs
predict_age <- 50 

# Compare prevalence by calendar year
by_year <- hf %>% group_by(type, year, sex) %>%
  summarise (hf = sum(hf), nohf = sum(nohf))
by_year$prev <- with (by_year, 100* hf / (hf + nohf))
by_year$type <- factor (by_year$type, levels = c("t2dm", "t1dm", "pop"),
                        labels = c("Type 2", "Type 1", "No diabetes"))
by_year$sex <- factor(by_year$sex, levels = 0:1, labels = c("Male", "Female"))
reshape2::dcast (by_year, sex + year ~ type, value.var = "prev")
prev_trend <- ggplot (by_year, aes(x = year, y = prev, colour = type)) + 
  geom_point() + geom_line() +
  facet_wrap (~sex) +
  scale_x_continuous("Year") +
  scale_y_continuous("Prevalence of heart failure (%)") +
  scale_color_discrete("")
prev_trend

# Remove all years except 2004 for prevalence estimation
hf <- filter(hf, year == 2004)

## Total number of events in each group ----
prev_smry <- hf %>% 
  group_by (type) %>%
  summarise (hf = sum(hf), nohf = sum(nohf))
prev_smry$hf_prop <- round (100 * prev_smry$hf / (prev_smry$hf + prev_smry$nohf),2)
prev_smry$n <- prev_smry$hf + prev_smry$nohf
prev_smry [,-1] <- lapply(prev_smry[,-1], format, big.mark = ",", trim = TRUE)
row.names(prev_smry) <- prev_smry$type

hf %>% group_by(sex) %>%  summarise (hf = sum(hf), nohf = sum(nohf))

## Regression models of associations between prevalence of diabetes after adjusting for age, sex and deprivation ----
mod_prev_base <- glm( cbind(hf, nohf) ~ age_ten*sex*dep_two - age_ten:sex:dep_two,
                      family = quasibinomial, data = hf)
summary (mod_prev_base)
# positive interactions all two-way interactions
mod_prev_type <- update(mod_prev_base, . ~  . + type + 
                          type:age_ten + type:sex + type:dep_two +
                          type:age_ten:sex + type:age_ten:dep_two + type:sex:dep_two)
summary (mod_prev_type)
# all two-way plus type interactions significant
# Create plots with modelled estimates, CIs and datapoints ----
hf_smry <- expand.grid(age_ten = seq(-3,4,0.1), sex = 0:1, dep_two = (1:10)/5, type = levels(hf$type))
pred_res <- MakePredict(mymodel = mod_prev_type, hf_smry, do_logit = FALSE)
hf_smry <- cbind(hf_smry, pred_res)
hf_smry <- UnTransform (hf_smry)

## Draw plot
plot1 <- ggplot(filter (hf_smry, dep_two == 5),
                aes(x = age_ten, y = est, colour = type,
                    group = interaction (type))) + 
  geom_line() +
  geom_ribbon(alpha = 0.1, mapping = aes(ymin = lower, ymax = upper, fill = type)) +
  facet_grid (~sex)  +
  scale_color_discrete("") + 
  scale_fill_discrete(("")) +
  scale_y_continuous("Heart failure prevalence", breaks = -(7:1), labels = round(100*plogis(-(7:1)),1)) +
  scale_x_continuous("Age (years)") +
  scale_size(guide = FALSE)
plot1
figure0_prev_diabetes_age_sex <- plot1
rm(plot1)

## Deprivation plot
plot1_dep <- ggplot(filter (hf_smry, age_ten == 50),
                aes(x = dep_two, y = est, colour = type,
                    group = interaction (type))) + 
  geom_line() +
  geom_ribbon(alpha = 0.1, mapping = aes(ymin = lower, ymax = upper, fill = type)) +
  facet_grid (~sex)  +
  scale_color_discrete("") + 
  scale_fill_discrete(("")) +
  scale_y_continuous("Heart failure prevalence", breaks = -(7:1), labels = round(100*plogis(-(7:1)),1)) +
  scale_x_continuous("Deprivation") +
  scale_size(guide = FALSE)
plot1_dep

## Statistics from model ouptuts ----
min_age <- min(hf$age)
# Prevalence young and old t2dm
# Women and men, t2dm ages 80 and 40
# Women aged 80
selectvar <- rep(0, length(coef(mod_prev_type)))
names(selectvar)<- rownames(vcov(mod_prev_type))
selectvar [ c('typet2dm', 'sex:typet2dm', 'dep_two:typet2dm', 'sex:dep_two:typet2dm')] <- 1
selectvar [ c('age_ten:typet2dm', 'age_ten:sex:typet2dm', 'age_ten:dep_two:typet2dm' )] <- (80-50)/10
women_80 <- CombineEstLgst (mod_prev_type, selectvar)
# women aged 40
selectvar [ c('age_ten:typet2dm', 'age_ten:sex:typet2dm', 'age_ten:dep_two:typet2dm' )] <- (40-50)/10
women_40 <- CombineEstLgst (mod_prev_type, selectvar)
# Men aged 80
selectvar <- rep(0, length(coef(mod_prev_type)))
names(selectvar)<- rownames(vcov(mod_prev_type))
selectvar [ c('typet2dm', 'dep_two:typet2dm')] <- 1
selectvar [ c('age_ten:typet2dm',  'age_ten:dep_two:typet2dm' )] <- (80-50)/10
men_80 <- CombineEstLgst (mod_prev_type, selectvar)
# Men aged 40
selectvar [ c('age_ten:typet2dm',  'age_ten:dep_two:typet2dm' )] <- (40-50)/10
men_40 <- CombineEstLgst (mod_prev_type, selectvar)
prev_out_t2dm <- rbind (men_40, men_80, women_40, women_80)
rm (men_40, men_80, women_40, women_80)

# Prevalence young and old t1dm
# Women and men, t2dm ages 60 and 40
# Women aged 60
selectvar <- rep(0, length(coef(mod_prev_type)))
names(selectvar)<- rownames(vcov(mod_prev_type))
selectvar [ c('typet1dm', 'sex:typet1dm', 'dep_two:typet1dm', 'sex:dep_two:typet1dm')] <- 1
selectvar [ c('age_ten:typet1dm', 'age_ten:sex:typet1dm', 'age_ten:dep_two:typet1dm' )] <- (60-50)/10
women_60 <- CombineEstLgst (mod_prev_type, selectvar)
# women aged 30
selectvar [ c('age_ten:typet1dm', 'age_ten:sex:typet1dm', 'age_ten:dep_two:typet1dm' )] <- (40-50)/10
women_40 <- CombineEstLgst (mod_prev_type, selectvar)
# Men aged 60
selectvar <- rep(0, length(coef(mod_prev_type)))
names(selectvar)<- rownames(vcov(mod_prev_type))
selectvar [ c('typet1dm', 'dep_two:typet1dm')] <- 1
selectvar [ c('age_ten:typet1dm',  'age_ten:dep_two:typet1dm' )] <- (60-50)/10
men_60 <- CombineEstLgst (mod_prev_type, selectvar)
# Men aged 30
selectvar [ c('age_ten:typet1dm',  'age_ten:dep_two:typet1dm' )] <- (40-50)/10
men_40 <- CombineEstLgst (mod_prev_type, selectvar)
prev_out_t1dm <- rbind (men_40, men_60, women_40, women_60)
rm (men_40, men_60, women_40, women_60)

# Prevalence deprived and not deprived and old t2dm
# Women and men, t2dm aged 50
selectvar <- rep(0, length(coef(mod_prev_type)))
names(selectvar)<- rownames(vcov(mod_prev_type))
selectvar [ c('typet2dm', 'sex:typet2dm')] <- 1
selectvar [ c('dep_two:typet2dm', 'sex:dep_two:typet2dm')] <- 1/10
women_dep1 <- CombineEstLgst (mod_prev_type, selectvar)
selectvar [ c('sex:typet2dm', 'sex:dep_two:typet2dm')] <- 0
men_dep1 <- CombineEstLgst (mod_prev_type, selectvar)
selectvar <- rep(0, length(coef(mod_prev_type)))
names(selectvar)<- rownames(vcov(mod_prev_type))
selectvar [ c('typet2dm', 'sex:typet2dm')] <- 1
selectvar [ c('dep_two:typet2dm', 'sex:dep_two:typet2dm')] <- 9/10
women_dep10 <- CombineEstLgst (mod_prev_type, selectvar)
selectvar [ c('sex:typet2dm', 'sex:dep_two:typet2dm')] <- 0
men_dep10 <- CombineEstLgst (mod_prev_type, selectvar)
prev_out_dep_t2dm <- rbind(women_dep1, men_dep1, women_dep10, men_dep10)
rm(women_dep1, men_dep1, women_dep10, men_dep10)

## Deprivation output for type 1 diabetes
selectvar <- rep(0, length(coef(mod_prev_type)))
names(selectvar)<- rownames(vcov(mod_prev_type))
selectvar [ c('typet1dm', 'sex:typet1dm')] <- 1
selectvar [ c('dep_two:typet1dm', 'sex:dep_two:typet1dm')] <- 1/10
women_dep1 <- CombineEstLgst (mod_prev_type, selectvar)
selectvar [ c('sex:typet1dm', 'sex:dep_two:typet1dm')] <- 0
men_dep1 <- CombineEstLgst (mod_prev_type, selectvar)
selectvar <- rep(0, length(coef(mod_prev_type)))
names(selectvar)<- rownames(vcov(mod_prev_type))
selectvar [ c('typet1dm', 'sex:typet1dm')] <- 1
selectvar [ c('dep_two:typet1dm', 'sex:dep_two:typet1dm')] <- 9/10
women_dep10 <- CombineEstLgst (mod_prev_type, selectvar)
selectvar [ c('sex:typet1dm', 'sex:dep_two:typet1dm')] <- 0
men_dep10 <- CombineEstLgst (mod_prev_type, selectvar)
prev_out_dep_t1dm <- rbind(women_dep1, men_dep1, women_dep10, men_dep10)
rm(women_dep1, men_dep1, women_dep10, men_dep10)

## Save models for supplement ----
summary_mod_prev_type <- summary(mod_prev_type)
vcov_mod_prev_type <- vcov(mod_prev_type)

## Save objects for use in paper and supplment via rmarkdown document ----
save(prev_smry, figure0_prev_diabetes_age_sex, plot1_dep, min_age, prev_out_t1dm, prev_out_t2dm,
     prev_out_dep_t1dm, prev_out_dep_t2dm,prev_trend, pt_age_2004,
     summary_mod_prev_type, vcov_mod_prev_type,
     file = paste0("data/HFprev", "_plots and data.Rdata"))