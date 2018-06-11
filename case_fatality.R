## Case fatality
# Load functions and packages and foldername
source("scripts/fnctns_pckgs_folder.r")

hf_cf <- read.csv (file = paste0(datastored, "/", "HFFatality", ".csv"), as.is = TRUE)

## Data processsing ----
names(hf_cf) <- tolower(names(hf_cf))

hf_cf <- hf_cf [ hf_cf$persontime != 0, ]
hf_cf <- hf_cf [ hf_cf$deprivation !=11, ]

# Sum admissions
sum(hf_cf$admission)
sum(hf_cf$death)

## As per JM, ignore heart failure deaths as PM diagnosis of heart failure considered unreliable
hf_cf$death <- 0

# Select only admitted or incident death patients, so looking at case-fatality among admitted patients
hf_cf <- filter (hf_cf, admission >=1 | death >= 1)
# Rescale and label variables for regression models
hf_cf <- filter(hf_cf, age  >=20)
hf_cf$sex <- hf_cf$sex - 1
hf_cf$type <- factor(hf_cf$type, levels = c(1,2,3,4), labels = c("other", "t1dm", "t2dm", "pop"))
hf_cf <- filter(hf_cf, type != "other")
hf_cf$type <- factor(hf_cf$type, levels = c("pop", "t1dm", "t2dm"))
hf_cf$age_ten <- hf_cf$age / 10
# Change ordering of deprivation so that 1 means low and 10 means high
hf_cf$deprivation <- 11 - hf_cf$deprivation
hf_cf$dep_two <- hf_cf$deprivation/5
## Centre year around mid-point (2008) for estimating associations
hf_cf$year_original <- hf_cf$year
hf_cf$year <- hf_cf$year - 2008
## Centre age around 50 for estimating associations 
hf_cf$age_ten <- hf_cf$age_ten - 5
## Set prediction age for graphs
predict_age <- 50 

## Summary statistics for case fatality ----
# By age, sex and type for table
cf_smr <- hf_cf %>%
  mutate (age_cut = cut2(age, cuts = c(30, 50, 70,90))) %>%
  group_by (age_cut, sex, type) %>%
  summarise (admission = sum(admission + death),
             ha30days = sum(ha30days + death),
             cf_adm30 = formatC(round(100* ha30days/admission,1), digits = 1, format = "f"),
             cf_res = paste0(cf_adm30, "% (", ha30days, "/", admission, ")")
             )
cf_age_sex_type <- reshape2::dcast (cf_smr, sex + age_cut ~ type, value.var = "cf_res")
# By sex, type
cf_smr <- hf_cf %>%
  mutate(age_cut = "overall") %>%
  group_by (age_cut, sex, type) %>%
  summarise (admission = sum(admission + death),
             ha30days = sum(ha30days + death),
             cf_adm30 = formatC(round(100* ha30days/admission,1), digits = 1, format = "f"),
             cf_res = paste0(cf_adm30, "% (", ha30days, "/", admission, ")")
  )
cf_sex_type <- reshape2::dcast (cf_smr, sex +age_cut ~ type, value.var = "cf_res")

cf_age_sex_type <- rbind(cf_age_sex_type, cf_sex_type)
cf_age_sex_type <- cf_age_sex_type[order(cf_age_sex_type$sex, cf_age_sex_type$age_cut),]
# By type for summary in text
cf_smr <- hf_cf %>%
  group_by (type) %>%
  summarise (admission = sum(admission + death),
             ha30days = sum(ha30days + death),
             cf_adm30 = formatC(round(100* ha30days/admission,1), digits = 1, format = "f"),
             cf_res = paste0(cf_adm30, "% (", ha30days, "/", admission, ")")
            )
cf_type <- cf_smr$cf_res
names(cf_type) <-cf_smr$type
cf_type

# Overall
cf_smr <- hf_cf %>%
  summarise (adm_death = sum(admission + death),
             death_or_30 = sum(ha30days + death),
             cf_adm30 = formatC(round(100* death_or_30/adm_death,1), digits = 1, format = "f"),
             cf_res = paste0(cf_adm30, "% (", death_or_30, "/", adm_death, ")")
  )
cf_overall <- cf_smr
## Modelling ----
hf_cf$success <- hf_cf$ha30days + hf_cf$death
hf_cf$failure <- pmax(0, hf_cf$admission + hf_cf$death - (hf_cf$ha30days + hf_cf$death))
mod_cf_final <- glm( cbind(success, failure)  ~ 
                       age_ten + I(age_ten^2) + sex + dep_two,
                     family = binomial, data = hf_cf)
summary (mod_cf_final) # all interactions weak

## Age, sex and deprivation interactions for diabetes tested
mod_cf_type <- update( mod_cf_final, . ~ . + type + type:sex)
summary (mod_cf_type) 
mod_cf_type_unad <- update( mod_cf_final, . ~  + type + type:sex)

mod_cf_type_age <- update( mod_cf_type, . ~ . + type:age_ten)
summary(mod_cf_type_age)

## Time trends
# No covariates
mod_cf_time_unad <- update (mod_cf_type, . ~  + year)
summary(mod_cf_time_unad) 
mod_cf_time_unad_type <- update(mod_cf_time_unad, . ~ . + type + year:type)
summary(mod_cf_time_unad_type)
# Statistically significant covariates
mod_cf_time_cov <- update(mod_cf_final, . ~ . + year)
summary(mod_cf_time_cov) 

mod_cf_time <- update(mod_cf_time_cov, . ~ . + type + year:type)
summary(mod_cf_time)

# smooth models
mod_cf_time_s <- gam (cbind(success, failure) ~  age_ten + I(age_ten^2) + sex + dep_two + 
                     s(year, by = type),
                   family = binomial, data = hf_cf)
summary(mod_cf_time_s)

## Plot estimates from models with data ----

## Rates by type, age and sex plot
# Create dataset for predictions
hf_cf_smry <- expand.grid(age_ten = seq(-3,4,0.1/10), sex = 0:1, dep_two = (1:10)/5, 
                       type = levels(hf_cf$type))
pred_res <- MakePredict(mymodel = mod_cf_type, hf_cf_smry, do_logit = TRUE)
hf_cf_smry <- cbind(hf_cf_smry, pred_res)
hf_cf_smry <- UnTransform (hf_cf_smry)

# Make datapoints for age plot
plot1_cf_points <- hf_cf %>% 
  group_by(cut2(age,  cuts = seq(40,80,10)), sex, type) %>%
  summarise (est = sum(success)/sum(success + failure),
             mysize = sum(success),
             age_ten = mean(age_ten))
plot1_cf_points <- UnTransform(plot1_cf_points, c("year", "dep_two"))
plot1_cf_points <- filter(plot1_cf_points, est != 0)

# Plot case fatality by age, sex and type
plot1 <- ggplot(filter (hf_cf_smry, dep_two == 5),
                aes(x = age_ten, y = est, colour = type,
                    group = interaction (type))) + 
  geom_line() +
  geom_ribbon(mapping = aes(ymin = lower, ymax = upper, fill = type), alpha = 0.1) +
  geom_point(data = plot1_cf_points, mapping = aes(size = mysize), alpha = 0.3) +
  facet_grid (~sex)  +
  scale_color_discrete("") + 
  scale_fill_discrete(("")) +
  scale_y_continuous("Case fatality (%)", labels = scales::percent, limits = c(0,NA)) +
  scale_x_continuous("Age (years)") +
  scale_size(guide = FALSE) 
plot1
figure_cf_diabetes_age_sex <- plot1
rm(plot1)

## Plot case fatality over time
mymodel <- "mod_cf_time"
hf_cf_smry_time <- expand.grid(age_ten = seq(-3,4,0.1), sex = 0:1, dep_two = (1:10)/5, 
                            type = levels(hf_cf$type), 
                            year = -4:5)
hf_cf_smry_time <- cbind (hf_cf_smry_time, 
                          MakePredict(mymodel = get(mymodel), hf_cf_smry_time, do_logit = TRUE))
hf_cf_smry_time <- UnTransform (hf_cf_smry_time, excludevars = "")

plot_cf_time <- ggplot(filter (hf_cf_smry_time, dep_two == 5 & age_ten == predict_age),
                    aes(x = year, y = est, colour = type, group = interaction(type))) + 
  geom_line() +
  geom_ribbon(alpha = 0.1, mapping = aes(ymin = lower, ymax = upper, fill = type)) +
  facet_grid (~sex)  +
  scale_color_discrete("") + 
  scale_fill_discrete(("")) +
  scale_x_continuous("Calendar years") +
  scale_y_continuous("Case fatality (%)", labels = scales::percent, limits = c(0,NA)) +
  scale_size(guide = FALSE) +
  scale_linetype(guide = FALSE)
plot_cf_time

## Process model estimates for paper and supplement ----
# Time trends estimates for text
cf_time_unad <- GlmToTable(mod_cf_time_unad, transformation = exp)
cf_time_diab <- GlmToTable(mod_cf_time, transformation = exp)

# Calculate estimates for women and men with t1 and t2
selectvar <- rep(0, length(coef(mod_cf_type)))
names(selectvar)<- rownames(vcov(mod_cf_type))
selectvar [ c("typet1dm",  "sex:typet1dm")] <- 1
cf_women_t1 <- CombineEst(mod_cf_type, selectvar)
selectvar["sex:typet1dm"] <- 0
cf_men_t1 <- CombineEst(mod_cf_type, selectvar)
selectvar[] <- 0
selectvar [ c("typet2dm",  "sex:typet2dm")] <- 1
cf_women_t2 <- CombineEst(mod_cf_type, selectvar)
selectvar["sex:typet2dm"] <- 0
cf_men_t2 <- CombineEst(mod_cf_type, selectvar)

## Models for supplement
summary_mod_cf_type <- summary(mod_cf_type)
vcov_mod_cf_type <- vcov(mod_cf_type)
summary_mod_cf_time <- summary(mod_cf_time)
vcov_mod_cf_time <- vcov(mod_cf_time)

# Odds ratios
mod_cf_type <- SimpleCI(mod_cf_type)
mod_cf_type_unad <- SimpleCI(mod_cf_type_unad)

# Compare with 5-year periods in Barasa et al
barasa <- tibble(yr_st = c(1987, 1992, 1997, 2002),
                 yr_end = c(1991, 1996, 2001, 2006),
                 hr = c(1, 0.69, 0.48, 0.40)) %>% 
  mutate(yr_mid = yr_st + (yr_st - yr_end),
         yr_diff = yr_mid - lag(yr_mid, 1),
         log_hr = log(hr),
         log_hr_diff = log(hr)- lag(log_hr),
         log_hr_diff_per_yr = log_hr_diff/yr_diff,
         hr_diff_per_yr = exp(log_hr_diff_per_yr) %>%  round(2)
  )

# Summary case fatality table ----
supp_table_cf <- hf_cf
supp_table_cf$age_bands <- cut2 (supp_table_cf$age, cuts = c(30, 50, 60, 70, 80))
lower_cut <- substr(levels(supp_table_cf$age_bands), 2, 3)
upper_cut <- as.numeric(substr(levels(supp_table_cf$age_bands), 5, 6)) -1
levels(supp_table_cf$age_bands) <- paste0(lower_cut, " to ", upper_cut)

supp_table_cf$dep_bands [supp_table_cf$deprivation %in% 1:5] <- "1-5"
supp_table_cf$dep_bands [supp_table_cf$deprivation %in% 6:10] <- "6-10"

supp_table_cf$year_cuts <- cut2 (supp_table_cf$year_original, cuts = seq(2006,2014,2))
supp_table_cf$sex <- factor(supp_table_cf$sex, levels = 0:1, labels = c("men", "women"))

supp_table_cf <- supp_table_cf %>% 
  group_by (age_bands, sex, dep_bands, type) %>%
  summarise(admission = sum(admission),
            death = sum(ha30days),
            persontime = round(sum(persontime),0))
supp_table_cf [ , c("admission", "death")] <- lapply(supp_table_cf [ , c("admission", "death")], 
                                                  function (x) ifelse(x <= 5, "<=5", x))
names(supp_table_cf) <- c("Age", "Sex", "Deprivation", "Diagnosis", "Admissions", "Deaths", "Persontime")
supp_table_cf <- supp_table_cf %>% 
  ungroup() %>% 
  select(-Persontime)

## Save objects for paper and supplement ----
save(cf_age_sex_type, cf_type, cf_overall, plot_cf_time, figure_cf_diabetes_age_sex,
     cf_time_unad, cf_time_diab,
     cf_women_t1, cf_men_t1,
     cf_women_t2, cf_men_t2,
     mod_cf_type, mod_cf_type_unad,
     summary_mod_cf_type, vcov_mod_cf_type,
     summary_mod_cf_time, vcov_mod_cf_time,
     supp_table_cf,
     file = "data/HF_cf_plots and data.Rdata")
