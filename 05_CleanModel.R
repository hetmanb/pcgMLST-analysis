#models.R
# Note: great tutorial available here: https://ladal.edu.au/regression.html#Mixed-Effects_Multinomial_Regression

## Load Package Libraries ##########

library(pacman)
pacman::p_load(here, 
               tidyverse,
               rio,
               clipr, 
               mclogit,
               tidymodels,
               sjPlot, 
               sjmisc,
               sjlabelled,
               sjstats,
               epiR,
               memisc)

# Note some bugfixes with mclogit requires me to use the dev version
#devtools::install_github("melff/mclogit",subdir="pkg")

## Prepare the modeling dataset ############

# model.data <- import(here("data", "model_snp_chicken.dta"))
model.data <- import(here("data","model_cgmlst99_chicken.tsv"))

model.data.final <- 
model.data %>% 
  dplyr::select(Isolate, incI1.cgMLST99, Serovar, Prov, Region, Year_Collected, Source, SourceLoc, Cluster) %>%
  mutate(incI1.cgMLST99 = as.factor(incI1.cgMLST99)) %>% 
  mutate(Cluster = as.numeric(Cluster)) %>%
  mutate(cat_group = case_when(incI1.cgMLST99 == "Group1" ~ 0,
                              incI1.cgMLST99 == "Group2" ~ 1, 
                              incI1.cgMLST99 == "Group3" ~ 2, 
                              incI1.cgMLST99 == "Group4" ~ 3, 
                              TRUE ~ NA), #other
         cat_group = factor(cat_group),
         cat_year = case_when(Year_Collected < 2013 ~ 0, 
                              Year_Collected == 2013 ~ 1, 
                              Year_Collected > 2013 ~ 2,
                              TRUE ~ NA),
         cat_region = case_when(Region == "ON" ~ 0,
                                Region == "QC" ~ 1,
                                Region == "WEST" ~ 2,
                                Region == "PRA" ~ 3,
                                Region == "ATL" ~ 4,
                                TRUE ~ NA),
         cat_region2 = case_when(Region == "ON" ~ 0,
                                 Region == "QC" ~ 1,
                                 Region == "ATL" ~ 1,
                                 Region == "WEST" ~ 2,
                                 Region == "PRA" ~ 2,
                                 TRUE ~ NA),
         # cat_source = case_when(Source == "Turkey" ~ 0,
         #                        Source == "Chicken" ~ 1, 
         #                        Source == "Human" ~ 2, 
         #                        TRUE ~ 3),
         cat_sourceloc = case_when(SourceLoc == "Hatchery" ~ 0,
                                   SourceLoc == "Farm" ~ 1,
                                   SourceLoc == "Abattoir" ~ 2,
                                   SourceLoc == "Diagnostic" ~ 3,
                                   SourceLoc == "Processing" ~ 4,
                                   SourceLoc == "Retail" ~ 5,
                                   TRUE ~ 6), #other
         cat_sourceloc2 = case_when(SourceLoc == "Hatchery" ~ 0,
                                    SourceLoc == "Farm" ~ 0,
                                    SourceLoc == "Diagnostic" ~ 0,
                                    SourceLoc == "Abattoir" ~ 1,
                                    SourceLoc == "Processing" ~ 1,
                                    SourceLoc == "Retail" ~ 2,
                                    TRUE ~ 3),
         cat_year = factor(cat_year), 
         cat_region = factor(cat_region), 
         cat_region2 = factor(cat_region2),
         cat_sourceloc = factor(cat_sourceloc), 
         cat_sourceloc2 = factor(cat_sourceloc2)) %>% 
  arrange(Cluster)

levels(model.data.final$cat_group) <- c("Group1", "Group2", "Group3", "Group4")
levels(model.data.final$cat_year) <- c("< 2013", "  2013", "> 2013")  

levels(model.data.final$cat_region) <- c("ON", "QC", "WEST", "PRA", "ATL")
levels(model.data.final$cat_region2) <- c("ON", "EAST", "WEST")

levels(model.data.final$cat_sourceloc) <- c("Hatchery", "Farm", "Abattoir", "Diagnostic", "Processing", "Retail", "Other")
levels(model.data.final$cat_sourceloc2) <- c("Pre-slaughter", "Abattoir", "Post-slaughter")


for(i in 1:nrow(model.data.final)){
  if(is.na(model.data.final$Cluster[i])){
    model.data.final$Cluster[i] <- model.data.final$Cluster[i-1] + 1 }
  else {
    model.data.final$Cluster[i] <- model.data.final$Cluster[i]
  }
}


## Estimate final models ####

# Set original referent categories for independent vars
model.data.final$cat_year <- relevel(model.data.final$cat_year, ref = '< 2013')
model.data.final$cat_region2 <- relevel(model.data.final$cat_region2, ref = "ON")
model.data.final$cat_sourceloc2 <- relevel(model.data.final$cat_sourceloc2, ref = "Pre-slaughter")


### Relevel the factor to change the referent to "Group1" ####
model.data.final$cat_group <- relevel(model.data.final$cat_group, ref = "Group1")
# Estimate the model 
final.model.group1 <- mblogit(cat_group ~ cat_year + cat_region2 + cat_sourceloc2, random = ~1|Cluster, data = model.data.final, method = "PQL")
# View the model results
tab_model(final.model.group1, show.reflvl = TRUE)


### Relevel for Group2 as base category ####
model.data.final$cat_group <- relevel(model.data.final$cat_group, ref = "Group2")
final.model.group2 <- mblogit(cat_group ~ cat_year + cat_region2 + cat_sourceloc2, random = ~1|Cluster, data = model.data.final, method = "PQL")
tab_model(final.model.group2, show.reflvl = TRUE)

### Relevel for Group3 as base category ####
model.data.final$cat_group <- relevel(model.data.final$cat_group, ref = "Group3")
final.model.group3 <- mblogit(cat_group ~ cat_year + cat_region2 + cat_sourceloc2, random = ~1|Cluster, data = model.data.final, method = "PQL")
tab_model(final.model.group3, show.reflvl = TRUE)


### Change the referent for Year, region and source######
# Note: once these have been changed, go back and run the above 3 models
# to obtain remaining estimates. 
model.data.final$cat_year <- relevel(model.data.final$cat_year, ref = '  2013')
model.data.final$cat_region2 <- relevel(model.data.final$cat_region2, ref = "EAST")
model.data.final$cat_sourceloc2 <- relevel(model.data.final$cat_sourceloc2, ref = "Abattoir")

## Test final model blups ####

# 1. Homogeneity of blups/random effects: 
# Extract the xB values from the model (non-conditional on random effects)
xb <- predict(object = final.model.group1, type = 'link', conditional = FALSE)

# Extract the linear predictors conditional on random effects 
xb.c <- predict(object = final.model.group1, type = 'link', conditional = TRUE)

# Extract the random intercepts by subtracting xb from xb.c 
blups <- xb.c - xb

# Plot the random effects as a function of xb
plot.new()
plot(x = xb[, 1], y = blups[, 1], type = "p", col = "blue")
points(x = xb[, 2], y = blups[, 2], col ="green") 
points(x= xb[, 3], y = blups[,3], col ="red")


## Test fit of random effects and model ####


model.data.final$cat_group <- relevel(model.data.final$cat_group, ref = "Group1")
# Estimate the model 
final.model <- mblogit(cat_group ~ cat_year + cat_region2 + cat_sourceloc2, random = ~1|Cluster, data = model.data.final, method = "PQL")
final.model.flat <- mblogit(cat_group ~ cat_year + cat_region2 + cat_sourceloc2, data = model.data.final, method = "PQL")

AIC(final.model)
AIC(final.model.flat)

model.0 <- mblogit(cat_group ~ 1, random = ~1|Cluster, data = model.data.final, method = "PQL")
model.1 <- mblogit(cat_group ~ cat_year, random = ~1|Cluster, data = model.data.final, method = "PQL")
model.2 <- mblogit(cat_group ~ cat_region2, random = ~1|Cluster, data = model.data.final, method = "PQL")
model.3 <- mblogit(cat_group ~ cat_sourceloc2, random = ~1|Cluster, data = model.data.final, method = "PQL")
model.1.2 <- mblogit(cat_group ~ cat_year + cat_region2, random = ~1|Cluster, data = model.data.final, method = "PQL")
model.1.3 <- mblogit(cat_group ~ cat_year + cat_sourceloc2, random = ~1|Cluster, data = model.data.final, method = "PQL")
model.2.3 <- mblogit(cat_group ~ cat_region2 + cat_sourceloc2, random = ~1|Cluster, data = model.data.final, method = "PQL")


model.1 <- mblogit(cat_group ~ cat_year, random = ~1|Cluster, data = model.data.final, method = "PQL")


model.fit <- data.frame(model = c("model.0", 
                                  "model.1",
                                  "model.2",
                                  "model.3",
                                  "model.1.2",
                                  "model.1.3",
                                  "model.2.3",
                                  "model.1.2.3"),
                        AIC = c(AIC(model.0),
                                AIC(model.1),
                                AIC(model.2),
                                AIC(model.3),
                                AIC(model.1.2),
                                AIC(model.1.3),
                                AIC(model.2.3),
                                AIC(final.model)))
model.fit

## Test for confounding using mixed logit models, not multinomial
library(GLMMadaptive)
model.1 <- mblogit(cat_group ~ cat_year, random = ~1|Cluster, data = model.data.final, method = "PQL")


######## Impacts on Year

cmodel.1 <- mixed_model(fixed = cat_group ~ cat_year,
                       random = ~ 1 | Cluster,
                       data = model.data.final,
                       family = binomial(link="logit"))

cmodel.1.2 <- mixed_model(fixed = cat_group ~ cat_year + cat_region2,
                       random = ~ 1 | Cluster,
                       data = model.data.final,
                       family = binomial(link="logit"))

cmodel.1.3 <- mixed_model(fixed = cat_group ~ cat_year + cat_sourceloc2,
                         random = ~ 1 | Cluster,
                         data = model.data.final,
                         family = binomial(link="logit"))

####### Impacts on Region

cmodel.2 <- mixed_model(fixed = cat_group ~ cat_region2,
                        random = ~ 1 | Cluster,
                        data = model.data.final,
                        family = binomial(link="logit"))

cmodel.2.1 <- mixed_model(fixed = cat_group ~ cat_region2 + cat_year,
                        random = ~ 1 | Cluster,
                        data = model.data.final,
                        family = binomial(link="logit"))

cmodel.2.3 <- mixed_model(fixed = cat_group ~ cat_region2 + cat_sourceloc2,
                        random = ~ 1 | Cluster,
                        data = model.data.final,
                        family = binomial(link="logit"))


######## Impacts on Source 

cmodel.3 <- mixed_model(fixed = cat_group ~ cat_sourceloc2,
                          random = ~ 1 | Cluster,
                          data = model.data.final,
                          family = binomial(link="logit"))

cmodel.3.1 <- mixed_model(fixed = cat_group ~ cat_sourceloc2 + cat_year,
                          random = ~ 1 | Cluster,
                          data = model.data.final,
                          family = binomial(link="logit"))

cmodel.3.2 <- mixed_model(fixed = cat_group ~ cat_sourceloc2 + cat_region2,
                          random = ~ 1 | Cluster,
                          data = model.data.final,
                          family = binomial(link="logit"))

######



## View model predictions and covariate patterns. ####

# Generate predictions from the final model using group1 as the base category
# Note:Response type returns predicted results for the model instead of the linear predictors. 
model_predictions <- predict(final.model.group1, type = "response", se.fit = TRUE, conditional = FALSE)

fit <- as.matrix(model_predictions[[1]])
se <- as.matrix(model_predictions[[2]])

# Generate upper and lower 95% CI's
ci.lower <- fit - 1.96*se
ci.upper <- fit + 1.96*se

fit.df <- as.data.frame(fit) %>% 
  distinct() %>% 
  mutate(cp = 1:n())

se.df <- as.data.frame(se) %>% 
  distinct() %>% 
  mutate(cp = 1:n())


# Gather the covariates from the final model dataset
covariate_df <- final.model.group1$data[, c("cat_year", "cat_region2", "cat_sourceloc2")] 
# Generate a running count of unique covariate patterns 
covariate_df$cp <- epi.cp(covariate_df[1:3])$id
# Remove duplicate covariate patterns
cp.df <- covariate_df %>% distinct()


# Join the covariate patterns with their predicted results
pred.df <- left_join(cp.df, fit.df, by = 'cp') %>% 
  left_join(se.df, by = 'cp') %>% 
  rename("Group1.x" = "Group1.fit", 
         "Group2.x" = "Group2.fit",
         "Group3.x" = "Group3.fit",
         "Group4.x" = "Group4.fit",
         "Group1.y" = "Group1.se",
         "Group2.y" = "Group2.se",
         "Group3.y" = "Group3.se",
         "Group4.y" = "Group4.se") %>% 
  arrange(-Group1.fit) %>% 
  mutate(Group1.rank = 1:n()) %>% 
  arrange(-Group2.fit) %>% 
  mutate(Group2.rank = 1:n()) %>% 
  arrange(-Group3.fit) %>% 
  mutate(Group3.rank = 1:n()) %>% 
  arrange(-Group4.fit) %>% 
  mutate(Group4.rank = 1:n()) %>%  
  write_clip()

    


