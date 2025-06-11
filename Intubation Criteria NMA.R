########## Packages and working drive ##########

# Packages
library(tidyverse)
library(googlesheets4)
library(UpSetR)
library(brms)
library(igraph)

# Working drive
setwd("~/Research/SHN Research/Criteria AHRF SR")

# Colours

c_navy <- "#001A49"
c_blue = "cornflower blue"
c_ice  <- "#AAC0CE"
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")



RR <- function(events_1, n_1, events_2, n_2){events_2*n_1/(events_1*n_2)}
se_logRR <- function(events_1, n_1, events_2, n_2){
  sqrt(1/events_1 + 1/events_2 - 1/n_1 - 1/n_2)
}
se_logit <- function(p,n){sqrt(1/(n*p) + 1/(n*(1-p)))}

logit <- function(p){log(p/(1-p))}
inv_logit <- function(x){exp(x)/(1+exp(x))}

mean_CrI <- function(x){
  c( mean = mean(x),
     lb95 = quantile(x, 0.025),
     ub95 = quantile(x, 0.975))
} 

########### Load criteria data #################

criteria_data <- 
  readxl::read_excel("Criteria_DataTable.xlsx") %>%
  fill(Study)

# how many studies use each theme / subtheme
criteria_data %>%
  count(Study, Theme) %>%
  pivot_wider(names_from = Theme,
              values_from = n,
              values_fill = 0) %>%
  pivot_longer(-Study, 
               names_to = "Theme") %>%
  group_by(Theme) %>%
  summarise(n = sum(value > 0)) %>%
  arrange(desc(n)) %>%
  print(n = 50)
  
# Redone to split out the respiratory subthemes
criteria_data %>%
  mutate(Theme = ifelse(
    Theme == "Respiratory" & (!is.na(Subtheme)),
    Subtheme, Theme)) %>%
  count(Study, Theme) %>%
  pivot_wider(names_from = Theme,
              values_from = n,
              values_fill = 0) %>%
  pivot_longer(-Study, 
               names_to = "Theme") %>%
  group_by(Theme) %>%
  summarise(n = sum(value > 0)) %>%
  arrange(desc(n)) %>%
  print(n = 50)

# upset plot figure 

upset_dat <-
  criteria_data %>%
    mutate(Theme = ifelse(
      Theme == "Respiratory" & (!is.na(Subtheme)),
      Subtheme, Theme)) %>%
    count(Study, Theme) %>%
    mutate( n = ifelse(n>0,1,0)) %>%
    pivot_wider(names_from = Theme,
                values_from = n,
                values_fill = FALSE) %>%
    relocate(
      Study,
      Oxygenation,
      Ventilation,
      `Respiratory arrest`,
      `Airway patency`,
      `Dyspnea/work of breathing`,
      `Device-specific criteria`,
      Neurologic,
      Hemodynamic,
      Trajectory,
      Radiologic
    )

sets <- names(upset_dat)[-1]
data_t <- as.data.frame(upset_dat)[,-1]

upset(data_t, sets = sets, 
      main.bar.color = "grey50", 
      matrix.color = "black",
      order.by = "freq", # Order by frequency of intersections
      mb.ratio = c(0.6,0.4),
      keep.order = FALSE,
      show.numbers = FALSE,
      sets.x.label = "Trials with criteria",
      mainbar.y.label = "Trials with each criteria combination")    
  

########## Network Meta-analysis with Study-level Covariates ################

########## Load outcome data and join criteria data ######################### 

gs4_deauth()
sheet_url <- "https://docs.google.com/spreadsheets/d/18ODjbVkKJwn0yn1NBx8c93q2SLTfiAdWdgbmwi_X9uA/edit?gid=1921780000#gid=1921780000"
general_data <- read_sheet(sheet_url, 
                           sheet = "General Data Extraction - Aggregated")

outcome_data <- read_sheet(sheet_url, 
                           sheet = "Outcomes - Aggregated",
                           col_types = "cccddccddccccccccccccccccccccc") %>%
  filter(!is.na(Study)) %>%
  filter(!(Study == "Added studies")) %>%
  # filter this until we sort out the usual care intervention classification
  filter(!(Study == "Arabi, 2022")) %>%
  # same with this study
  filter(!(Study == "Elkheshen, 2024"))

# need to make hierarchy of interventions to ensure same direction for all comparisons
interventions <- factor(outcome_data$Intervention)
unique(interventions)
order_interventions <- rev(c("H Bilevel",
                         "FM Bilevel",
                         "H CPAP",
                         "FM CPAP",
                         "HFNC",
                         "SOT"))

interventions_simplified <- factor(outcome_data$`Intervention (simplified)`)
unique(interventions_simplified)
order_interventions_simplified <- rev(c("Bilevel NIPPV",
                         "CPAP",
                         "HFNC",
                         "SOT"))

outcome_data$Intervention <- factor(outcome_data$Intervention,
                                    levels = order_interventions,
                                    labels = order_interventions,
                                    ordered = T)

outcome_data$`Intervention (simplified)` <- factor(outcome_data$`Intervention (simplified)`,
                                    levels = order_interventions_simplified,
                                    labels = order_interventions_simplified,
                                    ordered = T)
outcome_data <- 
  outcome_data %>%
    arrange(Study, Intervention)


# join criteria data

criteria_by_study <- 
  criteria_data %>%
  mutate(Theme = ifelse(
    Theme == "Respiratory" & (!is.na(Subtheme)),
    Subtheme, Theme)) %>%
  count(Study, Theme) %>%
  mutate( n = ifelse(n>0,1,0)) %>%
  pivot_wider(names_from = Theme,
              values_from = n,
              values_fill = 0)


outcome_data$Study[!(outcome_data$Study %in% criteria_by_study$Study)]
# labels match if this vector is empty:
criteria_by_study$Study[!(criteria_by_study$Study %in% outcome_data$Study)]

outcome_data_criteria <-
  outcome_data %>%
    mutate(Criteria = ifelse(Study %in% criteria_by_study$Study,
                             1, 0)) %>%
    left_join(criteria_by_study,
              by = "Study") %>%
    mutate_at(vars(Oxygenation:`Respiratory arrest`), ~replace_na(., 0))
  
################ Network plots ###########################

# mortality with each comparison on a line (good for network plots or some analyses)
data_mortality <- 
  outcome_data_criteria %>%
    filter(!is.na(`Mortality (n)`)) %>%
    select(Study:`Mortality (n)`, 
           Criteria:`Respiratory arrest`) %>%
  rename(events = `Mortality (n)`,
         total = `Number of participants`) %>%
  group_by(Study) %>%
  summarise(
    comparison = combn(`Intervention (simplified)`, 2, simplify = FALSE),
    events1 = combn(events, 2, FUN = function(x) x[1], simplify = TRUE),
    events2 = combn(events, 2, FUN = function(x) x[2], simplify = TRUE),
    total1 = combn(total, 2, FUN = function(x) x[1], simplify = TRUE),
    total2 = combn(total, 2, FUN = function(x) x[2], simplify = TRUE),
    Criteria = max(Criteria)
  ) %>%
  mutate(
    intervention1 = map_chr(comparison, ~ .x[1]),  # Extract first intervention in pair
    intervention2 = map_chr(comparison, ~ .x[2])
  ) %>%
  ungroup() %>%
  select(-comparison) %>%
  filter(intervention1 != intervention2) %>%
  relocate(Study, intervention1, intervention2) %>%
  mutate(intervention1 = factor(intervention1, 
                                levels = 1:4,
                                labels = order_interventions_simplified),
         intervention2 = factor(intervention2, 
                                levels = 1:4,
                                labels = order_interventions_simplified))

# intubation with each comparison on a line (good for network plots or some analyses)

data_intubation <-
  outcome_data_criteria %>%
    filter(!is.na(`Intubation (n)`)) %>%
    select(Study, `Intervention (simplified)`, `Number of participants - intubation outcome`,
           `Intubation (n)`, Criteria:`Respiratory arrest`) %>%
  rename(events = `Intubation (n)`,
         total = `Number of participants - intubation outcome`) %>%
  group_by(Study) %>%
  summarise(
    comparison = combn(`Intervention (simplified)`, 2, simplify = FALSE),
    events1 = combn(events, 2, FUN = function(x) x[1], simplify = TRUE),
    events2 = combn(events, 2, FUN = function(x) x[2], simplify = TRUE),
    total1 = combn(total, 2, FUN = function(x) x[1], simplify = TRUE),
    total2 = combn(total, 2, FUN = function(x) x[2], simplify = TRUE),
    Criteria = max(Criteria)
  ) %>%
  mutate(
    intervention1 = map_chr(comparison, ~ .x[1]),  # Extract first intervention in pair
    intervention2 = map_chr(comparison, ~ .x[2])
  ) %>%
  ungroup() %>%
  select(-comparison) %>%
  filter(intervention1 != intervention2) %>%
  relocate(Study, intervention1, intervention2) %>%
  mutate(intervention1 = factor(intervention1, 
                                levels = 1:4,
                                labels = order_interventions_simplified),
         intervention2 = factor(intervention2, 
                                levels = 1:4,
                                labels = order_interventions_simplified))



## Full network

# make list of edges

plot_network <- function(df_pairs){

# Count occurrences of each edge (pair)
edge_counts <- df_pairs %>%
  count(intervention1, intervention2, Criteria)

# Create the graph object
g <- graph_from_data_frame(edge_counts, directed = FALSE)

# Add weights as edge attributes
E(g)$weight <- edge_counts$n
#E(g)$weight <- 1:10

E(g)$color <- ifelse(edge_counts$Criteria == 1,
                     c_ice,
                     c_light)

E(g)$label.x <- c(0.5,
                  0.3,
                  -0.5,
                  -0.5,
                  0.3,
                  0.5,
                  -0.5,
                  -0.1,
                  0.1)  # Replace with specific x-coordinates
E(g)$label.y <- c(0.7,0.5,
                  0.2, -0.2,
                  -0.5, -0.7, 
                  0.5, -0.5, 
                  -0.5)  # Replace with specific y-coordinates

# Plot the network
plot(
  g,
  layout = layout_in_circle,         # Circular layout for rounded arrangement
  vertex.color = "grey90",             # Grey vertices
  vertex.size = 40,
  vertex.label.color = "black",      # Black vertex labels
  vertex.label.cex = 1.2,            # Adjust vertex label size
  edge.width = 5+E(g)$weight,          # Edge width proportional to count
  edge.label = E(g)$weight,          # Show weights on edges
  edge.alpha = 0.5,
  edge.label.color = "black",        # Black edge labels
  edge.label.cex = 1,                # Adjust edge label size
)
}

plot_network(data_mortality)
plot_network(data_intubation)

# for dyspnea criteria:

# intubation with each comparison on a line (good for network plots or some analyses)

data_intubation2 <-
  outcome_data_criteria %>%
  filter(!is.na(`Intubation (n)`)) %>%
  select(Study, `Intervention (simplified)`, `Number of participants - intubation outcome`,
         `Intubation (n)`, Criteria:`Respiratory arrest`) %>%
  rename(events = `Intubation (n)`,
         total = `Number of participants - intubation outcome`) %>%
  group_by(Study) %>%
  summarise(
    comparison = combn(`Intervention (simplified)`, 2, simplify = FALSE),
    events1 = combn(events, 2, FUN = function(x) x[1], simplify = TRUE),
    events2 = combn(events, 2, FUN = function(x) x[2], simplify = TRUE),
    total1 = combn(total, 2, FUN = function(x) x[1], simplify = TRUE),
    total2 = combn(total, 2, FUN = function(x) x[2], simplify = TRUE),
    Criteria = max(`Dyspnea/work of breathing`)
  ) %>%
  mutate(
    intervention1 = map_chr(comparison, ~ .x[1]),  # Extract first intervention in pair
    intervention2 = map_chr(comparison, ~ .x[2])
  ) %>%
  ungroup() %>%
  select(-comparison) %>%
  filter(intervention1 != intervention2) %>%
  relocate(Study, intervention1, intervention2) %>%
  mutate(intervention1 = factor(intervention1, 
                                levels = 1:4,
                                labels = order_interventions_simplified),
         intervention2 = factor(intervention2, 
                                levels = 1:4,
                                labels = order_interventions_simplified))

# Count occurrences of each edge (pair)
edge_counts <- data_intubation2 %>%
  count(intervention1, intervention2, Criteria)

# Create the graph object
g <- graph_from_data_frame(edge_counts, directed = FALSE)

# Add weights as edge attributes
E(g)$weight <- edge_counts$n
#E(g)$weight <- 1:10

E(g)$color <- ifelse(edge_counts$Criteria == 1,
                     c_ice,
                     c_light)

E(g)$label.x <- c(0.5,
                  0.3,
                  -0.5,
                  -0.5,
                  0.3,
                  0.5,
                  -0.5,
                  -0.1,
                  0.1,
                  -0.5)  # Replace with specific x-coordinates
E(g)$label.y <- c(0.7,0.5,
                  0.1, -0.1,
                  -0.5, -0.7, 
                  0.7, -0.5, 
                  -0.5,
                  0.3)  # Replace with specific y-coordinates

# Plot the network
plot(
  g,
  layout = layout_in_circle,         # Circular layout for rounded arrangement
  vertex.color = "grey90",             # Grey vertices
  vertex.size = 40,
  vertex.label.color = "black",      # Black vertex labels
  vertex.label.cex = 1.2,            # Adjust vertex label size
  edge.width = 5+E(g)$weight,          # Edge width proportional to count
  edge.label = E(g)$weight,          # Show weights on edges
  edge.alpha = 0.5,
  edge.label.color = "black",        # Black edge labels
  edge.label.cex = 1,                # Adjust edge label size
)


#################### Bayesian Network Meta-Analysis ######################

#################### With covariate for criteria (any) #####################

# Model - Intubation

df_ett <-
  outcome_data_criteria %>%
  filter(!is.na(`Intubation (n)`)) %>%
  select(Study, `Intubation (n)`,
         `Number of participants - intubation outcome`,
         `Intervention (simplified)`,
         Criteria:`Respiratory arrest`) %>%
  rename(events = `Intubation (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants - intubation outcome`) %>%
  select(Study, events, sample_size, treatment, Criteria)

# Set a reference treatment 
df_ett$treatment <- relevel(as.factor(as.character(df_ett$treatment)), ref = "SOT")

# Define the Bayesian NMA model
nma_model_ett <- brm(
  events | trials(sample_size) ~ 0 + 
    treatment*Criteria + 
    (1 | Study) + 
    # treat Study as fixed effect ie not hierarchical
    # this model assumes that between-study variance is same no matter which options are being compared
    # check model in gemtc to make sure similar
    # plot random effects
    # consistency of the network, checking that direct and indirect estimates match (node splitting approach)
    # applying the grades of evidence to the comparisons
    (0 + treatment | Study),  # No intercept, treatment effects
  family = binomial(link = "logit"),                           # Logistic regression for binary data
  data = df_ett,
  prior = c(
    prior(normal(0, 1), class = "b"),  # Priors for treatment effects
    prior(exponential(1), class = "sd")  # Prior for between-study heterogeneity
  ),
  iter = 4000, warmup = 1000, chains = 4, cores = 4
)

# Summarize results
summary(nma_model_ett)

# Plot posterior distributions of treatment effects
plot(nma_model_ett)

post_ett <- as_draws_df(nma_model_ett)

exp(mean_CrI(post_ett$b_Criteria))
exp(mean_CrI(post_ett$`b_treatmentBilevelNIPPV:Criteria`))
exp(mean_CrI(post_ett$`b_treatmentCPAP:Criteria`))
exp(mean_CrI(post_ett$`b_treatmentHFNC:Criteria`))

df <-
  outcome_data_criteria %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study:`Mortality (n)`, 
         Criteria:`Respiratory arrest`) %>%
  rename(events = `Mortality (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants`) %>%
  select(Study, events, sample_size, treatment, Criteria)

# Set a reference treatment 
df$treatment <- relevel(as.factor(as.character(df$treatment)), ref = "SOT")

# Define the Bayesian NMA model
nma_model <- brm(
  events | trials(sample_size) ~ 0 + 
    treatment*Criteria +
    (1 | Study) + 
    (0 + treatment | Study),  # No intercept, treatment effects
  family = binomial(link = "logit"),                           # Logistic regression for binary data
  data = df,
  prior = c(
    prior(normal(0, 1), class = "b"),  # Priors for treatment effects
    prior(exponential(1), class = "sd")  # Prior for between-study heterogeneity
  ),
  iter = 4000, warmup = 1000, chains = 4, cores = 4
)

# Summarize results
summary(nma_model)
post_mort <- as_draws_df(nma_model)

exp(mean_CrI(post_mort$b_Criteria))
exp(mean_CrI(post_mort$`b_treatmentBilevelNIPPV:Criteria`))
exp(mean_CrI(post_mort$`b_treatmentCPAP:Criteria`))
exp(mean_CrI(post_mort$`b_treatmentHFNC:Criteria`))

# Dyspnea / work of breathing

# Model - Intubation

df_ett <-
  outcome_data_criteria %>%
  filter(!is.na(`Intubation (n)`)) %>%
  select(Study, `Intubation (n)`,
         `Number of participants - intubation outcome`,
         `Intervention (simplified)`,
         `Dyspnea/work of breathing`) %>%
  rename(events = `Intubation (n)`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants - intubation outcome`,
         Criteria = `Dyspnea/work of breathing`) %>%
  select(Study, events, sample_size, treatment, Criteria)

# Set a reference treatment 
df_ett$treatment <- relevel(as.factor(as.character(df_ett$treatment)), ref = "SOT")

# Define the Bayesian NMA model
nma_model_ett <- brm(
  events | trials(sample_size) ~ 0 + 
    treatment*Criteria + 
    (1 | Study) + 
    (0 + treatment | Study),  # No intercept, treatment effects
  family = binomial(link = "logit"),                           # Logistic regression for binary data
  data = df_ett,
  prior = c(
    prior(normal(0, 1), class = "b"),  # Priors for treatment effects
    prior(exponential(1), class = "sd")  # Prior for between-study heterogeneity
  ),
  iter = 4000, warmup = 1000, chains = 4, cores = 4
)

# Summarize results
summary(nma_model_ett)

# Plot posterior distributions of treatment effects
plot(nma_model_ett)

post_ett <- as_draws_df(nma_model_ett)

exp(mean_CrI(post_ett$b_Criteria))
exp(mean_CrI(post_ett$`b_treatmentBilevelNIPPV:Criteria`))
exp(mean_CrI(post_ett$`b_treatmentCPAP:Criteria`))
exp(mean_CrI(post_ett$`b_treatmentHFNC:Criteria`))

df <-
  outcome_data_criteria %>%
  filter(!is.na(`Mortality (n)`)) %>%
  select(Study:`Mortality (n)`, 
         `Dyspnea/work of breathing`) %>%
  rename(events = `Mortality (n)`,
         Criteria = `Dyspnea/work of breathing`,
         treatment = `Intervention (simplified)`,
         sample_size = `Number of participants`) %>%
  select(Study, events, sample_size, treatment, Criteria)

# Set a reference treatment 
df$treatment <- relevel(as.factor(as.character(df$treatment)), ref = "SOT")

# Define the Bayesian NMA model
nma_model <- brm(
  events | trials(sample_size) ~ 0 + 
    treatment*Criteria +
    (1 | Study) + 
    (0 + treatment | Study),  # No intercept, treatment effects
  family = binomial(link = "logit"),                           # Logistic regression for binary data
  data = df,
  prior = c(
    prior(normal(0, 1), class = "b"),  # Priors for treatment effects
    prior(exponential(1), class = "sd")  # Prior for between-study heterogeneity
  ),
  iter = 4000, warmup = 1000, chains = 4, cores = 4
)

# Summarize results
summary(nma_model)
post_mort <- as_draws_df(nma_model)

exp(mean_CrI(post_mort$b_Criteria))
exp(mean_CrI(post_mort$`b_treatmentBilevelNIPPV:Criteria`))
exp(mean_CrI(post_mort$`b_treatmentCPAP:Criteria`))
exp(mean_CrI(post_mort$`b_treatmentHFNC:Criteria`))




# Oxygenation dependent on FiO2 vs not dependent on O2




# Oxygen saturation vs PF ratio
# Accessory muscle use or work of breathing vs not
# Neurologic criteria with GCS and without GCS
# respiratory rate (?dichotomized)
# pH (?dichotomized)
# Hemodynamics using vasopressors or not using vasopressors?
# Criteria with and without secretions



