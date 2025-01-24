library(gRain)
library(dagitty)
library(ucimlrepo)
library(fastDummies)
library(CCP)
library(nnet)
library(bnlearn)
library(graph)
library(pcalg)
library(MatchIt)
library(cobalt)
library(igraph)
library(ggraph)
library(gridExtra)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("graph")

iris_by_name <- fetch_ucirepo(name = "Heart Disease")
data = iris_by_name[["data"]][["original"]]

df <- as.data.frame(data)

# Specify data types
df$sex <- factor(df$sex, 
                 levels = c(0,1), 
                 labels = c("female", "male"))

df$cp <- factor(df$cp, levels = c(1, 2, 3, 4), labels = c("typical angina", "atypical angina", "non-anginal pain", "asymptomatic"))

df$fbs <- factor(df$fbs, 
                 levels = c(0,1), 
                 labels = c("fbs <= 120", "fbs > 120"))

df$restecg <- factor(df$restecg, 
                     levels = c(0,1,2), 
                     labels = c("normal", "ST-T wave abnormality", "left ventricular hypertrophy"))

df$exang <- factor(df$exang, 
                   levels = c(0,1), 
                   labels = c("no", "yes"))

df$slope <- factor(df$slope, 
                   levels = c(1,2,3), 
                   labels = c("upsloping", "flat", "downsloping"))

df$thal <- factor(df$thal, 
                  levels = c(3,6,7), 
                  labels = c("normal", "fixed defect", "reversible defect"))

df$num <- ifelse(df$num > 0, 1, 0)

df$num <- factor(df$num, 
                 levels = c(0, 1), 
                 labels = c("No Heart Disease", "Heart Disease"))

df$ca <- factor(df$ca, 
                levels = c(0,1,2,3), 
                labels = c("0", "1", "2", "3"))

df$age <- as.numeric(df$age)
df$trestbps <- as.numeric(df$trestbps)
df$chol <- as.numeric(df$chol)
df$thalach <- as.numeric(df$thalach)

#changed oldpeak from numerical to categorical because it is highly scewed 
df$oldpeak_ord <- cut(df$oldpeak,
                       breaks = c(-Inf, 1, 2, 3, Inf),
                       labels = c("No depression", "Low depression", "Moderate depression", "High depression"),
                      right = TRUE)
#
# # Convert the new variable to a factor
df$oldpeak_ord <- factor(df$oldpeak_ord,
                          levels = c("No depression", "Low depression", "Moderate depression", "High depression"),
                          ordered = TRUE)

# Keep oldpeak intact without introducing NA
df <- df[, !(names(df) %in% c("oldpeak"))]

# Check the cleaned dataset
nrow(df)  # Check the new size of the dataset
df <- na.omit(df)
#
variables <- colnames(df)
from_vars <- setdiff(variables, c("sex", "age"))

#whitelist
whitelist_combined <- data.frame(
  from = c("num", "num"),
  to = c("thalach", "cp")
)

no_incoming_edges <- data.frame(
  from = rep(variables, each = 2),  # All other variables
  to = c("sex", "age")              # No incoming edges to sex and age
)
no_incoming_edges <- no_incoming_edges[no_incoming_edges$from != no_incoming_edges$to, ]
test_vars <- c("thal", "restecg", "ca", "thalach", "oldpeak_ord")  # Test variables
no_outgoing_edges <- data.frame(
  from = rep(test_vars, each = length(variables)),
  to = variables
)
no_outgoing_edges <- no_outgoing_edges[no_outgoing_edges$from != no_outgoing_edges$to, ]
blacklist <- rbind(no_incoming_edges, no_outgoing_edges)

hc_learned <- pc.stable(df, whitelist = whitelist_combined, blacklist = blacklist, alpha = 0.05, undirected = FALSE)

plot(hc_learned)
ig <- as.igraph(hc_learned)

ggraph(ig, layout = "fr") +  
  geom_edge_link(arrow = arrow(length = unit(5, "mm"), type = "closed"), 
                 end_cap = circle(4, "mm"), 
                 color = "grey") +
  geom_node_point(size = 5, color = "skyblue") +
  geom_node_label(aes(label = name), repel = FALSE) +
  theme_void() +
  ggtitle("")

ci.test("num", "cp", data=df, test="x2")

#check assumptions (normally distributed) for all continuous vars on all levels of parent cat. or ord. var
# Plot histogram for 'thalach' conditioned on 'num'
ggplot(df, aes(x = thalach, fill = num)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~ num) +
  theme_minimal() +
  labs(title = "",
       x = "thalach",
       y = "Frequency")

#same for chol given sex 
ggplot(df, aes(x = chol, fill = sex)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~ sex) +
  theme_minimal() +
  labs(title = "",
       x = "chol",
       y = "Frequency")

#for continuous vars without cat parents:
continuous_vars <- c("thalach", "age", "trestbps", "chol")
hist_plots <- lapply(continuous_vars, function(var) {
  ggplot(df, aes_string(x = var)) +
    geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(title = paste("Histogram of", var), x = var, y = "Frequency")
})

# Combine all histograms into one plot
grid.arrange(grobs = hist_plots, ncol = 2)

#note: fbs <= 120 = 244 entries vs 38 for >= 120
fbs_frequency <- table(df$fbs)
restecg_frequency <- table(df$restecg)

#compute adjustment sets
adjustmentSets(hc_learned,"num","cp")

#linear regression
model_lr <- lm(thalach ~ num, data = df)  # Adjust as per identified confounders
summary(model_lr)

# Unadjusted linear regression model BASELINE
baseline_model <- lm(thalach ~ num, data = df)
summary(baseline_model)

#with confounders LR (Hand-crafted model)
model_lr <- lm(thalach ~ num + age + chol, data = df)  
summary(model_lr)

#without confounders LR (Structure Learned)
model_lr <- lm(thalach ~ num, data = df)  
summary(model_lr)

#propensity score matching for data given graph from assignment 1 
#Without Confounders PSM (Hand-crafted model)
ps_model <- glm(num ~ thalach, data = df, family = binomial(link = "logit"))
df$propensity_score <- predict(ps_model, type = "response")
matched_data <- matchit(num ~ age + chol, data = df, method = "nearest", distance = "logit")
matched_data <- match.data(matched_data)
lm_model_matched <- lm(thalach ~ num, data = matched_data)
summary(lm_model_matched)

#With Confounders PSM (Structure Learned)
ps_model <- glm(num ~ age + chol, data = df, family = binomial(link = "logit"))
df$propensity_score <- predict(ps_model, type = "response")
matched_data <- matchit(num ~ age + chol, data = df, method = "nearest", distance = "logit")
matched_data <- match.data(matched_data)
lm_model_matched <- lm(thalach ~ num + age + chol, data = matched_data)
summary(lm_model_matched)




