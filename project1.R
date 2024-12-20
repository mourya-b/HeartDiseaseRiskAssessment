library(gRain)
library(dagitty)
library(ucimlrepo)
library(fastDummies)
library(CCP)
library(nnet)
#cp, restecg, thal


iris_by_name <- fetch_ucirepo(name = "Heart Disease")
data = iris_by_name[["data"]][["original"]]

df <- as.data.frame(data)
# Specify data types
df$sex <- factor(df$sex, 
                levels = c(0,1), 
                labels = c("female", "male"))
##maybe ordinal
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

# Optional: Label the binary variable for clarity
df$num <- factor(df$num, 
                 levels = c(0, 1), 
                 labels = c("No Heart Disease", "Heart Disease"))

#df$num <- ordered(df$num, 
#                  levels = c(0,1,2,3,4), 
#                  labels = c("No Heart Disease", "Heart Disease - Level 1", 
#                                      "Heart Disease - Level 2", "Heart Disease - Level 3", 
#                                      "Heart Disease - Level 4"))


##remove rows containing missing values
df <- na.omit(df)

g <- dagitty('
dag {
bb="0,0,1,1"
age [pos="0.7,0.1"]         
sex [pos="0.2,0.15"]          
cp [pos="0.55,0.5"]          
trestbps [pos="0.4,0.8"]     
chol [pos="0.2,0.2"]         
fbs [pos="0.35,0.4"]        
restecg [pos="0.45,0.3"]     
thalach [pos="0.85,0.5"]    
exang [pos="0.75,0.35"]      
oldpeak [pos="0.5,0.4"]     
slope [pos="0.5,0.55"]       
ca [pos="0.4,0.2"]           
thal [pos="0.9,0.4"]         
num [pos="0.7,0.7"]          
age -> { cp trestbps chol fbs thalach ca}
sex -> { chol fbs restecg }
trestbps -> { cp restecg }
chol -> { restecg thalach trestbps num}
fbs -> { slope oldpeak trestbps}
slope -> { cp exang oldpeak }
ca -> { exang oldpeak }
thal -> { num }
thalach -> {restecg oldpeak}
exang -> {cp}
num -> {ca cp restecg thalach}
}
')
plot(g)

addEdge <- function(g, from, to) {
  # Extract the current graph definition
  graph_string <- as.character(g)
  updated_graph_string <- sub("\\}", paste0(from, " -> ", to, "\n}"), graph_string)
  g = dagitty(updated_graph_string)
  return(g)
}

g <- graphLayout(g)
impliedConditionalIndependencies(g)
r <- localTests(g,df, type='cis.pillai')
plotLocalTestResults(r)
print(mean(abs(r$estimate)))

#thalach -> slope
g <- addEdge( g, "thalach", "slope")
g <- graphLayout(g)
impliedConditionalIndependencies(g)
r <- localTests(g,df, type='cis.pillai')
plotLocalTestResults(r)
print(mean(abs(r$estimate)))

if (isAcyclic(g)) {
  cat("The graph is acyclic (DAG). No cycles detected.\n")
} else {
  cat("The graph contains cycles. It is not a valid DAG.\n")
}

#sex -> thal
g <- addEdge(g, "sex", "thal")
g <- graphLayout(g)
impliedConditionalIndependencies(g)
r <- localTests(g,df, type='cis.pillai')
plotLocalTestResults(r)
print(mean(abs(r$estimate)))

if (isAcyclic(g)) {
  cat("The graph is acyclic (DAG). No cycles detected.\n")
} else {
  cat("The graph contains cycles. It is not a valid DAG.\n")
}

#num -> exang
g <- addEdge(g, "num", "exang")
g <- graphLayout(g)
impliedConditionalIndependencies(g)
r <- localTests(g,df, type='cis.pillai')
plotLocalTestResults(r)
print(mean(abs(r$estimate)))

if (isAcyclic(g)) {
  cat("The graph is acyclic (DAG). No cycles detected.\n")
} else {
  cat("The graph contains cycles. It is not a valid DAG.\n")
}

# num -> oldpeak
g <- addEdge(g, "num", "oldpeak")
g <- graphLayout(g)
impliedConditionalIndependencies(g)
r <- localTests(g,df, type='cis.pillai')
plotLocalTestResults(r)
print(mean(abs(r$estimate)))

if (isAcyclic(g)) {
  cat("The graph is acyclic (DAG). No cycles detected.\n")
} else {
  cat("The graph contains cycles. It is not a valid DAG.\n")
}

# age -> num
g <- addEdge(g, "age", "num")
g <- graphLayout(g)
impliedConditionalIndependencies(g)
r <- localTests(g,df, type='cis.pillai')
plotLocalTestResults(r)
print(mean(abs(r$estimate)))

if (isAcyclic(g)) {
  cat("The graph is acyclic (DAG). No cycles detected.\n")
} else {
  cat("The graph contains cycles. It is not a valid DAG.\n")
}

##cannonical correlations
r <- c()
for( n in names(g) ){
  for( p in dagitty::parents(g,n) ){
    otherparents <- setdiff( dagitty::parents(g,n), p )
    tst <- ciTest( X=n, Y=p, Z=otherparents, df,
                   type="cis.pillai" )
    r <- rbind( r, data.frame(list(X=p,A="->",Y=n,
                                   cor=tst[,"estimate"],p=tst[,"p.value"])) )
  }
}

r.dagitty <- paste(r$X, r$A, r$Y, "[beta=",signif(r$cor,2),"] ", collapse="\n")
g_with_coefficients <- dagitty( r.dagitty )
coordinates( g_with_coefficients ) <- coordinates( g )
plot( g_with_coefficients, show.coefficients=TRUE )

#filtered_r <- r[r$p <= 0.05 & r$cor > 0.01,]
#r_del <- r[r$p > 0.05 | r$cor <= 0.01,]

#filtered_r_p <- r[r$p <= 0.05,]
#r_del_p <- r[r$p > 0.05,]
#length(filtered_r_p$X)
#length(r_del_p$X)

filtered_r <- r[abs(r$cor) > 0.1,]
r_del <- r[abs(r$cor) <= 0.1,]

#filtered_r_p <- filtered_r[filtered_r$p <= 0.05,]
del_r_p <- filtered_r[filtered_r$p > 0.05,]

filtered_r <- subset(filtered_r, !(X == "slope" & Y == "cp"))
filtered_r <- subset(filtered_r, !(X == "slope" & Y == "exang"))
filtered_r <- subset(filtered_r, !(X == "trestbps" & Y == "restecg"))
filtered_r <- subset(filtered_r, !(X == "fbs" & Y == "slope"))

length(filtered_r$X)
length(r_del$X)

filtered_r.dagitty <- paste(filtered_r$X, filtered_r$A, filtered_r$Y, "[beta=",signif(filtered_r$cor,2),"] ", collapse="\n")
g_with_coefficients <- dagitty( filtered_r.dagitty )
coordinates(g_with_coefficients ) <- coordinates( g )
plot(g_with_coefficients, show.coefficients=TRUE )

##sex -> num hypothesis
adjustmentSets(g_with_coefficients,"sex","num")
#coef(glm(num ~ sex, df, family="binomial" ))
model <- glm(num ~ sex, data = df, family = "binomial")
summary(model)
exp(coef(model))


adjustmentSets(g_with_coefficients,"num","cp")
#model2 <- glm(cp ~ num + trestbps, data = df, family = "multinom")
model2 <- multinom(cp ~ num + trestbps, data = df)
summary(model2)
exp(coef(model2))

model_cp_trestbps <- multinom(cp ~ num + trestbps, data = df)
model_cp_age <- multinom(cp ~ num + age, data = df)


predicted_probs_trestbps <- predict(model_cp_trestbps, type = "probs")

#asymptimatic
df_new = df
df_new$predicted_asymptomatic <- predicted_probs_trestbps[, "asymptomatic"]
df_hd <- df_new[df_new$num == "Heart Disease", ]

mean_asymptomatic_prob <- mean(df_hd$predicted_asymptomatic)
cat("P(cp = asymptomatic | num = Heart Disease, adjusted for trestbps):", mean_asymptomatic_prob)

#atypical anginna
df_new = df
df_new$predicted_atypical <- predicted_probs_trestbps[, "atypical angina"]
df_hd <- df_new[df_new$num == "Heart Disease", ]

mean_atypical_prob <- mean(df_hd$predicted_atypical)
cat("P(cp = atypical | num = Heart Disease, adjusted for trestbps):", mean_atypical_prob)


#typical angina
df_new = df
df_new$predicted_typical <- predicted_probs_trestbps[, "typical angina"]
df_hd <- df_new[df_new$num == "Heart Disease", ]

mean_typical_prob <- mean(df_hd$predicted_typical)
cat("P(cp = typical | num = Heart Disease, adjusted for trestbps):", mean_typical_prob)

#non-anginal pain
df_new = df
df_new$predicted_non_anginal <- predicted_probs_trestbps[, "non-anginal pain"]
df_hd <- df_new[df_new$num == "Heart Disease", ]

mean_non_anginal_prob <- mean(df_hd$predicted_non_anginal)
cat("P(cp = non_anginal | num = Heart Disease, adjusted for trestbps):", mean_non_anginal_prob)

mean_non_anginal_prob + mean_typical_prob + mean_atypical_prob + mean_asymptomatic_prob


