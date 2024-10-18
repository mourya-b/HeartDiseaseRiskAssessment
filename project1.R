library( gRain )
library(dagitty)
library(ucimlrepo)
library(fastDummies)
library(CCP)
iris_by_name <- fetch_ucirepo(name = "Heart Disease")
data = iris_by_name[["data"]][["original"]]

df <- as.data.frame(data)
# Specify data types
#categorical 
df$sex <- factor(df$sex, 
                levels = c(0,1), 
                labels = c("female", "male"))
##maybe ordinal
df$cp <- factor(df$cp, 
               levels = c(1, 2, 3, 4), 
               labels = c("typical angina", "atypical angina", "non-anginal pain", "asymptomatic"))

df$fbs <- factor(df$fbs, 
                 levels = c(0,1), 
                 labels = c("fbs <= 120", "fbs > 120"))

df$restecg <- factor(df$restecg, 
                 levels = c(0,1,2), 
                 labels = c("normal", "ST-T wave abnormality", "left ventricular hypertrophy"))

df$exang <- factor(df$exang, 
                     levels = c(0,1), 
                     labels = c("no", "yes"))

df$slope <- ordered(df$slope, 
                     levels = c(1,2,3), 
                     labels = c("upsloping", "flat", "downsloping"))

df$thal <- factor(df$thal, 
                     levels = c(3,6,7), 
                     labels = c("normal", "fixed defect", "reversible defect"))

df$num <- ordered(df$num, 
                  levels = c(0,1,2,3,4), 
                  labels = c("No Heart Disease", "Heart Disease - Level 1", 
                                      "Heart Disease - Level 2", "Heart Disease - Level 3", 
                                      "Heart Disease - Level 4"))

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
age -> { cp trestbps chol fbs thalach }
sex -> { trestbps chol fbs restecg thalach exang}
cp -> { num }
trestbps -> { cp restecg exang }
chol -> { restecg cp thalach trestbps }
fbs -> { chol slope exang oldpeak trestbps}
restecg -> { oldpeak thal num}
slope -> { cp exang oldpeak num }
ca -> { exang oldpeak thal num }
thal -> { exang num }
thalach -> {restecg oldpeak}
}
')
plot(g)

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

filtered_r <- r[r$p <= 0.05 & r$cor > 0.1, ]
length(filtered_r$X)

g <- dagitty('
dag {
bb="0,0,1,1"
age [pos="0.7,0.1"]         
cp [pos="0.6,0.55"]          
trestbps [pos="0.4,0.8"]     
chol [pos="0.2,0.2"]         
fbs [pos="0.35,0.4"]        
restecg [pos="0.35,0.3"]     
exang [pos="0.75,0.25"]      
oldpeak [pos="0.5,0.4"]     
slope [pos="0.5,0.55"]       
ca [pos="0.4,0.1"]           
thal [pos="0.9,0.4"]         
num [pos="0.5,0.7"]          
age -> { chol fbs trestbps }
slope -> { cp exang num oldpeak }
thal -> { exang num}
ca -> { num oldpeak thal }
cp -> { num }
chol -> {restecg}
fbs -> { trestbps }
}
')

g <- dagitty('
dag {
bb="0,0,1,1"
age [pos="0.7,0.1"]         
cp [pos="0.4,0.55"]          
trestbps [pos="0.4,0.8"]     
chol [pos="0.2,0.2"]         
fbs [pos="0.35,0.4"]        
restecg [pos="0.35,0.3"]     
exang [pos="0.75,0.25"]      
oldpeak [pos="0.5,0.7"]       
slope [pos="0.5,0.55"]       
ca [pos="0.8,0.65"]           
thal [pos="0.9,0.4"]         
num [pos="0.5,0.4"]        
age -> { chol fbs trestbps }
slope -> { cp exang num oldpeak }
thal -> { exang num}
ca -> { num oldpeak thal }
cp -> { num }
chol -> {restecg}
fbs -> { trestbps }
}
')

filtered_r.dagitty <- paste(filtered_r$X, filtered_r$A, filtered_r$Y, "[beta=",signif(filtered_r$cor,2),"] ", collapse="\n")
g_with_coefficients <- dagitty( filtered_r.dagitty )
coordinates( g_with_coefficients ) <- coordinates( g )
plot( g_with_coefficients, show.coefficients=TRUE )



