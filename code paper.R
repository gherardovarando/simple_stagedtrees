library(bnlearn)
library(magrittr)
library(arules)
library(deal)

### FUNCTIONS ###

source("functions.R")
## CORONARY ##
colnames(coronary) <- c("S","MW","PW","PRE","PRO","F")
dag_cor <- hc(coronary)
graphviz.plot(dag_cor)

simple_cor <- simple_marginal(coronary, order = c("S","PW","PRE","MW","F","PRO"))
bhc_cor <- stages_bhc(full(coronary,order = c("S","PW","PRE","MW","F","PRO"), join_unobserved = F,lambda = 1 ))
simplified_cor <- simplify(bhc_cor)
total_cor <- simple_total_bhc(full(coronary, order = c("S","PW","PRE","MW","F","PRO"), join_unobserved = F, lambda =1))


bic <- data.frame(BN = -2*BIC(dag_cor,coronary),Marginal = BIC(simple_cor),Simplified = BIC(simplified_cor), Total = BIC(total_cor), BHC = BIC(bhc_cor))
pos <- data.frame(Marginal = num_positions(simple_cor), Simplified = num_positions(simplified_cor), Total = num_positions(total_cor), BHC = num_positions(bhc_cor), Ratio = num_stages(bhc_cor)/num_positions(bhc_cor))


### CHDS ###

# READ FILE 
load("chds.RData")
dag_chds <- hc(chds)
graphviz.plot(dag_chds)

simple_chds <- simple_marginal(chds,  order = c("Social","Economic","Events","Admission"))
bhc_chds <- stages_bhc(full(chds, order = c("Social","Economic","Events","Admission"), join_unobserved = F,lambda = 1 ))
simplified_chds <- simplify(bhc_chds)
total_chds <- simple_total_bhc(full(chds, order = c("Social","Economic","Events","Admission"), join_unobserved = F, lambda =1))

bic <- rbind(bic,c(-2*BIC(dag_chds,chds),BIC(simple_chds),BIC(simplified_chds),BIC(total_chds),BIC(bhc_chds)))
pos <- rbind(pos,c(num_positions(simple_chds),num_positions(simplified_chds),num_positions(total_chds),num_positions(bhc_chds),num_stages(bhc_chds)/num_positions(bhc_chds)))


### PHD ###
dag_phd <- hc(PhDArticles)
graphviz.plot(dag_phd)

simple_phd <- simple_marginal(PhDArticles,  order = c("Articles","Mentor","Prestige","Gender", "Kids","Married"))
bhc_phd <- stages_bhc(full(PhDArticles,  order = c("Articles","Mentor","Prestige","Gender", "Kids","Married"), join_unobserved = F,lambda = 1 ))
simplified_phd <- simplify(bhc_phd)
total_phd <- simple_total_bhc(full(PhDArticles,  order = c("Articles","Mentor","Prestige","Gender", "Kids","Married"), join_unobserved = F, lambda =1))

bic <- rbind(bic,c(-2*BIC(dag_phd,PhDArticles),BIC(simple_phd),BIC(simplified_phd),BIC(total_phd),BIC(bhc_phd)))
pos <- rbind(pos,c(num_positions(simple_phd),num_positions(simplified_phd),num_positions(total_phd),num_positions(bhc_phd),num_stages(bhc_phd)/num_positions(bhc_phd)))


### Cachexia ###
# READ FILE
data <- cachexia_data
data <- cbind(discretizeDF(data[,-7]),MC = as.factor(data$MC))
dag_cac <- hc(data)
graphviz.plot(dag_cac)

simple_cac <- simple_marginal(data,  order = c("F","GC","MC","V", "GM","A","B"))
bhc_cac <- stages_bhc(full(data,  order = c("F","GC","MC","V", "GM","A","B"), join_unobserved = F,lambda = 1 ))
simplified_cac <- simplify(bhc_cac)
total_cac <- simple_total_bhc(full(data,  order = c("F","GC","MC","V", "GM","A","B"), join_unobserved = F, lambda =1))

bic <- rbind(bic,c(-2*BIC(dag_cac,data),BIC(simple_cac),BIC(simplified_cac),BIC(total_cac),BIC(bhc_cac)))
pos <- rbind(pos,c(num_positions(simple_cac),num_positions(simplified_cac),num_positions(total_cac),num_positions(bhc_cac),num_stages(bhc_cac)/num_positions(bhc_cac)))

## MathMarks
# Load Data
mathmarks <- discretizeDF(mathmarks)
dag_math <- hc(mathmarks)
graphviz.plot(dag_math)
order <- c("mechanics","vectors","algebra","analysis","statistics")

simple_math <- simple_marginal(mathmarks,  order = order)
bhc_math <- stages_bhc(full(mathmarks,  order = order, join_unobserved = F,lambda = 1 ))
simplified_math <- simplify(bhc_math)
total_math <- simple_total_bhc(full(mathmarks,  order = order, join_unobserved = F, lambda =1))

bic <- rbind(bic,c(-2*BIC(dag_math,mathmarks),BIC(simple_math),BIC(simplified_math),BIC(total_math),BIC(bhc_math)))
pos <- rbind(pos,c(num_positions(simple_math),num_positions(simplified_math),num_positions(total_math),num_positions(bhc_math),num_stages(bhc_math)/num_positions(bhc_math)))


## Titanic

titanic.df <- as.data.frame(Titanic)
titanic.df <- titanic.df[rep(row.names(titanic.df), titanic.df$Freq), 1:4]
dag_tit <- hc(titanic.df)
graphviz.plot(dag_tit)

simple_tit <- simple_marginal(titanic.df,  order = c("Class","Sex","Survived","Age"))
bhc_tit <- stages_bhc(full(titanic.df,  order = c("Class","Sex","Survived","Age"), join_unobserved = F,lambda = 1 ))
simplified_tit <- simplify(bhc_tit)
total_tit <- simple_total_bhc(full(titanic.df,  order = c("Class","Sex","Survived","Age"), join_unobserved = F, lambda =1))

bic <- rbind(bic,c(-2*BIC(dag_tit,titanic.df),BIC(simple_tit),BIC(simplified_tit),BIC(total_tit),BIC(bhc_tit)))
pos <- rbind(pos,c(num_positions(simple_tit),num_positions(simplified_tit),num_positions(total_tit),num_positions(bhc_tit),num_stages(bhc_tit)/num_positions(bhc_tit)))

## Asia
dag_asia <- hc(asia)
graphviz.plot(dag_asia)

order <- c("A","T","S","L","E","B","X","D")
simple_asia <- simple_marginal(asia,  order = order)
bhc_asia <- stages_bhc(full(asia,  order = order, join_unobserved = F,lambda = 1 ))
simplified_asia <- simplify(bhc_asia)
total_asia <- simple_total_bhc(full(asia,  order =order, join_unobserved = F, lambda =1))

bic <- rbind(bic,c(-2*BIC(dag_asia,asia),BIC(simple_asia),BIC(simplified_asia),BIC(total_asia),BIC(bhc_asia)))
pos <- rbind(pos,c(num_positions(simple_asia),num_positions(simplified_asia),num_positions(total_asia),num_positions(bhc_asia),num_stages(bhc_asia)/num_positions(bhc_asia)))

## FallEld
fall <- readRDS("FallEld.rds")

dag_fall <- hc(fall)
graphviz.plot(dag_fall)
order <- c("A","T","answer","R")

simple_fall <- simple_marginal(fall,  order = order)
bhc_fall <- stages_bhc(full(fall,  order = order, join_unobserved = F,lambda = 1 ))
simplified_fall <- simplify(bhc_fall)
total_fall <- simple_total_bhc(full(fall,  order =order, join_unobserved = F, lambda =1))

bic <- rbind(bic,c(-2*BIC(dag_fall,fall),BIC(simple_fall),BIC(simplified_fall),BIC(total_fall),BIC(bhc_fall)))
pos <- rbind(pos,c(num_positions(simple_fall),num_positions(simplified_fall),num_positions(total_fall),num_positions(bhc_fall),num_stages(bhc_fall)/num_positions(bhc_fall)))

## KSL ##
data(ksl)
ksl[,1] <- discretize(ksl[,1],breaks=2)
ksl[,2] <- discretize(ksl[,2],breaks=2)
ksl[,3] <- as.factor(ksl[,3])
ksl[,4] <- discretize(ksl[,4],breaks=2)

dag_ksl <- hc(ksl)
graphviz.plot(dag_ksl)
order <- c("Smok","FEV","Sex","Work", "Year","Alc","Kol","Hyp","logBMI")

simple_ksl <- simple_marginal(ksl,  order = order)
bhc_ksl <- stages_bhc(full(ksl,  order = order, join_unobserved = F,lambda = 1 ))
simplified_ksl <- simplify(bhc_ksl)
total_ksl <- simple_total_bhc(full(ksl,  order =order, join_unobserved = F, lambda =1))

bic <- rbind(bic,c(-2*BIC(dag_ksl,ksl),BIC(simple_ksl),BIC(simplified_ksl),BIC(total_ksl),BIC(bhc_ksl)))
pos <- rbind(pos,c(num_positions(simple_ksl),num_positions(simplified_ksl),num_positions(total_ksl),num_positions(bhc_ksl),num_stages(bhc_ksl)/num_positions(bhc_ksl)))

## TABLE CONSTRUCTION ##
rownames(bic) <- rownames(pos) <- c("coronary","chds","phd","cachexia","mathmarks","titanic","asia","falleld","ksl")
bic
pos

## CHDS PLOTS
plot(dag_chds)
plot(total_chds)
ceg.plot(total_chds)

## CORONARY PLOTS
plot(dag_cor)
plot(total_cor)
ceg.plot(total_cor)