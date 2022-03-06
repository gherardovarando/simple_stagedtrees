library("ggplot2")
library("reshape2")

qs <- c(0.2, 0.5, 0.7)
Ns <- c(25, 50, seq(from = 100, to = 1000, by = 100), 
        1250, 1500, 1750, 2000)
ns <- c(6)
comb <- expand.grid(ns,  Ns, qs)
res <- lapply(1:nrow(comb), function(i){
  q <- comb[i,3]
  n <- comb[i,1]
  N <- comb[i,2]
  fnam <- paste0("sim_results/", n, "_", N, "_", q , ".rds")
  print(fnam)
  tmp <- as.data.frame(readRDS(fnam))
  #tmp <- as.data.frame(t(colMeans(tmp)))
  tmp$q <- q
  tmp$n <- n 
  tmp$N <- N
  tmp
})

cols <- palette.colors(5, palette = "Classic Tableau")
names(cols) <- NULL
zalpha <- qnorm(0.975) # for 0.95 CI
results <- aggregate(value ~ q + n + N + variable,
                     data = melt(res, id.vars = c("q", "n", "N")), 
                     FUN = function(x) c(mean = mean(x), sd = sd(x)))
results <- do.call(data.frame, results)

pl <- ggplot(results, 
             aes(x = N, y = value.mean,
                 col = variable, 
                 fill = variable,
                 #linetype = variable,
                 group = variable,
                 ymin = value.mean - zalpha*value.sd/sqrt(100),
                 ymax = value.mean + zalpha*value.sd/sqrt(100))) + 
  geom_line() + 
  geom_ribbon(alpha = 0.4) + 
  facet_grid(rows = vars(q), scales = "free") +  
  ylab("normalized stage hamming distance") + theme_bw() +
  scale_y_log10()+
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  scale_color_manual(values =  cols,
                     labels = c("marginal", "total", "BHC", "simplified"))+ 
  scale_fill_manual(values =  cols,
                     labels = c("marginal", "total", "BHC", "simplified"))
  
pl
ggsave("plot_sim_order.pdf", plot = pl, width = 4, height = 5)

