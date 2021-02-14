## change to path of FallEld.rds 
fall_eld <- readRDS(file = "FallEld.rds")
order <- c("A", "R", "T", "answer")

model1 <- simple_marginal(data = fall_eld, lambda = 1, order = order, search = "bhc") %>% stndnaming(uniq = TRUE) 
plot(model1)

ceg1 <- ceg(model1)
Adj1 <- ceg2adjmat(ceg1)
g1 <- graph_from_adjacency_matrix(Adj1)
plot(g1, layout = layout_as_tree, edge.arrow.size = 0.2)


model2 <- simple_total_bhc(full(fall_eld, order = order))

plot(model2)

BIC(model1, model2)

ceg2 <- ceg(model2)
Adj2 <- ceg2adjmat(ceg2)
g2 <- graph_from_adjacency_matrix(Adj2)
plot(g2, layout = layout_as_tree, edge.arrow.size = 0.2)


