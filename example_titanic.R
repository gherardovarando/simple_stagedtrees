### Titanic

model1 <- stages_bhc(full(Titanic))
simplified1 <- simplify(model1)
model1
plot(model1)

model2 <- simple_marginal(Titanic, search = "hclust", k = 2)
model2
plot(model2)

ceg2 <- ceg(model2)
Adj2 <- ceg2adjmat(ceg2)
g2 <- graph_from_adjacency_matrix(Adj2)
plot(g2, layout = layout_as_tree, edge.arrow.size = 0.2)



model3 <- simple_total_bhc(full(Titanic))
plot(model3)
ceg3 <- ceg(model3)
Adj3 <- ceg2adjmat(ceg3)
g3 <- graph_from_adjacency_matrix(Adj3)
plot(g3, layout = layout_as_tree, edge.arrow.size = 0.2)
