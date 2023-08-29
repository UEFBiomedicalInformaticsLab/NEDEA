library(fgsea)

data(examplePathways)
data(exampleRanks)
set.seed(42)

examplePathways <- examplePathways[lengths(examplePathways) > 500]

fgseaRes_15_100 <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 100)

fgseaRes_15_500 <- fgsea(pathways = examplePathways, 
                         stats    = exampleRanks,
                         minSize  = 15,
                         maxSize  = 500)

fgseaRes_15_1000 <- fgsea(pathways = examplePathways, 
                         stats    = exampleRanks,
                         minSize  = 15,
                         maxSize  = 1000)


fgseaRes_15_Inf <- fgsea(pathways = examplePathways, 
                          stats    = exampleRanks,
                          minSize  = 15,
                          maxSize  = Inf)


names(examplePathways)[!names(examplePathways) %in% fgseaRes_15_1000$pathway]

lapply(examplePathways, function(x){length(intersect(x, names(exampleRanks)))})
