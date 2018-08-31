library(destiny) 
data(guo)

dm <- DiffusionMap(guo)
dpt <- DPT(dm)
example_obj <- dpt
save(example_obj, file="/home/rstudio/shared/example_obj.rda")

plot(dpt)
branches <- dpt@branch
dpt@branch
apply(dpt@tips,2,function(x){sum(na.omit(x))})
names(dpt)
dpt$Tspan8
length(dpt$Actb)
cell_names <- row.names(branches)

dpt <- branch_divide(dpt,c(1))
max(dpt@branch, na.rm = T)
hasStrDPT <- function(x) {
  !grepl("DPT", x, fixed = T)
}
dpt$Branch
notDPT <- names(dpt)[!grepl("DPT",names(dpt),fixed=T)]
cell_names <- notDPT[notDPT != "Cell" & notDPT != 'num_cells']
hasStrDPT("DPT1")
strsplit("DPT1", "DPT1")
if (is.null(cell_names)){
  cell_names <- 1:dim(branches)[1]
  row.names(branches) <- cell_names
}

nodes <- c()
g <- igraph::make_empty_graph(directed=F)
stem_name <- "stem"
nodes <- c(nodes, stem_name)
g <- igraph::add_vertices(g,1, name=stem_name)

for (branch in colnames(branches)){
  node_name <- paste(c("branch_", branch), collapse = "")
  nodes <- c(nodes, stem_name)
  g <- igraph::add_vertices(g, 1, name=node_name)
  g <- igraph::add_edges(g, c(stem_name, node_name))
}

edges <- igraph::as_edgelist(g)
edgeIds <- c()
nodeIds1 <- c()
nodeIds2 <- c()
# CM stands for cell mapping, its field in the common json format.
cellIdCM <- c()
edgeIdCM <- c()
pseudotimeCM <- c()
for (row in 1:(dim(edges)[1])){
  node1 <- edges[row,1]
  nodeIds1 <- c(nodeIds1, node1)
  node2 <- edges[row,2]
  nodeIds2 <- c(nodeIds2, node1)
  edge_id <- paste(c(node1, node2), collapse="_")
  edgeIds <- c(edgeIds, edge_id)
  branchId <- as.numeric(tail(strsplit(node2, "Branch")[[1]],1))
  branchId <- paste(c("Branch", branchId), collapse="")
  n_branch_cells <- sum(!is.na(branches[,branchId]))
  edgeIdCM <- c(edgeIdCM, rep(edge_id, n_branch_cells))
  branch_v <- branches[,branchId][!is.na(branches[,branchId])]
  cellIdCM <- c(cellIdCM, names(branch_v))
  pseudotimeCM <- c(pseudotimeCM, branch_v)
  
}

output <- list(
  nodes= list(nodeId=nodes),
  egdes= list(
    edgeId=edgeIds,
    nodeId1= edges[,1],
    nodeId2= edges[,2]
  ),
  cellMapping=list(
    cellId= cellIdCM,
    edgeId=edgeIdCM,
    psuedotime=pseudotimeCM
  )
)


plot(g)

colSums(is.na(branches))

row.names(branches)
###############################################333
source("/home/rstudio/shared/dpt_convert.r")
load("/home/rstudio/shared/example_obj.rda")
cxb <- to_cell_x_branch(dpt)
common_list <- to_common_list(dpt)
