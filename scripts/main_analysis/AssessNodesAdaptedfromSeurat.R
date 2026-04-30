library(MazamaCoreUtils)
library(pbapply)

# Set a default value if an object is null
#
# @param x An object to set if it's null
# @param default The value to provide if x is null
#
# @return default if x is null, else x
#
SetIfNull <- function(x, default) {
  if(is.null(x = x)){
    return(default)
  } else {
    return(x)
  }
}

DFT <- function(
  tree,
  node,
  path = NULL,
  include.children = FALSE,
  only.children = FALSE
) {
  if (only.children) {
    include.children = TRUE
  }
  children <- which(x = tree$edge[, 1] == node)
  child1 <- tree$edge[children[1], 2]
  child2 <- tree$edge[children[2], 2]
  if (child1 %in% tree$edge[, 1]) {
    if (!only.children) {
      path <- c(path, child1)
    }
    path <- DFT(
      tree = tree,
      node = child1,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <- c(path, child1)
    }
  }
  if (child2 %in% tree$edge[, 1]) {
    if (!only.children) {
      path <- c(path, child2)
    }
    path <- DFT(
      tree = tree,
      node = child2,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <- c(path, child2)
    }
  }
  return(path)
}

GetAllInternalNodes <- function(tree) {
  return(c(tree$edge[1, 1], DFT(tree = tree, node = tree$edge[1, 1])))
}

#' Build Random Forest Classifier
#'
#' Train the random forest classifier
#'
#'
#' @param object Seurat object on which to train the classifier
#' @param training.genes Vector of genes to build the classifier on
#' @param training.classes Vector of classes to build the classifier on
#' @param verbose Additional progress print statements
#' @param ... additional parameters passed to ranger
#'
#' @return Returns the random forest classifier
#'
#' @import Matrix
#'
#' @export
#'
#' @examples
#' pbmc_small
#' # Builds the random forest classifier to be used with ClassifyCells
#' # Useful if you want to use the same classifier with several sets of new data
#' classifier <- BuildRFClassifier(pbmc_small, training.classes = pbmc_small@ident)
#'
BuildRFClassifier <- function(
  object,
  training.genes = NULL,
  training.classes = NULL,
  verbose = TRUE,
  ...
) {
  #check('ranger')
  training.classes <- as.vector(x = training.classes)
  training.genes <- SetIfNull(x = training.genes, default = rownames(x = object@assays$RNA@data))
  training.data <- as.data.frame(
    x = as.matrix(
      x = t(
        x = object@assays$RNA@data[training.genes, ]
      )
    )
  )
  training.data$class <- factor(x = training.classes)
  if (verbose) {
    message("Training Classifier ...")
  }
  classifier <- ranger::ranger(
    data = training.data,
    dependent.variable.name = "class",
    classification = TRUE,
    write.forest = TRUE,
    ...
  )
  return(classifier)
}

#' Assess Cluster Split
#'
#' Method for determining confidence in specific bifurcations in
#' the cluster tree. Use the Out of Bag (OOB) error of a random
#' forest classifier to judge confidence.
#'
#' @param object Seurat object
#' @param node Node in the cluster tree in question
#' @param cluster1 First cluster to compare
#' @param cluster2 Second cluster to compare
#' @param genes.training A vector of genes to use to train the classifier,
#' defaults to \code{rownames(x = object@data)}
#' @param print.output Print the OOB error for the classifier
#' @inheritDotParams BuildRFClassifier -object
#' @return Returns the Out of Bag error for a random forest classifier
#' trained on the split from the given node
#' @export
#'
#' @examples
#' pbmc_small
#' pbmc_small <- FindClusters(object = pbmc_small, reduction.type = "pca",
#'                            dims.use = 1:10, resolution = 1.1, save.SNN = TRUE)
#' pbmc_small <- BuildClusterTree(pbmc_small, reorder.numeric = TRUE, do.reorder = TRUE)
#' # Assess based on a given node
#' AssessSplit(pbmc_small, node = 11)
#' # Asses based on two given clusters (or vectors of clusters)
#' AssessSplit(pbmc_small, cluster1 = 5, cluster2 = 6)
#' AssessSplit(pbmc_small, cluster1 = 4, cluster2 = c(5, 6))
#'
AssessSplit <- function(
  object,
  node,
  cluster1,
  cluster2,
  genes.training = NULL,
  print.output = TRUE,
  ...
) {
  genes.training <- SetIfNull(x = genes.training, default = rownames(x = object@assays$RNA@data))
  genes.training <- intersect(x = genes.training, rownames(x = object@assays$RNA@data))
  if (!length(x = genes.training)) {
    stop("None of the genes provided are in the data")
  }
  tree <- object@tools$BuildClusterTree
  if (!missing(x = node)) {
    if (!missing(x = cluster1) || !missing(x = cluster2)) {
      warning("Both node and cluster IDs provided. Defaulting to using node ID")
    }
    possible.nodes <- c(
      DFT(tree = tree, node = tree$edge[,1][1]),
      tree$edge[,1][1]
    )
    if (!node %in% possible.nodes) {
      stop("Not a valid node")
    }
    split <- tree$edge[which(x = tree$edge[,1] == node), ][,2]
    group1 <- DFT(tree = tree, node = split[1], only.children = TRUE)
    group2 <- DFT(tree = tree, node = split[2], only.children = TRUE)
    if (any(is.na(x = group1))) {
      group1 <- split[1]
    }
    if (any(is.na(x = group2))) {
      group2 <- split[2]
    }
  } else {
    group1 <- cluster1
    group2 <- cluster2
  }
  group1.cells <- WhichCells(object = object, idents = group1)
  group2.cells <- WhichCells(object = object, idents = group2)
  assess.data <- SubsetData(
    object = object,
    cells = c(group1.cells, group2.cells)
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells = group1.cells,
    value = "g1"
  )
  assess.data <- SetIdent(
    object = assess.data,
    cells = group2.cells,
    value = "g2"
  )
  rfc <- BuildRFClassifier(
    object = assess.data,
    # training.genes = assess.data@var.genes,
    training.genes = genes.training,
    training.classes = Idents(object = assess.data),
    ...
  )
  oobe <- rfc$prediction.error
  if (print.output) {
    message(paste0("Out of Bag Error: ", round(x = oobe, digits = 4) * 100, "%"))
  }
  return(oobe)
}

#' Assess Internal Nodes
#'
#' Method for automating assessment of tree splits over all internal nodes,
#' or a provided list of internal nodes. Uses AssessSplit() for calculation
#' of Out of Bag error (proxy for confidence in split).
#'
#' @param object Seurat object
#' @param node.list List of internal nodes to assess and return
#' @param all.below If single node provided in node.list, assess all splits
#' below (and including) provided node
#' @param genes.training A vector of genes to use to train the classifier,
#' defaults to \code{rownames(x = object@data)}
#' .
#' @return Returns the Out of Bag error for a random forest classifiers trained on
#' each internal node split or each split provided in the node list.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small
#' pbmc_small <- FindClusters(object = pbmc_small, reduction.type = "pca",
#'                            dims.use = 1:10, resolution = 1.1, save.SNN = TRUE)
#' pbmc_small <- BuildClusterTree(pbmc_small, reorder.numeric = TRUE, do.reorder = TRUE)
#' AssessNodes(pbmc_small)
#' }
#'
AssessNodes <- function(
  object,
  node.list,
  all.below = FALSE,
  genes.training = NULL
) {
  genes.training <- SetIfNull(x = genes.training, default = rownames(x = object@assays$RNA@data))
  genes.training <- intersect(x = genes.training, rownames(x = object@assays$RNA@data))
  if (!length(x = genes.training)) {
    stop("None of the genes provided are in the data")
  }
  tree <- object@tools$BuildClusterTree
  if (missing(x = node.list)) {
    node.list <- GetAllInternalNodes(tree = tree)
  } else {
    possible.nodes <- GetAllInternalNodes(tree = tree)
    if (any(!node.list %in% possible.nodes)) {
      stop(paste(
        node.list[!(node.list %in% possible.nodes)],
        "not valid internal nodes"
      ))
    }
    if (length(x = node.list == 1) && all.below) {
      node.list <- c(node.list, DFT(tree = tree, node = node.list))
    }
  }
  oobe <- pbsapply(
    X = node.list,
    FUN = function(x) {
      return(AssessSplit(
        object = object,
        node = x,
        genes.training = genes.training,
        print.output = FALSE,
        verbose = FALSE
      ))
    }
  )
  return(data.frame(node = node.list, oobe))
}