TMS_remerge <- merge(Tongue, c(Heart_and_Aorta, Marrow, Mammary_Gland, Fat, Kidney, Liver, Lung, Limb_Muscle, 
                               Pancreas, Spleen, Thymus, Bladder, Skin, Large_Intestine, Trachea))
TMS_remerge <- RunTSNE()


node.scores <- AssessNodes(cd4t0n6)
node.scores[order(node.scores$oobe,decreasing = TRUE),] -> node.scores
nodes.merge <- node.scores[which(node.scores$oobe > 0.04),]
nodes.to.merge <-sort(nodes.merge$node)