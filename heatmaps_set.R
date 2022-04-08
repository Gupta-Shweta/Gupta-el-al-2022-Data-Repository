# plot heatmap
library(ComplexHeatmap)

heatmaps_set <- function(input,selected_metabolite,selected_sample, color, title){
  selected_sample = selected_sample[order(selected_sample$comb),]
  heatmap.exp = input[row.names(input) %in% selected_metabolite, colnames(input) %in% selected_sample$Sample.ID]
  heatmap.exp = heatmap.exp[, match(selected_sample$Sample.ID,colnames(heatmap.exp))]
  heatmap.exp = t(scale(t(heatmap.exp),center = T, scale = T))
  stopifnot(all.equal(colnames(heatmap.exp), selected_sample$Sample.ID))
  colnames(heatmap.exp) = selected_sample$comb
  stopifnot(dim(heatmap.exp)[1]==length(selected_metabolite))
  #print(floor(min(heatmap.exp)))
  #print(ceiling(max(heatmap.exp)))
  
  heatmap = Heatmap(heatmap.exp,cluster_columns = TRUE,
                    column_title = paste(dim(heatmap.exp)[1]," significant metabolites have trait x treatment interactions"),
                    column_title_side = "top", 
                    column_dend_side = "top",
                    show_row_names = TRUE, 
                    name = "z-score")
  
  pdf(paste0("plots/Heatmap_",title,".pdf"),15,10)
  print(heatmap)
  dev.off()
}
