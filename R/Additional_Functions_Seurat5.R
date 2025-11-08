
add_annotation <- function(seu, annotations, to_label, annotation_name = "CellType"){
  Idents(seu) <- to_label
  names(annotations) <- levels(seu)
  seu <- RenameIdents(seu, annotations)
  seu <- AddMetaData(object=seu, metadata=Idents(seu), col.name = annotation_name)
  
}


plot_Vln_by_cluster <- function(seu, output_pathway, 
                                input_assay = "ALRA", 
                                cluster_cols_vector,
                                feature_vector = c("TH","CD56"),
                                plot_name = "Title", 
                                obj_name = "sample_name"){
  DefaultAssay(seu) <- input_assay # this is where to get the expression values from 
  for (cluster_col in cluster_cols_vector){
    Idents(seu) <- cluster_col
    cat("plotting ",cluster_col, "\n")
    p <- VlnPlot(
      seu,
      features = feature_vector,
      combine = FALSE,
      stack = TRUE,
      flip = TRUE,
      fill.by = "ident"
    ) + 
      RotatedAxis() + 
      ggtitle(plot_name)
    # Print or save
    n_rows = length(feature_vector)
    h = 120 * n_rows + 100
    n_col = length(unique(seu@meta.data[[cluster_col]]))
    w = n_col * 80 + 200
    png(paste0(output_pathway, obj_name,plot_name,"VlnPlot.png"), width = w, height = h)
    print(p)
    dev.off()
  }
}




run_umaps_clusters <- function(seu, kn_values = c(10, 20, 30), res_values = c(0.2, 0.5, 1.0), input_assay = "RNA", features_vector = c("TH","CD56"), output_pathway) {
  # Loop through kN values
  for (kN in kn_values) {
    cat("running using pca inputs \n")
    input_name = paste0("pca_", input_assay)
    umap_name <- paste0("kn", kN,"_",input_name, "_umap")
    graph_name <- paste0("kn", kN,"_",input_name,"_snn")  # Ensuring distinct graph names
    cat("Running umap: ",umap_name,"\n")
    # Run UMAP and FindNeighbors with UMAP input
    seu <- RunUMAP(seu, reduction = input_name, dims = 1:length(features_vector) -1, n.neighbors = kN, reduction.name = umap_name, spread = 0.6, min.dist = 0.2)
    seu <- FindNeighbors(seu, dims = 1:2, k.param = kN, reduction = umap_name, graph.name = graph_name)
    # Run FindClusters for each resolution
    for (res in res_values) {
      res_col_name <- paste0(graph_name,"_res.", res)
      cat("Finding clusters for kn",res_col_name,"\n")
      seu <- FindClusters(seu, resolution = res, graph.name = graph_name)
      cat("Plotting: ",res_col_name,"\n")
      filename <- paste0(output_pathway,umap_name,"_",res_col_name,".png")
      png(filename, width = 700, height = 400)
      print(DimPlot(seu, reduction = umap_name, label = TRUE, repel = TRUE, group.by = res_col_name, raster = FALSE))
      dev.off()
      plot_Vln_by_cluster(seu, output_pathway, 
                          input_assay = input_assay, 
                          cluster_cols_vector = res_col_name,
                          feature_vector = features_vector,
                          plot_name = umap_name, 
                          obj_name = input_assay)
    }
    cat("running using using direct feature inputs \n")    
    input_name = paste0("ft_", input_assay)
    umap_name <- paste0("kn", kN,"_",input_name, "_umap")
    graph_name <- paste0("kn", kN,"_",input_name,"_snn")  # Ensuring distinct graph names
    cat("Running umap: ",umap_name,"\n")
    # Run UMAP and FindNeighbors with UMAP input
    seu <- RunUMAP(seu, features = features_vector, n.neighbors = kN, reduction.name = umap_name, spread = 0.6, min.dist = 0.2)
    seu <- FindNeighbors(seu, dims = 1:2, k.param = kN, reduction = umap_name, graph.name = graph_name)
    # Run FindClusters for each resolution
    for (res in res_values) {
      res_col_name <- paste0(graph_name,"_res.", res)
      cat("Finding clusters for kn",res_col_name,"\n")
      seu <- FindClusters(seu, resolution = res, graph.name = graph_name)
      cat("Plotting: ",res_col_name,"\n")
      filename <- paste0(output_pathway,umap_name,"_",res_col_name,".png")
      png(filename, width = 700, height = 400)
      print(DimPlot(seu, reduction = umap_name, label = TRUE, repel = TRUE, group.by = res_col_name, raster = FALSE))
      dev.off()
      
      plot_Vln_by_cluster(seu, output_pathway, 
                          input_assay = input_assay, 
                          cluster_cols_vector = res_col_name,
                          feature_vector = features_vector,
                          plot_name = umap_name, 
                          obj_name = input_assay)
    }
  }
  return(seu)
}
