make.graphs <- function(data,i,outputFolder,projectName,graph_layout)
{
  tryCatch(
    {
      data <- data[Associations != '',list(Linked_SNP, Associations)]

      if(nrow(data) == 0)
        return(NULL)

      data = data[, .(Associations = unlist(strsplit(Associations, ","))), by = "Linked_SNP"]

      # remove unspecified associations
      data=data[!grepl(x=Associations, pattern = 'phenotype not specified',ignore.case = TRUE),]

      # trim the associations
      data[,Associations := trimws(Associations)]

      # make the graph nodes and edges
      nodes <- data.table(name = c(data$Linked_SNP,data$Associations), type= c(rep(c('SNP','Trait'),each=nrow(data))))
      nodes <- nodes[!duplicated(nodes)]

      edges <- data[,list(Linked_SNP,Associations)]

      net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T)


      fileName = sprintf('%s/%s_gSNP_%s.png',outputFolder,projectName,i)
      png(filename = fileName,width = 10,height = 8,units = 'in',res = 300)

      plot(net,
           edge.arrow.size=.05,
           vertex.label.color="black",
           vertex.label.dist=1.5,
           vertex.size=7,
           layout = graph_layout,
           vertex.color=c( "pink", "skyblue")[1+(V(net)$type=="SNP")] )

      dev.off()

      print_and_log(sprintf('Plot saved: %s',fileName))
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning')

    },
    warning = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning')

    })
}

make.eQTL.graphs.gtex <- function(data,i,outputFolder,projectName,graph_layout)
{
  tryCatch(
    {
      nodes <- data.table(name = c(data$Linked_SNP,data$geneSymbol), type= c(rep(c('SNP','Gene'),each=nrow(data))))
      nodes <- nodes[!duplicated(nodes)]

      edges <- data[,list(Linked_SNP,geneSymbol)]

      net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T)


      fileName = sprintf('%s/%s_eQTL_GTEx_gSNP_%s.png',outputFolder,projectName,i)
      png(filename = fileName,width = 10,height = 8,units = 'in',res = 300)

      plot(net,
           edge.arrow.size=.05,
           vertex.label.color="black",
           vertex.label.dist=1.5,
           vertex.size=7,
           layout=graph_layout,
           vertex.color=c( "pink", "skyblue")[1+(V(net)$type=="SNP")] )

      mtext('eQTL association plot [GTEx]',col = 'red',side = 3, adj = 0, line = 1)
      mtext(gsub(x = data[1,]$tissueSiteDetailId,pattern = '_',replacement = ' '),col = 'red',side = 3, adj = 0, line = 0)
      dev.off()

      print_and_log(sprintf('Plot saved: %s',fileName))
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning')

    },
    warning = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning')

    })
}

#deprecated
make.eQTL.graphs.ebi <- function(data,i,outputFolder,projectName,graph_layout)
{
  tryCatch(
    {
      nodes <- data.table(name = c(data$rsId,data$Gene), type= c(rep(c('SNP','Gene'),each=nrow(data))))
      nodes <- nodes[!duplicated(nodes)]

      edges <- data[,list(rsId,Gene)]

      net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T)


      fileName = sprintf('%s/%s_eQTL_gSNP_%s.png',outputFolder,projectName,i)
      png(filename = fileName,width = 10,height = 8,units = 'in',res = 300)

      plot(net,
           edge.arrow.size=.05,
           vertex.label.color="black",
           vertex.label.dist=1.5,
           vertex.size=7,
           layout=graph_layout,
           vertex.color=c( "pink", "skyblue")[1+(V(net)$type=="SNP")] )

      mtext('eQTL association plot',col = 'red',side = 3, adj = 0, line = 1)
      mtext(gsub(x = data[1,]$eQTL_group,pattern = '_',replacement = ' '),
            col = 'red',side = 3, adj = 0, line = 0)

      dev.off()

      print_and_log(sprintf('Plot saved: %s',fileName),display=FALSE)
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning',display=FALSE)

    },
    warning = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning',display=FALSE)

    })
}

#deprecated
make.eQTL.graphs.ebi.for.alldataset_noClustering <- function(data,outputFolder,projectName,graph_layout)
{
  tryCatch(
    {
      nodes <- data.table(name = c(data$rsId,data$Gene), type= c(rep(c('SNP','Gene'),each=nrow(data))))
      nodes <- nodes[!duplicated(nodes)]

      edges <- data[,list(rsId,Gene)]
      edges <- edges[!duplicated(edges)]

      net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T)


      fileName = sprintf('%s/%s_eQTL.png',outputFolder,projectName)
      png(filename = fileName,width = 10,height = 8,units = 'in',res = 300)

      plot(net,
           edge.arrow.size=.05,
           vertex.label.color="black",
           vertex.label.dist=1.5,
           vertex.size=7,
           layout=graph_layout,
           vertex.color=c( "pink", "skyblue")[1+(V(net)$type=="SNP")] )

      mtext('eQTL association plot',col = 'red',side = 3, adj = 0, line = 1)
      # mtext(gsub(x = data[1,]$eQTL_group,pattern = '_',replacement = ' '),
      #       col = 'red',side = 3, adj = 0, line = 0)

      dev.off()

      print_and_log(sprintf('Plot saved: %s',fileName),display=FALSE)
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning',display=FALSE)

    },
    warning = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning',display=FALSE)

    })
}

#deprecated
make.eQTL.graphs.ebi.for.alldataset <- function(data,outputFolder,projectName,graph_layout)
{
  tryCatch(
    {
      nodes <- data.table(name = c(data$rsId,data$Gene), type= c(rep(c('SNP','Gene'),each=nrow(data))))
      nodes <- nodes[!duplicated(nodes)]

      edges <- data[,list(rsId,Gene)]
      edges <- edges[!duplicated(edges)]

      net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T)
      fileName = sprintf('%s/%s_eQTL.png',outputFolder,projectName)

      # convert to undirected for clustering
      undirected_net <- as.undirected(net, mode = "collapse")

      # Apply the Louvain method
      clusters <- cluster_louvain(undirected_net)
      cluster_tbl <- data.table('item' = V(net)$name,
                                'cluster' = clusters$membership)

      cluster_tbl <- cluster_tbl[order(cluster),]
      cluster_tbl <- cluster_tbl[!(startsWith(item,prefix = 'rs'))]

      # plot
      eqtl_plot <- ggraph(undirected_net, layout = "fr") +
        geom_edge_link() +
        geom_node_point(aes(color = factor(membership(clusters))), size = 5) +
        geom_node_text(aes(label = name), repel = TRUE, max.overlaps = Inf) +
        theme_void()+
        theme(legend.position = "none")

      ggsave(eqtl_plot,filename =fileName,width = 20,height = 15,units = 'in',limitsize = FALSE)


      print_and_log(sprintf('Plot saved: %s',fileName),display=FALSE)
      return(cluster_tbl)
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning',display=FALSE)
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning',display=FALSE)
      return(NULL)

    })
}

#deprecated
make.eQTL.graphs.ebi.for.alldataset.clumped.with.categories <- function(data,output.index.table,outputFolder,projectName,graph_layout)
{
  tryCatch(
    {

      #data <- data[pValue<5e-8,c('#gSNP','Gene','eQTL_group')]
      data <- data[,c('#gSNP','Gene','eQTL_group')]

      data <- data[!is.na(Gene) & Gene !='',]
      data <- data[!duplicated(data)]

      data <- merge(x=data,y=output.index.table,by='#gSNP',all.x = TRUE,sort = FALSE)
      setDT(data)
      data$`#gSNP` <- NULL
      setcolorder(data,'gSNP')

      nodes <- data.table(name = c(data$gSNP,data$Gene,data$eQTL_group),
                          type= c(rep(c('SNP','Gene','Group'),each=nrow(data))))
      nodes <- nodes[!duplicated(nodes)]

      edges_var_gene <- data[, .(from = gSNP, to = Gene)]
      edges_gene_group <- data[, .(from = Gene, to = eQTL_group)]
      edges <- rbindlist(list(edges_var_gene, edges_gene_group))
      edges <- edges[!duplicated(edges),]

      net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T)
      fileName = sprintf('%s/%s_eQTL.png',outputFolder,projectName)

      # convert to undirected for clustering
      undirected_net <- as.undirected(net, mode = "collapse")

      # Apply the Louvain method
      clusters <- cluster_louvain(undirected_net)
      cluster_tbl <- data.table('item' = V(net)$name,
                                'type' = V(net)$type,
                                'cluster' = clusters$membership)

      cluster_tbl <- cluster_tbl[order(cluster),]
      cluster_tbl <- cluster_tbl[!(startsWith(item,prefix = 'rs'))]

      # plot
      eqtl_plot <- ggraph(undirected_net, layout = "fr") +
        geom_edge_link() +
        geom_node_point(aes(color = factor(membership(clusters))), size = 5) +
        geom_node_text(aes(label = name), repel = TRUE, max.overlaps = Inf) +
        theme_void()+
        theme(legend.position = "none")

      ggsave(eqtl_plot,filename =fileName,width = 20,height = 15,units = 'in',limitsize = FALSE)


      print_and_log(sprintf('Plot saved: %s',fileName),display=FALSE)
      return(cluster_tbl)
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in eqtl catalog plotting:', cond$message),level='warning',display=FALSE)
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Error occured in eqtl catalog plotting:', cond$message),level='warning',display=FALSE)
      return(NULL)

    })
}

make.eQTL.graphs.ebi.for.alldataset.clumped <- function(data,output.index.table,graph_layout,fileName)
{
  tryCatch(
    {

      #data <- data[pvalue<5e-8,c('#gSNP','Gene')]
      data <- data[,c('#gSNP','Gene')]
      data <- data[!is.na(Gene) & Gene !='',]
      data <- data[!duplicated(data)]

      data <- merge(x=data,y=output.index.table,by='#gSNP',all.x = TRUE,sort = FALSE)
      setDT(data)
      data$`#gSNP` <- NULL
      setcolorder(data,'gSNP')

      nodes <- data.table(name = c(data$gSNP,data$Gene),
                          type= c(rep(c('SNP','Gene'),each=nrow(data))))
      nodes <- nodes[!duplicated(nodes)]

      edges <- data[, .(from = gSNP, to = Gene)]

      net <- graph_from_data_frame(d=edges, vertices=nodes, directed=T)

      # convert to undirected for clustering
      undirected_net <- as.undirected(net, mode = "collapse")

      # Apply the Louvain method
      clusters <- cluster_louvain(undirected_net)
      cluster_tbl <- data.table('item' = V(net)$name,
                                'type' = V(net)$type,
                                'cluster' = clusters$membership)

      cluster_tbl <- cluster_tbl[order(cluster),]
      cluster_tbl <- cluster_tbl[!(startsWith(item,prefix = 'rs'))]

      # plot
      eqtl_plot <- ggraph(undirected_net, layout = "fr") +
        geom_edge_link() +
        geom_node_point(aes(color = factor(membership(clusters))), size = 5) +
        geom_node_text(aes(label = name), repel = TRUE, max.overlaps = Inf) +
        theme_void()+
        theme(legend.position = "none")

      ggsave(eqtl_plot,filename =fileName,width = 20,height = 15,units = 'in',limitsize = FALSE)


      print_and_log(sprintf('Plot saved: %s',fileName),display=FALSE)
      return(cluster_tbl)
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in eqtl catalog plotting:', cond$message),level='warning',display=FALSE)
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Error occured in eqtl catalog plotting:', cond$message),level='warning',display=FALSE)
      return(NULL)

    })
}
