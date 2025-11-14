find.associations.offline <- function(graph,data)
{
  out <- NULL

  out <- tryCatch(
    {
      vertices <- getAllVertices(graph,data$Linked_SNP,1)

      if(length(vertices) == 0)
        return(NULL)

      net=induced.subgraph(graph,c(vertices))
      out=as.data.table(get.data.frame(net))
      out <- out[startsWith(from,'rs'),]
      names(out) <- c('SNP','Phenotype')
      return(out)
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning')
      return(NULL)
    },
    warning = function(cond) {
      print_and_log(paste('Warning occured in plotting:',cond$message),level='warning')
      return(NULL)
    })

  return(out)
}

#deprecated
make.graphs.offline <- function(graph,outputFolder,projectName,i, data, graph.layout,level)
{
  tryCatch(
    {
      vertices <- getAllVertices(graph,data$Linked_SNP,level)

      if(length(vertices) == 0)
        return(NULL)

      net=induced.subgraph(graph,c(vertices))

      fileName = sprintf('%s/%s_gSNP_%s.png',outputFolder,projectName,i)
      png(filename = fileName,width = 10,height = 8,units = 'in',res = 300)

      plot(net,
           edge.arrow.size=.4,
           vertex.label.color="black",
           vertex.label.dist=1.5,
           vertex.size=7,
           layout=graph.layout,
           vertex.color=sapply(V(net)$type, function(x) ifelse(x=='Trait','red',ifelse(x=='SNP','blue','green'))))

      dev.off()

      print_and_log(sprintf('Plot saved: %s',fileName),display=FALSE)

      data=as.data.table(get.data.frame(net))
      data <- data[startsWith(from,'rs'),]

      return(invisible(data))
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning',display=FALSE)
    },
    warning = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message) ,level='warning',display=FALSE)
    })
}

#deprecated
make.graph.for.alldataset_noClustering <- function(graph,outputFolder,projectName, data, graph.layout,level)
{
  tryCatch(
    {
      vertices <- getAllVertices(graph,data$Linked_SNP,level)

      if(length(vertices) == 0)
        return(NULL)

      net=induced.subgraph(graph,c(vertices))

      fileName = sprintf('%s/%s_trait_graph.png',outputFolder,projectName)
      png(filename = fileName,width = 20,height = 16,units = 'in',res = 300)

      # add layout to decrease congestion
      layout <- layout_with_fr(net)

      plot(net,
           edge.arrow.size=.4,
           vertex.label.color="black",
           vertex.label.dist=1.5,
           vertex.size=7,
           layout=layout,
           #layout=graph.layout,
           vertex.color=sapply(V(net)$type, function(x) ifelse(x=='Trait','red',ifelse(x=='SNP','blue','green'))))

      mtext('Phenotype association plot',col = 'red',side = 3, adj = 0, line = 1)

      dev.off()

      print_and_log(sprintf('Plot saved: %s',fileName),display=FALSE)

      data=as.data.table(get.data.frame(net))
      data <- data[startsWith(from,'rs'),]

      return(invisible(data))
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning',display=FALSE)
    },
    warning = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message) ,level='warning',display=FALSE)
    })
}

#deprecated
make.graph.for.alldataset <- function(graph,outputFolder,projectName, data, graph.layout,level)
{
  tryCatch(
    {
      vertices <- getAllVertices(graph,data$Linked_SNP,level)

      if(length(vertices) == 0)
        return(NULL)

      net=induced.subgraph(graph,c(vertices))

      fileName = sprintf('%s/%s_trait_graph.png',outputFolder,projectName)


      # convert to undirected for clustering
      undirected_net <- as.undirected(net, mode = "collapse")

      # Apply the Louvain method
      clusters <- cluster_louvain(undirected_net)
      cluster_tbl <- data.table('item' = V(net)$name,
                                'cluster' = clusters$membership)

      cluster_tbl <- cluster_tbl[order(cluster),]
      cluster_tbl <- cluster_tbl[!(startsWith(item,prefix = 'rs'))]

      # plot
      trait_assoc_plot <- ggraph(undirected_net, layout = "fr") +
        geom_edge_link() +
        geom_node_point(aes(color = factor(membership(clusters))), size = 5) +
        geom_node_text(aes(label = name), repel = TRUE, max.overlaps = Inf) +
        theme_void()+
        theme(legend.position = "none")

      ggsave(trait_assoc_plot,filename=fileName,
             width = 20,height = 15,units = 'in',limitsize = FALSE)


      #mtext('Phenotype association plot',col = 'red',side = 3, adj = 0, line = 1)

      print_and_log(sprintf('Plot saved: %s',fileName),display=FALSE)

      data=as.data.table(get.data.frame(net))
      data <- data[startsWith(from,'rs'),]


      return(list("graph.data"=data,
                  "graph.analysis.output"=cluster_tbl))
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message),level='warning',display=FALSE)
    },
    warning = function(cond) {
      print_and_log(paste('Error occured in plotting:', cond$message) ,level='warning',display=FALSE)
    })
}

make.graph.for.alldataset.clumped <- function(data, graph.layout,fileName)
{
  tryCatch(
    {
      data <- data[,c('gSNP','Phenotype')]
      data <- data[!is.na(Phenotype) & Phenotype !='',]
      data <- data[, .(Phenotype = unlist(tstrsplit(Phenotype, "; ", type.convert = TRUE))), by = "gSNP"]
      data <- data[!is.na(Phenotype),]
      data <- data[!duplicated(data),]

      catalog_nodes=unique(c(data$gSNP,data$Phenotype))
      catalog_nodes=data.table(catalog_nodes)

      catalog_nodes[,type := ifelse(grepl(x = catalog_nodes,pattern = '^rs\\d'),'SNP','Trait')]

      catalog_edges = data
      names(catalog_edges)=c('from','to')


      net <- graph_from_data_frame(d=catalog_edges, vertices=catalog_nodes, directed=T)


      # convert to undirected for clustering
      undirected_net <- as.undirected(net, mode = "collapse")

      # Apply the Louvain method
      clusters <- cluster_louvain(undirected_net)
      cluster_tbl <- data.table('item' = V(net)$name,
                                'cluster' = clusters$membership)

      cluster_tbl <- cluster_tbl[order(cluster),]
      cluster_tbl <- cluster_tbl[!(startsWith(item,prefix = 'rs'))]

      # plot
      trait_assoc_plot <- ggraph(undirected_net, layout = "fr") +
        geom_edge_link() +
        geom_node_point(aes(color = factor(membership(clusters))), size = 5) +
        geom_node_text(aes(label = name), repel = TRUE, max.overlaps = Inf) +
        theme_void()+
        theme(legend.position = "none")

      ggsave(trait_assoc_plot,filename=fileName,
             width = 20,height = 15,units = 'in',limitsize = FALSE)


      #mtext('Phenotype association plot',col = 'red',side = 3, adj = 0, line = 1)

      print_and_log(sprintf('Plot saved: %s',fileName),display=FALSE)

      data=as.data.table(get.data.frame(net))
      data <- data[startsWith(from,'rs'),]


      return(list("graph.data"=data,
                  "graph.analysis.output"=cluster_tbl))
    }
    , error = function(cond) {
      print_and_log(paste('Error occured in GWAS catalog plotting:', cond$message),level='warning',display=FALSE)
    },
    warning = function(cond) {
      print_and_log(paste('Error occured in GWAS catalog plotting:', cond$message) ,level='warning',display=FALSE)
    })
}

getVertices.byName <- function(net,ID)
{
  if(!is.element(ID,attr(V(net),'name')))
    return(NULL)

  v1=V(net)[name==ID]
  v1s=neighbors(net,v1,mode = 'all')
  return(c(v1,v1s))
}

getVertices.byNumber <- function(net,number)
{
  v1=V(net)[number]
  v1s=neighbors(net,v1,mode = 'out')
  return(c(v1,v1s))
}


getAllVertices <- function(net,IDs,level=1)
{
  out.level1 <- c()
  out.level2 <- c()
  out.level3 <- c()

  for(ID in IDs)
  {
    ver <- getVertices.byName(net,ID)
    if(!is.null(ver))
      out.level1 <- union(out.level1,ver)
  }

  if(level==3)
  {

    for(out in out.level1)
    {
      ver <- getVertices.byNumber(net,out)
      if(!is.null(ver))
        out.level2 <- union(out.level2,ver)
    }

    for(out in out.level2)
    {
      ver <- getVertices.byNumber(net,out)
      if(!is.null(ver))
        out.level3 <- union(out.level3,ver)
    }

    return(out.level3)

  }else if(level==2){
    for(out in out.level1)
    {
      ver <- getVertices.byNumber(net,out)
      if(!is.null(ver))
        out.level2 <- union(out.level2,ver)
    }
    return(out.level2)

  }else{
    return(out.level1)
  }

}


