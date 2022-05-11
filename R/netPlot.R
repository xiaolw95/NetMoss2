#' @title Network plot.
#' @description This is a function to plot combined networks.
#'
#' @param result the result list of NetMoss.
#' @param num.top numeric. The number of top taxon to be highlighted in the paired networks.
#' @param num.score numeric. The number of taxon to be plotted in the NetMoss score barplot and point plot.
#' @param e.th numeric. The threshold of microbial correlations in the networks. Edges greater than this threshold are shown in the networks.
#' @param nodeSize numeric. The size of nodes in the networks.
#' @param nodeDize numeric. The size of nodes in the networks. The larger the value, the greater the difference in the size of the nodes.
#' @param edgeSize numeric. The size of edges in the networks.
#' @param edgeDize numeric. The size of edges in the networks. The larger the value, the greater the difference in the size of the edges
#' @param arrowSize numeric. The size of arrow in the networks.
#' @param my.layout the layout of network.
#' @param my.label logic. If TRUE then the label of the nodes will be plotted.
#'
#' @return case and control networks
#' @export
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics legend par plot title
#'
#' @examples
#' data(testData)
#' nodes_result = NetMoss(case_dir = mydata[[1]],
#'      control_dir = mydata[[2]],
#'      net_case_dir = mydata[[3]],
#'      net_control_dir = mydata[[4]])
#' netPlot(nodes_result)
#'
#'
netPlot <- function(result,
                    num.top = 5,
                    num.score = 30,
                    e.th = 0.4,
                    nodeSize = 2,
                    nodeDize = 1.2,
                    edgeSize = 0.2,
                    edgeDize = 0.3,
                    arrowSize = 0,
                    my.layout = layout_as_star,
                    my.label = TRUE)
{
  ############################################
  my.wd = getwd()
  nodes_result = result[[1]]
  nodes_result = nodes_result[order(nodes_result$NetMoss_Score,decreasing = T),]
  top.tax = as.character(nodes_result[1:num.top,1])

  ##abundance fold change
  case.all.sample = result[[4]]
  control.all.sample = result[[5]]
  case.all.sample = case.all.sample[order(as.character(case.all.sample$X)),]
  control.all.sample = control.all.sample[order(as.character(control.all.sample$X)),]
  all.sample.data = data.frame()
  for (k in 1:nrow(case.all.sample))
  {
    all.sample.data[k,1] = case.all.sample[k,1]
    all.sample.data[k,2] = mean(as.numeric(as.character(case.all.sample[k,-1])))
    all.sample.data[k,3] = mean(as.numeric(as.character(control.all.sample[k,-1])))
    all.sample.data[k,4] = log2((all.sample.data[k,2]+1)/(all.sample.data[k,3]+1))
  }
  colnames(all.sample.data) = c("genus","case","control","FC")
  rownames(all.sample.data) = all.sample.data$genus

  ####plot NetMoss score
  nodes_plot = nodes_result[1:num.score,]
  nodes_plot$taxon_names = factor(nodes_plot$taxon_names,levels = rev(nodes_plot$taxon_names))
  nodes_plot$FC = all.sample.data[as.character(nodes_plot$taxon_names),"FC"]
  nodes_plot[which(nodes_plot$FC < 0),"type"] = "Enriched in control"
  nodes_plot[which(nodes_plot$FC >= 0),"type"] = "Enriched in case"

  ##point
  p1 = ggplot(nodes_plot,aes(taxon_names,NetMoss_Score, color = p.adj))+
    scale_colour_distiller(palette = "Spectral")+
    geom_bar(stat = 'identity',width = 0.03, fill = "black")+
    geom_point(size = 4)+
    coord_flip()+
    theme_bw()+
    theme(panel.grid = element_blank())+
    labs(y = "NetMoss score", x = paste0("Top ",num.score," bacteria"))+
    theme(legend.position = "top")

  #fold change
  p2 = ggplot(nodes_plot,aes(taxon_names,FC,color = type))+
    geom_point(size = 4)+coord_flip()+
    theme_bw()+theme(panel.grid = element_blank())+
    geom_hline(yintercept = 0)+
    scale_color_brewer(palette = "Pastel1", name = "")+
    labs(x = "", y = "Log2(FC)")+
    theme(legend.position = "top")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  pp = ggarrange(p1,p2,NULL,ncol = 3,nrow = 1,widths = c(1,0.3,0.1),align = "h")


    #####highlight top taxon
  ##case
  e1.case = data.frame(result[[2]])
  e1.case[upper.tri(e1.case)] <- 0
  e1.case$genus = rownames(e1.case)
  e2.case = melt(e1.case, id.vars = "genus")
  e2.case = e2.case[which(e2.case$value != 0), ]
  ##
  e2.case$genus = as.character(e2.case$genus)
  e2.case$genus = gsub("\\.", "-", e2.case$genus)
  e2.case$genus = gsub(" ", "-", e2.case$genus)
  e2.case$genus = gsub('\\[', '', e2.case$genus)
  e2.case$genus = gsub('\\]', '', e2.case$genus)
  e2.case$genus = gsub("\\.", "-", e2.case$genus)
  e2.case$variable = as.character(e2.case$variable)
  e2.case$variable = gsub("\\.", "-", e2.case$variable)
  e2.case$variable = gsub(" ", "-", e2.case$variable)
  e2.case$variable = gsub('\\[', '', e2.case$variable)
  e2.case$variable = gsub('\\]', '', e2.case$variable)
  e2.case$variable = gsub("\\.", "-", e2.case$variable)

  edge.case = e2.case[which(e2.case$genus != e2.case$variable), ]
  edge.case2 = edge.case[which(abs(edge.case$value) > e.th), ]

  node.case = data.frame(unique(c(
    as.character(edge.case2$genus),
    as.character(edge.case2$variable)
  )))
  colnames(node.case) = "node"
  rownames(node.case) = node.case$node
  rownames(nodes_result) = nodes_result$taxon_names

  top.tax2 = intersect(top.tax,node.case$node)
  if (length(top.tax2) == 0)
  {
    warning("the threshold is too high to find the target taxon in the case networks!")
  }

  inter.genus = intersect(as.character(node.case$node),as.character(nodes_result$taxon_names))
  node.case[inter.genus,"score"] = nodes_result[inter.genus,"NetMoss_Score"]
  node.case2 = node.case
  if (ncol(node.case2[which(is.na(node.case2$score)),]) != 0)
  {
    node.case2[which(is.na(node.case2$score)),"score"] = 0
  }else
  {
    node.case2 = node.case2
  }

  node.case2$weight = abs(node.case2$score)
  node.case2$class = 1   ### gray
  node.case2[top.tax2, "class"] = 2  ###red

  edge.case2$weight = 0.1
  for (kk in 1:length(top.tax2))
  {
    for (mm in 1:nrow(edge.case2))
    {
      if (as.character(edge.case2[mm,1]) == top.tax2[kk] || as.character(edge.case2[mm,2]) == top.tax2[kk])
      {
        edge.case2[mm,"weight"] = 1
      }else
      {
        next
      }
    }
  }

  edge.case2$class = 1  ##positive red
  edge.case2[which(edge.case2$value < 0), "class"] = 2   #negative blue


  ###plot
  g1 <- graph.empty()
  g1 <- graph_from_data_frame(edge.case2, vertices = node.case2)
  #
  nodeSize <- nodeSize
  nodeDize <- nodeDize
  edgeSize <- edgeSize
  edgeDize <- edgeDize
  VColor <- c("#8f8f8f", "#e57265")
  EColor <- c("#e57265", "#78a0c4")
  VText <- c(0.2,1)
  V(g1)$size <-
    nodeSize + nodeDize * 10 * as.numeric(as.vector(node.case2$weight))
  V(g1)$color <- VColor[node.case2$class]
  V(g1)$label.cex <- VText[node.case2$class]
  V(g1)$frame.color <- "black"
  E(g1)$width <-
    edgeSize + (edgeDize * abs(3 * as.numeric(as.vector(
      edge.case2$weight
    ))))
  E(g1)$color <- EColor[edge.case2$class]
  E(g1)$arrow.size <- arrowSize


  ##control
  e1.control = data.frame(result[[3]])
  e1.control[upper.tri(e1.control)] <- 0
  e1.control$genus = rownames(e1.control)
  e2.control = melt(e1.control, id.vars = "genus")
  e2.control = e2.control[which(e2.control$value != 0), ]
  ##
  e2.control$genus = as.character(e2.control$genus)
  e2.control$genus = gsub("\\.", "-", e2.control$genus)
  e2.control$genus = gsub(" ", "-", e2.control$genus)
  e2.control$genus = gsub('\\[', '', e2.control$genus)
  e2.control$genus = gsub('\\]', '', e2.control$genus)
  e2.control$genus = gsub("\\.", "-", e2.control$genus)
  e2.control$variable = as.character(e2.control$variable)
  e2.control$variable = gsub("\\.", "-", e2.control$variable)
  e2.control$variable = gsub(" ", "-", e2.control$variable)
  e2.control$variable = gsub('\\[', '', e2.control$variable)
  e2.control$variable = gsub('\\]', '', e2.control$variable)
  e2.control$variable = gsub("\\.", "-", e2.control$variable)

  edge.control = e2.control[which(e2.control$genus != e2.control$variable), ]
  edge.control2 = edge.control[which(abs(edge.control$value) > e.th), ]

  node.control = data.frame(unique(c(
    as.character(edge.control2$genus),
    as.character(edge.control2$variable)
  )))
  colnames(node.control) = "node"
  rownames(node.control) = node.control$node
  rownames(nodes_result) = nodes_result$taxon_names

  top.tax2 = intersect(top.tax,node.control$node)
  if (length(top.tax2) == 0)
  {
    warning("the threshold is too high to find the target taxon in the control networks!")
  }

  inter.genus = intersect(as.character(node.control$node),as.character(nodes_result$taxon_names))
  node.control[inter.genus,"score"] = nodes_result[inter.genus,"NetMoss_Score"]
  node.control2 = node.control
  if (ncol(node.control2[which(is.na(node.control2$score)),]) != 0)
  {
    node.control2[which(is.na(node.control2$score)),"score"] = 0
  }else
  {
    node.control2 = node.control2
  }

  node.control2$weight = abs(node.control2$score)
  node.control2$class = 1   ### gray
  node.control2[top.tax2, "class"] = 2  ###red

  edge.control2$weight = 0.1
  for (kk in 1:length(top.tax2))
  {
    for (mm in 1:nrow(edge.control2))
    {
      if (as.character(edge.control2[mm,1]) == top.tax2[kk] || as.character(edge.control2[mm,2]) == top.tax2[kk])
      {
        edge.control2[mm,"weight"] = 1
      }else
      {
        next
      }
    }
  }

  edge.control2$class = 1  ##positive red
  edge.control2[which(edge.control2$value < 0), "class"] = 2   #negative blue


  ###plot
  g2 <- graph.empty()
  g2 <- graph_from_data_frame(edge.control2, vertices = node.control2)
  #
  nodeSize <- nodeSize
  nodeDize <- nodeDize
  edgeSize <- edgeSize
  edgeDize <- edgeDize
  VColor <- c("#8f8f8f", "#e57265")
  EColor <- c("#e57265", "#78a0c4")
  VText <- c(0.2,1)
  V(g2)$size <-
    nodeSize + nodeDize * 10 * as.numeric(as.vector(node.control2$weight))
  V(g2)$color <- VColor[node.control2$class]
  V(g2)$label.cex <- VText[node.control2$class]
  V(g2)$frame.color <- "black"
  E(g2)$width <-
    edgeSize + (edgeDize * abs(3 * as.numeric(as.vector(
      edge.control2$weight
    ))))
  E(g2)$color <- EColor[edge.control2$class]
  E(g2)$arrow.size <- arrowSize

  ########
  if (!my.label)
  {
    par(mfrow = c(1, 1))
    ##half circle
    p11 = ggplot() +
      geom_ellipse(aes(x0 = 0, y0 = 0, a = 10, b = 6, angle = 0),
                   fill = "#f6a5c0", alpha = 0.1, color = NA) +
      scale_y_continuous(limits = c(-6,0),expand = c(0,0))+coord_flip()+
      theme_void()
    ##half circle
    p12 = ggplot() +
      geom_ellipse(aes(x0 = 0, y0 = 0, a = 10, b = 6, angle = 0),
                   fill = "#93cf96", alpha = 0.1, color = NA) +
      scale_y_continuous(limits = c(0,6),expand = c(0,0))+coord_flip()+
      theme_void()

    multiplot(p11,p12,cols = 2)

    ##case
    par(fig=c(0,0.5,0.1,0.9),new=TRUE)
    plot(
      g1,
      layout = my.layout,
      vertex.frame.color = "NA",
      vertex.shapes = "none",
      vertex.label = ""
    )

    ##title
    par(fig=c(0.3,0.4,0.8,0.95),new=TRUE)
    title("case network")

    ##legend
    par(fig=c(0.1,0.4,0.1,0.3),new=TRUE)
    legend(
      x = -1,
      y = -1.5,
      bty = "n",
      ####edges color
      c("positive correlation", "negative correlation"),
      lty = 1,
      lwd = 2,
      col = EColor
    )

    ##control
    par(fig=c(0.5,1,0.1,0.9),new=TRUE)
    plot(
      g2,
      layout = my.layout,
      vertex.frame.color = "NA",
      vertex.shapes = "none",
      vertex.label = ""
    )

    ##title
    par(fig=c(0.5,0.7,0.8,0.95),new=TRUE)
    title("control network")

    ##legend
    par(fig=c(0.5,0.8,0.1,0.3),new=TRUE)
    legend(
      x = -1.5,
      y = -1.5,
      bty = "n",
      ####nodes color
      c("others", "key taxon"),
      pch = 21,
      pt.bg = VColor
    )
  }else
  {
    par(mfrow = c(1, 1))
    ##half circle
    p11 = ggplot() +
      geom_ellipse(aes(x0 = 0, y0 = 0, a = 10, b = 6, angle = 0),
                   fill = "#f6a5c0", alpha = 0.1, color = NA) +
      scale_y_continuous(limits = c(-6,0),expand = c(0,0))+coord_flip()+
      theme_void()
    ##half circle
    p12 = ggplot() +
      geom_ellipse(aes(x0 = 0, y0 = 0, a = 10, b = 6, angle = 0),
                   fill = "#93cf96", alpha = 0.1, color = NA) +
      scale_y_continuous(limits = c(0,6),expand = c(0,0))+coord_flip()+
      theme_void()

    multiplot(p11,p12,cols = 2)

    ##case
    par(fig=c(0,0.5,0.1,0.9),new=TRUE)
    plot(
      g1,
      layout = my.layout,
      vertex.frame.color = "NA",
      vertex.shapes = "none"
    )

    ##title
    par(fig=c(0.3,0.4,0.8,0.95),new=TRUE)
    title("case network")

    ##legend
    par(fig=c(0.1,0.4,0.1,0.3),new=TRUE)
    legend(
      x = -1,
      y = -1.5,
      bty = "n",
      ####edges color
      c("positive correlation", "negative correlation"),
      lty = 1,
      lwd = 2,
      col = EColor
    )

    ##control
    par(fig=c(0.5,1,0.1,0.9),new=TRUE)
    plot(
      g2,
      layout = my.layout,
      vertex.frame.color = "NA",
      vertex.shapes = "none"
    )

    ##title
    par(fig=c(0.5,0.7,0.8,0.95),new=TRUE)
    title("control network")

    ##legend
    par(fig=c(0.5,0.8,0.1,0.3),new=TRUE)
    legend(
      x = -1.5,
      y = -1.5,
      bty = "n",
      ####nodes color
      c("others", "key taxon"),
      pch = 21,
      pt.bg = VColor
    )
  }

    setwd(my.wd)

    ###score
    ggsave("NetMoss_score.pdf", pp)

    print (paste0("the NetMoss score saved in ", my.wd))
  }
