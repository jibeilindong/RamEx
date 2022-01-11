#################################################################
# Function:  Ramanome IRCA util  
# Last update: 2021-09-22, Jing Gongchao; He Yuehui; Shi Huang
#################################################################
# install necessary libraries
p <- c("optparse","RColorBrewer","igraph", "circlize", "permute", "ggplot2", "reshape2", "proxy", "pheatmap", "Cairo", "grid", "plyr")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
    suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
#--------------------------------------------------
CleanData <-function(bdata, removeNA=T, removeNeg=T){
    if(sum(bdata==Inf)>0){
        inx <- bdata == Inf;
        bdata[inx] <- NA;
        bdata[inx] <- max(bdata, na.rm=T)*2
    }
    if(sum(bdata==-Inf)>0){
        inx <- bdata == -Inf;
        bdata[inx] <- NA;
        bdata[inx] <- min(bdata, na.rm=T)/2
    }
    if(removeNA){
        if(sum(is.na(bdata))>0){
            bdata[is.na(bdata)] <- min(bdata, na.rm=T)/2
        }
    }
    if(removeNeg){
        if(sum(bdata<=0) > 0){
            inx <- bdata <= 0;
            bdata[inx] <- NA; 
            bdata[inx] <- min(bdata, na.rm=T)/2
        }
    }
    bdata;
}

JSD<-function(mat, eps=10^-4, overlap=TRUE,...)
{
    
    if(!is.numeric(mat))
        stop("mat must be a numeric matrix\n")
    mat<-t(mat)
    z <- matrix(NA, nrow=ncol(mat), ncol=ncol(mat))
    colnames(z) <- rownames(z) <- colnames(mat)
    
    w <- mat < eps
    if (any(w)) mat[w] <- eps
## If you takes as input a matrix of density values 
## with one row per observation and one column per 
## distribution, add following statement below.
   mat <- sweep(mat, 2, colSums(mat) , "/")
    
    for(k in seq_len(ncol(mat)-1)){
      for(l in 2:ncol(mat)){
        ok <- (mat[, k] > eps) & (mat[, l] > eps)
        if (!overlap | any(ok)) {
          m=0.5*(mat[,k]+mat[,l])
          z[k,l] <- sqrt(0.5*sum(mat[,k] *(log(mat[,k]) - log(m)))+0.5*sum(mat[,l] *(log(mat[,l]) - log(m))))
          z[l,k] <- sqrt(0.5*sum(mat[,l] *(log(mat[,l]) - log(m)))+0.5*sum(mat[,k] *(log(mat[,k]) - log(m))))
        }
      }
    }
    diag(z)<-0
    z
}


vectorize_dm<-function(dm, group=NULL, duplicate=TRUE)
{
        if(ncol(dm)!=nrow(dm) & any(is.na(dm))==TRUE)
        stop('The distance matrix is not squared')
        dm<-data.matrix(dm)
        require("reshape2")
    if(!is.null(group)){
        if( length(unique(group))==1)
        stop('At least two levels for a given sample category in your metadata file are required.')
        if( length(group)!=nrow(dm))
        stop('The number of rows in metadata and distance matrix are not equal')
        if(is.factor(group) & nlevels(group)>length(group)*0.9)
        stop('The number of levels in a certain category can not exceed 90% of total number of samples')
        
        colnames(dm)<-rownames(dm)<-paste(rownames(dm), group, sep="____")
        if(duplicate){ 
        melt_dm<-subset(melt(dm), value!=0 & value!=1) 
        }else{ 
        dm[lower.tri(dm)]<-NA
        melt_dm<-subset(melt(dm), ! is.na(value) & value!=0 & value!=1) 
        }
        
        Row_Info<-data.frame(do.call(rbind,strsplit(as.character(melt_dm[,1]),"____",fixed=TRUE)))
        Col_Info<-data.frame(do.call(rbind,strsplit(as.character(melt_dm[,2]),"____",fixed=TRUE)))
        VS<-paste(Row_Info[,2],"_VS._",Col_Info[,2],sep="")
        dm_value<-data.frame(VS,Row_Info,Col_Info,d=melt_dm$value)
        
        colnames(dm_value)<-c("GroupPair","Sample_1","Group_1","Sample_2","Group_2","value")
        if(is.factor(group)){
        DistType<-as.factor(dm_value$Group_1==dm_value$Group_2)
        DistType<-factor(DistType,levels=levels(DistType),labels=c("AllBetween","AllWithin"))
        dm_value<-data.frame(DistType,dm_value)}
    }else{
        if (duplicate) { 
        dm_value<-subset(melt(dm), value!=0 & value!=1)  
        }else{ 
        dm[lower.tri(dm)] <- NA
        dm_value <- subset(melt(dm), ! is.na(value) & value!=0 & value!=1)  }
        colnames(dm_value)<-c("Sample_1","Sample_2","value")
        }
   
dm_value

}

getPermuteMatrix <- function(perm, N,  strata = NULL) {
    ## 'perm' is either a single number, a how() structure or a
    ## permutation matrix
    if (length(perm) == 1) {
        perm <- how(nperm = perm) 
    }
    ## apply 'strata', but only if possible: ignore silently other cases
    if (!missing(strata) && !is.null(strata)) {
        if (inherits(perm, "how") && is.null(getBlocks(perm)))
            setBlocks(perm) <- strata
    }
    ## now 'perm' is either a how() or a matrix
    if (inherits(perm, "how"))
        perm <- shuffleSet(N, control = perm)
    ## now 'perm' is a matrix (or always was). If it is a plain
    ## matrix, set minimal attributes for printing. This is a dirty
    ## kluge: should be handled more cleanly.
    if (is.null(attr(perm, "control")))
        attr(perm, "control") <-
            structure(list(within=list(type="supplied matrix"),
                           nperm = nrow(perm)), class = "how")
    perm
}
#-------------------------------
log_mat<-function(mat, base=2, peudozero=1e-6){
         mat[mat==0] <- peudozero
         log_mat <- log(mat, base)
         return(log_mat)
         }
#-------------------------------
PAM.best<-function(matrix,dm){
          if(!is.numeric(matrix))
                 stop("matrix must be a numeric matrix\n")
          if(!is.numeric(dist.mat) && class(dm)=="dist")
                 stop("dist.mat must be numeric distance matrix\n")
          require("cluster")
          require("fpc")
          require("clusterSim")
          
          min_nc=2
          if(nrow(matrix)>20){
          max_nc=20} else {
          max_nc=nrow(matrix)-1}
          res <- array(0,c(max_nc-min_nc+1, 2))
          res[,1] <- min_nc:max_nc
          siavgs <- array(0,c(max_nc-min_nc+1, 2))
          siavgs[,1] <- min_nc:max_nc
          clusters <- NULL
          for (nc in min_nc:max_nc)
          {
          cl <- pam(dm, nc, diss=TRUE)
          res[nc-min_nc+1,2] <- CH <- index.G1(matrix,cl$cluster,d=dm,centrotypes="medoids")
          siavgs[nc-1,2]<-cl$silinfo$avg.width
          clusters <- rbind(clusters, cl$cluster)
          }
          CH_nc<-(min_nc:max_nc)[which.max(res[,2])]
          Si_nc<-(min_nc:max_nc)[which.max(siavgs[,2])]
          print(paste("max CH for",CH_nc,"clusters=",max(res[,2])))
          print(paste("max Si for",Si_nc,"clusters=",max(siavgs[,2])))
          
          CH_cluster<-clusters[which.max(res[,2]),]
          Si_cluster<-clusters[which.max(siavgs[,2]),]
          objectList      <- list()
             objectList$min_nc <- min_nc
             objectList$max_nc <- max_nc
             objectList$CH     <- CH_cluster
             objectList$Si     <- Si_cluster
             objectList$CH_nc  <- CH_nc
             objectList$Si_nc  <- Si_nc
             objectList$res    <- res
             objectList$siavgs <- siavgs
             return(objectList)
}


Plot_network_graph<-function(Corr_mat, Pos_Edge=TRUE, Neg_Edge=TRUE, Threshold=0.6, node_size="degree", node_color="degree", layout="layout.fruchterman.reingold", outdir="./out.pdf" ){
                    #init graph
                    if (Pos_Edge & Neg_Edge){
                       g <- graph.adjacency(abs(Corr_mat) > Threshold & abs(Corr_mat)<1,mode="undirected")
                    }else if (Pos_Edge){
                         g <- graph.adjacency(Corr_mat > Threshold & Corr_mat < 1,mode="undirected")
                    }else if (Neg_Edge){
                         g <- graph.adjacency(Corr_mat < -1 * Threshold & Corr_mat > -1 ,mode="undirected")
                    }else stop('Please check the pos/neg edges parameters')
                    
                    E(g)$color <- 'grey60'
                    V(g)$size <- 5
                    V(g)$color <- 'cornflowerblue'
                    E(g)$label <- 1
                    
                    Pos_Edge_N <- 0
                    Neg_Edge_N <- 0
                    
                    for(i in 2:nrow(Corr_mat)){
                        for(j in 1:(i-1)){
                                    if ((Neg_Edge) & (Corr_mat[i,j] < -1 * Threshold) & (Corr_mat[i,j] > -1)){
                                            E(g, path=c(rownames(Corr_mat)[i],colnames(Corr_mat)[j]))$color <- 'green'
                                            E(g, path=c(rownames(Corr_mat)[i],colnames(Corr_mat)[j]))$weight <- round(abs(Corr_mat[i,j]),3)
                                            Neg_Edge_N <- Neg_Edge_N + 1
                                    }
                                    if ((Pos_Edge) & (Corr_mat[i,j] > Threshold) & (Corr_mat[i,j] < 1)){
                                            E(g, path=c(rownames(Corr_mat)[i],colnames(Corr_mat)[j]))$color <- 'red'
                                            E(g, path=c(rownames(Corr_mat)[i],colnames(Corr_mat)[j]))$weight <- round(abs(Corr_mat[i,j]),3)
                                            Pos_Edge_N <- Pos_Edge_N + 1
                                    }       
                            }
                    }
                    
                    map <- function(x, range = c(0,1), from.range=NA) {
                        if(any(is.na(from.range))) from.range <- range(x, na.rm=TRUE)
                        
                        ## check if all values are the same
                        if(!diff(from.range)) return(
                                matrix(mean(range), ncol=ncol(x), nrow=nrow(x), 
                                        dimnames = dimnames(x)))
                        
                        ## map to [0,1]
                        x <- (x-from.range[1])
                        x <- x/diff(from.range)
                        ## handle single values
                        if(diff(from.range) == 0) x <- 0 
                        
                        ## map from [0,1] to [range]
                        if (range[1]>range[2]) x <- 1-x
                        x <- x*(abs(diff(range))) + min(range)
                        
                        x[x<min(range) | x>max(range)] <- NA
                        
                        x
                    }
                    
                    
                    Node_num<-length(V(g))
                    Edge_num<-length(E(g))
                    Ave_degree<-mean(igraph::degree(g))
                    Degree_assortativity<-assortativity_degree(g)
                    Ave_path_length<-average.path.length(g)
                    Density<-graph.density(g)
                    Diameter<-diameter(g)
                    Radius<-radius(g)
                    Cluster_num<-clusters(g)$no
                    Modularity<-modularity(g, clusters(g)$membership)
                    Transitivity<-transitivity(g)
                    Closeness_centralization<-centralization.closeness(g)$centralization
                    Betweenness_centralization<-centralization.betweenness(g)$centralization
                    Degree_centralization<-centralization.degree(g)$centralization
                    Cluster_membership<-clusters(g)$membership
                    
                    
                    if(node_size=="degree" & var(igraph::degree(g))!=0) {vertex_size<-map(igraph::degree(g),c(1,10))
                    }else{
                    vertex_size<-node_size
                    }
                    
                    if(node_color=="degree" & var(igraph::degree(g))!=0){
                    vertex_color<-try(map(igraph::degree(g),c(1,10))) 
                    }else if (node_color=="Cluster_membership"){
                    vertex_color<-Cluster_membership
                    }else{
                    vertex_color<-node_color
                    }
                    
                    if(!is.null(outdir)){
                    png(outdir,width=4000,height=4000)
                    op<-par(mar=c(2,2,6,2))
                    str = paste(
                    "  Threshold: ",Threshold,
                    "  Nodes: ", Node_num,
                    "  Edges: ", Edge_num,
                    "  # of Pos edge: ",Pos_Edge_N,  
                    "  # of Neg edge: ",Neg_Edge_N,
                    "\nAverage_degree: ", Ave_degree,
                    "  Degree_assortativity: ", Degree_assortativity,
                    "  Clusters: ", Cluster_num, 
                    "  Modularity: ", Modularity,
                    "  Average_path_length: ",Ave_path_length,
                    "\nDensity: ", Density,
                    "  Diameter: ", Diameter,
                    "  Radius: ",Radius,  
                    "  Transitivity: ",Transitivity,
                    "\nCloseness_centralization:", Closeness_centralization,
                    "  Betweenness_centralization:", Betweenness_centralization,
                    "  Degree_centralization:", Degree_centralization
                    )
                    l<-get(layout)(g)
                    plot(g, layout=l,  main=str, edge.label=E(g)$weight, edge.label.cex = 0.25, vertex.frame.color="white", vertex.size=vertex_size, vertex.color=vertex_color, vertex.label.cex=0.5) 
                    if(Ave_degree!=0){
                    legend("bottomright", c("pos","neg"), lty=c(1,1), col=c("red","green") )
                    }
                    par(op)
                    invisible(dev.off())
                    }
                    objectList      <- list()
                    objectList$Node_num                      <- Node_num
                    objectList$Edge_num                      <- Edge_num 
                    objectList$Pos_Edge_num                  <- Pos_Edge_N 
                    objectList$Neg_Edge_num                  <- Neg_Edge_N 
                    objectList$Ave_degree                    <- Ave_degree
                    objectList$Degree_assortativity          <- Degree_assortativity 
                    objectList$Ave_path_length               <- Ave_path_length
                    objectList$Density                       <- Density
                    objectList$Diameter                      <- Diameter
                    objectList$Radius                        <- Radius
                    objectList$Cluster_num                   <- Cluster_num
                    objectList$Modularity                    <- Modularity
                    objectList$Transitivity                  <- Transitivity
                    objectList$Closeness_centralization      <- Closeness_centralization
                    objectList$Betweenness_centralization    <- Betweenness_centralization
                    objectList$Degree_centralization         <- Degree_centralization
                    objectList$g                             <- g
                    objectList$Cluster_membership            <- Cluster_membership
                    
                    invisible(objectList)
                               
}                              
                               
                               
Plot_local_chordDiagram<-function(Corr_mat, main_tem=main_tem,Pos_Edge=TRUE, Neg_Edge=TRUE, Threshold=0.6, bands_ann=NULL){

                   require("RColorBrewer")
                   require("circlize")
                   colours <- c(brewer.pal(9, "Set1"),brewer.pal(9, "Pastel1"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),brewer.pal(8, "Accent"),brewer.pal(8, "Dark2"), brewer.pal(8, "Pastel2"))

                   if (Pos_Edge & Neg_Edge){
                       col_mat=colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
                    }else if (Pos_Edge){
                       col_mat=colorRamp2(c(-1, -1e-10+Threshold, Threshold, 1), c("#FFFFFFFF", "#FFFFFFFF", "white", "red"))
                    }else if (Neg_Edge){
                       col_mat=colorRamp2(c(-1, -Threshold, -Threshold+1e-10, 1), c("blue", "white", "#FFFFFFFF", "#FFFFFFFF"))
                    }else stop('Please check the pos/neg edges parameters')

                   grid.col = colours[unclass=bands_ann$Group]
                   
                   op<-par(mar=c(0,0,2,0))
                   circos.par(gap.degree = c(rep(1,6),8,rep(1,6),8,rep(1,2),8)) # circle gap widths
                   chordDiagram(Corr_mat, symmetric = TRUE, transparency=0, col = col_mat, grid.col = grid.col, 
                                grid.border = "black", annotationTrack = "grid", preAllocateTracks = list(list(track.height = 0.02)))
                   circos.trackPlotRegion(track.index = 2, 
                                          panel.fun = function(x, y) {
                                          xlim = get.cell.meta.data("xlim")
                                          ylim = get.cell.meta.data("ylim")
                                          sector.index = get.cell.meta.data("sector.index")
                                          circos.text(mean(xlim), mean(ylim), sector.index, col = "black", cex = 0.5, 
                                              facing = "bending.inside", niceFacing = TRUE)
                                          }, 
                                          bg.border = NA)
                   invisible(lapply(levels(bands_ann$Group), function(x){
                                                 b<-paste("B", as.character(subset(bands_ann, Group==x)$Wave_num), sep="")
                                                 highlight.sector(sector.index = b, track.index = 1, col = "grey10", text = x, text.vjust = -1, niceFacing = TRUE, cex=0.8)
                                                 }
                                                 ))
                   title(main=main_tem,cex=2)
				   circos.clear()
                   par(op)
}


Trim_v_corr_mat<-function(v_Corr_mat, Pos_Edge=TRUE, Neg_Edge=TRUE, Threshold=0.6){
                 if (Pos_Edge & Neg_Edge){
                               v_Corr_mat <- subset(v_Corr_mat, value > Threshold | value < (-Threshold))
                           }else if (Pos_Edge){
                               v_Corr_mat <- subset(v_Corr_mat, value > Threshold & value < 1)
                           }else if (Neg_Edge){
                               v_Corr_mat <- subset(v_Corr_mat, value > -1 & value < (-Threshold))
                           }else stop('Please check the pos/neg edges parameters')
                  v_Corr_mat 
                  }
Plot_global_chordDiagram<-function(Corr_mat, main_tem=main_tem, Pos_Edge=TRUE, Neg_Edge=TRUE, Threshold=0.6){
                          require("circlize")
                          w_num<-sort(unique(as.numeric(gsub("[A-Z]", "", rownames(Corr_mat)))))
                          if(length(unique(w_num)) < length(w_num)) 
                          stop("Row/column names of input correlation matrix should be unique. Please check!")
                          bands_ord<-data.frame(Band=rownames(Corr_mat), ord=order(w_num))
                          v_Corr_mat<-vectorize_dm(Corr_mat, group=bands_ord$ord, duplicate=FALSE) 
                                                 
                          if (Pos_Edge & Neg_Edge){
                                                 v_Corr_mat <- subset(v_Corr_mat, value > Threshold | value < (-Threshold))
                                              }else if (Pos_Edge){
                                                 v_Corr_mat <- subset(v_Corr_mat, value > Threshold & value < 1)
                                              }else if (Neg_Edge){
                                                 v_Corr_mat <- subset(v_Corr_mat, value > -1 & value < (-Threshold))
                                              }else stop('Please check the pos/neg edges parameters')
                                                                     
                          op<-par(mar=c(0,0,2,0))
                          circos.par(cell.padding = c(0.02, 1.00, 0.02, 1.00))
                          circos.initialize("a", xlim = c(1,nrow(Corr_mat)))
                          circos.trackPlotRegion(ylim = c(0,1), bg.border = NA, track.height = 0.05)
                          circos.rect(seq(0, nrow(Corr_mat)-1), rep(0, nrow(Corr_mat)), 1:nrow(Corr_mat), rep(1, nrow(Corr_mat)), col = "#53A1DEFF", border=NA)
                          
                          breaks = c(1, seq(from=100, to=nrow(Corr_mat), by=100), nrow(Corr_mat))
                          #breaks = c(w_num[1], seq(from=round(w_num[1]+100, digit=-2), to=max(w_num), by=100), max(w_num))
                          circos.axis(h = "top", major.at = breaks, labels = paste0("B", breaks), major.tick.percentage = 0.5, labels.cex = 1, labels.away.percentage = 0.5)
                          col_mat=colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
                          if(nrow(v_Corr_mat)>1){
                          apply(v_Corr_mat, 1, function(i) { circos.link("a", as.numeric(i[3]), "a", as.numeric(i[5]), col = col_mat(as.numeric(i[6]))) })
                          }
						  title(main=main_tem,cex=2)
                          circos.clear()
                          par(op)

                          }
                         

init.graph <- function(data, dir=F) {
              labels<-union(unique(data[,1]), unique(data[,2]))
              ids<-1:length(labels); names(ids)<-labels
              from<-as.character(data[,1]); to<-as.character(data[,2]);
              edges<-matrix(c(ids[from], ids[to]), nc=2)
              g<-graph.empty(direct=dir)
              g<-add.vertices(g, length(labels))
              V(g)$label<-labels
              g<-add.edges(g, edges)
              
              g
}

#E:\RWAS\Network_20200426\Plot_local_chordDiagram_new_20200426.r
Plot_local_chordDiagram_new<-function(Corr_mat, a_degree,
                                      main_title,
                                      Pos_Edge=TRUE, Neg_Edge=TRUE, 
                                      Threshold=0.6, bands_ann=NULL){
  require("circlize")
  w_num<-sort(unique(as.numeric(gsub("[A-Z]", "", rownames(Corr_mat)))))
  if(length(unique(w_num)) < length(w_num)) 
    stop("Row/column names of input correlation matrix should be unique. Please check!")
  bands_ord<-data.frame(Band=rownames(Corr_mat), ord=order(w_num))
  v_Corr_mat<-vectorize_dm(Corr_mat, group=bands_ord$ord, duplicate=FALSE) 
  #v_Corr_mat_tem<-v_Corr_mat[1:20,]
  
  if (Pos_Edge & Neg_Edge){
    v_Corr_mat <- subset(v_Corr_mat, value > Threshold | value < (-Threshold))
  }else if (Pos_Edge){
    v_Corr_mat <- subset(v_Corr_mat, value > Threshold & value < 1)
  }else if (Neg_Edge){
    v_Corr_mat <- subset(v_Corr_mat, value > -1 & value < (-Threshold))
  }else stop('Please check the pos/neg edges parameters')
  
  
#  a<-data.frame(factors = rep("a",nrow(Corr_mat)),
#                x=as.numeric(gsub("[A-Z]", "", row.names(degree[,1]))),
#                y=degree[,group_name_degree]/max(degree))
  
  v_Corr_mat$Group_1<-as.numeric(gsub("[A-Z]", "", v_Corr_mat$Sample_1))
  v_Corr_mat$Group_2<-as.numeric(gsub("[A-Z]", "", v_Corr_mat$Sample_2))
  
  #pdf(paste(output,group_name,"+++.pdf",sep=''),width = 10,height = 10)
  #circos.clear()
  #par(mar = c(0, 0, 0, 0))
  par(mar = c(1, 1, 1, 1))
  circos.par("start.degree" = 90)
  raw_wn <- as.numeric(gsub("[A-Z]", "", colnames(Corr_mat)))
  circos.initialize("a", xlim = c(raw_wn[1]-1,raw_wn[length(raw_wn)]+0.1*(raw_wn[length(raw_wn)]-raw_wn[1]))) # 'a` just means there is one sector
  col_line= "grey20"#"#EAA041"#橙色
  col_mat=colorRamp2(c(-1, -Threshold, -Threshold+1e-10, 1), c("#2E3192", "#56B3DE", "#FFFFFFFF", "#FFFFFFFF"))
  circos.trackPlotRegion(ylim = c(0, 1),  track.height = 0.3,
                         #ylim = c(0, 0.1),  track.height = 0.05,
                         bg.border = NA, panel.fun = function(x, y) {
                           xlim = get.cell.meta.data("xlim") # in fact, it is c(0, 100)
                           
                           #circos.lines(a$x, a$y, type = "l",area = T,sector.index = "a",track.index = 1,
                           #             col = col_line,border = col_line)
                           #circos.polygon(a$x, a$y,sector.index = "a",track.index = 1,
                           #             col = col_line,border = col_line)
                           #circos.points(a$x, a$y, type = "l",area = T,sector.index = "a",track.index = 1,
                           #             col = col_line,border = col_line)
                           circos.trackLines(a_degree$factor, a_degree$x, a_degree$y, type = "h",area = T,
                                             col = col_line,border = col_line)
                           circos.text(xlim[1]-10, 0.3, expression('D'*'e'*'g'*'r'*'e'*'e'),
                                       facing = "downward", adj = c(1, 1),col=col_line,cex = 1)
                           
                           circos.lines(c(raw_wn[1],raw_wn[length(raw_wn)]), c(0, 0), col = "#CCCCCC")
                           circos.text(xlim[1]-10, 0.1, expression('S'*'h'*'i'*'f'*'t'~'('*'c'*'m'^-1*')'),
                                       facing = "downward", adj = c(1, 1),col="grey20",cex = 1)
                           
                           breaks = round(c(raw_wn[1],c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)*raw_wn[length(raw_wn)]),0)
                           circos.axis(h=0, major.at = breaks, labels =breaks,
                                       col="grey80",
                                       direction="inside",
                                       labels.facing="inside",
                                       major.tick.length=0.02,
                                       #lwd=1,
                                       #major.tick.percentage = 0.5, 
                                       #labels.away.percentage = 0,
                                       labels.col ="grey80",
                                       labels.cex = 0.7)
                           if(nrow(v_Corr_mat)>1){
                             apply(v_Corr_mat, 1, function(i) 
                             { circos.link("a", as.numeric(i[3]), 
                                           "a", as.numeric(i[5]), 
                                           col = col_mat(as.numeric(i[6]))) })
                           }
                         })
  #title(main=main_title,cex=2)
  circos.clear()
#  par(op)
  #dev.off()
}

#require(devtools)
#install_version("data.table", version = "1.11.8", repos = "http://cran.us.r-project.org")
#install.packages("/mnt/data6/heyh/Network_20200514/Hmisc_4.4-0.tar.gz")

#---------------------------------------------------------------------------------
## rcorr_df() replace cor() for calculating both r and p-value
rcorr_df<-function(df){
  #Installing older versions of packages
  # for data.table pakage, only 1.11.8 version  works well.
  #require(devtools)
  #install_version("data.table", version = "1.11.8", repos = "http://cran.us.r-project.org")
  #install.packages("/mnt/data6/heyh/Network_20200514/Hmisc_4.4-0.tar.gz")
  #require(Hmisc)
  library(Hmisc)
  df_mat<-as.matrix(df)
  return(Hmisc::rcorr(df_mat))
}
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
## vectorize_dm both r and p-value
vectorize_dm_rcorr<-function(dm, group=NULL, duplicate=TRUE)
{
  #dm<-cor_mats_tem[[1]]
  dm_r<-data.matrix(dm$r)
  dm_pvalue<-data.matrix(dm$P)
  
  dm_r_value<-vectorize_dm(dm_r, group=NULL, duplicate=FALSE)
  dm_p_value<-vectorize_dm_pvalue(dm_pvalue, group=NULL, duplicate=FALSE)
  dm_value<-cbind(dm_r_value,dm_p_value[,"value"])
  colnames(dm_value)[ncol(dm_value)]<-"p_value"
  dm_value
}
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
## vectorize_dm both r and p-value
vectorize_dm_rcorr_new<-function(dm, group=NULL, duplicate=TRUE)
{
  #dm<-cor_mats_tem[[1]]
  dm_r<-data.matrix(dm[[1]])#dm$r
  dm_pvalue<-data.matrix(dm[[3]])#dm$P
  
  dm_r_value<-vectorize_dm(dm_r, group=NULL, duplicate=FALSE)
  dm_p_value<-vectorize_dm_pvalue(dm_pvalue, group=NULL, duplicate=FALSE)
  dm_value<-cbind(dm_r_value,dm_p_value[,"value"])
  colnames(dm_value)[ncol(dm_value)]<-"p_value"
  dm_value
}
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
## vectorize_dm_pvalue 
vectorize_dm_pvalue<-function(dm, group=NULL, duplicate=TRUE)
{
  #if(ncol(dm)!=nrow(dm) & any(is.na(dm))==TRUE)
  if(ncol(dm)!=nrow(dm))
    stop('The distance matrix is not squared')
  dm<-data.matrix(dm)
  require("reshape2")
  if(!is.null(group)){
    if( length(unique(group))==1)
      stop('At least two levels for a given sample category in your metadata file are required.')
    if( length(group)!=nrow(dm))
      stop('The number of rows in metadata and distance matrix are not equal')
    if(is.factor(group) & nlevels(group)>length(group)*0.9)
      stop('The number of levels in a certain category can not exceed 90% of total number of samples')
    
    colnames(dm)<-rownames(dm)<-paste(rownames(dm), group, sep="____")
    if(duplicate){ 
      #melt_dm<-subset(melt(dm), value!=0 & value!=1)
      melt_dm<-melt(dm) 
    }else{ 
      dm[lower.tri(dm)]<-NA
      #melt_dm<-subset(melt(dm), ! is.na(value) & value!=0 & value!=1) 
      melt_dm<-subset(melt(dm), ! is.na(value))
    }
    
    Row_Info<-data.frame(do.call(rbind,strsplit(as.character(melt_dm[,1]),"____",fixed=TRUE)))
    Col_Info<-data.frame(do.call(rbind,strsplit(as.character(melt_dm[,2]),"____",fixed=TRUE)))
    VS<-paste(Row_Info[,2],"_VS._",Col_Info[,2],sep="")
    dm_value<-data.frame(VS,Row_Info,Col_Info,d=melt_dm$value)
    
    colnames(dm_value)<-c("GroupPair","Sample_1","Group_1","Sample_2","Group_2","value")
    if(is.factor(group)){
      DistType<-as.factor(dm_value$Group_1==dm_value$Group_2)
      DistType<-factor(DistType,levels=levels(DistType),labels=c("AllBetween","AllWithin"))
      dm_value<-data.frame(DistType,dm_value)}
  }else{
    if (duplicate) { 
      #dm_value<-subset(melt(dm), value!=0 & value!=1) 
      dm_value<-melt(dm)
    }else{ 
      dm[lower.tri(dm)] <- NA
      #dm_value <- subset(melt(dm), ! is.na(value) & value!=0 & value!=1)  }
      dm_value <- subset(melt(dm), ! is.na(value))  }
    colnames(dm_value)<-c("Sample_1","Sample_2","value")
  }
  dm_value
}
#---------------------------------------------------------------------------------

#---------------------------------------------------
# plot_local_chorddiagram  (connectedness all & p.value<=0.05)
#---------------------------------------------------
Plot_local_chordDiagram_rpvalue_new<-function(x, Pos_Edge=TRUE, Neg_Edge=TRUE, Threshold=0.6, bands_ann=NULL){
  # x<-"PMMA"
  # Pos_Edge=F
  # Neg_Edge=T
  # Threshold=0.6
  # bands_ann=local_bands_ann_new
  
  #Corr_mats_rpvalue<-bands_cor_mats_rpvalue[1]
  #Corr_mat<-Corr_mats_rpvalue[[1]][[1]] # r
  #Corr_mat_pvalue<-Corr_mats_rpvalue[[1]][[2]] # p.value
  #x<-"Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817"
  #x<-"Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20151101"
  #x<-"PMMA"
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  #Bands_label<-round(as.numeric(gsub("B","",colnames(Corr_mat))),0)
  Bands_label<-bands_ann$Wave_num
  colnames(Corr_mat)<-rownames(Corr_mat)<-Bands_label
  colnames(Corr_mat_pvalue)<-rownames(Corr_mat_pvalue)<-Bands_label
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-1.1 # p.value<0.05
  #Corr_mat[which(Corr_mat>=Threshold)]<-1.1 # r cutoff
  #Corr_mat[which(Corr_mat==0)]<-1.1
  
  require("RColorBrewer")
  require("circlize")
  colours <- c(brewer.pal(9, "Set1"),brewer.pal(9, "Pastel1"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"),brewer.pal(8, "Accent"),brewer.pal(8, "Dark2"), brewer.pal(8, "Pastel2"))
  #colours <- c("#339900","#9900FF","#FF6666","#CC6600")#green,purple,pink,yellow
  
  #if (Pos_Edge & Neg_Edge){
  #  col_mat=colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  #}else if (Pos_Edge){
  #  col_mat=colorRamp2(c(-1, -1e-10+Threshold, Threshold, 1), c("#FFFFFFFF", "#FFFFFFFF", "white", "red"))
  #}else if (Neg_Edge){
  #  col_mat=colorRamp2(c(-1, -Threshold, -Threshold+1e-10, 1), c("blue", "white", "#FFFFFFFF", "#FFFFFFFF"))
  #}else stop('Please check the pos/neg edges parameters')
  
  if (Pos_Edge & Neg_Edge){
    col_mat=colorRamp2(c(-1, 0, 1), c("Blue4", "White", "Red2"))
    Corr_mat[which(Corr_mat<Threshold&Corr_mat>-(Threshold))]<-1.1 # r cutoff
    
  }else if (Pos_Edge){
    col_mat=colorRamp2(c(-1, -1e-10+Threshold, Threshold, 1), c("#FFFFFF00", "#FFFFFF00", "Red", "DarkRed"))
    Corr_mat[which(Corr_mat<Threshold)]<-1.1 # r cutoff
    
  }else if (Neg_Edge){
    col_mat=colorRamp2(c(-1, -Threshold, -Threshold+1e-10, 1), c("DarkBlue", "Blue", "#FFFFFF00", "#FFFFFF00"))
    Corr_mat[which(Corr_mat>(-Threshold))]<-1.1 # r cutoff
    
  }else stop('Please check the pos/neg edges parameters')
  col_link<-col_mat(Corr_mat)
  col_link[which(Corr_mat==1.1)]<-"#FFFFFF00"
  
  bands_ann<-bands_ann[order(bands_ann$Group),]
  grid.col = colours[unclass=bands_ann$Group]
  
  op<-par(mar=c(0,0,0,0))
  
  n_bands_Group<-as.numeric(lapply(levels(bands_ann$Group), function(x){nrow(bands_ann[which(bands_ann$Group==x),])}))-1
  x_vec<-c()
  for (i in n_bands_Group) {
    tem<-rep(1,i)
    x_vec<-c(x_vec,tem,8)
  }
  #device="png"
  #require(Cairo)
  #Cairo(file = paste(outpath_png, "Lobal_chordDiagram_", x, ".",device , sep=""), 
  #      unit="in", dpi=300, width=6, height=6, type=device, bg="white")
  circos.clear()
  circos.par(gap.degree = x_vec) # circle gap widths
  
  chordDiagram(Corr_mat, order=as.character(bands_ann$Wave_num),symmetric = T, 
               #  chordDiagram(Corr_mat, order=paste("B", bands_ann$Wave_num,sep = ""),symmetric = T, 
               link.visible=T,
               keep.diagonal = TRUE,
               #link.border=2,
               #link.lwd=4,
               #grid.border=10,
               #scale=T,
               transparency=0, 
               #col = col_mat, 
               col = col_link, 
               grid.col = grid.col, 
               #grid.col = "White", 
               grid.border = NULL, annotationTrack = "grid",
               #grid.border = "black", annotationTrack = "grid")
               preAllocateTracks = list(list(track.height = 0.02)))
  #preAllocateTracks = list(list(track.height = 0.5)))
  circos.trackPlotRegion(ylim = c(0, 1),track.index = 2, track.height = 0.3,
                         #bg.border = NA,
                         #bg.col="red",
                         panel.fun = function(x, y) {
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           sector.index = get.cell.meta.data("sector.index")
                           circos.text(mean(xlim), mean(ylim), sector.index, col = "White", font=2, cex = 0.4, 
                                       facing = "bending.inside", niceFacing = TRUE)
                         }, 
                         bg.border = NA)
  invisible(lapply(levels(bands_ann$Group), function(x){
    #b<-paste("B", as.character(subset(bands_ann, Group==x)$Wave_num), sep="")
    b<-as.character(subset(bands_ann, Group==x)$Wave_num)
    #highlight.sector(sector.index = b, track.index = 1, col = "grey10", text = x, 
    #                 text.vjust = -1, niceFacing = TRUE, cex=0.8)
    #circos.trackPlotRegion(ylim = c(0, 1),track.index = 1, track.height = 0.01,bg.col = NA, bg.border = NA )
    highlight.sector(sector.index = b, 
                     track.index = 1, 
                     #col = "grey10",
                     col = "grey80",
                     #col = grid.col,					 
                     lwd = 0.001,
                     #lwd = 10,
                     text = x,
                     text.col="grey30",					 
                     text.vjust = -0.5,
                     facing = "bending.inside", 
                     niceFacing = TRUE, cex=0.6)
  }
  ))
  circos.clear()
  par(op)
  #dev.off()
  #dev.off()
}
#
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------
# plot_global_chorddiagram  (connectedness all & p.value<=0.05)
#---------------------------------------------------
#E:\RWAS\Network_20200426\Plot_local_chordDiagram_new_20200426.r
#E:\RWAS\Fig_Cosmetic\Fig6
Plot_global_chordDiagram_rpvalue<-function(Corr_mat, 
                                           a_degree=NULL,
                                           a_SumCorr=NULL,
                                           a_MeanCorr=NULL,
                                           main_title=NULL,
                                           Pos_Edge=FALSE, 
                                           Neg_Edge=FALSE, 
                                           Threshold=0.6, 
                                           bands_ann=NULL){
  require("circlize")
  w_num<-sort(unique(as.numeric(gsub("[A-Z]", "", rownames(Corr_mat)))))
  if(length(unique(w_num)) < length(w_num)) 
    stop("Row/column names of input correlation matrix should be unique. Please check!")
  bands_ord<-data.frame(Band=rownames(Corr_mat), ord=order(w_num))
  v_Corr_mat<-vectorize_dm(Corr_mat, group=bands_ord$ord, duplicate=FALSE) 
  #v_Corr_mat_tem<-v_Corr_mat[1:20,]
  
  if (Pos_Edge & Neg_Edge){
    v_Corr_mat <- subset(v_Corr_mat, value > Threshold | value < (-Threshold))
  }else if (Pos_Edge){
    v_Corr_mat <- subset(v_Corr_mat, value > Threshold & value < 1)
  }else if (Neg_Edge){
    v_Corr_mat <- subset(v_Corr_mat, value > -1 & value < (-Threshold))
  }else stop('Please check the pos/neg edges parameters')
 
  
  v_Corr_mat$Group_1<-as.numeric(gsub("[A-Z]", "", v_Corr_mat$Sample_1))
  v_Corr_mat$Group_2<-as.numeric(gsub("[A-Z]", "", v_Corr_mat$Sample_2))
  
  #pdf(paste(output,group_name,"+++.pdf",sep=''),width = 10,height = 10)
  #circos.clear()
  #par(mar = c(0, 0, 0, 0))
  par(mar = c(1, 1, 1, 1))
  circos.clear()
  circos.par("start.degree" = 90)
  raw_wn <- as.numeric(gsub("[A-Z]", "", colnames(Corr_mat)))
  circos.initialize("a", xlim = c(raw_wn[1]-1,raw_wn[length(raw_wn)]+0.1*(raw_wn[length(raw_wn)]-raw_wn[1]))) # 'a` just means there is one sector
   if (Pos_Edge & Neg_Edge){
    col_line="grey20"
  }else if (Pos_Edge){
  col_line="#FF6666"
  }else if (Neg_Edge){
  col_line="#66CCFF"
  }else stop('Please check the pos/neg edges parameters')
  #col_line= "grey20"#"#EAA041"#橙色
  col_mat=colorRamp2(c(-1, -Threshold, -Threshold+1e-10,Threshold-1e-10,Threshold, 1), c("#0033FF", "#66CCFF","#FFFFFF00", "#FFFFFF00", "#FF6666", "#FF3300"))
  #col_mat=colorRamp2(c(-1, -Threshold, -Threshold+1e-10,Threshold-1e-10,Threshold, 1), c("DarkBlue", "Blue", "#FFFFFF00", "#FFFFFF00", "Red", "DarkRed"))
  circos.trackPlotRegion(ylim = c(0, 1),  track.height = 0.3,
                         #ylim = c(0, 0.1),  track.height = 0.05,
                         bg.border = NA, panel.fun = function(x, y) {
                           xlim = get.cell.meta.data("xlim") # in fact, it is c(0, 100)
                           
                           #circos.lines(a$x, a$y, type = "l",area = T,sector.index = "a",track.index = 1,
                           #             col = col_line,border = col_line)
                           #circos.polygon(a$x, a$y,sector.index = "a",track.index = 1,
                           #             col = col_line,border = col_line)
                           #circos.points(a$x, a$y, type = "l",area = T,sector.index = "a",track.index = 1,
                           if (!is.null(a_degree)){
                             circos.trackLines(a_degree$factor, a_degree$x, a_degree$y, type = "h",area = T,
                                             col = col_line,border = col_line)
                             circos.text(xlim[1]-10, 0.3, expression('D'*'e'*'g'*'r'*'e'*'e'),
                                       facing = "downward", adj = c(1, 1),col=col_line,cex = 1)
                           }
                           if (!is.null(a_SumCorr)){
                              circos.trackLines(a_SumCorr$factor, a_SumCorr$x, a_SumCorr$y, type = "h",area = T,
                                             col = col_line,border = col_line)
                              circos.text(xlim[1]-10, 0.3, expression('S'*'u'*'m'*'C'*'o'*'r'*'r'),
                                       facing = "downward", adj = c(1, 1),col=col_line,cex = 1)
                           } 
                           if (!is.null(a_MeanCorr)){
                              circos.trackLines(a_MeanCorr$factor, a_MeanCorr$x, a_MeanCorr$y, type = "h",area = T,
                                             col = col_line,border = col_line)
                              circos.text(xlim[1]-10, 0.3, expression('M'*'e'*'a'*'n'*'C'*'o'*'r'*'r'),
                                       facing = "downward", adj = c(1, 1),col=col_line,cex = 1)
                           }
						   
                           circos.lines(c(raw_wn[1],raw_wn[length(raw_wn)]), c(0, 0), col = "#CCCCCC")
                           circos.text(xlim[1]-10, 0.1, expression('S'*'h'*'i'*'f'*'t'~'('*'c'*'m'^-1*')'),
                                       facing = "downward", adj = c(1, 1),col="grey60",cex = 0.8)
                           
                           breaks = round(c(raw_wn[1],c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)*raw_wn[length(raw_wn)]),0)
                           circos.axis(h=0, major.at = breaks, labels =breaks,
                                       col="grey50",
                                       #col="grey80",
                                       direction="inside",
                                       labels.facing="inside",
                                       major.tick.length=0.02,
                                       #lwd=1,
                                       #major.tick.percentage = 0.5, 
                                       #labels.away.percentage = 0,
                                       labels.col ="grey50",
                                       #labels.col ="grey80",
                                       labels.cex = 0.7)
                           if(nrow(v_Corr_mat)>1){
                             apply(v_Corr_mat, 1, function(i) 
                             { circos.link("a", as.numeric(i[3]), 
                                           "a", as.numeric(i[5]), 
                                           col = col_mat(as.numeric(i[6]))) })
                           }
                         })
  #title(main=main_title,cex=2)
  circos.clear()
  #par(op)
  #dev.off()
}
#---------------------------------------------------------------------------------------------------------








