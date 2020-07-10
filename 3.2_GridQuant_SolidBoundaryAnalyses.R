
# 3.2_GridQuant_SolidBoundaryAnalyses.R
# Analyses of immune cell infiltration in solid pattern as a function of the distance from boundary

# Author: Daniele Tavernari
# Please cite: Tavernari et al., 2021 (see full citation on GitHub repository's README)

library(reshape2)
library(ggplot2)
library(ggpubr)

########## Input
InDir_matrices = "Grid_matrices/"
InDir_boundaries = "All_boundaries/"
OutDir = "SolidBoundaryAnalyses/"
s = "s8B"
grid_spacing = 100
cellTypes = c("Bcell", "CD4", "CD8", "macrophages", "Ki67")
####################

########## Other parameters
cellTypes_colors_df = data.frame(row.names = cellTypes, cellType = cellTypes, color = c("cyan3", "blue4", "goldenrod1", "red3","green3"), stringsAsFactors = F )
normal_color = "gray"
lepidic_color = "dodgerblue4"
acinar_color = "orange"
papillary_color = "lightseagreen"
solid_color = "red"
PatternColorsDf = data.frame(row.names = c("normal","lepidic", "papillary", "acinar", "solid"), Pattern = c("normal","lepidic", "papillary", "acinar", "solid"), color = c(normal_color,lepidic_color,papillary_color,acinar_color,solid_color))
crown_steps_inside = c(0,500,1000)
crown_steps_outside = c(0,500)
####################

################## FUNCTIONS ##################
compute_shortest_distance = function(meaMat, col, row, all_x_coos, all_y_coos){

   xp = all_x_coos[col]
   yp = all_y_coos[row]
   nas = data.frame(which(is.na(meaMat), arr.ind=T))
   finalDist = Inf
   for (rn in rownames(nas))
   {
      xna = all_x_coos[nas[rn,"col"]]
      yna = all_x_coos[nas[rn,"row"]]
      this_dist = sqrt((xp-xna)^2+(yp-yna)^2)
      if (this_dist<finalDist)
      {
         finalDist = this_dist
      }
   }
   return( finalDist )
}

correct_erroneous_padding_from_top = function(meaMat_is_in_boundary){

  for ( cn in 1:ncol(meaMat_is_in_boundary) )
  {
    for (rn in 1:nrow(meaMat_is_in_boundary))
    {
      if (is.na(meaMat_is_in_boundary[rn,cn]) ) 
      { 
        break
      } else {
        meaMat_is_in_boundary[rn,cn] = 0
      }
    }
    
  }

  return(meaMat_is_in_boundary)

}

compute_distances_remaining_markers = function(distdf, meaMat, ct){

   for (rn in rownames(distdf))
   {
      distdf[rn,ct] = meaMat[distdf[rn,"y"],distdf[rn,"x"]]
   }
   return(distdf)
}

valid_changePoints = function( meaMat, row, col, mode = "ClosedBoundary" ){
   v = 1*is.na(meaMat[row,col:ncol(meaMat)])
   cp_after = (sum(v[c(TRUE, !v[-length(v)] == v[-1])]) %% 2) == 1
   v = 1*is.na(meaMat[row,1:col])
   cp_before = (sum(v[c(TRUE, !v[-length(v)] == v[-1])]) %% 2) == 1
   return (cp_before & cp_after)
}

pad_and_correct = function(meaMat_is_in_boundary){

   meaMat_is_in_boundary[1,!is.na(meaMat_is_in_boundary[1,])] = 0
   meaMat_is_in_boundary[!is.na(meaMat_is_in_boundary[,1]),1] = 0
   meaMat_is_in_boundary[nrow(meaMat_is_in_boundary),!is.na(meaMat_is_in_boundary[nrow(meaMat_is_in_boundary),])] = 0
   meaMat_is_in_boundary[!is.na(meaMat_is_in_boundary[,ncol(meaMat_is_in_boundary)]),ncol(meaMat_is_in_boundary)] = 0
   thism = meaMat_is_in_boundary
   maxRow = (nrow(thism)-1)
   maxCol = (ncol(thism)-1)
   for ( row in 2:maxRow )
   {
      for ( col in 2:maxCol )
      {
         if (((thism[row,col]==0) %in% c(T)) & (( sum(thism[row-1,(col-1):(col+1)],na.rm=T)>0   ) )) # & ( sum(thism[row+1,(col-1):(col+1)],na.rm=T)>0   )
         { 
            thism[row,col] = 1
         }
         if (((thism[row,col]==1) %in% c(T)) & (( sum(thism[row-1,(col-1):(col+1)],na.rm=T)==0   ) & ( sum(thism[row+1,(col-1):(col+1)],na.rm=T)==0   ))) { thism[row,col] = 0 }
      }
   }
   for ( row in maxRow:2 )
   {
      for ( col in maxCol:2 )
      {
         if (((thism[row,col]==0) %in% c(T)) & (( sum(thism[row+1,(col-1):(col+1)],na.rm=T)>0   ) )) # & ( sum(thism[row+1,(col-1):(col+1)],na.rm=T)>0   )
         { 
            thism[row,col] = 1
         }
         if (((thism[row,col]==1) %in% c(T)) & (( sum(thism[row-1,(col-1):(col+1)],na.rm=T)==0   ) & ( sum(thism[row+1,(col-1):(col+1)],na.rm=T)==0   ))) { thism[row,col] = 0 }
      }
   }
   return (thism)
}

compute_distances_ClosedBoundary_first_marker = function(meaMat, meaMat_uncensored, all_x_coos, all_y_coos, ct, polyline_file, fileRoot){

   last_row = nrow(meaMat)
   meaMat_is_in_boundary = matrix(data = 0, nrow = nrow(meaMat), ncol = ncol(meaMat))
   distdf = data.frame(matrix(nrow = 0, ncol = 8, dimnames = list(c(),c("x","y","dist","macrophages", "Ki67", "CD4", "CD8", "Bcell"))))
   for (row in 1:last_row)
   {
      # if (row==round(nrow(meaMat)/10) ) { cat('\n',"10% done",'\n')}
      start = which(is.na(meaMat[row,]))
      if (length(start)==0) { next }
      startx = start[1]+1

      for (col in startx:ncol(meaMat))
      {
         if (is.na(meaMat[row,col])) { next }
         if ( sum(is.na(meaMat[row:nrow(meaMat),col]))==0 ) { next }
         if ( !(valid_changePoints( meaMat, row, col, mode = "ClosedBoundary" )) ) { next }
         meaMat_is_in_boundary[row,col] = 1
      }
   }

   meaMat_is_in_boundary = polyLine_censoring(polyline_file, meaMat_is_in_boundary, all_x_coos, all_y_coos, do_close_boundary = T)
   meaMat_is_in_boundary = pad_and_correct(meaMat_is_in_boundary)
   
   meaMat_is_in_boundary = correct_erroneous_padding_from_top( meaMat_is_in_boundary ) # use only ad-hoc
   

   fileName = paste0(fileRoot,"_meaMat_is_in_boundary_after.pdf")
   pdf(fileName, 10, 10, useDingbats = F )
   df = melt(meaMat_is_in_boundary)
   df[(df$value<0) %in% c(T),"value"] = NA
   levelz = levels(factor(df$Var1))
   p=ggplot(df) + theme_bw() + 
     geom_tile(aes(Var2,ordered(Var1, levels=rev(levelz)),fill=value))+
     scale_fill_gradientn(colours=colorRampPalette(c("white","black" ))(n = 100))+
     theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0))+
     labs(title=paste0(""), x="",y="")+
     coord_fixed()
   print(p)
   dev.off()
   meaMat_is_in_boundary[is.na(meaMat_is_in_boundary)] = 1

   distdf = data.frame(matrix(nrow = sum(meaMat_is_in_boundary), ncol = 8, dimnames = list(c(1:sum(meaMat_is_in_boundary)),c("x","y","dist","macrophages", "Ki67", "CD4", "CD8", "Bcell"))))

   indices = which(meaMat_is_in_boundary == 1, arr.ind = TRUE)
   for (index in 1:nrow(indices))
   {
      distdf[index,"x"] = as.numeric(indices[index,"col"])
      distdf[index,"y"] = as.numeric(indices[index,"row"])
      # distdf[index,"dist"] = compute_shortest_distance(meaMat, distdf[index,"x"], distdf[index,"y"], all_x_coos, all_y_coos)
      distdf[index,ct] = meaMat_uncensored[distdf[index,"y"],distdf[index,"x"]]
   }

   distdf = distdf[!is.na(distdf[,ct]),]

   return( distdf )
}



intersect_segment_tile = function(px1,px2,py1,py2,xs1,xs2,ys1,ys2){
   # upper segment
   x_intup = (ys1-py1)*(px2-px1)/(py2-py1)+px1
   x_intbottom = (ys2-py1)*(px2-px1)/(py2-py1)+px1
   y_right = (xs2-px1)*(py2-py1)/(px2-px1)+py1
   if ( ((x_intup<xs2) & (x_intup>xs1)) %in% c(T) ){ 
      return(TRUE)
      } else if ( ((x_intbottom<xs2) & (x_intbottom>xs1)) %in% c(T) ) {
         return(TRUE)
      } else if ( ((y_right<ys2) & (y_right>ys1)) %in% c(T) ) {
         return(TRUE)
      } else {
         return(FALSE)
      }
}

polyLine_censoring = function(polyline_file, meaMat, all_x_coos, all_y_coos, do_close_boundary = FALSE){

   hspac = round((all_x_coos[2]-all_x_coos[1])/2)
   pl = read.table(file = polyline_file, sep = "-")
   pixelSize = as.numeric(substr(as.character(pl[2,"V1"]),10,nchar(as.character(pl[2,"V1"]))))
   pl_points = unlist(strsplit(as.character(pl[1,"V1"]), split=";"))
   plp_df = data.frame(x = rep(NA,length(pl_points)+1*do_close_boundary), y = rep(NA,length(pl_points)+1*do_close_boundary))
   index = 1
   for (point in pl_points)
   {
      poiSpl = unlist(strsplit(point, split = ", "))
      plp_df[index,"x"] = as.numeric(substr(poiSpl[1],8,nchar(poiSpl[1])))*pixelSize
      plp_df[index,"y"] = as.numeric((poiSpl[2]))*pixelSize
      index = index+1
   }

   if ( do_close_boundary )
   {
      plp_df[nrow(plp_df),"x"] = plp_df[1,"x"]
      plp_df[nrow(plp_df),"y"] = plp_df[1,"y"]
   }
   
   for ( segment in c(1:(nrow(plp_df)-1) ) )
   {
      px1 = plp_df[segment,"x"]
      py1 = plp_df[segment,"y"]
      px2 = plp_df[segment+1,"x"]
      py2 = plp_df[segment+1,"y"]
      xs = all_x_coos[( all_x_coos<(max(c(px1,px2))+hspac) ) & ( all_x_coos>(min(c(px1,px2))-hspac) )]
      ys = all_y_coos[( all_y_coos<(max(c(py1,py2))+hspac) ) & ( all_y_coos>(min(c(py1,py2))-hspac) )]
      for (this_xs in xs)
      {
         for (this_ys in ys)
         {
            if (( intersect_segment_tile(px1,px2,py1,py2,this_xs-hspac,this_xs+hspac,this_ys-hspac,this_ys+hspac) )) { meaMat[which(all_y_coos==this_ys),which(all_x_coos==this_xs)] = NA }
         }
      }
      
   }

   return( meaMat )

}

distance_squared = function(start, end){
   return((start[1]-end[1])^2+(start[2]-end[2])^2)
}

compute_shortest_distance_limited = function(xquery, yquery, plp_df){
   finalDist = Inf
   start = c(plp_df[1,"x"],plp_df[1,"y"])
   point = c(xquery,yquery)
   for (rn in rownames(plp_df)[2:nrow(plp_df)])
   {
      end = c(plp_df[rn,"x"],plp_df[rn,"y"])
      d2 = distance_squared(start, end)
      if (d2 == 0){
         this_dist = sqrt(distance_squared(point, start))
      } else {
         t = ((point[1] - start[1]) * (end[1] - start[1]) + (point[2] - start[2]) * (end[2] - start[2]))/d2
         if (t < 0){
            this_dist = sqrt(distance_squared(point, start))
         }
        else if (t > 1.0){
            this_dist = sqrt(distance_squared(point, end))
         }
        else {
             proj = c(start[1] + t * (end[1] - start[1]), start[2] + t * (end[2] - start[2]))
             this_dist = sqrt(distance_squared(point, proj))
        }   
      }
      start = end
      if (is.na(this_dist)) { this_dist = Inf }
      if (this_dist<finalDist)
      {
         finalDist = this_dist
      }
   }
   return( finalDist )
}

compute_distances_TrueBoundary = function(sdf, polyline_file, grid_spacing, closed_boundary = F){
   tbDist = sdf
   pl = read.table(file = polyline_file, sep = "-")
   pixelSize = as.numeric(substr(as.character(pl[2,"V1"]),10,nchar(as.character(pl[2,"V1"]))))
   pl_points = unlist(strsplit(as.character(pl[1,"V1"]), split=";"))
   if (closed_boundary) { pl_points = c(pl_points, pl_points[1]) }
   plp_df = data.frame(x = rep(NA,length(pl_points)), y = rep(NA,length(pl_points)), stringsAsFactors=F)
   index = 1
   for (point in pl_points)
   {
      poiSpl = unlist(strsplit(point, split = ", "))
      plp_df[index,"x"] = ceiling(as.numeric(substr(poiSpl[1],8,nchar(poiSpl[1])))*pixelSize/grid_spacing)
      plp_df[index,"y"] = ceiling(as.numeric((poiSpl[2]))*pixelSize/grid_spacing)
      index = index+1
   }
   for (rn in rownames(tbDist))
   {
      tbDist[rn,"dist"] = compute_shortest_distance_limited(xquery = tbDist[rn,"x"], yquery = tbDist[rn,"y"], plp_df = plp_df )*grid_spacing
   }
   return( tbDist )
}

crowns_analysis = function(crown_steps, tbDist, s, meaMat_censored, mode){
   crowndf = data.frame(matrix(nrow=0,ncol=3,dimnames=list(NULL,c("crown","cellType","density"))),stringsAsFactors=F)
   from = 0
   for (to in crown_steps[2:length(crown_steps)])
   {
      resetted_meaMat = meaMat_censored
      resetted_meaMat[!is.na(resetted_meaMat)] = 0
      this_tbDist = tbDist[(tbDist$dist>=from) & (tbDist$dist<to),]
      crowndf = rbind(crowndf,data.frame(crown=paste0(from,"_",to,"_",mode),cellType="macrophages",density=this_tbDist[,"macrophages"]/((grid_spacing/1000)^2),x=this_tbDist$x,y=this_tbDist$y,stringsAsFactors=F))
      crowndf = rbind(crowndf,data.frame(crown=paste0(from,"_",to,"_",mode),cellType="Bcell",density=this_tbDist[,"Bcell"]/((grid_spacing/1000)^2),x=this_tbDist$x,y=this_tbDist$y,stringsAsFactors=F))
      crowndf = rbind(crowndf,data.frame(crown=paste0(from,"_",to,"_",mode),cellType="CD4",density=this_tbDist[,"CD4"]/((grid_spacing/1000)^2),x=this_tbDist$x,y=this_tbDist$y,stringsAsFactors=F))
      crowndf = rbind(crowndf,data.frame(crown=paste0(from,"_",to,"_",mode),cellType="CD8",density=this_tbDist[,"CD8"]/((grid_spacing/1000)^2),x=this_tbDist$x,y=this_tbDist$y,stringsAsFactors=F))
      crowndf = rbind(crowndf,data.frame(crown=paste0(from,"_",to,"_",mode),cellType="Ki67",density=this_tbDist[,"Ki67"]/((grid_spacing/1000)^2),x=this_tbDist$x,y=this_tbDist$y,stringsAsFactors=F))
      resetted_meaMat[cbind(this_tbDist$y,this_tbDist$x)] = resetted_meaMat[cbind(this_tbDist$y,this_tbDist$x)]+1
      fileName = paste0(OutDir,"crowns_masks/",s,"_Crown_from",from,"_to",to,"_mask_gridSpacing",grid_spacing,"_",mode,".pdf")
      pdf(fileName, 10, 10, useDingbats = F )
      df = melt(resetted_meaMat)
      # df[(df$value<0) %in% c(T),"value"] = NA
      levelz = levels(factor(df$Var1))
      p=ggplot(df) + theme_bw() + 
        geom_tile(aes(Var2,ordered(Var1, levels=rev(levelz)),fill=value))+
        scale_fill_gradientn(colours=colorRampPalette(c("white","black" ))(n = 100))+
        theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0))+
        labs(title=paste0(""), x="",y="")+
        coord_fixed()
      print(p)
      dev.off()
      from = to
   }
   if (mode == "inside")
   {
      resetted_meaMat = meaMat_censored
      resetted_meaMat[!is.na(resetted_meaMat)] = 0
      this_tbDist = tbDist[(tbDist$dist>=to),]
      crowndf = rbind(crowndf,data.frame(crown=paste0(from,"_core","_",mode),cellType="macrophages",density=this_tbDist[,"macrophages"]/((grid_spacing/1000)^2),x=this_tbDist$x,y=this_tbDist$y,stringsAsFactors=F))
      crowndf = rbind(crowndf,data.frame(crown=paste0(from,"_core","_",mode),cellType="Bcell",density=this_tbDist[,"Bcell"]/((grid_spacing/1000)^2),x=this_tbDist$x,y=this_tbDist$y,stringsAsFactors=F))
      crowndf = rbind(crowndf,data.frame(crown=paste0(from,"_core","_",mode),cellType="CD4",density=this_tbDist[,"CD4"]/((grid_spacing/1000)^2),x=this_tbDist$x,y=this_tbDist$y,stringsAsFactors=F))
      crowndf = rbind(crowndf,data.frame(crown=paste0(from,"_core","_",mode),cellType="CD8",density=this_tbDist[,"CD8"]/((grid_spacing/1000)^2),x=this_tbDist$x,y=this_tbDist$y,stringsAsFactors=F))
      crowndf = rbind(crowndf,data.frame(crown=paste0(from,"_core","_",mode),cellType="Ki67",density=this_tbDist[,"Ki67"]/((grid_spacing/1000)^2),x=this_tbDist$x,y=this_tbDist$y,stringsAsFactors=F))
      resetted_meaMat[cbind(this_tbDist$y,this_tbDist$x)] = resetted_meaMat[cbind(this_tbDist$y,this_tbDist$x)]+1
      fileName = paste0(OutDir,"crowns_masks/",s,"_Crown_from",from,"_Core_mask_gridSpacing",grid_spacing,"_",mode,".pdf")
      pdf(fileName, 10, 10, useDingbats = F )
      df = melt(resetted_meaMat)
      # df[(df$value<0) %in% c(T),"value"] = NA
      levelz = levels(factor(df$Var1))
      p=ggplot(df) + theme_bw() + 
        geom_tile(aes(Var2,ordered(Var1, levels=rev(levelz)),fill=value))+
        scale_fill_gradientn(colours=colorRampPalette(c("white","black" ))(n = 100))+
        theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0))+
        labs(title=paste0(""), x="",y="")+
        coord_fixed()
      print(p)
      dev.off()
   }
   return(crowndf)
}

crowns_plots = function( crowndf, fileNameRoot, include_Ki67 = TRUE ){

   fileName = paste0(fileNameRoot,"_boxplot.pdf")
   sorting_crown = data.frame(crown = unique(crowndf$crown), from = NA, inout = "autside", stringsAsFactors=F)
   sorting_crown[substr(sorting_crown$crown,nchar(sorting_crown$crown)-5,nchar(sorting_crown$crown))=="inside","inout"] = "inside"
   sorting_crown_out = sorting_crown[sorting_crown$inout=="autside",]
   sorting_crown_out$from = as.numeric(unlist(strsplit(sorting_crown_out$crown,split="_"))[seq(1,length(unlist(strsplit(sorting_crown_out$crown,split="_"))),by=3)])
   sorting_crown_out = sorting_crown_out[order(-sorting_crown_out$from),]
   
   sorting_crown_in = sorting_crown[sorting_crown$inout=="inside",]
   sorting_crown_in$from = as.numeric(unlist(strsplit(sorting_crown_in$crown,split="_"))[seq(1,length(unlist(strsplit(sorting_crown_in$crown,split="_"))),by=3)])
   sorting_crown_in = sorting_crown_in[order(sorting_crown_in$from),]

   sorting_crown = rbind(sorting_crown_out,sorting_crown_in)
   fill = factor(crowndf$crown,  levels = as.character(sorting_crown$crown))
   levels_fill_inside = levels(fill)[substr(levels(fill),nchar(levels(fill))-5,nchar(levels(fill)))=="inside"]
   levels_fill_outside = levels(fill)[substr(levels(fill),nchar(levels(fill))-6,nchar(levels(fill)))=="outside"]
   fillColors_out = rev(c("mediumorchid2","mediumorchid3","mediumorchid4","purple4")[1:length(levels_fill_outside)])
   fillColors_in = c("firebrick1","firebrick2","firebrick3","firebrick","darkred")[1:length(levels_fill_inside)]
   fillColors = c(fillColors_out,fillColors_in)
   x = crowndf$cellType
   y = crowndf$density
   pdf( fileName, width = 10, height = 7 , useDingbats = F )
   p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(color = "black", outlier.shape = NA) + scale_fill_manual(values = fillColors) + 
               ggtitle( "" ) + xlab( "" ) + ylab( paste0("Density, N/(mm^2)") ) + labs( fill = "Crown" ) + stat_compare_means(label = "p.signif",method = "kruskal") + ylim(0, 4000) +
               theme_bw() + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
   print(p)
   dev.off()

   crowndf_aggr = aggregate(density~crown+cellType, data = crowndf, FUN=mean)

   if (!(include_Ki67))
   {
      crowndf_aggr$cellType = as.character(crowndf_aggr$cellType)
      crowndf_aggr = crowndf_aggr[crowndf_aggr$cellType!="Ki67",]
   }
   
   fileName = paste0(fileNameRoot,"_Means_LinePlot.pdf")
   x = factor(crowndf_aggr$crown, levels = levels(fill))
   y = crowndf_aggr$density
   group = factor(crowndf_aggr$cellType)
   groupColors = as.character(cellTypes_colors_df[levels(group),"color"])
   pdf( fileName, width = 10, height = 7 , useDingbats = F )
   p = ggplot(mapping = aes(y = y, x = x, group = group)) + geom_line(aes(color=group), show.legend = T, size = 2) + geom_point(aes(color=group),size = 10, alpha = 0.5) + ggtitle( "" ) + xlab( "" ) + ylab( "Mean density, N/(mm^2)" ) + labs( color = "Cell type" ) + theme_bw()
   p = p + scale_color_manual( values=groupColors )
   print(p)
   dev.off()

}

###############################################

#################### Main #####################

dir.create(OutDir, showWarnings = F)
dir.create(paste0(OutDir,"distdf_expanded/"), showWarnings = F)
dir.create(paste0(OutDir,"crowns_masks/"), showWarnings = F)
dir.create(paste0(OutDir,"colorized_crowns/"), showWarnings = F)

### Computing distdf for very large areas
cat("\n",s,"\n")
clb = paste0(s,"_SBoundaryVeryLarge.txt")
polyline_file = paste0(InDir_boundaries,clb)
first = TRUE
for (ct in cellTypes)
{
  load(paste0(InDir_matrices,s,"_tileSizeMicrons",grid_spacing,"_",ct,"_meaMat_and_coos.RData"))
  meaMat_uncensored = meaMat
  meaMat_censored = polyLine_censoring(polyline_file, meaMat, all_x_coos, all_y_coos, do_close_boundary = T)
  
  if (first)
  {
     distdf = compute_distances_ClosedBoundary_first_marker(meaMat_censored, meaMat_uncensored, all_x_coos, all_y_coos, ct, polyline_file, fileRoot = paste0(OutDir,substr(clb,1,nchar(clb)-4),"_tileSizeMicrons",grid_spacing) )
     first = FALSE
  } else {
     distdf = compute_distances_remaining_markers(distdf, meaMat_uncensored, ct)
  }
}
save(distdf, file = paste0(OutDir,"distdf_expanded/",s,"_gridSpacing",grid_spacing,"_distdf_solid_VeryLarge.RData"))

### Computing distdf for exact areas
clb = paste0(s,"_closedBoundary_solid_right.txt")
polyline_file = paste0(InDir_boundaries,clb)
first = TRUE
for (ct in cellTypes)
{
  load(paste0(InDir_matrices,s,"_tileSizeMicrons",grid_spacing,"_",ct,"_meaMat_and_coos.RData"))
  meaMat_uncensored = meaMat
  meaMat_censored = polyLine_censoring(polyline_file, meaMat, all_x_coos, all_y_coos, do_close_boundary = T)
  
  if (first)
  {
     distdf = compute_distances_ClosedBoundary_first_marker(meaMat_censored, meaMat_uncensored, all_x_coos, all_y_coos, ct, polyline_file, fileRoot = paste0(OutDir,substr(clb,1,nchar(clb)-4),"_tileSizeMicrons",grid_spacing) )
     first = FALSE
  } else {
     distdf = compute_distances_remaining_markers(distdf, meaMat_uncensored, ct)
  }
}
save(distdf, file = paste0(OutDir,"distdf_expanded/",s,"_gridSpacing",grid_spacing,"_distdf_solid_Exact.RData"))

load(file = paste0(OutDir,"distdf_expanded/",s,"_gridSpacing",grid_spacing,"_distdf_solid_Exact.RData"))
sdf = distdf # counts df restricted to the solid closed region
### removing tiles too close to BoundaryToExclude
clb = paste0(s,"_SBoundaryToExclude_1.txt")
polyline_file = paste0(InDir_boundaries,clb)
if (file.exists(polyline_file))
{
   tbDist = compute_distances_TrueBoundary(sdf, polyline_file, grid_spacing, FALSE)
   sdf = tbDist[tbDist$dist>=max(crown_steps_inside),]
   sdf$dist = NA
}
clb = paste0(s,"_RealSBoundary_1.txt")
polyline_file = paste0(InDir_boundaries,clb)
# using sdf, compute distances wrt true boundary
closed_boundary = FALSE
tbDist = compute_distances_TrueBoundary(sdf, polyline_file, grid_spacing, closed_boundary)
maxDist = max(tbDist$dist)
clb = paste0(s,"_closedBoundary_solid_right.txt")
polyline_file = paste0(InDir_boundaries,clb)
load(paste0(InDir_matrices,s,"_tileSizeMicrons",grid_spacing,"_Bcell_meaMat_and_coos.RData")) # meaMat, all_x_coos, all_y_coos
meaMat_censored = polyLine_censoring(polyline_file, meaMat, all_x_coos, all_y_coos, do_close_boundary = T)
crowndf_inside = crowns_analysis(crown_steps = crown_steps_inside, tbDist, s, meaMat_censored, mode = "inside")

distdf_exact = distdf
rownames(distdf_exact) = paste0(distdf_exact$x,"_",distdf_exact$y)
load(file = paste0(OutDir,"distdf_expanded/",s,"_gridSpacing",grid_spacing,"_distdf_solid_VeryLarge.RData"))
rownames(distdf) = paste0(distdf$x,"_",distdf$y)
sdf = distdf[rownames(distdf)[!(rownames(distdf) %in% rownames(distdf_exact))],] # counts df restricted to the solid closed region
rownames(sdf) = 1:nrow(sdf)
### removing tiles too close to BoundaryToExclude
clb = paste0(s,"_SBoundaryToExclude_1.txt")
polyline_file = paste0(InDir_boundaries,clb)
if (file.exists(polyline_file))
{
   tbDist = compute_distances_TrueBoundary(sdf, polyline_file, grid_spacing, FALSE)
   sdf = tbDist[tbDist$dist>=max(crown_steps_outside),]
   sdf$dist = NA
}
clb = paste0(s,"_RealSBoundary_1.txt")
polyline_file = paste0(InDir_boundaries,clb)
# using sdf, compute distances wrt true boundary
closed_boundary = FALSE
tbDist = compute_distances_TrueBoundary(sdf, polyline_file, grid_spacing, closed_boundary)
maxDist = max(tbDist$dist)
print(maxDist)
clb = paste0(s,"_closedBoundary_solid_right.txt")
polyline_file = paste0(InDir_boundaries,clb)
load(paste0(InDir_matrices,s,"_tileSizeMicrons",grid_spacing,"_Bcell_meaMat_and_coos.RData")) # meaMat, all_x_coos, all_y_coos
meaMat_censored = polyLine_censoring(polyline_file, meaMat, all_x_coos, all_y_coos, do_close_boundary = T)
crowndf_outside = crowns_analysis(crown_steps = crown_steps_outside, tbDist, s, meaMat_censored, mode = "outside")

crowndf = rbind(crowndf_inside,crowndf_outside)

save(crowndf,file = paste0(OutDir,s,"_CrownAnalysis_gridSpacing",grid_spacing,".RData"))

ResDir = OutDir
grid_spacing = 100

load(file = paste0(ResDir,s,"_CrownAnalysis_gridSpacing",grid_spacing,".RData"))
fileNameRoot = paste0(ResDir,s,"_CrownAnalysis_gridSpacing",grid_spacing)
crowns_plots(crowndf, fileNameRoot, include_Ki67 = FALSE)

load(file = paste0(OutDir,s,"_CrownAnalysis_gridSpacing",grid_spacing,".RData"))
crowndf_aggr = aggregate(density~crown+cellType, data = crowndf, FUN=mean)
clb = paste0(s,"_closedBoundary_solid_right.txt")
for (ct in as.character(cellTypes))
{
   polyline_file = paste0(InDir_boundaries,clb)
   load(paste0(InDir_matrices,s,"_tileSizeMicrons",grid_spacing,"_Bcell_meaMat_and_coos.RData")) # meaMat, all_x_coos, all_y_coos
   meaMat_censored = polyLine_censoring(polyline_file, meaMat, all_x_coos, all_y_coos, do_close_boundary = T)
   resetted_meaMat = meaMat_censored
   resetted_meaMat[!is.na(resetted_meaMat)] = 0
   for (cr in unique(crowndf$crown))
   {
      this_crowndf = crowndf[(crowndf$crown==cr) & (crowndf$cellType==ct),]
      resetted_meaMat[cbind(this_crowndf$y,this_crowndf$x)] = crowndf_aggr[(crowndf_aggr$crown==cr) & (crowndf_aggr$cellType==ct),"density"]
   }
   resetted_meaMat = polyLine_censoring(polyline_file, resetted_meaMat, all_x_coos, all_y_coos, do_close_boundary = T)
   rs = which((rowSums(resetted_meaMat)==0) %in% c(T))
   min_rs = rs[c(rs[2:length(rs)]-rs[1:(length(rs)-1)],1)>1]-1
   max_rs = rs[c(1,rs[2:length(rs)]-rs[1:(length(rs)-1)])>1]+1
   rs = which((colSums(resetted_meaMat)==0) %in% c(T))
   min_cs = rs[c(rs[2:length(rs)]-rs[1:(length(rs)-1)],1)>1]-1
   max_cs = rs[c(1,rs[2:length(rs)]-rs[1:(length(rs)-1)])>1]+1
   resetted_meaMat = resetted_meaMat[min_rs:max_rs,min_cs:max_cs]
   fileName = paste0(OutDir,"colorized_crowns/",s,"_",ct,"_colorizedCrowns_gridSpacing",grid_spacing,".pdf")
   pdf(fileName, 4, 4, useDingbats = F )
   dfmelt = melt(resetted_meaMat)
   dfmelt[(dfmelt$value<0) %in% c(T),"value"] = NA
   levelz = levels(factor(dfmelt$Var1))
   p=ggplot(dfmelt) +
    geom_tile(aes(Var2,ordered(Var1, levels=rev(levelz)),fill=value))+
    scale_fill_gradientn(colours=colorRampPalette(c("white",as.character(cellTypes_colors_df[ct,"color"]) ))(n = 100))+
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(title=paste0(ct), x="",y="",fill="Density")+
    coord_fixed()
   print(p)
   dev.off()
}

