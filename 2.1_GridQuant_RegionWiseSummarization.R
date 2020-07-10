
# 2.1_GridQuant_RegionWiseSummarization.R
# This is step 2.1 of GridQuant pipeline. In this step, grid values are summarized at the level of regions previously annotated on the image, and whose boundary coordinates were saved.

# --- Usage ---
# The script can be run directly after setting up inputs and other parameters. 

# --- Input ---
# InDir_boundaries: folder containing boundaries saved at previous step
# InDir_matrices: folder containing grid matrices
# OutDir_RegionWiseTables: folder where region wise tables should be saved
# s: sample name
# grid_spacing = grid spacing/tile size in microns
# cellTypes = vector of all cell types whose detections are to be processed

# --- Output ---
# Region-wise counts tables and heatmap masks (to check whether the region considered is as desired) 

# --- Other parameters ---
# Cell type colors

# --- Requires ---
# R, reshape2, ggplot2

# Author: Daniele Tavernari
# Please cite: Tavernari et al., 2021 (see full citation on GitHub repository's README)

library(reshape2)
library(ggplot2)

########## Input
InDir_matrices = "Grid_matrices/"
InDir_boundaries = "All_boundaries/"
OutDir_RegionWiseTables = "RegionWiseTables/"
s = "s8B"
grid_spacing = 100
cellTypes = c("Bcell", "CD4", "CD8", "macrophages", "Ki67")
####################

########## Other parameters
cellTypes_colors_df = data.frame(row.names = cellTypes, cellType = cellTypes, color = c("cyan3", "blue4", "goldenrod1", "red3","green3"), stringsAsFactors = F )
####################

########## Functions

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

compute_distances_ClosedBoundary_first_marker = function(meaMat, all_x_coos, all_y_coos, ct, polyline_file, fileRoot){

   last_row = nrow(meaMat)
   meaMat_is_in_boundary = matrix(data = 0, nrow = nrow(meaMat), ncol = ncol(meaMat))
   distdf = data.frame(matrix(nrow = 0, ncol = 8, dimnames = list(c(),c("x","y","dist","macrophages", "Ki67", "CD4", "CD8", "Bcell"))))
   for (row in 1:last_row)
   {
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
   meaMat_is_in_boundary = correct_erroneous_padding_from_top( meaMat_is_in_boundary )
   
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
   meaMat_is_in_boundary[is.na(meaMat_is_in_boundary)] = 0
   distdf = data.frame(matrix(nrow = sum(meaMat_is_in_boundary), ncol = 8, dimnames = list(c(1:sum(meaMat_is_in_boundary)),c("x","y","dist","macrophages", "Ki67", "CD4", "CD8", "Bcell"))))
   indices = which(meaMat_is_in_boundary == 1, arr.ind = TRUE)
   for (index in 1:nrow(indices))
   {
      distdf[index,"x"] = as.numeric(indices[index,"col"])
      distdf[index,"y"] = as.numeric(indices[index,"row"])
      distdf[index,ct] = meaMat[distdf[index,"y"],distdf[index,"x"]]
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

####################


########## Main
dir.create(OutDir_RegionWiseTables, showWarnings = F)

cat("\n", "Sample:", s, "\n")
AllClosedBoundaries = list.files(path = InDir_boundaries, pattern = paste0(s,'_closedBoundary_*'))
for (clb in AllClosedBoundaries[substr(AllClosedBoundaries,nchar(AllClosedBoundaries)-2,nchar(AllClosedBoundaries))=="txt" ])
{
  polyline_file = paste0(InDir_boundaries,clb)
  first = TRUE
  for (ct in cellTypes)
  {
    load(paste0(InDir_matrices,s,"_tileSizeMicrons",grid_spacing,"_",ct,"_meaMat_and_coos.RData"))
    meaMat_censored = polyLine_censoring(polyline_file, meaMat, all_x_coos, all_y_coos, do_close_boundary = T)
    if (first)
    {
       distdf = compute_distances_ClosedBoundary_first_marker(meaMat_censored, all_x_coos, all_y_coos, ct, polyline_file, fileRoot = paste0(OutDir_RegionWiseTables,substr(clb,1,nchar(clb)-4),"_tileSizeMicrons",grid_spacing) )
       first = FALSE
    } else {
       distdf = compute_distances_remaining_markers(distdf, meaMat_censored, ct)
    }
  }
  save(distdf,file=paste0(OutDir_RegionWiseTables,substr(clb,1,nchar(clb)-4),"_tileSizeMicrons",grid_spacing,"_counts_table.RData"))
}
####################
