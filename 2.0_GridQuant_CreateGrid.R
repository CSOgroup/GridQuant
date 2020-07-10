
# 2.0_GridQuant_CreateGrid.R
# This is step 2.0 of GridQuant pipeline. In this step, cell counts will be summarized as grid matrices with tunable tile size.

# --- Usage ---
# The script can be run directly after setting up inputs and other parameters. Importantly, you should change cell types according to the ones you assayed.

# --- Input ---
# InDir_detections: folder containing cell detection tables saved at previous step
# InDir_boundaries: folder containing boundaries saved at previous step
# OutDir_matrices: folder where grid matrices should be saved
# OutDir_heatmaps: folder where heatmaps should be saved
# s: sample name
# grid_spacing = grid spacing/tile size in microns
# cellTypes = vector of all cell types whose detections are to be processed


# --- Output ---
# Grid matrices with cell counts and heatmaps for all cell types

# --- Other parameters ---
# Cell type colors

# --- Requires ---
# R, reshape2, ggplot2

# Author: Daniele Tavernari
# Please cite: Tavernari et al., 2021 (see full citation on GitHub repository's README)

library(reshape2)
library(ggplot2)

########## Input
InDir_detections = "All_detections/"
InDir_boundaries = "All_boundaries/"
OutDir_matrices = "Grid_matrices/"
OutDir_heatmaps = "Grid_heatmaps/"
s = "s8B"
grid_spacing = 100
cellTypes = c("Bcell", "CD4", "CD8", "macrophages", "Ki67")
####################

########## Other parameters
cellTypes_colors_df = data.frame(row.names = cellTypes, cellType = cellTypes, color = c("cyan3", "blue4", "goldenrod1", "red3","green3"), stringsAsFactors = F )
####################

########## Functions

construct_tunable_df = function( all_detects, grid_spacing, maxX, maxY ){

  tadp = all_detects[,c("Class","Centroid.X.µm","Centroid.Y.µm")]
  xmax = ceiling(maxX/grid_spacing)
  ymax = ceiling(maxY/grid_spacing)
  tadp$x = factor(ceiling(tadp$Centroid.X.µm/grid_spacing), levels = c(1:xmax))
  tadp$y = factor(ceiling(tadp$Centroid.Y.µm/grid_spacing), levels = c(1:ymax))
  aa = t(table(tadp$x, tadp$y))
  return( aa )

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

####################

########## Main
dir.create(OutDir_heatmaps, showWarnings = F)
dir.create(OutDir_matrices, showWarnings = F)

cat("\n", "Sample:", s, "\n")
maxX = 0
maxY = 0
for (ct in as.character(cellTypes))
{
   df = read.table(file = paste0(InDir_detections,s,"_",ct,"_allDetections.txt") , header = T, stringsAsFactors = F, sep = "\t", quote = '')
   maxX = max(maxX,ceiling(max(df$Centroid.X.µm)))
   maxY = max(maxY,ceiling(max(df$Centroid.Y.µm)))
}

for (ct in as.character(cellTypes))
{
    cat("\n", "... cell type:", ct, "\n")
    all_detects = read.table(file = paste0(InDir_detections,s,"_",ct,"_allDetections.txt"), header = T, stringsAsFactors = F, sep = "\t", quote = '')
    meaMat = construct_tunable_df( all_detects, grid_spacing, maxX, maxY )
    all_x_coos = seq(grid_spacing/2, ncol(meaMat)*grid_spacing, by = grid_spacing )
    all_y_coos = seq(grid_spacing/2, nrow(meaMat)*grid_spacing, by = grid_spacing )
    if ( (length(all_y_coos)!=nrow(meaMat)) | (length(all_x_coos)!=ncol(meaMat)) ) { next }
    save(meaMat, all_x_coos, all_y_coos, file = paste0(OutDir_matrices,s,"_tileSizeMicrons",grid_spacing,"_",ct,"_meaMat_and_coos.RData"))
    
    fileName = paste0(OutDir_heatmaps,s,"_tileSizeMicrons",grid_spacing,"_",ct,"_Heatmap.pdf")
    pdf(fileName, 9, 6, useDingbats = F )
    dfmelt = melt(meaMat)
    dfmelt[(dfmelt$value<0) %in% c(T),"value"] = NA
    levelz = levels(factor(dfmelt$Var1))
    colnames(dfmelt) = c("Var1","Var2","Counts")
    p=ggplot(dfmelt) +
      geom_tile(aes(Var2,ordered(Var1, levels=rev(levelz)),fill=Counts))+
      scale_fill_gradientn(colours=colorRampPalette(c("white",as.character(cellTypes_colors_df[ct,"color"]) ))(n = 100))+
      theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
      labs(title=paste0(ct," cells, grid spacing = ",grid_spacing," microns"), x="",y="")+
      coord_fixed()
    print(p)
    dev.off()

    AllClosedBoundaries = list.files(path = paste0(InDir_boundaries), pattern = paste0(s,'_closedBoundary_*') )
    meaMat_censored = meaMat
    for (clb in AllClosedBoundaries[substr(AllClosedBoundaries,nchar(AllClosedBoundaries)-2,nchar(AllClosedBoundaries))=="txt" ])
    {
       polyline_file = paste0(InDir_boundaries,clb)
       meaMat_censored = polyLine_censoring(polyline_file, meaMat_censored, all_x_coos, all_y_coos, do_close_boundary = T)
    }
    
    fileName = paste0(OutDir_heatmaps,s,"_tileSizeMicrons",grid_spacing,"_",ct,"_Heatmap_withAnnotations.pdf")
    pdf(fileName, 9, 6, useDingbats = F )
    dfmelt = melt(meaMat_censored)
    dfmelt[(dfmelt$value<0) %in% c(T),"value"] = NA
    levelz = levels(factor(dfmelt$Var1))
    colnames(dfmelt) = c("Var1","Var2","Counts")
    p=ggplot(dfmelt) +
      geom_tile(aes(Var2,ordered(Var1, levels=rev(levelz)),fill=Counts))+
      scale_fill_gradientn(colours=colorRampPalette(c("white",as.character(cellTypes_colors_df[ct,"color"]) ))(n = 100))+
      theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.border=element_rect(fill = NA, colour='black',size=0.7)) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
      labs(title=paste0(ct," cells, grid spacing = ",grid_spacing," microns"), x="",y="")+
      coord_fixed()
    print(p)
    dev.off()
  
}

####################


