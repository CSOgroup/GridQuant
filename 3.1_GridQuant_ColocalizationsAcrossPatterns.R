
# 3.1_GridQuant_ColocalizationsAcrossPatterns.R
# Colocalizations computed as Spearman's correlation between Ki-67 and immune cells

# Author: Daniele Tavernari
# Please cite: Tavernari et al., 2021 (see full citation on GitHub repository's README)

library(reshape2)
library(ggplot2)

########## Input
InDir_RegionWiseTables = "RegionWiseTables/"
InDir_boundaries = "All_boundaries/"
OutDir = "Colocalization_acrossPatterns/"
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
####################

########## Functions
dan.boxplots = function( fileName, x, y, fill = NULL, xlab = "default", ylab = "default", filllab = "default", plotTitle = "", signifTest = "kruskal", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = "default", jitterColors = "black", labelJitteredPoints = NULL, fileWidth = 7, fileHeight = 7, hlines_coo = NULL, hlines_labels = NULL )
{
  
  if (xlab=="default") { xlab = deparse(substitute(x)) }
  if (ylab=="default") { ylab = deparse(substitute(y)) }
  if (filllab=="default") { filllab = deparse(substitute(fill)) }

  pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
  if (!is.null(fill))
  {
    if (!(fillColors=="default") )
    {
      p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(color = xColors, outlier.shape = NA,fatten=3) + scale_fill_manual(values = fillColors) + geom_point(pch = 16, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
        ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
        theme_bw() + theme(text = element_text(size=14),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
    } else {
      p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(color = xColors, outlier.shape = NA,fatten=3) + geom_point(pch = 16, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
        ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
        theme_bw() + theme(text = element_text(size=14),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
    }   
  } else {
    # cat(xlab,ylab) label = "p.format",
    p = ggplot(mapping = aes(y = y, x = x)) + geom_boxplot(color = xColors, outlier.shape = NA,fatten=3) +
      ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) +
      theme_bw() + theme(text = element_text(size=14),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    if (!is.null(labelJitteredPoints))
    {
      p = p + geom_text(aes(label=labelJitteredPoints),hjust=0, vjust=0, size = 3, fontface = "bold", position = position_jitter(width = 0.2, height = 0), color = jitterColors)
    } else {
      p = p + geom_jitter(width=.2, size = 4.5, alpha = 0.5, color = jitterColors)
    }
  }
  if (!(is.null(signifTest))) { p = p + stat_compare_means(comparisons = comparisons, method = signifTest, label.y = labelycoo)}
  if (!is.null(hlines_coo))
  {
    p = p + geom_hline( yintercept = hlines_coo, linetype="dashed", color = "gray44" ) + geom_text( aes(0.5, hlines_coo, label = hlines_labels, vjust = -1, hjust = -1), color = "gray44")
  }
  if (!is.null(ylimLeft))
  {
    p = p + ylim(ylimLeft, ylimRight)
  }
  print(p)
  dev.off()
  return()
}
####################


########## Main
dir.create(OutDir, showWarnings = F)

cat("\n", "Sample:", s, "\n")
colocdf_Ki67 = data.frame(matrix(nrow = 0, ncol = 9, dimnames = list(c(),c("Sample","Pattern","grid_spacing","closedBoundaryFile","macrophages", "Ki67", "CD4", "CD8", "Bcell"))))
countsdf_tot = data.frame(matrix(nrow = 0, ncol = 6, dimnames = list(NULL,c("macrophages", "Ki67", "CD4", "CD8", "Bcell","sumPN"))))
AllClosedBoundaries = list.files(path = InDir_boundaries, pattern = paste0(s,'_closedBoundary_*'))
for (clb in AllClosedBoundaries[substr(AllClosedBoundaries,nchar(AllClosedBoundaries)-2,nchar(AllClosedBoundaries))=="txt" ])
{
  load(file=paste0(InDir_RegionWiseTables,substr(clb,1,nchar(clb)-4),"_tileSizeMicrons",grid_spacing,"_counts_table.RData"))
  distdf$sumPN = distdf$macrophages+distdf$Ki67+distdf$CD4+distdf$CD8+distdf$Bcell
  distdf = distdf[,c("macrophages", "Ki67", "CD4", "CD8", "Bcell","sumPN")]
  countsdf_tot = rbind(countsdf_tot,distdf)

}
countsdf_tot$sumPN_density = countsdf_tot$sumPN/((grid_spacing/1000)^2)
### cutoff: max(c(200, 5%))
quantilez = quantile(countsdf_tot$sumPN_density,seq(0,1,0.05))
cutoff = max(c(200, as.numeric(quantilez["5%"])))
AllClosedBoundaries = list.files(path = InDir_boundaries, pattern = paste0(s,'_closedBoundary_*'))
for (clb in AllClosedBoundaries[substr(AllClosedBoundaries,nchar(AllClosedBoundaries)-2,nchar(AllClosedBoundaries))=="txt" ])
{
  load(file=paste0(InDir_RegionWiseTables,substr(clb,1,nchar(clb)-4),"_tileSizeMicrons",grid_spacing,"_counts_table.RData"))
  distdf = distdf[(distdf$macrophages+distdf$Ki67+distdf$CD4+distdf$CD8+distdf$Bcell)/((grid_spacing/1000)^2)>cutoff,]
  # Co-localization heatmap
  method = "spearman"
  cordf = matrix(nrow = length(cellTypes), ncol = length(cellTypes), dimnames = list(cellTypes,cellTypes))
  for (ct1 in cellTypes)
  {
     for (ct2 in cellTypes)
     {
        cordf[ct1,ct2] = signif(cor(distdf[,ct1],distdf[,ct2],method="spearman"),3)
     }
  }
  

  this_pattern = tolower(substr(clb,nchar(s)+17,nchar(clb)-4))
  this_pattern = sub("\\_.*", "", this_pattern)
  colocdf_Ki67 = rbind(colocdf_Ki67,data.frame(Sample=s,Pattern=this_pattern,grid_spacing=grid_spacing,closedBoundaryFile=clb,
   macrophages=cordf["Ki67","macrophages"], Ki67=cordf["Ki67","Ki67"], CD4=cordf["Ki67","CD4"], CD8=cordf["Ki67","CD8"], Bcell=cordf["Ki67","Bcell"]))
}

save(colocdf_Ki67, file = paste0(OutDir,"colocdf_Ki67.RData"))

load(file = paste0(OutDir,"colocdf_Ki67.RData"))
colocdf_Ki67$Pattern = factor(colocdf_Ki67$Pattern, levels = c("normal","lepidic", "papillary", "acinar", "solid"))

for (immune_cell in c("Bcell","CD4","CD8","macrophages"))
{
    mea = colocdf_Ki67[colocdf_Ki67$grid_spacing==grid_spacing,]
    fileName = paste0(OutDir,"/Ki67",immune_cell,"_Coloc_AcrossPatterns","_tileSizeMicrons",grid_spacing,"_Boxplot.pdf")
    y = mea[,immune_cell]
    x = factor(mea$Pattern,  levels = c("normal","lepidic", "papillary", "acinar", "solid"))
    xColors = as.character(PatternColorsDf[levels(x)[levels(x) %in% unique(x)],"color"])
    jitterColors = NULL#mea$color
    dan.boxplots( fileName, x, y, xlab = "", ylab = paste0("Colocalization (Spearman's r)"), signifTest = NULL, comparisons = NULL, labelycoo = NULL, xColors = xColors, jitterColors = "black", labelJitteredPoints = mea$Sample, fileWidth = 10, fileHeight = 7, hlines_coo = NULL, hlines_labels = NULL )
    save(mea, file = paste0(substr(fileName, 1,nchar(fileName)-4),".RData") )
}
####################