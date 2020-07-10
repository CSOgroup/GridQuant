
# 3.0_GridQuant_DensityAcrossPatterns.R
# Densities for all cell types are computed across patterns.

# Author: Daniele Tavernari
# Please cite: Tavernari et al., 2021 (see full citation on GitHub repository's README)

library(reshape2)
library(ggplot2)

########## Input
InDir_RegionWiseTables = "RegionWiseTables/"
OutDir_Densities = "Density_acrossPatterns/"
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

########## Main
dir.create(OutDir_Densities, showWarnings = F)
counts_tables = list.files(path = InDir_RegionWiseTables)
densdf = data.frame(matrix(nrow = 0, ncol = 5, dimnames = list(c(),c("Sample","Pattern","CellType","Counts","Density"))),stringsAsFactors=F)
these_counts_tables = counts_tables[grepl(paste0("tileSizeMicrons",grid_spacing,"_counts"),counts_tables)]
for (tabb in these_counts_tables)
{
   load(file=paste0(InDir_RegionWiseTables,tabb)) # distdf
   tabbs = unlist(strsplit(tabb,split="_"))
   distdf$Sample = tabbs[1]
   distdf$Pattern = tolower(tabbs[3])
   this_densdf = melt(distdf, id.vars = c("Sample","Pattern"), measure.vars = c(cellTypes))
   colnames(this_densdf) = c("Sample","Pattern","CellType","Counts")
   this_densdf$Density = this_densdf$Counts/((grid_spacing/1000)^2)
   densdf = rbind(densdf,this_densdf)
   all_ct = as.character(unique(this_densdf$CellType))
   this_densdf2 = data.frame(Sample = rep(as.character(unique(this_densdf$Sample)),length(all_ct)), 
      Pattern = rep(as.character(unique(this_densdf$Pattern)),length(all_ct)), 
      counts_table = rep(tabb,length(all_ct)), CellType = all_ct, Mean_density = NA )
   for (ct in all_ct)
   {
      this_densdf2[this_densdf2$CellType==ct,"Mean_density"] = mean(this_densdf[this_densdf$CellType==ct,"Density"])
   }
}

x = densdf$CellType
y = densdf[,"Density"]
fill = factor(densdf$Pattern,  levels = c("normal","lepidic", "papillary", "acinar", "solid"))
fillColors = as.character(PatternColorsDf[levels(fill)[levels(fill) %in% unique(fill)],"color"])
fileName = paste0(OutDir_Densities,"AllSamples_Density_AcrossPatterns_tileSizeMicrons",grid_spacing,"_Boxplot_noJitter.pdf")
pdf( fileName, width = 10, height = 7 , useDingbats = F )
p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(color = "black", outlier.shape = NA) + scale_fill_manual(values = fillColors) + 
            ggtitle( "" ) + xlab( "" ) + ylab( paste0("Density, N/(mm^2)") ) + labs( fill = "Pattern" ) + ylim(0, 4000) + 
            theme_bw() + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
print(p)
dev.off()
####################   

