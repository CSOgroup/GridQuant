# GridQuant
Gridding framework for quantification and statistical analyses of spatially-resolved protein expression from high resolution images

**Pipeline**

GridQuant involves the following steps, ranging from automated cell detection in QuPath to summarization of cell counts into discrete matrices of tunable grid size. 
The .groovy scripts are to be run from QuPath script editor with the image opened.
Detailed explanations of scripts usage are available in their initial lines.

_1.0_GridQuant_calibrate.groovy_

Calibration of parameters for automated cell detection.

_1.1_GridQuant_getAllDetections.groovy_

Cell detection and saving of tables with cell coordinates.

_1.2_GridQuant_getBoundaries.groovy_

Saving of the coordinates of annotations drawn on the image.

_2.0_GridQuant_CreateGrid.R_

Summarization of cell counts for each cell type into matrices of tunable grid size and plotting of counts heatmaps.

_2.1_GridQuant_RegionWiseSummarization.R, 3.0_GridQuant_DensityAcrossPatterns.R, 3.1_GridQuant_ColocalizationsAcrossPatterns.R, 3.2_GridQuant_SolidBoundaryAnalyses.R_

Various downstream analyses.

**Citation**

[Tavernari et al., Cancer Discovery 2021](https://cancerdiscovery.aacrjournals.org/content/11/6/1490)

**Contacts**

daniele.tavernari@unil.ch
