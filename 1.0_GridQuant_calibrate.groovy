
// 1.0_GridQuant_calibrate.groovy
// This is step 1.0 of GridQuant pipeline. It consists in calibrating the parameters for automated cell detection on a small representative region of interest (ROI), in order to allow visual inspection before performing cell detection on the entire image.

// --- Usage ---
// 1) Load the image in QuPath.
// 2) Draw a small representative ROI
// 3) Uncomment the line below corresponding to the channel of interest. You might need to change the channel name.
// 4) Tune the parameters for cell detection
// 5) Run the script
// 6) Check visually that the automated cell detection is reasonable. You might find it useful to switch to grayscale.
// 7) Repeat steps 4-6 until you are satisfied
// 8) Repeat steps 3-7 for all channels
// 9) Repeat steps 2-8 on a few more representative ROIs for validation
// 10) Copy the script lines with your final parameters and paste them in the right position in the next script of GridQuant, i.e. '1.1_GridQuant_getAllDetections.groovy'

// --- Input ---
// Rectangular ROIs drawn by the user on the image

// --- Output ---
// Cells detected for a given channel

// --- Other parameters ---
// Cell detection parameters

// --- Requires ---
// QuPath (tested on: v0.2.0-m4, on a Unix system)

// Author: Daniele Tavernari
// Please cite: Tavernari et al., 2021 (see full citation on GitHub repository's README)

println "Script starting..."

import qupath.lib.objects.PathAnnotationObject
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.roi.RectangleROI
import qupath.lib.gui.scripting.QPEx

// Input
//

clearDetections();
selectAnnotations();

def imageData = QPEx.getCurrentImageData()

// Uncomment the appropriate line

// macrophages
// runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "Alexa 594",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 600.0,  "watershedPostProcess": false,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
// CD8
// runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "RFP",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 5.0,  "maxAreaMicrons": 400.0,  "threshold": 400.0,  "watershedPostProcess": false,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
// CD4
// runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "CFP",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 600.0,  "watershedPostProcess": false,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
// Bcell
// runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "FITC",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 300.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
// Ki67
// runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "CY5",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 5.0,  "maxAreaMicrons": 400.0,  "threshold": 200.0,  "watershedPostProcess": false,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

println "Done"
