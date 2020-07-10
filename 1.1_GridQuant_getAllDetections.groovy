
// 1.1_GridQuant_getAllDetections.groovy
// This is step 1.1 of GridQuant pipeline. In this step, cell detection (with parameters tuned at the step before) for all channels is performed and results are saved.

// --- Usage ---
// 1) Load the image in QuPath.
// 2) Set up the inputs and possibly output file names (keep in mind however that downstream scripts will read these files in the format: sample_cellType_allDetections.txt )
// 3) Paste the cell detection lines of code tuned at the step before
// 4) Run the script

// --- Input ---
// OutDir: folder where results should be saved
// Sample: sample name

// --- Output ---
// Tables with cell detections and their coordinates

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
def OutDir = '/mnt/ndata/daniele/lung_multiregion/Reports/GridQuant_for_GitHub/Processed/All_detections/'
def Sample = 's8B'
//

clearDetections();
clearAnnotations();

def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()
int imageWidth = server.getWidth()
int imageHeight = server.getHeight()
roi = new RectangleROI(1, 1, imageWidth-2, imageHeight-2)
annotation = new PathAnnotationObject(roi, PathClassFactory.getPathClass("Region"))
imageData.getHierarchy().addPathObject(annotation, false)
selectAnnotations();



// macrophages
selectAnnotations();
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "Alexa 594",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 600.0,  "watershedPostProcess": false,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
saveDetectionMeasurements(buildFilePath(OutDir, Sample + '_macrophages_allDetections.txt'))
// CD8
selectAnnotations();
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "RFP",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 5.0,  "maxAreaMicrons": 400.0,  "threshold": 400.0,  "watershedPostProcess": false,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
saveDetectionMeasurements(buildFilePath(OutDir, Sample + '_CD8_allDetections.txt'))
// CD4
selectAnnotations();
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "CFP",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 600.0,  "watershedPostProcess": false,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
saveDetectionMeasurements(buildFilePath(OutDir, Sample + '_CD4_allDetections.txt'))
// Bcell
selectAnnotations();
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "FITC",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 300.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
saveDetectionMeasurements(buildFilePath(OutDir, Sample + '_Bcell_allDetections.txt'))
// Ki67
selectAnnotations();
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "CY5",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 5.0,  "maxAreaMicrons": 400.0,  "threshold": 200.0,  "watershedPostProcess": false,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
saveDetectionMeasurements(buildFilePath(OutDir, Sample + '_Ki67_allDetections.txt'))


clearDetections();
clearAnnotations();


println "Done"
