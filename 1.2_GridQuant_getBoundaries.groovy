
// 1.2_GridQuant_getBoundaries.groovy
// This is step 1.2 of GridQuant pipeline. In this step, coordinates of boundaries (open or closed) drawn on the image are saved.

// --- Usage ---
// 1) Load the image in QuPath.
// 2) Draw your boundaries as polygons or polylines
// 3) Set up the inputs
// 4) Select the polygon or polyline representing your boundary of interest
// 5) Run the script
// 6) Repeat points 3-5 for all of your boundaries

// --- Input ---
// OutDir: folder where results should be saved
// OutFile: name of the output file

// --- Output ---
// Text file with coordinates of boundary and pixel size

// --- Other parameters ---
// 

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
def OutDir = 'All_boundaries/'
def OutFile = "s8B_closedBoundary_papillary_top.txt"
//

def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()
double pixelDimMicrons = server.getPixelCalibration().getAveragedPixelSizeMicrons()
def polyPoints = getSelectedROI().getPolygonPoints()
print polyPoints
String result = polyPoints.join(";")

File file = new File(OutDir+OutFile)
file.write result+"\n"
file << "PixelDim="+pixelDimMicrons.toString()+"\n"

println file.text

println "Done"
