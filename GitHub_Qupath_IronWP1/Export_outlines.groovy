setImageType('FLUORESCENCE')
clearDetections()

// Get the current image name (without extension)
def imageName = getProjectEntry().getImageName().replaceFirst(/\.[^.]+$/, "")

// Define output folder and build file path
def outputDir = "\\\\vs03\\VS03-IMM-2\\Iron_project\\Iron_WP1_pathology\\Data_Analysis\\Outline_files\\HLA-Iba1"
def exportPath = outputDir + "\\" + imageName + ".geojson"

// Export outlines with dynamic file name
exportAllObjectsToGeoJson(exportPath, "EXCLUDE_MEASUREMENTS")