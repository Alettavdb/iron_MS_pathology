// ==========================
//  ImageJ Macro: Lysosomal Load Measurement
//  Works on currently opened image
//  Aletta, 14-08-2025
// ==========================

imp = getTitle();  // get base image title
lastDash = lastIndexOf(imp, " - ");
baseName = substring(imp, 0, lastDash);

// create max intensity projections
channels = newArray("DAPI", "LAMP1", "GAL3", "iba1");
for (c = 0; c < channels.length; c++) {
    selectWindow(baseName + " - " + channels[c] + ".Confocal");
    run("Z Project...", "projection=[Max Intensity]");
}

// Close original z-stacks
for (c = 0; c < channels.length; c++) {
    selectWindow(baseName + " - " + channels[c] + ".Confocal");
    close();
}

// set colors
selectWindow("MAX_" + baseName + " - iba1.Confocal"); run("Green");
selectWindow("MAX_" + baseName + " - GAL3.Confocal"); run("Yellow");
selectWindow("MAX_" + baseName + " - LAMP1.Confocal"); run("Magenta");
selectWindow("MAX_" + baseName + " - DAPI.Confocal"); run("Blue");

// enhance contrast for all
for (c = 0; c < channels.length; c++) {
    selectWindow("MAX_" + baseName + " - " + channels[c] + ".Confocal");
    run("Enhance Contrast...", "saturated=0.35");
}

// Create composites and save as .tiff and .jpg
// This allows to check results and see progress

selectWindow("MAX_" + baseName + " - iba1.Confocal");
run("Duplicate...", " ");
selectWindow("MAX_" + baseName + " - LAMP1.Confocal");
run("Duplicate...", " ");
selectWindow("MAX_" + baseName + " - GAL3.Confocal");
run("Duplicate...", " ");
selectWindow("MAX_" + baseName + " - DAPI.Confocal");
run("Duplicate...", " ");


cmd = "c2=[MAX_" + baseName + " - iba1.Confocal-1] " +
      "c3=[MAX_" + baseName + " - DAPI.Confocal-1] " +
      "c6=[MAX_" + baseName + " - LAMP1.Confocal-1] " +
      "c7=[MAX_" + baseName + " - GAL3.Confocal-1] create";
run("Merge Channels...", cmd);

saveAs("Tiff", "//vs03/VS03-IMM-2/Iron_project/Iron_WP1_pathology/Data_Analysis/LAMP1 + GAL3 + iba1/Composites/" + baseName + ".tif");
saveAs("Jpeg", "//vs03/VS03-IMM-2/Iron_project/Iron_WP1_pathology/Data_Analysis/LAMP1 + GAL3 + iba1/Composites/" + baseName + ".jpg");

close;


// merge DAPI and iba1 to be able to properly perform the outlines
run("Merge Channels...", "c2=[MAX_" + baseName + " - iba1.Confocal] c3=[MAX_" + baseName + " - DAPI.Confocal] create");

// Run ROI manager
run("ROI Manager...");
//setTool("freehand");
waitForUser("Outline the cells, add them to ROI Manager\n\nClick OK to continue.");

roiManager("Save", "//vs03/VS03-IMM-2/Iron_project/Iron_WP1_pathology/Data_Analysis/LAMP1 + GAL3 + iba1/ROIs/" + baseName + "_RoiSet.zip");


// Create a mask of LAMP1
selectWindow("MAX_" + baseName + " - LAMP1.Confocal");
setAutoThreshold("Li dark no-reset");
run("Convert to Mask");

// Create a mask of GAL3
selectWindow("MAX_" + baseName + " - GAL3.Confocal");
setAutoThreshold("MaxEntropy dark no-reset");
run("Convert to Mask");

// Create the overlapping mask of LAMP1 and GAL3
imageCalculator("AND create", "MAX_" + baseName + " - LAMP1.Confocal", "MAX_" + baseName + " - GAL3.Confocal");

// Final check for measurements
waitForUser("If all is good, we will measure the ROIs now\n\nClick OK to continue.");



run("Set Measurements...", "area mean limit display redirect=None decimal=2");

n = roiManager("count");       // get total number of ROIs
run("Clear Results");          // clear any previous results

for (i = 0; i < n; i++) {
    roiManager("Select", i);   // select current ROI

    // 1. Total ROI area on Composite (reference)
    selectImage("Composite");
    roiManager("Measure");

    // 2. LAMP1 positive area
    selectWindow("MAX_" + baseName + " - LAMP1.Confocal");
    setAutoThreshold("Li dark no-reset");
    roiManager("Measure");

    // 3. GAL3 & LAMP1 overlap
    selectImage("Result of MAX_" + baseName + " - LAMP1.Confocal");
    setAutoThreshold("MaxEntropy dark no-reset");
    roiManager("Measure");
}

saveAs("Results", "//vs03/VS03-IMM-2/Iron_project/Iron_WP1_pathology/Data_Analysis/LAMP1 + GAL3 + iba1/Measurements/Results" + baseName +".csv");

waitForUser("Done :) MAKE SURE TO COPY THE RESULTS!\n\nClick OK to continue.");

// Close all open image windows
while (nImages > 0) {
    selectImage(nImages);
    close();
}

// Empty the ROI Manager
if (roiManager("count") > 0) {
    roiManager("Deselect");
    roiManager("Delete");
}

// Clear the Results table
if (isOpen("Results")) {
    selectWindow("Results");
    run("Clear Results");
}