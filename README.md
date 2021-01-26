# seahorse-data-analysis-PIXI
pixi (‘R-integrated pixel intensity (PIXI) analysis’) 

Packages needed EBImage, stringr, tidyverse
please install other packages if needed as well.

### Tips:
* input filename of the TIF files is important. The input filename is now set to be for a default Cytation Gen5 file naming convention
That means it should be similar to "Bright Field_D7_1_001.tif"
More specifically: the second string between underscores is assumed to be the wellname (this is default for Cytation Gen5)
Even more specifically: the function "str_split" from the stringr package splits up the filename into strings with the
underscore as separating character. It is now set that the second string from the left is the wellname (tempString <- str_split(fileNameToconvert, fixed("_"), simplify = TRUE)[,2])

* input format of the TIF file is important, we noticed that some TIF files from cytation output have one channel (or frame) and others have 4 channels (or frames)
this impacts the calculation and writing of the newly generated tif files.

* The script outputs both the total intensity of all non background pixels ("total_intensity") and the
number of non-black pixels ("non_blackPixels"), that can both be used for normalizing Seahorse data.
The "total_intensity" was used in the Janssen et al 2021 Scientific reports paper.

### Publication
This script is associated with:
"Novel standardized method for extracellular flux analysis of oxidative and glycolytic metabolism in peripheral blood mononuclear cells"
Janssen et al. 2021 Scientific Reports 11:1662
https://rdcu.be/cebd6

