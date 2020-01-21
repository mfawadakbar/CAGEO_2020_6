
# CAGEO_2020_6 
This repository contains code used to generate all the results in the manuscript "GRANITE AND SHALE CLAY MAPPING IN SHEWA-SHAHBAZGARHI COMPLEX (PAKISTAN): APPLICATION OF DECORRELATION STRETCHING AND FEATURE ORIENTED PRINCIPAL COMPONENT ANALYSIS ON MULTISPECTRAL DATA".

The code in this repository can be compiled using Google Earth Engine Javascript Platform. You will need to learn how to run a code on Google Earth Engine: https://www.youtube.com/watch?v=2i6cw7nTbhI&t=28s
 
The code will automatically run on Google Earth Engine using this link but it needs will need a Google Account: https://code.earthengine.google.com/e350f0edb0af533acb7f9545171a359c



## Requirements
The code can be executed online by opening the link (given above) in any browser (Chrome Recommended). No specific system requirements needed. The code is executed on Google Earth Engine cloud computing platform. In our study, this code was developed and executed on a system with the following specifications, Intel® Core™ i5-8400 CPU @ (2.8 GHz and 2.81 GHz) and 8 GB RAM.

## Contact details
The authors of the manuscript are reserarchers at Intelligent Information Processing Lab, National Center of Artificial Intelligence University of Engineering and Technology Peshawar, Khyber Pakhtunkhwa, Pakistan;
* Emails:  mfawadakbar@uetpeshawar.edu.pk; khan.m@uetpeshswar.edu.pk; shahabuddin@uetpeshawar.edu.pk

## Code Description
1. Sentinel image was filtered and experimented with to obtain the most suitable image which we loaded through an image ID. The filtering code is commented in the code file

2. All the scatter plots/charts, bands TCCs, FCCs, and Individual PCs outputs are commented to avoid overloading the Google Earth Engine cloud processing platform. You can uncomment the 'print(...)' and Map.addLayer(...) to get these outputs.

3. The code is ordered as follows:
  
   - Filtering and Loading Images and Selecting bands
  
   - Defining DS and PCA main functions
  
   - Bands TCCs and FCCs (Outputs are commented)
  
   - Raw bands scatter plots (Outputs are commented)
  
   - Apply DS algorithm on all dataset to stretch bands data
  
   - Stretched bands scatter plots (Outputs are commented)
  
   - Apply PCA on stretched bands data of Landsat -8, ASTER, Sentinel - 2 (in this order) print results and plot PC1 vs PC2 and PC3 scatter plots for each dataset
  
   - Apply PCA on Raw data --> Obtain Scatter Plot of PCs ---> Stretch the PCs data ---> Obtain Scatter Plots of stretched PCs.
            Repeat this for Landsat - 8, ASTER and Sentinel - 2, respectively.
  
   - Center the map to the specified coordinates (The last line of code)

4. To visualize invidual PCs as individual layers copy and paste this code snippet below PCA application:

`
for (var i = 0; i < bandNames.length().getInfo(); i++) {
  var band = pcImage.bandNames().get(i).getInfo();
  Map.addLayer(pcImage.select([band]), {min: -2, max: 2}, band);
}
`
