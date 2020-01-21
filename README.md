
# CAGEO_2020_6
This repository contains code used to generate all the results in the manuscript "GRANITE AND SHALE CLAY MAPPING IN SHEWA-SHAHBAZGARHI COMPLEX (PAKISTAN): APPLICATION OF DECORRELATION STRETCHING AND FEATURE ORIENTED PRINCIPAL COMPONENT ANALYSIS ON MULTISPECTRAL DATA".

The code is this repository can be compiled using Google Earth Engine Javascript Platform.
 
Google Earth Engine Link of the Code: https://code.earthengine.google.com/ac66af0f58994f62485bee76f42e5527

Description: A readme section is available at the start of the code which describes all the necessary details.
Requirements: The code can be executed online by opening the link (given above) in any browser (Chrome Recommended). No specific system requirements needed. The code is executed on Google Earth Engine cloud computing platform. In our study, this code was developed and executed on a computer system with the following specifications, Intel® Core™ i5-8400 CPU @ (2.8 GHz and 2.81 GHz) and 8 GB RAM.

Contact details: Intelligent Information Processing Lab, National Center of Artificial Intelligence University of Engineering and Technology Peshawar; 

Emails:  khan.m@uetpeshswar.edu.pk; mfawadakbar@uetpeshawar.edu.pk 


1) Sentinel image was filtered and experimented with to obtain the most suitable image which we loaded through image ID. The filtering code is commented
2) All the scatter plots/charts, bands TCCs, FCCs, and Individual PCs outputs are commented to avoid overloading the system. You can uncomment the 'print(...)' and Map.addLayer(...) to get these outputs
3) The code is ordered as follows:
  a) Filtering and Loading Images, Selecting bands
  b) DS and PCA main function definition
  c) Bands TCCs and FCCs (Outputs are commented)
  d) Raw bands scatter plots (Outputs are commented)
  e) Apply DS algorithm on all dataset to stretch bands data
  f) Stretched bands scatter plots (Outputs are commented)
  g) Apply PCA of stretched bands data of Landsat -8, ASTER, Sentinel - 2 (in this order) output results and plot PC1 vs PC2 and PC3 scatter plots for each dataset
  f) Apply PCA on Raw data --> Obtain Scatter Plot of PCs ---> Stretche the PCs data ---> Obtain Scatter Plots of stretched PCs
            Done this for Landsat -8, ASTER and Sentinel - 2 in this order.
  h) Center the map to the specified coordinates (Last line of code)
4) To get invidual PCs copy and paste this code snippet below PCA application:
5) 

for (var i = 0; i < bandNames.length().getInfo(); i++) {
  var band = pcImage.bandNames().get(i).getInfo();
  //Map.addLayer(pcImage.select([band]), {min: -2, max: 2}, band);
}