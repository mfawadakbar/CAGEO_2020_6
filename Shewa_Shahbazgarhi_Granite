var aster = ee.ImageCollection("ASTER/AST_L1T_003"),
    sentinel = ee.ImageCollection("COPERNICUS/S2_SR"),
    landsat = ee.ImageCollection("LANDSAT/LC08/C01/T1_RT"),
    sent = ee.ImageCollection("COPERNICUS/S2_SR"),
    regoin = /* color: #d63000 */ee.Geometry.Point([72.27225926886842, 34.25461952433233]),
    sentroi = 
    /* color: #98ff00 */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[72.17590684462516, 34.562927543817764],
          [72.17590684462516, 34.32254586241205],
          [72.61707658339469, 34.32254586241205],
          [72.61707658339469, 34.562927543817764]]], null, false),
    DSregion = 
    /* color: #bf04c2 */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[72.16930936026165, 34.27247480092394],
          [72.16930936026165, 34.19385006923403],
          [72.35264371084759, 34.19385006923403],
          [72.35264371084759, 34.27247480092394]]], null, false),
    roi = /* color: #d63000 */ee.Geometry.MultiPoint(
        [[72.25484224373645, 34.26258130906891],
         [72.30668397957629, 34.25264966064848]]),
    sentroipoint = /* color: #d63000 */ee.Geometry.MultiPoint(
        [[72.43097463815593, 34.407355551433646],
         [72.5834099408903, 34.391775057471115],
         [72.4948326703825, 34.384408723298876],
         [72.44470754831218, 34.334243764850996]]);

//Filter and Load Sentinel 2 L2C Image
/*var sentImage1 = ee.Image(sent
.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10))
.filterDate("2018-01-02", "2020-12-30")
.filterBounds(sentroipoint)
.sort('CLOUDY_PIXEL_PERCENTAGE')
.first()
);*/

// Mosaic Image for validation of results in Ambela Complex
/*var sentImage2 = ee.Image(sent
.filterDate("2018-01-01", "2019-05-02")
.sort('CLOUD_COVER')
.filterBounds(sentroi)
.first());
var sentImage = ee.ImageCollection([sentImage1, sentImage2]).mosaic();
print("Sentinel 2 Scene", sentImage);
*/

//Load a selected sentinel image, we have filetered this image from the sentinel database
var sentImage = ee.Image('COPERNICUS/S2_SR/20190310T055641_20190310T060759_T43SBT');
//Select Sentinel - 2 Bands to use in the study
var sentbands = ['B2','B3','B4','B8','B11','B12'];

//Filter and load Landsat 8 Image
var landImage = ee.Image(landsat
.filterDate("2018-01-01", "2019-05-02")
.sort('CLOUD_COVER')
.filterBounds(roi)
.first());
print("Landsat 8 Scene", landImage);
//Select Landsat - 8 Bands to use in the study
var landbands = ['B5','B6','B7','B10','B11'];

//Filter and load ASTER Image
var asterImage = ee.Image(aster
.filterDate("2003-09-27", "2019-04-20")
.sort('CLOUDCOVER')
.filterBounds(roi)
.first());
print("ASTER Scene", asterImage);
//Select ASTER Bands to use in the study
var asterbands = ['B3N','B04', 'B05', 'B06', 'B07', 'B08', 'B09'];

//Ignore this, this is only used to avoid errors while concatinting eigen pairs
var eigenCollection = ee.Array([[0],[1],[2],[3],[4]]);


// Decorrelation Stretching Main Function     
function decStr(bandsImage, location, scale){
  var bandNames = bandsImage.bandNames();
  // Naming the axis for intuition
  var dataAxis = 0;
  var bandsAxis = 1;
  // Calculate the mean for each band image
  var meansAll = bandsImage.reduceRegion(ee.Reducer.mean(), location, scale);
  // Generate an array (1D Matrix) of mean of each band
  var arrayOfmeans = ee.Image(meansAll.toArray());
  // Collapse the bands data such that each pixel is a matrix of pixel values of each band
  var pixelArrays = bandsImage.toArray();
  // Use the means array and the collapsed band data to center each pixel of each band by subtracting its corresponding mean from it
  var meanCent = pixelArrays.subtract(arrayOfmeans);
  // Calculate the Covariance matrix for the bands data
  var covar = meanCent.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: location,
    scale: scale
  });
  
  // Get the covariance in array format which shows the band-band covarince of the data
  var covarArray = ee.Array(covar.get('array'));
  // Perform eigen decomposition of the covariance matrix to obtain eigen values and eigen vector pairs
  var eigenPairs = covarArray.eigen();
  var eigenValues = eigenPairs.slice(bandsAxis, 0, 1); // slice(axis, start, end, step)
  var eigenVectors = eigenPairs.slice(bandsAxis, 1);
  // Rotate by the eigenvectors, scale to a variance of 30, and rotate back.
  //Store a diagonal matrix in i
  var i = ee.Array.identity(bandNames.length()); // i will be used to isolate each band data and scale its variance e.g i = [1,0,0,0,0] = isolate first band from 5 bands
  // Calculate variance from the eigenvalues ---> variance = 1/sqrt(eigenvalues)
  // matrixToDiag = Computes a square diagonal matrix from a single column matrix for multiplication purposes
  var variance = eigenValues.sqrt().matrixToDiag();
  //Multiply diagonal matrix i by 30 and divide by vaiance to obtain scaling variance matrix
  var scaled = i.multiply(30).divide(variance); //Changed from 30 -> 50, It was observed that changing variance scale increases contrast. Best contrast obtained for 30
  // Calculate a rotation matrix ---> rotationMatrix =  Eigenvect.Transpose * ScaledVariance * Eigenvect
  var rotation = eigenVectors.transpose()
    .matrixMultiply(scaled)
    .matrixMultiply(eigenVectors);
  // Convert 1-D nomalized array image data to 2-D and transpose it so it can be multiplied with rotation matrix
  var transposed = meanCent.arrayRepeat(bandsAxis, 1).arrayTranspose();
  // Multiply the transposed data with the rotation matrix
  return transposed.matrixMultiply(ee.Image(rotation))
    .arrayProject([bandsAxis])   //This drop unecessary axis from the transposed data and only retains 2 axis
    .arrayFlatten([bandNames])  //Flatten collections of collections
    .add(127).byte(); // Conver pixel values to 127 means so it can be visualized between 0 - 255 range.
    
    // .byte is used to force element wise operation
}
// Principal Component Analysis Main Function
function PCA(meanCent, scale, location){
  // Flatten the band image data in from 2D to a 1D array
  var arrays = meanCent.toArray();
  print('PCA applying on', meanCent);
  // Calculate the covariance matrix for the bands data of the region
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: location,
    scale: scale,
    maxPixels: 1e9
  });
  // Get the band to band covariance of the region in 'array' format. Here .get('array') --> casts to an array
  var covarArray = ee.Array(covar.get('array'));
  // Perform an eigen analysis and slice apart the values and vectors.
  var eigenPairs = covarArray.eigen();
  // This is a P-length vector of Eigenvalues. Here P = number of PCs
  var eigenValues = eigenPairs.slice(1, 0, 1);
  // This is a PxP matrix with eigenvectors in rows.
  var eigenVectors = eigenPairs.slice(1, 1);
  //Print and store eigen pairs in eigenCollection variable and export to drive
  print('eigen Values', eigenValues);
  print('eigen Vector', eigenVectors);
    //Make feature collection out of eigenpairs so it can be exported to excel. From there we Convert it to a table using a python script
  eigenCollection = ee.Feature(null,{values:ee.Array.cat([eigenValues,eigenVectors],1)}); 
  print('Eigen Collection Length',eigenCollection);
    // Export the FeatureCollection to excel sheet in drive
  Export.table.toDrive({
  collection: ee.FeatureCollection([eigenCollection]),
  description: 'eigenAnalysis',
  fileFormat: 'CSV'
  });
  // Convert the 1D image array back to 2D matrix for multiplication
  var imageMat = arrays.toArray(1);
  // To obtain PC = EigenVectors * 2D Image Matrix
  var PCs = ee.Image(eigenVectors).matrixMultiply(imageMat);
  // Turn the square roots of the Eigenvalues into a P-band image.
  var sdImage = ee.Image(eigenValues.sqrt())
    .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);
  // Turn the PCs into a P-band image, normalized by SD.
  return PCs
    // Throw out an an unneeded dimension, [[]] -> [].
    .arrayProject([0])
    // Make the one band array image a multi-band image, [] -> image.
    .arrayFlatten([getNewBandNames('pc')])
    // Normalize the PCs by their SDs.
    .divide(sdImage);
}

          //TCCs and FCCs
  //ASTER L1T FCC
var trueColor = {
  bands: ["B05", "B04", "B3N"],
  min: 0,
  max: 300
};
//Map.addLayer(asterImage, trueColor, "ASTER False-color");

  //Sentinel - 2 L2A TCC
var trueColor = {
  bands: ["B4", "B3", "B2"],
  min: 0,
  max: 3000
};
//Map.addLayer(sentImage, trueColor, "Sentinel 2 True Color");

  //Sentinel - 2 L2A FCC
var trueColor = {
  bands: ["B11", "B12", "B8"],
  min: 0,
  max: 3000
};
//Map.addLayer(sentImage, trueColor, "Sentinel 2 False Color");

  //Landsat - 8 raw TCC
var trueColor = {
  bands: ["B4", "B3", "B2"],
  min: 0,
  max: 30000
};
//Map.addLayer(landImage, trueColor, "Landsat True Color");

  //Landsat - 8 raw FCC
var trueColor = {
  bands: ["B6", "B7", "B5"],
  min: 0,
  max: 30000
};
//Map.addLayer(landImage, trueColor, "Landsat False Color");

        //RAW BANDS CHARTS
  // ASTER raw bands scatter chart
// Get a dictionary with band names as keys, pixel lists as values.
var result = asterImage.reduceRegion(ee.Reducer.toList(), DSregion, 120);
// Convert the band data to plot on the y-axis to arrays.
var y1 = ee.Array(result.get('B04'));
var y2 = ee.Array(result.get('B05'));
var y3 = ee.Array(result.get('B06'));
var y4 = ee.Array(result.get('B07'));
var y5 = ee.Array(result.get('B08'));
var y6 = ee.Array(result.get('B09'));
// Concatenate the y-axis data by stacking the arrays on the 1-axis.
var yValues = ee.Array.cat([y1, y2, y3, y4, y5, y6], 1);
// The band data to plot on the x-axis is a List.
var xValues = ee.Array(result.get('B3N'));
// Make a band correlation chart.
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['B04', 'B05', 'B06', 'B07', 'B08', 'B09'])
    .setOptions({
      title: 'ASTER L1T Radiance',
      hAxis: {'title': 'B3N(NIR)'},
      vAxis: {'title': 'B04(SWIR4), B05(SWIR5), B06(SWIR6), B07(SWIR7), B08(SWIR8), B09(SWIR9)'},
      pointSize: 3,
});
// Print the chart.
//print(chart);

  // Landsat - 8 raw bands scatter chart
var result = landImage.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('B6'));
var y2 = ee.Array(result.get('B7'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('B5'));
print('xValues: ', xValues);
print('yValues: ', yValues);
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['B6', 'B7'])
    .setOptions({
      title: 'Landsat - 8 OLI Raw Data',
      hAxis: {'title': 'B5(NIR)'},
      vAxis: {'title': 'B6 (SWIR1) & B7 (SWIR2)'},
      pointSize: 1,
});
// Print the chart.
//print(chart);

  // Sentinel - 2 RAW BAND CHART
var result = sentImage.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('B11'));
var y2 = ee.Array(result.get('B12'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('B8'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['B11', 'B12'])
    .setOptions({
      title: 'Sentinel - 2 L2A Data',
      hAxis: {'title': 'B8(NIR)'},
      vAxis: {'title': 'B11 (SWIR1) & B12 (SWIR2)'},
      pointSize: 1,
});
//print(chart);


              //APPLYING DS ON ALL DATA
// Selecting bands to apply DS
var landBandsImage = landImage.select(landbands);
var asterBandsImage = asterImage.select(asterbands);
var sentBandsImage = sentImage.select(sentbands);

//Obtain DS Results for All Satelites using dcs function
var DSLand = decStr(landBandsImage, DSregion, 1000);
var DSaster = decStr(asterBandsImage, DSregion, 1000);
var DSsent = decStr(sentBandsImage, DSregion, 1000);

//FCC of 3 bands of DS results for all satelites
var selectBands = [0,1,2]; //B5,B6,B7
Map.addLayer(DSLand.select(selectBands), {}, 'DCS Landsat Image');
var selectBands = [3,4,5]; //B6,B7,B8
Map.addLayer(DSaster.select(selectBands), {}, 'DCS Aster Image');
var selectBands = [3,4,5]; //B8,B11,B12
Map.addLayer(DSsent.select(selectBands), {}, 'DCS Sentinel Image');


      //DS BANDS CHARTS
  //Landsat - 8 stretched data scatter chart
var result = DSLand.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('B6'));
var y2 = ee.Array(result.get('B7'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('B5'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['B6', 'B7'])
    .setOptions({
      title: 'Landsat - 8 Data after Decorrelation Stretching',
      hAxis: {'title': 'B5(NIR)'},
      vAxis: {'title': 'B6 (SWIR1) & B7 (SWIR2)'},
      pointSize: 2,
});
//print(chart);

  //ASTER stretched data scatter chart
var result = DSaster.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('B04'));
var y2 = ee.Array(result.get('B05'));
var y3 = ee.Array(result.get('B06'));
var y4 = ee.Array(result.get('B07'));
var y5 = ee.Array(result.get('B08'));
var y6 = ee.Array(result.get('B09'));
var yValues = ee.Array.cat([y1, y2, y3, y4, y5, y6], 1);
var xValues = ee.Array(result.get('B3N'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['B04', 'B05', 'B06', 'B07', 'B08', 'B09'])
    .setOptions({
      title: 'ASTER L1T Data after Decorrelation Stretching',
      hAxis: {'title': 'B3N(NIR)'},
      vAxis: {'title': 'B04(SWIR4), B05(SWIR5), B06(SWIR6), B07(SWIR7), B08(SWIR8), B09(SWIR9)'},
      pointSize: 3,
});
//print(chart);

  //Sentinel 2 stretched data scatter chart
var result = DSsent.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('B11'));
var y2 = ee.Array(result.get('B12'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('B8'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['B11', 'B12'])
    .setOptions({
      title: 'Sentinel - 2 Data after Decorrelation Stretch',
      hAxis: {'title': 'B8(NIR)'},
      vAxis: {'title': 'B11 (SWIR1) & B12 (SWIR2)'},
      pointSize: 2,
});
//print(chart);

            //PRINCIPAL COMPONENT ANALYSIS
      //Applying PCA on DS of Landsat 8
// Obtain the geometry of the region from the raw bands since stretched bands doesnt have geometry
var region = landImage.geometry();
var bands = [0,1,2,3,4];
var image =  DSLand.select(bands);
// Set some paramters for PCA application.
var scale = 30;
var bandNames = image.bandNames();
//Obtain means and stored it in a dictionary
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
// For faster covariance reduction we mean center the data
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);
// This function returns the names of bands in list form
var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  });
};
// Here wer used the centered imagery and the parameters defined to obtain PC of the region specified
var pcImage = PCA(centered, scale, region);
// Obtain FCC of the selected PCs based on Crosta Technique
Map.addLayer(pcImage, {bands: ['pc2', 'pc3', 'pc5'], min: -2, max: 2}, 'Landsat - 8 PCA of DS');


  //Landsat - 8 stretched data PCs scatter chart
var result = pcImage.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('pc2'));
var y2 = ee.Array(result.get('pc3'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('pc1')); 
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['pc2', 'pc3'])
    .setOptions({
      title: 'Landsat 8 PCA of DS',
      hAxis: {'title': 'PC1'},
      vAxis: {'title': 'PC2 and PC3'},
      pointSize: 2,
});
//print(chart);

//Applying PCA on DS of ASTER
var region = asterImage.geometry();
var bands = [0,1,2,3,4,5,6];
var image =  DSaster.select(bands);
var scale = 30;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);
var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  });
};
var pcImage = PCA(centered, scale, region);
Map.addLayer(pcImage, {bands: ['pc6', 'pc7', 'pc5'], min: -2, max: 2}, 'ASTER L1T - PCA of DS');


  //ASTER stretched data PCs scatter chart
var result = pcImage.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('pc2'));
var y2 = ee.Array(result.get('pc3'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('pc1'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['pc2', 'pc3'])
    .setOptions({
      title: 'ASTER L1T - PCA of DS',
      hAxis: {'title': 'PC1'},
      vAxis: {'title': 'PC2 and PC3'},
      pointSize: 2,
});
//print(chart);

  //Applying PCA on DS of Sentinel - 2
var region = sentImage.geometry();
var bands = [0,1,2,3,4,5];
var image =  DSsent.select(bands);
var scale = 30;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);
var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  });
};
var pcImage = PCA(centered, scale, region);
Map.addLayer(pcImage, {bands: ['pc1', 'pc3', 'pc4'], min: -2, max: 2}, 'Sentinel 2 L2C - PCA of DS'); //changed from PC1, PC2 and PC3



  //Sentinel - 2 stretched data PCs scatter chart
var result = pcImage.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('pc2'));
var y2 = ee.Array(result.get('pc3'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('pc1'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['pc2', 'pc3'])
    .setOptions({
      title: 'Sentinel 2 L2C - PCA of DS',
      hAxis: {'title': 'PC1'},
      vAxis: {'title': 'PC2 and PC3'},
      pointSize: 2,
});
//print(chart);

        // DECORRELATION STRETCHING of PRINCIPAL COMPONENT ANALYSIS OF RAW DATA
      //LANDSAT - 8 DS of PCA
  //PCA
var region = landImage.geometry();
var image =  landImage.select(landbands);
var scale = 30;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);
var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  });
};
var pcImage = PCA(centered, scale, region);
Map.addLayer(pcImage, {bands: ['pc2', 'pc3', 'pc1'], min: -2, max: 2}, 'Landsat 8 - PCA');



  //Landsat - 8 raw data PCs scatter chart
var result = pcImage.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('pc2'));
var y2 = ee.Array(result.get('pc3'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('pc1'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['pc2', 'pc3'])
    .setOptions({
      title: 'Landsat 8 - PCA',
      hAxis: {'title': 'PC1'},
      vAxis: {'title': 'PC2 and PC3'},
      pointSize: 2,
});
//print(chart);

  //Stretch PCs of Landsat - 8 raw data using DS
var DSofPCA = decStr(pcImage, DSregion, 1000);
selectBands = [1,2,0];
Map.addLayer(DSofPCA.select(selectBands), {}, 'Landsat 8 - DS of PCA');

  //Landsat - 8 raw data streched PCs scatter chart
var result = DSofPCA.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('pc2'));
var y2 = ee.Array(result.get('pc3'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('pc1'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['pc2', 'pc3'])
    .setOptions({
      title: 'Landsat 8 - DS of PCA',
      hAxis: {'title': 'PC1'},
      vAxis: {'title': 'PC2 and PC3'},
      pointSize: 2,
});
//print(chart);

      //ASTER L1T DS of PCA
  //PCA
var region = asterImage.geometry();
var image =  asterImage.select(asterbands);
var scale = 30;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);

var pcImage = PCA(centered, scale, region);
Map.addLayer(pcImage, {bands: ['pc5', 'pc6', 'pc7'], min: -2, max: 2}, 'ASTER L1T - PCA');

  // Plot each PC as a new layer
for (var i = 0; i < bandNames.length().getInfo(); i++) {
  var band = pcImage.bandNames().get(i).getInfo();
  //Map.addLayer(pcImage.select([band]), {min: -2, max: 2}, band);
}

  //ASTER raw data PCs scatter chart
var result = pcImage.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('pc2'));
var y2 = ee.Array(result.get('pc3'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('pc1'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['pc2', 'pc3'])
    .setOptions({
      title: 'ASTER L1T - PCA ',
      hAxis: {'title': 'PC1'},
      vAxis: {'title': 'PC2 and PC3'},
      pointSize: 2,
});
//print(chart);


  // Stretch PCs of ASTER raw data using DS
var DSofPCA = decStr(pcImage, DSregion, 1000);
selectBands = [0,1,2];
Map.addLayer(DSofPCA.select(selectBands), {}, 'ASTER L1T - DS of PCA');

  //ASTER raw data streched PCs scatter chart
var result = DSofPCA.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('pc2'));
var y2 = ee.Array(result.get('pc3'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('pc1'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['pc2', 'pc3'])
    .setOptions({
      title: 'ASTER L1T - DS of PCA',
      hAxis: {'title': 'PC1'},
      vAxis: {'title': 'PC2 and PC3'},
      pointSize: 2,
});
//print(chart);

      //Sentinel - 2 DS of PCA
  //PCA
var region = sentImage.geometry();
var image =  sentImage.select(sentbands);
var scale = 30;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);


var pcImage = PCA(centered, scale, region);
// Plot each PC as a new layer
Map.addLayer(pcImage, {bands: ['pc4', 'pc5', 'pc3'], min: -2, max: 2}, 'Sentinel 2 - PCA');

  // Plot each PC as a new layer
for (var i = 0; i < bandNames.length().getInfo(); i++) {
  var band = pcImage.bandNames().get(i).getInfo();
  //Map.addLayer(pcImage.select([band]), {min: -2, max: 2}, band);
}

  // Sentinel - 2 raw data PCs scatter chart
var result = pcImage.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('pc2'));
var y2 = ee.Array(result.get('pc3'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('pc1'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['pc2', 'pc3'])
    .setOptions({
      title: 'Sentinel 2 - PCA',
      hAxis: {'title': 'PC1'},
      vAxis: {'title': 'PC2 and PC3'},
      pointSize: 2,
});
//print(chart);


  // Stretch PCs of Sentinel - 2 raw data using DS
var DSofPCA = decStr(pcImage, DSregion, 1000);
selectBands = [2,3,4];
Map.addLayer(DSofPCA.select(selectBands), {}, 'Sentinel 2 - DS of PCA');

  //Sentinel - 2 raw data stretched PCs scatter chart
var result = DSofPCA.reduceRegion(ee.Reducer.toList(), DSregion, 120);
var y1 = ee.Array(result.get('pc2'));
var y2 = ee.Array(result.get('pc3'));
var yValues = ee.Array.cat([y1, y2], 1);
var xValues = ee.Array(result.get('pc1'));
var chart = ui.Chart.array.values(yValues, 0, xValues)
    .setSeriesNames(['pc2', 'pc3'])
    .setOptions({
      title: 'Sentinel 2 - DS of PCA',
      hAxis: {'title': 'PC1'},
      vAxis: {'title': 'PC2 and PC3'},
      pointSize: 2,
});
//print(chart);

Map.setCenter(72.2653928137903,34.25603836981476, 13);
