/*
The following Earth Engine (EE) script was created to demonstrate a pure-data approach to predictive modeling in archeology.
It compares various GIS datasets with the coordinates for a representative 
sample of sites, using T-score analysis to determine the influence of each geographic variable on archeological site location.
Then, it evaluates all points in the region to determine their correlation to the site collection.

This is intended to test the hypothesis that higher-conforming locations are 
more likely to play host to an undiscovered archeological site that is 
culturally consistent with the reference sample.

However, the number of calculations required to cover large areas quickly exceeds the memory available to EE users,
limiting the utility of this specific script to demonstrative purposes.
The following code can be used as a prototype to develop an independent algorithm that does not rely on EE memory.

Currently, the analysis performed is limited to simple derivations of 
elevation data, but it can be expanded to consider any number of datasets, limited only by availability and user computing power.
(For example, an implementation could include cost-distance analysis for nearby hydrological features, or other such resources.)
*/

/// -----------------------------------------------

// Define the map center & region of interest.
var mapCenter = [16.59, 20.45]; // TODO: [14, 20];
var roiDims = [0.003, 0.002]; // TODO: [30, 20];
var roi = ee.Geometry.Rectangle({
  coords: [
	mapCenter[0] - (roiDims[0] / 2),
	mapCenter[1] - (roiDims[1] / 2),
	mapCenter[0] + (roiDims[0] / 2),
	mapCenter[1] + (roiDims[1] / 2),
  ],
  geodesic: false,
});
Map.setCenter(mapCenter[0], mapCenter[1], 18);

// GIS images that comprise the model:
var srtm = ee.Image('USGS/SRTMGL1_003').select('elevation');
var slope = ee.Terrain.slope(srtm);

// Collection of coordinates describing the set of representative known sites.
// Here, they are randomly generated within ROI bounds.
var siteGeoPoints = ee.FeatureCollection.randomPoints(roi, 10);
print('Sample Sites', siteGeoPoints);

/// -----------------------------------------------

// Calculate a site's prominence relative to its immediate surroundings.
function relativeProminence(point) {
  var bufferPoint = point.buffer(90); // Considers a 90 m radius around the site.
  // Retrieve local min, max, and actual elevation.
  var elev = srtm.reduceRegion({reducer: ee.Reducer.mean(), geometry: point, scale: 30}).get('elevation');
  var min = srtm.reduceRegion({reducer: ee.Reducer.min(), geometry: bufferPoint, scale: 30}).get('elevation');
  var max = srtm.reduceRegion({reducer: ee.Reducer.max(), geometry: bufferPoint, scale: 30}).get('elevation');
  var dmax = ee.Number(max).subtract(elev);
  var dmin = ee.Number(elev).subtract(min);
  /*
  Represents the relative prominence as a ratio of nearby highs to lows.
  Larger values represent a site that is found near a local elevation minimum,
  whereas fractional values represent a site found near the local maximum.
  Does not account for the flatness or irregularity of the terrain;
  future implementations could also consider the elevation range as a factor.
  */
  return dmax.divide(dmin);
}

// Calculates, sets, and then extracts geographic data from a site collection.
function extractData(siteCollection) {
  // Reducer function for retrieving site data from a GIS image.
  var reduce = function(image, site) {
	return image.reduceRegion({ reducer: ee.Reducer.mean(), geometry: site.geometry(), scale: 30 });
  };
  // Call analysis functions in a map chain to set site properties.
  var siteDataCollection = siteCollection
	.map(function(site) { 	// Calculate the elevation.
  	var siteElev = reduce(srtm, site).get('elevation');
  	return site.set('elevation', siteElev);
	}).map(function(site) {   // Calculate the slope.
  	var siteGrade = reduce(slope, site).get('slope');
  	return site.set('grade', siteGrade);
	}).map(function(site) {   // Calculate the prominence.
  	var siteProm = relativeProminence(site.geometry());
  	return site.set('prominence', siteProm);
	});                   	// Chain additional calculations here.
 
  // Separate collection data into individual arrays, by field.
  var aggregatedElev = siteDataCollection.aggregate_array('elevation');
  var aggregatedGrade = siteDataCollection.aggregate_array('grade');
  var aggregatedProm = siteDataCollection.aggregate_array('prominence');

  return {
	elevation: aggregatedElev,
	grade: aggregatedGrade,
	prominence: aggregatedProm,
  };
}

// Store the aggregated data arrays in a more convenient way.
var siteData = extractData(siteGeoPoints);
var elevData = siteData.elevation;    
var slopeData = siteData.grade;  	 
var promData = siteData.prominence;   
print('Elevation data', elevData);  
print('Slope data', slopeData);  
print('Prominence data', promData);

/// -----------------------------------------------

// Calculate an array's mean, standard deviation, and coefficient of variance together.
function calculateMeanStdDev(dataArray) {
  var mean = ee.Number(dataArray.reduce(ee.Reducer.mean()));
  var sqDiff = dataArray.map(function (val) { return ee.Number(val).subtract(mean).pow(2); });
  var stdDev = ee.Number(sqDiff.reduce(ee.Reducer.mean())).sqrt();
  var cv = stdDev.divide(mean).multiply(100);

  return {
	mean: mean,
	stdDev: stdDev,
	cv: cv,
  };
}

// Calculate the t-score for a single value.
/* NOTE: We use t-score here for simplicity and flexibility, assuming a 
   roughly normal distribution, but accounting for possibly small or 
   non-normal distributions in the representative sample of sites.
   Similar methods of statistical analysis could be used in its place, 
   depending on the site data.  */
function calculateTScore(val, mean, stdDev) {
  return ee.Number(val).subtract(mean).divide(stdDev);
}

// Holds the mean and std deviation for cross reference.
var elevFit = calculateMeanStdDev(elevData);
var slopeFit = calculateMeanStdDev(slopeData);
var promFit = calculateMeanStdDev(promData);
print('Elevation Fit', elevFit);
print('Slope Fit', slopeFit);
print('Prominence Fit', promFit);

// Weigh scores based on the coefficient of variation's inverse.
// (Higher conformity is weighted more heavily.)
var elevWeight = ee.Number(1).divide(elevFit.cv);
var slopeWeight = ee.Number(1).divide(slopeFit.cv);
var promWeight = ee.Number(1).divide(promFit.cv);

/// -----------------------------------------------

// Analysis function for site candidate points, similar to extractData().
function analyzePoint(point) {
  // Reducer function for retrieving point data from a GIS image.
  var reduce = function(image) {
	return image.reduceRegion({ reducer: ee.Reducer.mean(), geometry: point, scale: 30 });
  };
  // Calculate all variables as you would for a site.
  var elevation = reduce(srtm).get('elevation');
  var grade = reduce(slope).get('slope');
  var prominence = relativeProminence(point);
  // Calculate the T-score for each variable, using the mean and standard of deviation of the site collection.
  /* NOTE: Another approach would be to recalculate the mean and stdDev, 
     including each point as a temp site.
     This would effectively "fill in" potential gaps in the site collection 
     using the site candidate, to test if the candidate point would increase 
     the site collection's coherency or diminish it.  */
  var elevScore = calculateTScore(elevation, elevFit.mean, elevFit.stdDev).abs().multiply(elevWeight);
  var slopeScore = calculateTScore(grade, slopeFit.mean, slopeFit.stdDev).abs().multiply(slopeWeight);
  var promScore = calculateTScore(prominence, promFit.mean, promFit.stdDev).abs().multiply(promWeight);
  /* NOTE: Another approach would be to set these scores as separate fields to 
     the point feature, for independent analysis at a later stage.
     For this example we can simply combine them now and return the overall 
     correlation.   */
  return elevScore.add(slopeScore).add(promScore);
}

/// -----------------------------------------------

// Create an image of raw longitude & latitude bands.
var lonLat = ee.Image.pixelLonLat();
// Get the corners of the ROI.
var roiCorners = ee.List(roi.coordinates().get(0));
// Create a vertical and a horizontal line along the boundaries of the ROI.
var roiHeight = ee.Geometry.LineString([roiCorners.get(1), roiCorners.get(2)]);
var roiWidth = ee.Geometry.LineString([roiCorners.get(0), roiCorners.get(1)]);
// Count the number of pixels intersecting the line, to get ROI dimensions in number of pixels.
var yPixels = ee.Number(lonLat.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: roiHeight,
  scale: 30,
  maxPixels: 2e10,
}).get('longitude'));
var xPixels = ee.Number(lonLat.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: roiWidth,
  scale: 30,
  maxPixels: 2e10,
}).get('latitude'));

// Creates a feature collection of points representing pixels within the ROI.
function pointMatrix() {
  var corner1 = roiCorners.get(0);
  var corner2 = roiCorners.get(1);
  var corner3 = roiCorners.get(3);
  // Calculate the spacing between pixels.
  var spacingX = ee.Number(ee.List(corner2).get(0)).subtract(ee.List(corner1).get(0)).divide(xPixels);
  var spacingY = ee.Number(ee.List(corner3).get(1)).subtract(ee.List(corner2).get(1)).divide(yPixels);
  // The last pixel to calculate.
  var endX = xPixels.subtract(1);
  var endY = yPixels.subtract(1);
 
  // Create a new feature with point geometry, from lon/lat coordinates.
  function newPointFeat(coords) {
	return ee.Feature(ee.Geometry.Point(coords));
  }

  // Iteratively generate points based on ROI pixel dimensions.
  var points = ee.List.sequence(0, endX).map(function(x) {
	return ee.List.sequence(0, endY).map(function(y) {
  	var xIncr = ee.Number(x).multiply(spacingX);
  	var xCoord = ee.Number(ee.List(corner1).get(0)).add(xIncr);
  	var yIncr = ee.Number(y).multiply(spacingY);
  	var yCoord = ee.Number(ee.List(corner1).get(0)).add(yIncr);
  	var pointFeat = newPointFeat([xCoord, yCoord]);
  	return ee.Feature(pointFeat);
	});
  }).flatten();

  return ee.FeatureCollection(points);
}

// Initialize the matrix of point features.
var pointMtx = pointMatrix();

/// -----------------------------------------------

// Analyze each point in the matrix, saving the result as part of the feature.
var analyzedFeats = pointMtx.map(function(feat) {
  return feat.set('tScore', analyzePoint(feat.geometry()));
});

// Finds the highest tScore (least suitable).
var leastFit = ee.Number(analyzedFeats.reduceColumns({ reducer: ee.Reducer.max(), selectors: ['tScore'] }).get('max'));
// Finds the lowest tScore (most suitable).
var mostFit = ee.Number(analyzedFeats.reduceColumns({ reducer: ee.Reducer.min(), selectors: ['tScore'] }).get('min'));

print('tScore range:', leastFit, mostFit);


// Color gradient from red to green.
var suitabilityGradient = ['red', 'green'];

// Determine relative fitness value for each feature point.
var suitability = analyzedFeats.map(function(feat) {
  var pointScore = ee.Number(feat.get('tScore'));
  var fitRange = leastFit.subtract(mostFit);
  // Normalize fitness score to the range of t-score values.
  var fitness = pointScore.subtract(mostFit).divide(fitRange);
  // A lower fitness score indicates higher suitability.
  return feat.set('normalizedTScore', fitness);
});
print('Suitability:', suitability);

// Create a new image displaying the matrix of color-coded points.
var suitabilityImage = ee.Image().byte().paint({
  featureCollection: suitability,
  color: 'normalizedTScore',
});

// DEBUG: I have no idea why this doesn't display anything:
//Map.addLayer(suitabilityImage, {palette: suitabilityGradient}, 'Suitability');

/// -----------------------------------------------

// Sort by fitness in ascending order, then take the top 10% of points.
/* NOTE: Another approach could be to set a threshold for the suitability 
   value, and only take the points that exceed that threshold.
   That would avoid limiting the number of sites returned in a region with 
   many highly-conforming  points, and would also avoid returning a number of 
   useless points for low-conformity areas.
   One possible way to determine this threshold would be to perform another
   round of t-score analysis, this time on the inverse of the normalizedTScore 
   values, and select the points based on the resulting standard of deviation, 
   for example, by selecting only sites that are two standards of deviation 
   greater than the mean, discarding the points that are close to the average 
   fitness value, or less. (Remember that lower normalizedTScore values are 
   indicative of superior fitness.)
   That said, the top 10% of the collection is suitable for the purposes of 
   this prototype.  */
var suitabilityAsc = suitability.sort('normalizedTScore', true);
var bestFitting = suitabilityAsc.limit(ee.Number(suitabilityAsc.size()).multiply(0.1).round());
print('Top 10% Most Suitable Points:', bestFitting);

var suitabilityDesc = suitability.sort('normalizedTScore', false);
var leastFitting = suitabilityDesc.limit(ee.Number(suitabilityDesc.size()).multiply(0.1).round());
print('Top 10% Least Suitable Points:', leastFitting);
