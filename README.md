# arth498_earth_engine
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
