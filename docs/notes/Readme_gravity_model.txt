The idea is to compute transport costs between regions in less than proportional way to the distance (e.g. log(1.05,x)), and as distance compute it from the center of gravity of forest resources of the export region to the center of gravity of the population for the importing region.


Does a package (julia preferibly) exist to compute the center of gravity by zones/regions, provided a region shapefile and a thematic raster? For example I want the center of gravity of forest recourses so I provide a raster of corine land cover for forest and a shapefiles of administrative borders for NUTS1 (states) and I got a gravity center for each state (either as a point shapefile or, perhaps better, as a table) ?

