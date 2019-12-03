## MeteoSwiss Data

The directory presents dataset organized from the meteorological measurements provided by [Swiss Federal Office of Meteorology and Climatology (MeteoSwiss)](https://www.meteoswiss.admin.ch/home/climate/swiss-climate-in-detail/climate-normals/normal-values-per-measured-parameter.html).

The dataset is a compilation of 17 types of measurements including temperature, snowfall, precipitation, humidity, sunshine duration, recorded in weather stations distributed over Switzerland. Monthly normals and yearly averages of the measurements calculated based on the time period 1981-2010 are available at 91 stations. For the stations, we are also provided geographical locations in GPS format and altitude values, i.e., meters above sea level.

The file `Meteo.mat` contains the following attributes:
* `records_81_10`: each cell corresponding to one station with 17 type monthly average measurements in the first 12 columns, and the last column belongs to the yearly average.
* `stations`:  name of the stations corresponding to the cell order `records_81_10`.
* `records_name`: name of the measurements corresponding to the order of the rows in one cell of `records_81_10`.
* `GPS`: GPS(lat,long) location of the stations ordered in `stations`.
* `altitude`: Altitude of the stations ordered in `stations`.


