%# Contents of "toolbox/geo/private"
%   
%   almanac    - Almanac data for the Solar System
%   angl2str   - Angle conversion to a string
%   angledim   - Converts angles from one unit system to another
%   aut2geod   - Converts from authalic latitude to geodetic latitude
%   axes2ecc   - Computes eccentricity given semimajor and semiminor axes
%   azimuth    - Calculates azimuth between points on a geoid
%   cen2geod   - Converts from geocentric latitude to geodetic latitude
%   clipdata   - Clip map data at the -pi to pi border of a display
%   cnf2geod   - Converts from conformal latitude to geodetic latitude
%   deg2dm     - Converts angles from degrees to deg:min vector format
%   deg2dms    - Converts angles from degrees to deg:min:sec vector format
%   deg2km     - Converts distances from degrees to kilometers
%   deg2nm     - Converts distances from degrees to nautical miles
%   deg2rad    - Converts angles from degrees to radians
%   deg2sm     - Converts distances from degrees to statute miles
%   dist2str   - Distance conversion to a string
%   distance   - Calculates distances between points on a geoid
%   distdim    - Converts distances from one unit system to another
%   dms2deg    - Converts angles from deg:min:sec to degrees
%   dms2dm     - Converts angles from deg:min:sec to deg:min vector format
%   dms2mat    - Converts a dms vector format to a [deg min sec] matrix
%   dms2rad    - Converts angles from deg:min:sec to radians
%   ecc2flat   - Computes the flattening of an ellipse given an eccentricity
%   ecc2n      - Computes the parameter n of an ellipse given an eccentricity
%   elpcalc    - Computes volume and surface area for an oblate spheroid
%   epsm       - Calculate the accuracy of the map computations
%   eqa2grn    - Equal area coordinates to Greenwich coordinates
%   eqaazim    - Lambert Equal Area Azimuthal Projection
%   eqaconic   - Albers Equal Area Conic Projection
%   eqacylin   - Equal Area Cylindrical Projection
%   eqdazim    - Equidistant Azimuthal Projection
%   eqdconic   - Equidistant Conic Projection
%   eqdcylin   - Equidistant Cylindrical Projection
%   flat2ecc   - Computes the eccentricity of an ellipse given a flattening
%   gc2sc      - Converts a great circle definition to a center and radius
%   geod2aut   - Converts from geodetic latitude to authalic latitude
%   geod2cen   - Converts from geodetic latitude to geocentric latitude
%   geod2cnf   - Converts from geodetic latitude to conformal latitude
%   geod2iso   - Comutes the isometric latitude given the geodetic latitude
%   geod2par   - Converts from geodetic latitude to parametric latitude
%   geod2rec   - Converts from rectifying latitude to geodetic latitude
%   geoidtst   - Tests for a valid geoid vector
%   grn2eqa    - Greenwich coordinates to equal area cartesian coordinates
%   hms2hm     - Converts time from hrs:min:sec to hr:min vector format
%   hms2hr     - Converts time from hrs:min:sec to hours
%   hms2mat    - Converts a hms vector format to a [hrs min sec] matrix
%   hms2sec    - Converts time from hrs:min:sec to seconds
%   hr2hm      - Converts time from hours to hrs:min format
%   hr2hms     - Converts time from hours to hrs:min:sec vector format
%   hr2sec     - Converts time from hours to seconds
%   interpm    - Interpolates vector data to a specified data separation
%   isieee     - Returns ONE
%   iso2geod   - Computes the geodetic latitude given the isometric latitude
%   km2deg     - Converts distances from kilometers to degrees
%   km2nm      - Converts distances from kilometers to nautical miles
%   km2rad     - Converts distances from kilometers to radians
%   km2sm      - Converts distances from kilometers to statute miles
%   lambcyln   - LAMBCYLIN  Lambert Equal Area Cylindrical Projection
%   lambert    - Lambert Conformal Conic Projection
%   ltln2val   - Returns map code value associated with positions
%   majaxis    - Computes semimajor axis given a semiminor axis and an eccentricity
%   mat2dms    - Converts a [deg min sec] matrix to vector format
%   mat2hms    - Converts a [hrs min sec] matrix to vector format
%   meanm      - Mean for geographic data
%   merccalc   - Transformation of data to and from a Mercator space
%   minaxis    - Computes semiminor axis given a semimajor axis and an eccentricity
%   n2ecc      - Computes the eccentricity of an ellipse given the parameter n
%   nm2deg     - Converts distances from nautical miles to degrees
%   nm2km      - Converts distances from nautical miles to kilometers
%   nm2rad     - Converts distances from nautical miles to radians
%   nm2sm      - Converts distances from nautical miles to statute miles
%   npi2pi     - Truncates angles into the -180 deg to 180 deg range
%   org2pol    - Computes the location of the north pole in a transformed map
%   par2geod   - Converts from parametric latitude to geodetic latitude
%   rad2deg    - Converts angles from radians to degrees
%   rad2dm     - Converts angles from radians to deg:min vector format
%   rad2dms    - Converts angles from radians to deg:min:sec vector format
%   rad2km     - Converts distances from radians to kilometers
%   rad2nm     - Converts distances from radians to nautical miles
%   rad2sm     - Converts distances from radians to statute miles
%   rcurve     - Computes various radii of curvature for an ellipsoid
%   rec2geod   - Converts from rectifying latitude to geodetic latitude
%   reckon     - Computes points at a specified azimuth and range
%   rotatem    - Rotate map data for specified origin and orientation
%   roundn     - Rounds input data at specified power of 10
%   rsphere    - Computes radii for a auxiliary spheres
%   scircle2   - Small circle defined by its center and perimeter
%   sec2hm     - Converts time from seconds to hrs:min vector format
%   sec2hms    - Converts time from seconds to hrs:min:sec vector format
%   sec2hr     - Converts time from seconds to hours
%   sm2deg     - Converts distances from statute miles to degrees
%   sm2km      - Converts distances from statute miles to kilometers
%   sm2nm      - Converts distances from statute miles to nautical miles
%   sm2rad     - Converts distances from statute miles to radians
%   sphcalc    - Computes volume and surface area for a sphere
%   time2str   - Time conversion to a string
%   track      - Connects navigational waypoints with track segments
%   track1     - Track lines defined by starting point, azimuth and range
%   track2     - Track lines defined by a starting and ending point
%   trimdata   - Trim map data exceeding projection limits
%   undoclip   - Removes object clips introduced by CLIPDATA
%   undotrim   - Removes object trims introduced by TRIMDATA
%   unitstr    - Tests for valid unit string or abbreviations
%   utm        - Universal Transverse Mercator Cylindrical Projection
%   utmgeoid   - Recommended UTM geoids for zone
%   utmzone    - Universal Transverse Mercator zone
%   utmzoneui  - UTM zone picker
%   vec2mtx    - Regular matrix map from vector data
%   zero22pi   - Truncates angles into the 0 deg to 360 deg range
%   zerom      - Constructs a regular matrix map of all zeros
%   
