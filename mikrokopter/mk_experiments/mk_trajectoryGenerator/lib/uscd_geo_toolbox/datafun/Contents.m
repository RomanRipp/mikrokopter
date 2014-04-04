% Data analysis and Modifications # Contents of "toolbox/datafun"
%   
% Functions take care about NaN's:  COVNAN, GRADNAN, SUMNAN
%                                   MEDNAN, MEANNAN, STDNAN
% 
% Interpolation: AKIMA, INTERP0, INTERPG, INTERPC, INTERPM
%                REPL_NAN, FILLNAN
% 
% Filtering: MEANIND1, MEANIND2, MFILTER, GRPMEAN, NOISE
% 
% Indexing: IND2GRP, GRP2IND, INSERT, REMFOLLW
% 
% Mapping: OBJMAP, OBANA, KRIGING
% 
% Griding: MONOGRID, MKGRID
% 
% other: INELLIP, ONLINE, PSPIKE, ZEROSIZE
%        HEPIRO, JOINVEC, GET_EOF, ISOAREA
%   
%---------------------------------------------------------------------------------
%   
%   obana/  - MEX-Files for OBANA2
%   
%   akima     - Akima spline interpolation
%   bfilter   - Butterworth digital and analog filter (BUTTER)
%   center    - Transforms Data into Intervall
%   chk_lon   - Check if Longitude in correct Period for Limits
%   cor4flat  - Correct Matrices to use SURFACE with Flat shading
%   cos_line  - 1-D Interpolation via Cosine
%   cos_surf  - 2-D Interpolation via Cosine
%   covnan    - Covariance matrix and Correlation Coefficients, ignores NaN
%   cyclic    - Creates a monotonic Vector of cyclic Values
%   fillnan   - Interpolate over NaN's in Matrice
%   findextr  - Find Extrema of Vector
%   gauss     - Gaussian function
%   get_eof   - Calculates eigen-vectors, -values and the covariance matrix
%   grad      - Gradient in 1. and 2. Order
%   gradnan   - Approximate gradient
%   grid2tri  - Triangulation of Grid
%   grp2ind   - Built IndexVector from StartIndex and Length
%   grpmean   - Mean Value at specified Groups
%   hepiro    - Rotates XYZ-Coordinates by Heading, Pitch and Roll
%   ind2grp   - Built StartIndex and Length from IndexVector
%   inellip   - True for points inside elliptical region
%   inpoly    - True for points inside polygonal region
%   insert    - Inserts Values into Array
%   insert00  - INSERT Inserts Values into Array
%   interp0   - Interpolate linear over gabs
%   interp21  - Interpolates a Single Point into a 2-D function
%   interpc   - Interpolates Values of Matrice to Edges of specified Column
%   interpg   - Interpolate coordinates of polygon on edges of grid
%   interpm   - N-D interpolation (table lookup)
%   ismesh    - True if the inputs should be automatically meshgridded
%   isoarea   - Returns Area of IsoPatches
%   joinvec   - Returns continuous Lines of NaN-Row separated Segments
%   joinvek   - Returns continuous Liness of NaN-Row separated Segments
%   kriging   - Interpolation from irregular points by Kriging
%   mat2ind   - Converts Matrice to Indexed
%   maxnan    - MAXMISS Column maximum with missing data
%   meanind1  - 1-Dimensional Running Mean
%   meanind2  - 2-Dimensional Running Mean of Matrice
%   meannan   - Column mean with missing data
%   meanval1  - Running Value-Mean along 1. Dimension
%   meanval2  - Running Value-Mean along 1. and 2. Dimension
%   medmean   - 1-Dimensional Median Filter
%   mednan    - Median value, ignores NaN
%   medsamp   - Lowsampling of TimeSeries using median filter
%   mfilter   - FIR digital filter
%   minnan    - Column minimum with missing data
%   mkgrid    - Generation of common Grid from Inputs
%   modline   - Returns modified Line with fixed derivatives at Start and End
%   mono      - Monotonize vector
%   monogrid  - Make grid Monotonic in X and Y
%   noise     - Generate noises with periodicity
%   obana     - Simple objective mapping
%   obana2    - OBANA3 fast data gridding using Gaussian weights.
%   obana3    - Fast data gridding using Gaussian weights.
%   objmap    - Objective mapping interpolation
%   online    - True for Points near a Curve
%   paxes     - Principle Axes of 2-Component TimeSeries
%   pspike    - Markiere einzelne Spikes und interpoliere ggf
%   remfollw  - Remove Following equal Elements from Matrice
%   repl_nan  - Replace NaN-Values by mean of surrounding non-NaN
%   sampline  - Smooth and lowsample Vectors
%   sep_line  - Transforms and Separates Data by NaN
%   sort_grp  - Sort a Matrix by Rows and in Order of specified Columns
%   splinen   - Calculates Spline
%   stdnan    - Column standard deviation with missing data
%   str2int   - Converts bitwise between Strings and Vectors of Integers
%   sumnan    - Sum of elements with missing data
%   surfarea  - Surface Area
%   uniqued   - Returns a Matrice with no repetitions along a specified Dimension
%   xcorr6    - Cross-correlation function estimates, Matlab-6.5-Version
%   xy_rot    - Rotates Matrix in X-Y-LAYER
%   zerosize  - Returns Depth of ZeroFields in Matrice
%   
