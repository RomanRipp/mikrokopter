function flat = ecc2flat(parm)

%ECC2FLAT  Computes the flattening of an ellipse given an eccentricity
%
%  f = ECC2FLAT(mat) computes the flattening of an ellipse (or
%  ellipsoid of revolution) given the eccentricity.  If the
%  input is a column vector, then each input is assumed to be an
%  eccentricity.  If the input has two columns, then the second
%  column is assumed to be the eccentricity.  This allows geoid
%  vectors from ALMANAC to be used as inputs.  If the input
%  is an n x m matrix, where m,n > 1, then each element is assumed
%  to be an eccentricity.
%
%  See also:  FLAT2ECC, ECC2N, MAJAXIS, MINAXIS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.7 $    $Date: 1998/08/10 17:47:39 $

if nargin == 0
    error('Incorrect number of arguments')
end

%  Dimension tests

if min(size(parm)) == 1 & ndims(parm) <= 2

	col = min(size(parm,2), 2);   % Select first or second column
	eccent = parm(:,col);         % First col if single vector input
	                              % Second col if two column inputs (eg. geoid vecs)
else
    eccent = parm;        %  General matrix input
end

%  Ensure real inputs

if ~isreal(eccent)
   warning('Imaginary parts of complex argument ignored')
   eccent = real(eccent);
end

%  Compute the flattening

flat = 1 - sqrt(1 - eccent.^2);