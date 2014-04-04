function tf = ismesh(varargin)
%ISMESH True if the inputs should be automatically meshgridded.
%    ISMESH(X,Y) returns true if X and Y are vectors of
%    different orientations.
%
%    ISMESH(X,Y,Z) returns true if X,Y,Z are vectors of
%    different orientations.
%
%    ISMESH(...) returns true if all the inputs are vectors of
%    different orientations.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%    $Revision: 1.4 $ $Date: 1997/11/21 23:40:30 $

for i=1:length(varargin)
  ns{i} = size(varargin{i})~=1; % Location of non-singleton dimensions
  isvec(i) = sum(ns{i})<=1;     % Is vector.
  nd(i) = ndims(varargin{i});    % Number of dimensions.
end

% True if inputs are 2-D, all vectors, and their non-singleton
% dimensions aren't along the same dimension.
tf = all(isvec) & ~isequal(ns{:});
