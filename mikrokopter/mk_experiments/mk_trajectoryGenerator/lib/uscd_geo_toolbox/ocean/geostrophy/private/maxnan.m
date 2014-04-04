function [y,i] = maxnan(x,dim)
%MAXMISS Column maximum with missing data.
%  Y = MAXNAN(X,DIM) returns the largest element in Dimension DIM  of X.
%  Missing data values must be encoded as NaNs. For vectors, MAXMISS(X)
%  returns the largest value of the elements in X.
%  [Y,I] = MAXNAN(Y,DIM) stores the indices of the maximum values in vector I.
%
%  See also MEANNAN, MINNAN.

%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1995/03/08 14:27:31 $
%
% added ~isempty	G.Krahmann, IfM Kiel, Oct 1995
% added help function	G.Krahmann, IfM Kiel, Jun 1996

if nargin==0
  help maxnan
  return
end

% check for version
v=version;
if strcmp(v(1),'5')
  if nargin==1
    dim = 1;
    if size(x,1)==1
      dim = 2;
    end
  end
else
  dim=0;
end


if ~isempty(x)
  notvalid = find(isnan(x));
  x(notvalid) = -Inf*notvalid;
  if dim == 0
   [y,i] = max(x);
  else
   [y,i] = max(x,[],dim);
  end
else
  y=nan;
end
