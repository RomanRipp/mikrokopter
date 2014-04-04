function [y,i] = minnan(x,dim)
%MINNAN Column minimum with missing data.
%  Y = MINNAN(X,DIM) returns the smallest element in the Dimension DIM of X.
%  Missing data values must be encoded as NaNs. For vectors, MINNAN(X)
%  returns the smallest value of the elements in X.
%  [Y,I] = MINNAN(X,...) stores the indices of the minimum values in vector I.
%
%  See also MAXNAN, MEANNAN.

%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1995/03/08 14:27:31 $
%
% added ~isempty	G.Krahmann, IfM Kiel, Oct 1995
% added help function	G.Krahmann, IfM Kiel, Jun 1996

if nargin==0
  help minnan
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
  x(notvalid) = Inf*notvalid;
  if dim == 0
   [y,i] = min(x);
  else
   [y,i] = min(x,[],dim);
 end
else
  y=nan;
end
