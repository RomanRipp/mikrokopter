function y = sumnan(x,dim)
%SUMNAN Sum of elements with missing data.
%  Y = SUMNAN(X) returns the sum of each column of X as a row vector
%  where missing data values are encoded as NaNs. For vectors, SUMNAN(X)
%  returns the sum of the elements in X.

%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1995/09/27 14:45:28 $
%
% added ~isempty	G.Krahmann, IfM Kiel, Oct 1995

if nargin==0
  help sumnan
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
  x(notvalid) = 0*notvalid;
  if dim == 0
   y = sum(x);
  else
   y = sum(x,dim);
  end
else
  y=nan;
end

