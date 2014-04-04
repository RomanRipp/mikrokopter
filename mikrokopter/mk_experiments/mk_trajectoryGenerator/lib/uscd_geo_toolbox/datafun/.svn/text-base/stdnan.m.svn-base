function y = stdnan(x)

% STDNAN   Column standard deviation with missing data.
%
%  STDNAN(X) returns the standard deviation of each column of X as a row
%  vector where missing data values are encoded as NaNs. For vectors,
%  STDNAN(X) returns the standard deviation of the elements in X.

%  This code has been suggested by Douglas M. Schwarz (schwarz@kodak.com)
%  in the news group comp.soft-sys.matlab.
%
%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1995/03/08 14:27:31 $
%
% added ~isempty	G.Krahmann, IfM Kiel, Oct 1995
% added n==0 check	G.Krahmann, IfM Kiel, Jun 1997

if ~isempty(x)
  notvalid = isnan(x);
  [m,n] = size(x(notvalid));
  x(notvalid) = zeros(m,n);
  n = sum(1 - notvalid);
  bad = find( n==0 );
  if ~isempty( bad )
    n(bad) = n(bad)*nan;
  end
  y = sqrt((sum(x.^2) - (sum(x).^2) ./ n) ./ (n-1));
else
  y=nan;
end

