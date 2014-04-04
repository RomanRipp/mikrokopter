function [ok,v,n] = chkcvec(v,opt)

% CHKCVEC  Checks Input for CellArray of numeric Vectors
%
%  [ok,C,N] = chkcvec(V,Option)
%
%  V CellArray of Vectors or 2-dimensional numeric Array
%
%  Option ~= 0 ==> Numeric Arrays not allowed,
%
%   default: Option == 0   ==>  Numeric Arrays --> CellArray
%
%  C CellArray of Vectors (single Row or single Column)
%  N Numeric Array with a Vector of C in each Row
%
%
 
if nargin < 2
   opt = 0;
end

n = [];

if isnumeric(v)
   n = v;
   v = {};
   ok = ( isequal(opt,0) & ( ndims(v) == 2 ) );
   if ok & ~isempty(n)
      m = size(n,1);
      v = cell(m,1);
      for ii = 1 : m
          v(ii) = { n(ii,:) };
      end
   end
   return
end

ok = iscell(v);
if ~ok
   return
end

for d = [ 2  1 ]
    ok = 1;
    try
      n = cat(d,v{:});
    catch
      ok = 0;
    end
    if ok
       s = size(n); p = prod(s);
       ok = ( ( p == s(d) ) & isnumeric(n) );
    end
    if ok
       break
    end
end

if ~ok
    return
end

m = prod(size(v));

n = zeros( m , min(1,p) );

if isempty(n)
   return
end

for ii = 1 : m 
    w = v{ii};
    if isempty(w)
       n(ii,:) = NaN;
    else
       w = w(:)';
       s = size(n,2);
       t = size(w,2);
       if s < t
          n = cat( 2 , n , NaN * zeros(m,t-s) );
       end
       n(ii,:) = NaN;
       n(ii,1:t) = w;
    end
end 
