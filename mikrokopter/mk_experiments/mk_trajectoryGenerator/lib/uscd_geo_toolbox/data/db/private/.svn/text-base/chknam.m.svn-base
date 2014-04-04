function ini = chknam(ini,mode);

% CHKNAM  Checks VariableIdentifer by DBIDENT for duplicate Names
%
% Duplicate Names will removed
%

if isempty(ini)
   return
end

if isempty(ini.ID)
   return
end

if nargin < 2
   mode = 1;
end

%--------------------------------------------------------
% Check for Duplicate Names

[h,si] = sort(ini.Name);

h = diff(double(char(h)),1,1);     % Diff of sorted Strings

h = ( sum(h==0,2) == size(h,2) );  % equal Strings

h = cat( 1 , 0 , h );              % 1. String unique

if ~any(h)
    return
end

%--------------------------------------------------------
% Check InputTypes for Header & Data

if isequal(mode,1)

   %--------------------------------------------------
   % Find Groups of equal Names,
   % compare MinTyp and MaxTyp 

   g  = double(~h);
   i0 = find(g);
   g  = cumsum(g,1);    % Groups of equal Names

   ts = ini.Type(si);

   t0 = ts;
   [hh,ii] = sort(t0);
   [hh,jj] = sort(g(ii));
   t0 = t0(ii(jj));
   t0 = t0(i0);

   t1 = ts;
   [hh,ii] = sort(-t1);
   [hh,jj] = sort(g(ii));
   t1 = t1(ii(jj));
   t1 = t1(i0);

   ini.Type(si) = 3 + ( ts - 3 ) .* double( t0(g) == t1(g) );

end

%--------------------------------------------------------
% Remove Duplicate Names

 h  = si(find(h));

fld = fieldnames(ini);

for f = fld(:)'
    v      = getfield(ini,f{1});
    v(h,:) = [];
    ini    = setfield(ini,f{1},v);
end
