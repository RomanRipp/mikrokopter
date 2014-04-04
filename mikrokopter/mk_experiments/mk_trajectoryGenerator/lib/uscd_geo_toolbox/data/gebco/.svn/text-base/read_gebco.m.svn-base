function [xy,i0,lg,lv,mm,ok] = read_gebco(file,area,lev);

% READS IsoLines from GEBCO-Ascii-Files
%
% [ XY , I0 , LG , LV , MXY , OK ] = READ_GEBCO( File )
%
% File:  FileName(s) to read, CharacterArray or CellStringArray
%
%   XY = [ Lon  Lat ], NaN-separated Segments
%   I0 = StartIndex of Segment in XY
%   LG = Length of Segment
%   LV = Level of Segment
%  MXY = Center of Segment [ MX  MY ]
%
%   OK = Ok-Flag for File(s)
%
% READ_GEBCO( File , Area , Level ) read Data matching specified Area:
%
%  Area:   [ LonMin LonMax  LatMin  LatMax ], [ 1 by 4 ]
%  Level:  Vector for Levels to read
%
%

Nin = nargin;

msg = cell(0,1);

if Nin < 2
   area = [];
end

if Nin < 3
   lev = [];
end

%---------------------------------------------

[ok,file] = chkcstr(file);

if ~ok
    msg = cat(1,msg,{'FileName(s) must be Strings.'});
end

ok = isnumeric(area);
if ~ok
    msg = cat(1,msg,{'Area must be numeric.'});
elseif ~isempty(area)
    ok = isequal(size(area),[1 4]);
    if ~ok
        ok = ( Nin < 3 );
        if ok
           lev  = area;
           area = [];
        end
    end
    if ok
       ok = all(isfinite(area));
    end
    if ~ok
        msg = cat(1,msg,{'Area must have 4 finite Elements.'});
    else
        area([1 2]) = area([1 2]+[1 -1]*(area(1)>area(2)));
        area([3 4]) = area([3 4]+[1 -1]*(area(3)>area(4)));
    end
end

ok = isnumeric(lev);
if ~ok
    msg = cat(1,msg,{'Level must be numeric.'});
elseif ~isempty(lev)
    lev = lev(:);
    ok = all(isfinite(lev));
    if ~ok
        msg = cat(1,msg,{'Level must have finite Elements.'});
    end
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    error(msg);
end

%---------------------------------------------

file = file(:);
nf   = size(file,1);

if nf == 1
   [msg,xy,i0,lg,lv,mm] = load_gebco(file{1},area,lev);
   if ~isempty(msg)
       warning(msg);
   end
   return
end

msg = cell(nf,1);

xy = cell(nf,1);
i0 = cell(nf,1);
lg = cell(nf,1);
lv = cell(nf,1);
mm = cell(nf,1);

z0 = 0;

for ii = 1 : nf

    if exist(file{ii},'file') == 2

       [msg{ii},xy{ii},i0{ii},lg{ii},lv{ii},mm{ii}] = load_gebco(file{ii},area,lev);

       if isempty(msg{ii})

          i0{ii} = i0{ii} + z0;

          z0 = z0 + size(xy{ii},1);

       end

    else

        msg{ii} = sprintf('File "%s" doesn''t exist.',file{ii});

    end

    if ~isempty(msg{ii})
        warning(msg{ii});
    end

end

ok = strcmp(msg,'');

if ~any(ok)
    xy = zeros(0,2);
    i0 = zeros(0,1);
    lg = zeros(0,1);
    lv = zeros(0,1);
    mm = zeros(0,2);
    return
end

ii = find(ok);
   
xy = cat(1,xy{ii});
i0 = cat(1,i0{ii});
lg = cat(1,lg{ii});
lv = cat(1,lv{ii});
mm = cat(1,mm{ii});
    
%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,xy,i0,lg,lv,mm] = load_gebco(file,area,lev)

msg = '';

xy = zeros(0,2);
i0 = zeros(0,1);
lg = zeros(0,1);
lv = zeros(0,1);
mm = zeros(0,2);
 
try
   xy = load(file);
catch
   msg = sprintf('Can''t load Data from File "%s".\n%s',file,lasterr);
end

if isempty(msg)
   if ~isnumeric(xy)
       msg = 'Data must be numeric.';
   elseif isempty(xy)
       msg = 'Empty Data.';
   elseif ~( ( ndims(xy) == 2 ) & ( size(xy,2) == 2 ) )
       msg = 'Data must be a 2-Column Matrice.';
   end
   if ~isempty(msg)
       msg = sprintf('Invalid Data in File "%s".\n%s',msg);
   end
end

if ~isempty(msg)
    xy = zeros(0,2);
    return
end

nn = size(xy,1);

z = 1;

while z < nn

      n    = xy(z,2);  % Counts

      if ~( ( mod(n,1) == 0 ) & ( n >= 0 ) )
          msg = sprintf('Counts in Line %.0f of File "%s" must be Positive Integers.',z,file);
          xy = zeros(0,2);
          return
      end

      xy(z,2) = NaN;

      ok = isempty(lev);
      if ~ok
          ok = any( xy(z,1) == lev );
      end
      if ~ok
          xy(z+(0:n),1) = NaN;
      elseif ~isempty(area)
          jj = ( 1 : n ) + z;
          ok = ( ( area(1) <= xy(jj,2) ) & ( xy(jj,2) <= area(2) ) & ...
                 ( area(3) <= xy(jj,1) ) & ( xy(jj,1) <= area(4) )       );
          if ~all(ok)
              if ~any(ok)
                  xy(z+(0:n),1) = NaN;
              else
                  ok = find(~ok);
                  xy(jj(ok),1) = NaN;
              end
          end
      end
   
      z    = z + n + 1;
end

i0 = isnan(xy(:,1));
if any(i0)
   if all(i0)
      xy = zeros(0,2);
      return
   end
   i0 = find(~i0);
   xy = xy(i0,:);
   nn = size(xy,1);
end

i0 = find(isnan(xy(:,2)));
lg = diff(cat(1,i0,nn+1));
lv = xy(i0,1);

xy = xy(:,[2 1]);  % [ Lon Lat ]

xy(i0,:) = 0;

mm = cumsum(xy,1);

% One Record is NaN-Line: "lg-1"
mm = ( mm(i0+lg-1,:) - mm(i0,:) ) ./ (lg(:,[1 1])-1);

xy(i0,:) = NaN;


%******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,str] = chkcstr(str,opt)

% CHKCSTR  Checks Input for CellString, contains Strings !
%
%  [ok,str] = chkcstr(str,Option)
%
%  Option ~= 0 ==> CharacterArrays not allowed,
%
%   default: Option == 0   ==>  CharacterArrays --> CellString
%
 
if nargin < 2
   opt = 0;
end

if strcmp(class(str),'char') & isequal(opt,0)
   n = size(str,1);
   if n == 1
      str = strrep(str,char(32),char(1));
   end
   str = cellstr(str);
   if n == 1
      str = strrep(str,char(1),char(32));
   end
end

ok = iscellstr(str);
if ~ok
   return
end

try
  s = cat(2,str{:});
catch
  ok = 0;
  return
end
 
ok = ( strcmp(class(s),'char')  &  ( prod(size(s)) == size(s,2) ) );


