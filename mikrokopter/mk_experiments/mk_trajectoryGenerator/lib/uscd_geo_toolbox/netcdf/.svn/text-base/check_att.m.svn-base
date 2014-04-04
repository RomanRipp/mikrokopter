function [val,par,ind] = check_att(val,att,mfv);

% CHECK_ATT  Checks Value with NetCDF-Attributes for Read/Write
%
% VAL = CHECK_ATT(ATT,VAL,[MissFillValue])  prepare to write
% VAL = CHECK_ATT(VAL,ATT,[MissFillValue])  read Value
%
% ATT = { AttName [AttType] AttValue ... }
%

Nin = nargin;

if Nin < 1
   val = [];
end

if Nin < 2
   att = {};
end

rev = iscell(val);   % True for Reverse

if rev
   tmp = val;
   val = att;
   att = tmp;
end

%-------------------------------------------------------------------

msg = cell(0,1);

is_char = ischar(val);
if ~( isnumeric(val) | is_char | isempty(val) )
     msg = cat(1,msg,{'Value must be a numeric or character-array.'});
end

%--------------------------------------------------------------

s = size(att);
ok = ( iscell(att) & ( size(att,2) >= 2 ) & ( ndims(att) == 2 ) );
if ok & ( s(1) > 0 )
   ok = chkcstr(att(:,1));
   if ok & ( s(2) > 2 ) 
      att = att(:,[ 1  3 ]);
   end
else
   ok = isempty(att);
end


if ~ok
    msg = cat(1,msg,{'ATT must be a CellArray with 2 or 3 Columns: { AttName [AttType] AttValue }'});
end

%--------------------------------------------------------------

if Nin < 3
   mfv = [];
end

if isempty(mfv)
   mfv = [NaN NaN];
else
   p = prod(size(mfv));
   ok = ( isnumeric(mfv) & ( p <= 2 ) );
   if ok & ( p == 1 )
      mfv = mfv([1 1]);
   elseif ~ok
      msg = cat(1,msg,{'MissFillValue must be a max. 2-elemnt numeric.'});
   end
end

%--------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    error(msg)
end

%***************************************************************

par = { 'MissVal' {'missing_value'} {[  NaN  mfv(1) ]} 
        'FillVal' {'_FillValue'}    {[  NaN  mfv(2) ]}
        'Minimum' {'valid_min'}     {[ -Inf ]}
        'Maximum' {'valid_max'}     {[ +Inf ]}
        'Range'   {'valid_range'}   {[ -Inf +Inf ]}
        'Offset'  {'add_offset'}    {[  0 ]}
        'Scale'   {'scale_factor'}  {[  1 ]}            };

par = permute(par,[2 1]);

fld = par(1,:);

ind = par([1 2],:);
ind(2,:) = { {[]} };

nom = par([1 2],:);
par = par([1 3],:);

par = struct(par{:});
ind = struct(ind{:});
nom = struct(nom{:});

if ~isempty(att)
    for f = fld
        ii = strcmp( att(:,1) , getfield(nom,f{1}) );
        if any(ii)
           ii = find(ii);
           at = att{ii,2};
           if isequal(at,'\0') 
              at = char(0);
           end
           na = 1 + strcmp(f{1},'Range');
           if prod(size(at)) == na
              ind  = setfield( ind , f{1} , ii );
              v    = getfield( par , f{1} ); 
              v(1:na) = at(:)';
              par  = setfield( par , f{1} , v );
           end
        end
    end
end

if isempty(val) | isempty(att)
   return
end

%***************************************************************
% Check for Miss/Fill-Values first
%***************************************************************

rpl = struct( 'MissVal' , {[]} , ...
              'FillVal' , {[]}  );

rr = 1 + rev;

if ~isempty(att)
    for f = fieldnames(rpl)'
        if ~isempty(getfield(ind,f{1}))
            v = getfield(par,f{1});
            v = v(rr);
            if isnan(v)
               ii = isnan(val);
            else
               ii = ( val == v );
               if rev & strcmp(f{1},'MissVal')
                  ii = ( ii | isnan(val) );
               end               
            end
            rpl = setfield(rpl,f{1},find(ii));
        end
    end
end


%***************************************************************
if ~rev
%***************************************************************
% Check Array which was Read from NetCDF-File with Attributes
%***************************************************************

    %---------------------------------------
    % Scale and Add Offset

    if ~isequal(par.Scale,1)
        val = double(val) * par.Scale;
    end

    if ~isequal(par.Offset,0)
        val = double(val) + par.Offset;
    end

 
%***************************************************************
else
%***************************************************************
% Prepare Array to write
%***************************************************************

    %------------------------------------------------
    % Check for Valid_Min, Valid_Max ofr ValidRange,
    %   Set to FillValue or MissingValue

    is_out = zeros(size(val)); % Values outside valid Range

    if ~isempty(ind.Minimum)
        is_out = ( is_out | ( val < par.Minimum(1) ) );
    end

    if ~isempty(ind.Maximum)
        is_out = ( is_out | ( val > par.Maximum(1) ) );
    end

    if ~isempty(ind.Range)
        is_out = ( is_out | ( val < min(par.Range) ) ...
                          | ( max(par.Range) < val )     );
    end

    is_out = find(is_out);

    %------------------------------------------------
    % Check for Scale and Offset

    if ~isequal(par.Offset(1),0)
        val = val - par.Offset(1);
    end

    if ~isequal(par.Scale(1),1)
        val = val / par.Scale(1);
    end

    %------------------------------------------------

    if ~isempty(is_out) 
        for f = { 'FillVal'  'MissVal' }
            if ~isempty(getfield(ind,f{1}))
                v = getfield(par,f{1});
                val(is_out) = v(1);
                break
            end
        end
    end

end

%***************************************************************
% Replace Miss and Fill
%***************************************************************

rr = 3 - rr;

for f = fieldnames(rpl)'
    ii = getfield(rpl,f{1});
    if ~isempty(ii) & ~isempty(getfield(ind,f{1}))
        v = getfield(par,f{1});
        v = v(rr);
        if ischar(val) & isnan(v)
           v = ' ';
        end
        val(ii) = v;
    end
end

%**********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );

%**********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

