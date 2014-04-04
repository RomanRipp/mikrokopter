function [hl,ht] = lineaxis(axe,varargin);

% LINEAXIS Creates an AxisLine inside an Axes
%
%  LINEAXIS( AxesHandle , 'Type'   , 'x' | 'y'        , ...
%                      'Location'  ,  yc |  xc | NaN  , ...
%                      'TickAlign' , -1 | 0 | 1       , ...
%                      'TextAlign' , -1 | 0 | 1       , ...
%                      AxesProperty1 , AxesPropertyValue1  , ...  ) 
%
%  The imaginary part defines the Alignment of the TickLabels in the
%   Direction of Type.
% 
% [ LineHandle , TextHandle ] = LINEAXIS( ... ) returns the Handles
%  of the AxisLine and TextLabels
%
% Example:
%
%  figure;
%  
%  axe = axes('xlim',[0 10],'ylim',[0 10]);
%
%  lineaxis(axe,'Type','x','Location',NaN,'TickAlign', 0,'TextAlign',-1);
%  lineaxis(axe,'Type','y','Location',8.4,'TickAlign',-1,'TextAlign', 1-i);
%  lineaxis(axe,'Type','y','Location',6.8,'TickAlign',-1,'TextAlign',-1);
%  lineaxis(axe,'Type','y','Location',2.5,'TickAlign', 0,'TextAlign', 0+i);
%
%  % set( axe , 'visible' , 'off' )  % looks better 
%
             
Msg = '';
nl  = char(10);


cnf = struct( 'Type'       , { 'y' } , ...
              'Location'   , { NaN } , ...
              'TickAlign'  , {  0  } , ...
              'TextAlign'  , { -1  }       );


%---------------------------------------------------------------------------

VarIn = varargin;
VarIn = VarIn(:);

Nin = size(VarIn,1);

%---------------------------------------------------------------------------

if nargin < 1
  error('Not enough Input Arguments.')
end

%---------------------------------------------------------------------------
% Check Axes

ok = ( isnumeric(axe) & ( prod(size(axe)) == 1 ) );
if ok
 ok = ishandle(axe);
 if ok
    ok = strcmp(get(axe,'type'),'axes');
 end
end

if ~ok
  error('First Input must be a AxesHandle.');
end

%---------------------------------------------------------------------------
% Check VarIn

if ~isempty(VarIn)

  if ( mod(Nin,2) ~= 0 )
    error( 'Additional Inputs must contain Property-Value-Pairs.' );
  end
 
  Nin   = Nin / 2;
  VarIn = reshape( VarIn , 2 , Nin );

  if ~iscellstr(VarIn(1,:))
    error([ 'Additional Inputs must contain Property-Value-Pairs.' nl ...
           'Properties must be Strings.' ]);
  end 

  %---------------------------------------------------------------------------
  % Check for Configuration


  field = fieldnames(cnf);

  for ii = 1 : size(field,1)

    jj = find( strcmp( lower(VarIn(1,:)) , lower(field{ii}) ) );
 
    if ~isempty(jj)

      val = VarIn{2,jj(end)};

      VarIn(:,jj) = [];

      ok = ( prod(size(val)) == 1 );

      if ok

        switch field{ii}

         case 'Type'

           ok = ischar(val);
           if ok
              val = lower(val);
              ok = any(strcmp(val,{'x' 'y'}));
           end
 
         otherwise
 
           ok = ( isnumeric(val) & ( isfinite(val) | ...
                   ( isnan(val) & strcmp(field{ii},'Location') ) ) );
         end

      end

      if ok
        
         cnf = setfield(cnf,field{ii},val);

      else

         Msg = [ Msg nl(1:(end*(~isempty(Msg)))) ...
                 'Invalid Value for "'  field{ii} '".' ];
 
      end   

    end

  end

  if ~isempty(Msg)
     error(Msg)
  end

end


%---------------------------------------------------------------------------
% Check AxesProperties

if ~isempty(VarIn)

  try
    set( axe , VarIn{:} );
  catch
    error(['Invalid Additional Inputs.' nl lasterr ]);
  end

end
%---------------------------------------------------------------------------
% Prepare CenterAxes

 a = lower( cnf.Type );

xy = 'xy';

jj = sum(cumprod(double(~(a==xy)))) + 1;

 b = xy(3-jj);

xy_ind = [ jj  3-jj ];


%----------------------------------------------------------

alim = get( axe , [ a 'lim' ] );
blim = get( axe , [ b 'lim' ] );


if isnan( cnf.Location )
   cnf.Location = mean(blim);
end

if ( cnf.Location < blim(1) ) | ...
   ( cnf.Location > blim(2) )
   fprintf([ nl 'Warning: Location out of AxesLimits.' nl ]);
end

%----------------------------------------------------------

at = get( axe , [ a 'tick'      ] );
al = get( axe , [ a 'ticklabel' ] );

at = at(:)';

st = size(at,2);

%--------------------------------
% Fill up Labels
if ~isempty(al)

  if ischar(al)
    al = cellstr(al); 
  end

  al = al(:)';

  sl = size(al,2);
  if sl < st
    ind = ( 1 : st );
    ind = ind - sl * floor(ind/sl);
    ind = ind + sl * ( ind == 0 );      
    al  = al(ind);
 end

end

%--------------------------------
% Remove Values outside

jj = find( ( at < alim(1) ) | ...
           ( at > alim(2) )       );

at(jj) = [];
if ~isempty(al)
  al(jj) = [];
end


%--------------------------------
% Modify Axes

set( axe , [ a 'lim' ]      , alim     , ...
           [ b 'lim' ]      , blim     , ...
           [ a 'limmode' ]  , 'manual' , ...
           [ b 'limmode' ]  , 'manual' , ...
           'nextplot'       , 'add'          );

%%%        [ a 'tick' ]     , []     , ...
%%%        [ a 'tickmode' ] , 'manual' , ...

%--------------------------------
% Line

st = size(at,2);
 
 n = st + 1;

adata = NaN*ones(3,n);
bdata = NaN*ones(3,n);

adata([1 2],1)   = alim(:);         % Line
adata([1 2],2:n) = at([1 1],:);  % Ticks

tl = get( axe , 'ticklength' );
tl = tl(1) * diff(blim);

% TickAlign: -1    -1 .. 0
%             0    -1 .. 1
%             1     0 .. 1 

f = cat( 1 , -1*(cnf.TickAlign<=0) , 1*(cnf.TickAlign>=0) );

bdata([1 2],1)   = cnf.Location;
bdata([1 2],2:n) = cnf.Location + f(:,ones(1,st))*tl;

acolor = get( axe , [ a 'color' ] );

hl = line( 'parent'     , axe , ...
           [ a 'data' ] , adata(:) , ...
           [ b 'data' ] , bdata(:) , ...
           'marker'     , 'none'   , ...
           'linestyle'  , '-'      , ...
           'linewidth'  , get( axe , 'linewidth' ) , ...
           'color'      , acolor     );

%--------------------------------
% Text

if isempty(al)
   ht = [];
   return
end

 
fn = get( axe , 'fontname'    );
fu = get( axe , 'fontunits'   );
fs = get( axe , 'fontsize'    );
fw = get( axe , 'fontweight'  );


aalign = { 'top'       'right' 
           'bottom'    'left'   };

balign = { 'right'     'top'
           'center'    'middle'
           'left'      'bottom' };

talign = { 'vertical'  'horizontal' };


% TextAlign: -1    left          1
%             0    left  right   1  2
%             1          right      2

ip = [ -1  1 ];

ih = sign(real(cnf.TextAlign));
iv = sign(imag(cnf.TextAlign));

ht = NaN * ones( st , 1+(ih==0) );

ih = ( 1+(ih>0) : 2-(ih<0) );
iv = iv + 2;

for ii = 1 : st

  for jj = ih

      apos = adata( 1,ii+1);
      bpos = bdata(jj,ii+1) + ip(jj)*tl;

      pos = zeros(1,3);

      pos(xy_ind) = [ apos  bpos ];

      ht(ii,jj) = text( 'parent'     , axe    , ...
                        'position'   , pos    , ...
                        'string'     , al{ii} , ...
                        'fontname'   , fn     , ...
                        'fontunits'  , fu     , ...
                        'fontsize'   , fs     , ...
                        'fontweight' , fw     , ...
                        'color'      , acolor     , ...
       [  talign{xy_ind(1)} 'alignment' ] , aalign{jj,xy_ind(1)} , ...
       [  talign{xy_ind(2)} 'alignment' ] , balign{iv,xy_ind(1)}        );    

  end

end