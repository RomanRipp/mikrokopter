function [Msg,hf,ht] = textlogo(fig,txt,pos,varargin);

% TEXTLOGO  Creates a LogoFrame, contains Text
%
% [ Msg, FrameHandle, TextHandle ] = ...
%      TEXTLOGO( FigureHandle , LogoText , Position )
%
%-------------------------------------------------------------------
% Inputs:
%
%   FigureHandle specifies the Figure, in which the Logo will created
%
%   LogoText     Text of the Logo, a CharArray or CellStringArray
%                 the Elements in Text will displayed as multiple Lines
%
%   Position     PixelPosition of the Logo in the Figure
%                  [ DX  DY  Width  Height ]
%                  DX > 0   Left;    DX < 0  Right
%                  DY > 0   Bottom;  DY < 0  Top
%                  Width, Height == NaN  ==>  automaticly detecting   
%
%  The FontSize of the Text will automaticly fitted into the Position.
%
%-------------------------------------------------------------------
% OutPuts:
%
%   Msg          contains ErrorMessages 
%   FrameHandle  Handle  of FrameUIControl
%   TextHandle   Handles of TextUIControls
%          
%  Note:  The     Units of the UIControls are set to 'pixels'
%         The FontUnits of the UIControls are set to 'points'
% 
%-------------------------------------------------------------------
% Additional Inputs:
%
%  [ ... ] = TEXTLOGO( ... , 'Offset'  , PixelOffset   , ...
%                            'FAspect' , FontAspect    , ...
%                            'HAspect' , HeightAspect  , ...
%                            'FName'   , FontName      , ...
%                            'FWeight' , FontWeight    , ...
%                            'FAngle'  , FAngle        , ...      );
%
%  PixelOffset   Offset of the Text in the Frame: [ Horiz Vert ]
%
%  FontAspect    AspectRatio in FontSize for LogoText, must have the 
%                 same Number of Elements like the LogoText
%
%  HeightAspect  Factor for the Height of the TextLine
%                 The default Height of the TextLine will determined
%                  by getting the Extension of an TextObject with the
%                  same FontProperties.
%
%  FontName      FontNames for LogoText, String or CellArray of Strings,
%                 corresponding to LogoText
%
%
%  FontWeight    FontWeight for LogoText, String or CellArray of Strings,
%                 corresponding to LogoText
%
%
%  FontAngle     FontAngle for LogoText, String or CellArray of Strings,
%                 corresponding to LogoText
%
%
% More Inputs can be UIControl-Properties-PropertyValue-Pairs, like:
%   BackGroundColor, ForeGroundColor etc.
%
%
%-------------------------------------------------------------------
% example:
%
%  % Create a Figure
%  fig = figure('position',[ 400 200 300 300 ], ...
%               'color'   , [ 1  1  1 ]); 
%
%  txt = { 'O'  'OceanData'  'ToolBox' };  % 3-Line Text
%
%  [Msg,hf1,ht1] = textlogo(fig,txt,[10 10 80 80]);
%
%  % Check the ErrorMessage
%  if ~isempty(Msg)
%     fprintf([ 'Error using TEXTLOGO: '  Msg ]);
%     return
%  end
% 
%  [Msg,hf2,ht2] = textlogo(fig,txt,[110 10 80 80], ...
%                            'Offset'  , [ 0  0 ] , ...
%                            'FAspect' , [ 3  1.5  1 ] ); 
%
%  [Msg,hf3,ht3] = textlogo(fig,txt,[210 10 80 80], ...
%                            'Offset'  , [ 0  0 ] , ...
%                            'FAspect' , [ 3    1    1 ]   , ...
%                            'HAspect' , [ 5/6  5/6  1 ]   , ...
%                            'backgroundcolor',[ 1 1   1 ] , ...
%                            'foregroundcolor',[ 0 0.3 1 ]        ); 
%
%  set(ht3(1),'fontweight','bold');
%
%  % Test automatic Detection of Position
%
%  txt = { 'E'  'ExperimentSetup'  'ToolBox' };
%
%  % Automatic Detection of Width
%
%  [Msg,hf4,ht4] = textlogo(fig,txt,[ 10 110 NaN 80 ], ...
%                            'Offset'  , [ 0  0 ] , ...
%                            'FAspect' , [ 3    1    1 ]   , ...
%                            'HAspect' , [ 5/6  5/6  1 ]  );
%
%  set(ht4(1),'fontweight','bold');
%
%
%  % No Width and High given, use DefaultFontSize for Minimum (FAspect)
%
%  [Msg,hf5,ht5] = textlogo(fig,txt,[-10 110 NaN NaN ], ...
%                            'Offset'  , [ 0  0 ] , ...
%                            'FAspect' , [ 3    1    1 ]   , ...
%                            'HAspect' , [ 5/6  5/6  1 ]  );
%
%  % No Width and High given, specify FontSize for Minimum (FAspect)
%
%  [Msg,hf6,ht6] = textlogo(fig,txt,[-10 -10 NaN NaN ], ...
%                            'Offset'  , [ 0  0 ] , ...
%                            'FAspect' , [ 3    1    1 ]   , ...
%                            'HAspect' , [ 5/6  5/6  1 ]   , ...
%                            'FontSize' , 12                     );
%
%

nl = char(10);

Msg  = '';
hf   = [];
ht   = [];

Msg0 = 'TEXTLOGO: ';

nm0 = size(Msg0,2);

nl0 = char([ 10 32*ones(1,nm0+0) ]);


FrameWidth  = 2;

FontUnits = 'points';


% Offset for correct Width of TextUIControl

spp = get(0,'screenpixelsperinch');

OffsetChar = 'H';

OffsetLen = [ 10  20  40        % TextWidth
               2   3   5  ];    % Lenght of OffsetString

OffsetLen(2,:) = OffsetLen(2,:) + ( spp > 80 );


%*********************************************************
% Check Inputs

Nin = nargin;

if Nin < 3
  Msg = [ Msg0  'Not enough InputArguments.' ];
  return
end

%-----------------------------------------------------------
% Figure

ok = ( isnumeric(fig) & ( prod(size(fig)) == 1 ) );
if ok
   ok = ishandle(fig);
   if ok
      ok = strcmp(get(fig,'type'),'figure');
   end
end

if ~ok
 Msg = [  Msg  nl0(1:(end*(~isempty(Msg)))) ...
         'First Input must be a FigureHandle.' ];
end

%-----------------------------------------------------------
% LogoText

ok = ~isempty(txt);
if ok
   [ok,txt] = chkcstr(txt,0);
end

if ~ok
   Msg = [  Msg  nl0(1:(end*(~isempty(Msg)))) ...
          'LogoText must be a nonempty CharArray or CellStringArray.' ];
end

txt = txt(:);

Nt = size(txt,1);

%-----------------------------------------------------------
% Position

pos = pos(:)';
ok = ( isnumeric(pos)       &  ...
      any( size(pos,2) == [ 2 4 ] ) );
if ok
  pos = cat(2,pos,NaN*ones(1,4-size(pos,2)));
  ok  = all( isfinite(pos([1 2])) );
  ok1 = isnan(pos([3 4]));
  ok2 = ( isfinite(pos([3 4]))  &  ( pos([3 4]) > 0 ) );
  ok  = ( ok  &  all( ok1 | ok2 ) );
end 

if ~ok
 Msg = [  Msg  nl0(1:(end*(~isempty(Msg)))) ...
         'Position must be a numeric Vector: ' ...
         '[DX DY Width Height]. ' ...
         nl0 'Width and Height must larger ZERO or NaN.' ];
end

  
%-----------------------------------------------------------
% Other

[m,cnf,prop] = checkin(Nt,varargin);

% cnf: FAspect, HAspect, Offset, FName, FWeight, FAngle
 
if ~isempty(m)
    frm = sprintf('%s%s','%s',nl0);
    Msg = [ Msg  nl0(1:(end*(~isempty(Msg)))) ...
       sprintf(frm,m{:}) ];
end

if ~isempty(Msg)
  Msg = [ Msg0   Msg ];
  return
end 

%*********************************************************

% Check Widht and Height

if ~all( ( pos([3 4]) > 2*cnf.Offset ) | isnan(pos([3 4])) );
  Msg = [  Msg0 ...
          'Position too small, Width and Height must be larger then '  ...
           sprintf('[ %.0f  %.0f ]', 2*cnf.Offset) ];
  return
end
  

pos_nan = find(isnan(pos));
pos_neg = find( pos < 0 );


cnf.Offset = cnf.Offset + FrameWidth;

cnf.FAspect = cnf.FAspect / min(cnf.FAspect);


  % Space for Text

  Space    = pos([3 4]) - 2*cnf.Offset;

  Space(pos_nan-2) = inf;   % !!!


  % FigurePosition in Pixels

  figuni = get(fig,'units');
           set(fig,'units','pixels');
  figpos = get(fig,'position');
           set(fig,'units',figuni);


%****************************************************************
%  FramePosition

  pos1          = pos;

  % Width and Height if NaN

  pos1(pos_nan) = 2*cnf.Offset(pos_nan-2) + 10;


  % Left and Bottom if Negative

  pos1(pos_neg) = figpos(pos_neg+2) + pos(pos_neg) - pos1(pos_neg+2) + 1;


%****************************************************************
% Create Frame

  try

   hf = uicontrol( prop{:} ,  ...
                   'parent'   , fig     , ...
                   'style'    , 'frame' , ...
                   'tag'      , 'LOGO_FRAME'     );

  catch

    Msg = lasterr;

  end


  if ~isempty(Msg)
    Msg = [ Msg0  'Error create Frame.' nl Msg ];
    return
  end

  set( hf , 'units'     , 'pixels' , ...
            'position'  ,  pos1    , ...
            'fontunits' , FontUnits       );

  fn = get(hf,'FontName');
  fw = get(hf,'FontWeight');
  fa = get(hf,'FontAngle');

%****************************************************************
%  Fit Text into Frame, determine FontSize

  % DummyAxes 
  axe = axes('parent' , fig , ...
             'units'    , 'pixels'  , ...
             'position' , pos1      , ...
             'visible'  , 'off'            );

  % DummyText
  ht  = text('parent' , axe , ...
          'units' , 'pixels' , ...
          'position',[ 1 1 0 ] , ...
          'string'  , '' , ... 
          'interpreter' , 'none' , ...
          'fontunits'  , get(hf,'FontUnits') , ... 
          'fontsize'   , get(hf,'FontSize')  , ... 
          'fontname'   , fn     , ... 
          'fontangle'  , fa     , ... 
          'fontweight' , fw     , ...
          'visible'    , 'off'  , ...
          'interpreter' , 'none'     );
   
  %-------------------------------
  % OffsetString

   c = polyfit( OffsetLen(1,:) , OffsetLen(2,:) , 2 );
 
  spc  = [ 0  0 ];
  ext  = zeros(Nt,4);

  ini = { 'FName'   'FontName'    fn
          'FWeight' 'FontWeight'  fw 
          'FAngle'  'FontAngle'   fa  };

  is_nan = all(isnan(pos([3 4])));

  if is_nan
     fs = get(hf,'fontsize'); % Use DefaultFontSize for Minimum FontSize
  else
     fs = 5;
  end

  zz = 0;  % Counter

  while ( all( spc < Space )  |  ( zz == 0 ) ) & ( fs < 100 )

     spc3 = spc;
     ext3 = ext;

     if ~is_nan
         fs = fs + 1;
     end

     for ii = 1 : Nt

         for jj = 1 : size(ini,1)
             v = getfield(cnf,ini{jj,1});
             if ~isempty(v{ii} )
                 set( ht , ini{jj,2} , v{ii} );
             end
         end

         set(ht,'FontSize' , round( fs*cnf.FAspect(ii) ) );

         set(ht,'String',OffsetChar);
         ext0 = get(ht,'extent');

         set(ht,'String',[ txt{ii} 'gjqp' ]);
         ext1 = get(ht,'extent');

         fa = get(ht,'fontangle');
         bl = ~isempty(txt{ii});
         if bl
            bl = strcmp(txt{ii}(end),' ');
         end
         if strcmp(fa,'italic') & ~bl
            txt{ii} = cat(2,txt{ii},' ');
         end

         set(ht,'String',txt{ii});

         ext(ii,:) = get(ht,'extent');

         ext(ii,3) = ceil( ext(ii,3) + ext0(3)*polyval(c,size(txt{ii},2)) );
         ext(ii,4) = ceil( ext1(4) );

         p = permute(ini(:,[2 3]),[2 1]);

         set( ht , p{:} );

     end
     % ii

     spc = ceil( [ max(ext(:,3))  sum(ceil(ext(:,4).*cnf.HAspect)) ] );

     if is_nan
        break
     end

     if any( ext(:,4) < ext3(:,4)-1 )
        spc = spc3;
        ext = ext3;
        fs  = fs - ~any( spc > Space );
        break
     end

     zz = zz + 1;

  end
  % while

  delete(ht)
  delete(axe)

  if ~is_nan & any( spc > Space )
      spc = spc3;
      ext = ext3;
       fs = fs - 1;
  end

%*****************************************************
% Set FramePosition

  % Width and Height if NaN

  if ~isempty(pos_nan)

      Space(pos_nan-2) = spc(pos_nan-2);

      pos(pos_nan) =  Space(pos_nan-2) + 2*cnf.Offset(pos_nan-2);

  end

  % Left and Bottom if Negative

  if ~isempty(pos_neg)

      pos(pos_neg) = figpos(pos_neg+2) + pos(pos_neg) - pos(pos_neg+2) + 1;

  end

  set(hf,'units'    , 'pixels' , ...
         'position' ,  pos          );


%*****************************************************
% Create Text

 ht = zeros(Nt,1);

 hh = ceil( ext(:,4) .* cnf.HAspect );


 vpos = zeros(Nt,1);

 vpos(2:Nt) = hh(Nt:-1:2);

 vpos = cumsum(vpos);
 vpos = vpos(Nt:-1:1) + pos(2) + cnf.Offset(2);

 % Move Vertcial if HAspect > 1
 
 dh =  floor( ( hh - ext(:,4) ) / 2 ) .* ( hh > ext(:,4) );

 hh = min( hh , ext(:,4) );

 vpos = vpos + dh;

 for ii = 1 : Nt

     ht(ii) = uicontrol( prop{:}               , ...
                        'parent'     , fig        , ...
                        'style'      , 'text'  );

         for jj = 1 : size(ini,1)
             v = getfield(cnf,ini{jj,1});
             if ~isempty(v{ii} )
                 set( ht(ii) , ini{jj,2} , v{ii} );
             end
         end

     set( ht(ii) , 'units'      , 'pixels'   , ...
                   'position'   , [ pos(1)+cnf.Offset(1)  vpos(ii) Space(1) hh(ii) ] , ...
                   'FontUnits'  , get(hf,'FontUnits')     , ...
                   'FontSize'   , round( fs*cnf.FAspect(ii) ) , ...
                   'string'     ,  txt{ii}              );

 end

%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,cnf,prop] = checkin(n,vin)

msg = cell(0,1);

cnf = struct( 'FAspect' , {  ones(n,1) } , ...
              'HAspect' , {  ones(n,1) } , ...
              'Offset'  , { zeros(1,2) } , ...
              'FName'   , {  cell(n,1) } , ...
              'FWeight' , {  cell(n,1) } , ...
              'FAngle'  , {  cell(n,1) }        );

prop = {};


if isempty(vin)
   return
end

nin = prod(size(vin));

ok = ( mod(nin,2) == 0 );
if ok
   nin = nin/2;
   vin = reshape(vin,2,nin);
   ok = chkcstr(vin(1,:));
end

if ~ok
    msg = cat(1,msg,{'Additional Inputs must be Property-Value-Pairs.'}, ...
                    {'Properties must be Strings.'} );
    return
end

ok = zeros(1,nin);

fn = fieldnames(cnf);
fl = lower(fn);

for ii = 1 : nin

    p = lower(vin{1,ii});
    v = vin{2,ii};

    jj = strcmp( fl , p );

    ok(ii) = ~any(jj);

    if ~ok(ii)

        jj = find(jj);
 
        str = '';

       if isempty(v)

          str = 'Empty Value for Property "%s".';

       else
 
          if ~ischar(v)
              v = v(:);
          end

          m = size(v,1);

          switch p

          %--------------------------------------------
          case { 'faspect'  'haspect' }

              vok = isnumeric(v);
              if vok
                 vok = all( v > 0 );
              end

              if ~vok
                  str = 'Value for Property "%s" must be positive numeric.';
              elseif m == 1
                 v = v(ones(n,1));
              elseif ( m < n )
                 str = sprintf('Value for Property "%s" must have 1 or %.0f Elements.', ...
                               '%s',n);
              elseif ( m > n )
                 v = v(1:n);
              end

          %--------------------------------------------
          case   'offset'

              if m == 1
                 v = v([1 1]);
              else
                 v = v(1:2);
              end
              vok = isnumeric(v);
              if vok
                 vok = all( mod(v,1) == 0 );
              end
              if vok
                 v = v(:)';
              else
                 str = 'Value for Property "%s" must be Integer.';
              end
  
          %--------------------------------------------
          case { 'fweight'  'fname'  'fangle' }

              [vok,v] = chkcstr(v);
              if ~vok
                  str = 'Value for Property "%s" must be a StringArray.';
              else
                  v = v(:);
                  m = size(v,1);
                  if m == 1
                     v = v(ones(n,1));
                  elseif ( m < n )
                     str = sprintf('Value for Property "%s" must have 1 or %.0f Elements.', ...
                                    '%s',n);
                  elseif ( m > n )
                     v = v(1:n);
                  end
               end

          end

       end

       if isempty(str)
          cnf = setfield(cnf,fn{jj},v);
       else
          msg = cat(1,msg,{strrep(str,'%s',fn{jj})});
       end

    end

end

if any(ok)
   ok = find(ok);
   prop = vin(:,ok);
end

%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


