function  [Msg,out] = axe_zoom(action,axe,varargin);

% AXE_ZOOM  Zoom for Axes
%
%------------------------------------------------------------------
% 1. Create new Zoom
%
% [ Msg , Lim ] = axe_zoom( 'New' , AxeHandle , UpFcn );
%
%  Lim = [ XLim  YLim ];
%
% The current AxesLimits [ XLim  YLim ] can''t be exeeded !!!
%
% UpFcn will evaluated after successfull Zoom-Action, using FEVAL:
%
%   Msg = FEVAL( UpFcn{:} , AxeHandle , Lim );
%
%  UpFcn is a CellArray, with a MatlabFunction as String in the 1. Element.
%
%  UpFcn = { MatlabFunction , [VarArgIn] }, see also below !!!
%
%------------------------------------------------------------------
% 2. Activate Zoom
%
% Msg = axe_zoom( 'On'  , AxeHandle );
%
%------------------------------------------------------------------
% 3. Disactivate Zoom
%
% Msg = axe_zoom( 'Off' , AxeHandle );
%
%------------------------------------------------------------------
% 4. Zoom In
%
% [ Msg , Lim ] = axe_zoom( 'In' , AxeHandle , Lim , CallUpFcn );
%
%  CallUpFcn =  0  |  1   default: 0  no call of UpFcn
%                          intern: 1     call of UpFcn
%
%------------------------------------------------------------------
% 5. Zoom Out
% 
% [ Msg , Lim ] = axe_zoom( 'Out' , AxeHandle , N , CallUpFcn ); 
%
% Zooms N-times out,  N == 0 zooms to start
%
%------------------------------------------------------------------
% Special Statements for UpFcn
%
% If the 1. Element (MatlabFunction) of UpFcn is a String 
%    which start and ends with "@", it will evaluate by using EVAL:
%
%    eval( UpFcn{1}( 2 : (end-1) ) ); 
%
%   In this case the AxeHandle has the VariableName "axe", 
%    which can be used in UpFcn{2} to assign the AxeHandle.
%
% If the 2. Element (1. Input for MatlabFunction) of UpFcn is a String 
%    which start and ends with "@", it will evaluate by using EVAL before: 
%
%    UpFcn{2} = eval( UpFcn{2}( 2 : (end-1) ) );  
%
%   In this case the AxeHandle has the VariableName "axe", 
%    which can be used in UpFcn{2} to assign the AxeHandle.
%   That gives more flexibility to define the first Input
%    of the MatlabFunction.
%
%

Msg = '';
out = [];

nl = char(10);

app = mfilename;
fsp = ( app == filesep );
if any(fsp)
    app = app(max(find(fsp))+1:end);
end
app = upper(app);


Msg0 = sprintf('%s: ',app);

Nin = nargin;

%-----------------------------------------------------------

if Nin < 1
   axe    = gca;
   action = 'New';
end

if Nin < 2
   if isnumeric(action)
      axe    = action;
      action = 'New';
   else
      axe = gca;
   end
end

%-----------------------------------------------------------
% Check Action

if ~chkstr(action,1)
    Msg = sprintf('%sAction must be a String.',Msg0);
    return
end

Msg0 = sprintf('%s(%s): ',app,action);
 
action = upper(action);

%------------------------------------------------------
% Check Handle

ok = ( isnumeric(axe) & ( prod(size(axe)) == 1 ) );
if ok
   ok = ishandle(axe);
   if ok
      typ = get(axe,'type');
      ok = strcmp( typ , 'axes' );
      if ~ok
          ok = strcmp(typ,'figure');
          if ok
             axe = get(axe,'currentaxes');
             if isempty(axe)
                return
             end
          else
             [m,p] = recpar(axe,'axes');
             ok = isempty(m);
             if ok
                axe = p(1);
             end
          end
      end
      if ok &  ~strcmp(action,'NEW')
         ok = isappdata(axe,app);
      end
   end
end

if ~ok
    Msg = sprintf('%sInput H must be an %s-Handle or Children of it.',Msg0,app);
    return
end

%------------------------------------------------------

if ~strcmp(action,'NEW')
    ud = getappdata(axe,app);
end

%------------------------------------------------------

switch upper(action)

%******************************************
case 'NEW'

  evl  = zeros(1,2);
  narg = zeros(1,2);

  if Nin < 3

    UpFcn = {};
    
  else

    UpFcn = varargin{1};

    if ischar(UpFcn)
       UpFcn = cellstr(UpFcn);
    end

    ok = iscell(UpFcn);
    if ok & ~isempty(UpFcn)
       ok = chkstr(UpFcn{1},1);
    end

    if ~ok
       Msg = [  Msg0  'Value for UpFcn must be a CharacterArray' ...
                      ' or CellArray with a String in the 1. Element.' ];
       return
    end

    UpFcn = UpFcn(:)';

    % Check for EVAL of Elements

    for jj = 1 : 2
        if size(UpFcn,2) >= jj
           if chkstr(UpFcn{jj}) & ( size(UpFcn{jj},2) >= 3 )
              evl(jj) = strcmp(UpFcn{jj}([1 end]),'@@');
              if evl(jj)
                 UpFcn{jj} = UpFcn{jj}( 2 : (end-1) );
              end
           end
        end
    end
    
    if ~evl(1)
       try
          narg =  [ nargin(UpFcn{1}) nargout(UpFcn{1}) ];
       end
    end

    
  end

  %------------------------------------------------------------

  xl = get(axe,'xlim');
  yl = get(axe,'ylim');

  lim = [ xl  yl ];

  DownFcn = [ 'axe_zoom(''Down'',gcbo,1);' ];

  ud = struct( 'History' , { lim        } , ...
               'DownFcn' , { DownFcn    } , ...
               'UpFcn'   , { UpFcn      } , ...
               'UpEval'  , { evl        } , ...
               'UpNarg'  , { narg       } , ...
               'Pointer' , { pointer('zoom')    }      );

  setappdata(axe,app,ud);

  set( axe , 'xlim' , xl , ...
             'ylim' , yl , ...
    'ButtonDownFcn' , DownFcn   );


  fig = get(axe,'parent');

  if ~strcmp( get(fig,'pointer') , 'custom')
      set(fig,pointer('default'));
  end


  out = lim;

%******************************************
case 'ON'

  set( axe , 'ButtonDownFcn' , ud.DownFcn );

%******************************************
case 'OFF'

  set( axe , 'ButtonDownFcn' , '' );

%******************************************
case 'DOWN'
 
  if Nin < 3
     clb = 0;
  else
     clb = varargin{1};
  end

  %---------------------------------------------
  % View X-Y only

  v = get(axe,'view');
  if ~( mod(v(2)/90,2) == 1 )
     return
  end

  fig = get(axe,'parent');

  %---------------------------------------------

  sel = get(fig,'selectiontype');
  sel = sel(1:3);
  
  %---------------------------------------------
  % Check for inside Axe

  cp = get( axe , 'currentpoint' );

  cp = cp( 1 , [ 1  2 ] );

  op = [ -1   1 ];

  if ~all( ( cp(1)*op <= op.*get(axe,'xlim') )  &  ...
           ( cp(2)*op <= op.*get(axe,'ylim') )  )  & ...
      strcmp(get(axe,'visible'),'on') & strcmp(sel,'nor')
      return
  end

  %---------------------------------------------

  is_cont = ~isempty(get(axe,'uicontextmenu'));

  switch( sel )

   %-----------------------------------  
   case { 'ope'  'alt'  'ext' }

     if ( size(ud.History,1) == 1 ) | ...
        (  is_cont  & strcmp(sel,'alt') )  |  ...
        ( ~is_cont  & strcmp(sel,'ext') )
 
        return

     end

     act = 'Out';
     in  = any(strcmp(sel,{'alt' 'ext'}));  % Zoom-Out-Number
 
   %-----------------------------------  
   case 'nor'

    c = ud.Pointer;
    p = c;

    f = fieldnames(p);
    for ff = f(:)'
        p = setfield( p , ff{1} , get(fig,ff{1}) );
    end

    set(fig,c);

    rbbox;

    cp1 = get(axe,'currentpoint');

    set(fig,p);

    act = 'In';

    in  = cat( 2 , cp(1) , cp1(1,1) , ...
                   cp(2) , cp1(1,2) );

    if any( abs(in([1 3])-in([3 4])) < 1e3*eps )
       return 
    end

   %-----------------------------------  
   otherwise
 
     return

  end

  %------------------------------------------------

  [Msg,out] = axe_zoom(act,axe,in,clb);

  if ~isempty(Msg)
      Msg = sprintf('%sError call %s.\n%s',Msg0,app,Msg);
      return
  end

 
%******************************************
case 'IN'

  if Nin < 3
      out = ud.History(end,:);
      return
  end

  if Nin < 4
     clb = 0;
  else
     clb = varargin{2};
  end

  lim = varargin{1};

  ok = ( isnumeric(lim)  &  (  prod(size(lim)) >= 4 ) );
  if ok
     lim = lim(1:4);
     lim = lim(:)';
     ok = all(isfinite(lim));
  end

  if ~ok
     Msg = [ Msg0  'Input must be 4-Element RowVector: [ XLim  YLim ].' ];
     return
  end

  lim0 = ud.History(  1  , : );
  lim1 = ud.History( end , : );

  op  = [ 1  -1   1  -1 ];  % Operator

  il  = ( 1 : 4 );

  %------------------------------------
  % Check for Min <=  Max

  ok = ( lim(il).*op <= lim(il+op).*op );

  lim = lim(il) .* ok + lim(il+op) .* (~ok);

  %------------------------------------
  % Check for Min ~= Max

  ok  = ~( lim(il) == lim(il+op) );

  lim = lim .* ok + lim1 .* (~ok);

  %------------------------------------
  % Check for Limit lim0

  ok  = ( ( lim0.*op - lim.*op ) < 0 );

  lim = lim .* ok  + lim0 .* (~ok);

  %------------------------------------

  ud.History = cat( 1 , ud.History , lim );

  set( axe , 'xlim' , lim([1 2]) , ...
             'ylim' , lim([3 4])       );

  setappdata(axe,app,ud);

  out = lim;

  
  %------------------------------------------------
  % CallBack

  if ~isequal(clb,1)  |  isempty(ud.UpFcn)
    return
  end

  %------------------------------------------------

   fcn = cat( 2 , ud.UpFcn , { axe lim } );

   MsgF = evalfcn(axe,fcn,ud.UpEval,ud.UpNarg,'Up');

   if ~isempty(MsgF)
       fprintf(1,'%s %s\n',Msg0,MsgF);
   end

%******************************************
case 'OUT'

  n = size(ud.History,1);

  if n == 1
     return
  end

  if Nin < 3
     nr = 1;
  else
     nr = varargin{1};
  end

  if Nin < 4
     clb = 0;
  else
     clb = varargin{2};
  end

  ok = ( isnumeric(nr)  &  ( prod(size(nr)) == 1 ) );
  if ok
     ok = ( isfinite(nr)  &  ( nr >= 0 )  &  ( mod(nr,1) == 0 ) );
  end

  if ~ok
    Msg = [ Msg0  'Input must be a positive Integer.' ];
  end



  n = n - min( nr+(n-1)*(nr==0) , n-1 );


  ud.History = ud.History(1:n,:);

  lim = ud.History( n , : );

  set( axe , 'xlim' , lim([1 2]) , ...
             'ylim' , lim([3 4])  );

  setappdata(axe,app,ud);
     
  out = lim;

  %------------------------------------------------
  % CallBack

  if ~isequal(clb,1)  |  isempty(ud.UpFcn)
      return
  end

  %------------------------------------------------

   fcn = cat( 2 , ud.UpFcn , { axe lim } );

   MsgF = evalfcn(axe,fcn,ud.UpEval,ud.UpNarg,'Up');

   if ~isempty(MsgF)
       fprintf(1,'%s %s\n',Msg0,MsgF);
   end

end

%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = evalfcn(axe,fcn,evl,narg,name)

% EVALFCN  Evaluates Expression
%

msg = '';
   
if evl(1)
   try
       eval(fcn{1});
   catch
       msg = sprintf('Error Evaluate 1. User%sExpression: %s\n%s', ...
                      name,fcn{1},lasterr);
   end
   return
end

if evl(2)
   try
       fcn{2} = eval(fcn{2});
   catch
       msg = sprintf('Error Evaluate 2. User%sExpression: %s\n%s', ...
                      name,fcn{2},lasterr);
       return
   end
end

nf  = size(fcn,2);
nf  = min( nf , 1+narg(1) + (nf-narg(1)) * (narg(1)<0) );

fcn = fcn(1:nf);

try
  if narg(2) == 0
     feval(fcn{:});
  else
     msg = feval(fcn{:});
  end
catch
  msg = lasterr;
end

if chkstr(msg,1)
   msg = sprintf('Error call User%sFcn: %s\n%s',name,fcn{1},msg);
else
   msg = '';
end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Msg,c,t] = recpar(h,typ);

% RECPAR  returns ParentHistory of Handle
%
% [Msg,HandleHist,TypeHist] = RECPAR( Handle , StopType )
%
%  recurse trough the Parents unto ( ParentType == StopType )
%
%  StopType starts with UpperCase, returns History excluding
%     Handle with ( ParentType == StopType )
%
%  default: StopType = 'Root'  (recurse unto 'figure')
%
%  HandleHist(end) == Handle
%    TypeHist(end) == HandleType
%

Msg = '';
 c  = zeros(0,1);
 t  =  cell(0,1);

if nargin < 1
   Msg = 'Input Handle is missing.';
   return
end

if nargin < 2
   typ = 'Root';
end

%-----------------------------------------------

if isempty(h)
   return
end

ok = ( isnumeric(h) &  ( prod(size(h)) == 1 ) );
if ok
   ok = ishandle(h);
end

if ~ok
   Msg = 'First Input must be a Single Handle.';
   return
end

if ~( ischar(typ) & ~isempty(typ) & ...
      ( prod(size(typ)) == size(typ,2) ) )
   Msg = 'Type must be a String';
   return
end

%-----------------------------------------------

c = h;
t = { get(h,'type') };

z = 1;

t0 = lower(typ);

while ~( c(1) == 0 )  &  ( ~strcmp(t{1},t0) | ( z == 1 ) )

   z = z + 1;

   c = cat( 1 ,         get(c(1),'parent')  , c );

   t = cat( 1 , { lower(get(c(1),'type')) } , t );

end


if strcmp( t{1} , t0 )

  n = 1 + strcmp( typ(1) , upper(typ(1)) );

  c = c(n:z);
  t = t(n:z);

else

   Msg = [ 'Handle has no Parents with Type '''  typ '''.' ]; 

end

%-----------------------------------------------
   

%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = pointer(mode);

%------------------------------------
switch mode

%------------------------------------
case 'default'
%------------------------------------

c = [ ...
'ooo......#......'
'o##oo...###.....'
'o####oo..#......'
'.o#####oo...###.'
'.o#######oo.....'
'..o########oo...'
'..o##########o..'
'...o#####oooo...'
'...o#####o......'
'....o##oo#o.....'
'....o##o.o#o....'
'.....o#o..o#o...'
'.....o#o...o#o..'
'......o.....o#o.'
'.............o#o'
'..............oo'  ];

focus = [ 1  1 ];

%------------------------------------
case 'zoom'
%------------------------------------

c = [ ...
'....###.........'
'..##o.o##.......'
'.#.o...o.#......'
'.#o..#..o#......'
'#o...#...o#.....'
'#..#####..#.....'
'#o...#...o#.....'
'.#o..#..o#......'
'.#.o...o.#o.....'
'..##o.o####o....'
'....###.o###o...'
'.........o###o..'
'..........o###o.'
'...........o###.'
'............o#o.'
'................' ];

focus = [ 6  6 ];

%------------------------------------
end
%------------------------------------

c = double(c);

c(find(c==double('.'))) = NaN;
c(find(c==double('#'))) = 1;
c(find(c==double('o'))) = 2;


c = struct( 'Pointer'             , { 'custom' }  , ...
             'PointerShapeCdata'  , { c }         , ...
            'PointerShapeHotSpot' , { focus    }       );

