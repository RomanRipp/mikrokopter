function [Msg,h] =  objmove(varargin)

% OBJMOVE  moves a GraphicObject
%
% [Msg, Handle] =  OBJMOVE( 'New' , Handle , ... )
%
% following Inputs can be ObjectProperty-PropertyValue-Pairs
%
% For action 'new'  more Options are:
% 
%   OBJMOVE( ... , 'MotionFcn' , MotionFCN        , ...
%                      'UpFcn' ,     UpFCN        , ...
%                   'Children' , AssosiatedHandles      )
%
%  The MotionFCN and UpFCN will evaluate like:
%
%   Msg = FEVAL( FCN{:} , Handle , AxesCurrentPoint , FigureCurrentPoint );
%
%  FCN is a CellArray, with a MatlabFunction as String in the 1. Element.
%
%  FCN = { MatlabFunction , [VarArgIn] }
%
% If the 1. Element (MatlabFunction) of FCN is a String 
%    which start and ends with "@", it will evaluate by using EVAL:
%
%    eval( FCN{1}( 2 : (end-1) ) ); 
%
%   In this case the ObjectHandle has the VariableName "h", 
%    which can be used in FCN{2} to assign the Handle.
%
% If the 2. Element (1. Input for MatlabFunction) of FCN is a String 
%    which start and ends with "@", it will evaluate by using EVAL before: 
%
%    FCN{2} = eval( FCN{2}( 2 : (end-1) ) );  
%
%   In this case the ObjectHandle has the VariableName "h", 
%    which can be used in FCN{2} to assign the Handle.
%   That gives more flexibility to define the first Input
%    of the MatlabFunction.
%
% see also: TEXTMOVE, AXEMOVE
%

Msg = '';
 h  = [];

nl = char(10);


fcn = mfilename;
fsp = ( fcn == filesep );
if any(fsp)
   fcn = fcn(max(find(fsp))+1:end);
end
app = upper(fcn);

Msg0 = sprintf('%s: ',app);

%------------------------------------------------------
% Check for Handle

[h,Vin] = hndlvarg(varargin,2,'');

uic = 0;  % True for UIControl

ok = ~isempty(h);
if ok
   axe = get(h,'parent');
    ok = strcmp( get(axe,'type') , 'axes' );
   if ~ok
       ok = ( strcmp( get(axe,'type') , 'figure' ) & ...
              strcmp( get( h ,'type') , 'uicontrol' )     );
       uic = ok;
   end
end
if ~ok
    Msg = [ Msg0 'A valid ObjectHandle of an Axes is missing.'];
    return
end


%-----------------------------------------------------------
% Check Action

if isempty(Vin)

   action = 'NEW';

else

   action = Vin{1};

   if ~( ischar(action) &  ~isempty(action) & ...
         ( prod(size(action)) == size(action,2) )   );
       Msg = [ Msg0 'First and Second Inputs must be a Handle' nl ...
                    '  or a String for Action or Property.' ];
       return
   end

   action = upper(action);

   if any( strcmp( action , { 'NEW' 'DOWN' 'MOVE' 'UP' } ) )
      Vin = Vin(2:end);
   else
      action = 'NEW';
   end

end

Msg0 = sprintf('%s(%s): ',app,action);

%------------------------------------------------------
% Check Handle

if ~strcmp(action,'NEW')

    if ~isappdata(h,app)
       Msg = [ Msg0 'Invalid Handle.'];
       return
    end

    tud = getappdata(h,app);

    if isequal(fieldnames(tud),{'Parent'})
         h = tud.Parent;
        ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );
       if ok
          ok = ishandle(h);
          if ok
             ok = isappdata(h,app);
          end
       end
       if ~ok
           Msg = [ Msg0  'Invalid ParentHandle.' ];
           return
       end
       tud = getappdata(h,app);
    end

    hdl = tud.Children;

    if isempty(hdl)
       hdl = h;
    else
        ok = ishandle(hdl);
       if ~any(ok)
           hdl = h;
       else
           if ~all(ok)
               ok = find(ok);
              hdl = hdl(ok);
           end
           hdl = cat( 1 , h , hdl(:) );
       end
    end
 
end

%------------------------------------------------------

DownFcn = sprintf('%s(''Down'',gcbo);',fcn); 

%------------------------------------------------------

switch action

%******************************************
case 'NEW'

   %------------------------------------------------------
   % Check following Inputs

   ok = isempty(Vin);
   if ~ok
     nv = prod(size(Vin));
     ok = ( mod(nv,2) == 0 );
     if ok & ~isempty(Vin)
        Vin = reshape( Vin , 2 , nv/2 );
        ok = iscellstr( Vin(1,:) );
     end
   end

   if ~ok
       Msg = [ Msg0 'Following Inputs must be Property-Value-Pairs.' ...
                    ' Properties must be Strings.' ];
       return
   end

   %------------------------------------------------------

   tud = struct( 'Selection'      , { '' } , ...  % SelectionType
                 'ActPointer'     , {pointer('def')} , ...  % Pointer of Figure
                 'MovPointer'     , {pointer('mov')} , ...
                 'RotPointer'     , {pointer('rot')} , ...
                 'Axes'           , { [] } , ...  % AxesHandle
                 'StartPoint'     , { [] } , ...  % StartPoint   [pixel]
                 'AxesPoint'      , { [] } , ...  % CurrentPoint [data]
                 'EraseMode'      , { '' } , ...  % Original EraseMode
                 'Children'       , { [] } , ...  % Assosiated Handles
                 'MotionFcn'      , { {} } , ...
                 'MotionEval'     , { [0 0] } , ...  % True for Eval Elements
                 'MotionNarg'     , { [0 0] } , ...  % [ Nin Nout ]
                 'UpFcn'          , { {} } , ...
                 'UpEval'         , { [0 0] } , ...
                 'UpNarg'         , { [0 0] } , ...
     'WindowButtonMotionFcn'      , { '' } , ...
     'WindowButtonUpFcn'          , { '' }        );

  
   %------------------------------------------------------------
   % Search for MotionFcn, UpFcn and EditAble

   sets = { 'on'  'off' };

   if ~isempty(Vin)

     vok = ones(1,size(Vin,2));

     for ff = { 'MotionFcn'  'UpFcn'  'Children' }

       ii = find( strcmp( Vin(1,:) , ff{1} ) );

       if ~isempty(ii)

          vok(ii) = 0;

          val = Vin{2,ii(end)};

          switch ff{1}

            %-------------------------------------------------------------
            case 'Children'
            %-------------------------------------------------------------

             ok = ( strcmp(class(val),'double') | isempty(val) );
             if ok & ~isempty(val)
                val = val(:);
                ok  = all(ishandle(val));
                if ok
                   for hh = val(:)'
                       if uic
                          ok = strcmp(get(hh,'type'),'uicontrol');
                       else
                          ok = strcmp(get(get(hh,'parent'),'type'),'axes');
                       end
                       if ~ok
                           break
                       end
                   end
                end
             end

             if ~ok
                 if uic
                    Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
                           'Value for '  ff{1} ' must be UIControls.' ];
                 else
                    Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
                           'Value for '  ff{1} ' must be Children of Axes.' ];
                 end
             end

            %-------------------------------------------------------------
            case { 'MotionFcn'  'UpFcn' }
            %-------------------------------------------------------------

             if isempty(val)
                val = {};
             elseif ischar(val)
                val = cellstr(val);
             end

             ok = iscell(val);
             if ok & ~isempty(val)
                ok = chkstr(val{1},1);
             end
 
             if ~ok

                 Msg = [  Msg  nl(1:(end*(~isempty(Msg)))) ...
                         'Value for '  ff{1} ' must be a CharacterArray' ...
                         ' or CellArray with a String in the 1. Element.' ];

             else

                 val = val(:)';

                 % Check for EVAL of Elements

                 evl = getfield( tud , [ff{1}(1:(end-3)) 'Eval'] );

                 for jj = 1 : 2
                     if size(val,2) >= jj
                        if chkstr(val{jj}) & ( size(val{jj},2) >= 3 )
                           evl(jj) = strcmp(val{jj}([1 end]),'@@');
                           if evl(jj)
                              val{jj} = val{jj}( 2 : (end-1) );
                           end
                        end
                     end
                 end

                 tud = setfield( tud , [ff{1}(1:(end-3)) 'Eval'] , evl );
                 
                 if ~evl(1)
                    try
                       tud = setfield( tud , [ff{1}(1:(end-3)) 'Narg'] , ...
                                       [ nargin(val{1}) nargout(val{1}) ] );
                    end
                 end

             end

          %-------------------------------------------------------------
          end
          %-------------------------------------------------------------

          if ok
             tud = setfield( tud , ff{1} , val );
          end

       end

     end
     % ff

     if ~isempty(Msg)
         Msg = [ Msg0 Msg ];
         return
     end
  
     if ~all(vok)
         Vin = Vin(:,find(vok));
     end

   end

   %------------------------------------------------------------
   % Create / Set Handle

   try

      if ~isempty(Vin)
         set( h , Vin{:} );
      end

   catch

      Msg = lasterr;
 
   end

   if ~isempty(Msg)
       Msg = [ Msg0 'Invalid Property-Values-Pairs.' nl Msg ];
       return
   end
     
   %------------------------------------------------------------
   % Prepare Handle
   
   set( h , 'ButtonDownFcn' , DownFcn );

   setappdata(h,app,tud);

   if ~isempty(tud.Children)
       set( tud.Children , 'ButtonDownFcn' , DownFcn );
       cud = struct('Parent',{h});
       for hh = tud.Children(:)'
           setappdata(hh,app,cud);
       end
   end

%******************************************
case 'DOWN'

   cp = get(axe,'CurrentPoint');
   cp = cp(1,[1 2]);

   %***************************************
   % Check SelectionType

   if uic
      fig = axe;
   else
      fig = get( axe , 'Parent' );
   end

   sel = get( fig , 'SelectionType' );
   hc  = get(  h  , 'uicontextmenu' );

   if       strcmp(sel,'normal')
            select = 'move';
   elseif   strcmp(sel,'alt') & isempty(hc) & uic
            cb = get(h,'CallBack');
            ok = ( size(cb,2) < size(DownFcn,2) );
            if ~ok
                ok = isempty(findstr(DownFcn,cb));
            end
            if ok
               select = 'move';
            else
               return
            end
   elseif   strcmp(sel,'alt') & ~isempty(hc)
            if isappdata(hc,'CurrentPosition')
               cp    = get(axe,'CurrentPoint');
               cp    = cp(1,[2 1]);
               try
                  ext = mercator(axe,'info');
               catch
                  ext = [];
               end
               if ~isempty(ext)
                   cp(1) = mercator({cp(1)},0,ext);
               end
               setappdata(hc,'CurrentPosition',cp);
            end
            if  isappdata(hc,'CurrentHandle')
               setappdata(hc,'CurrentHandle',h);
            end
            if  isappdata(hc,'CurrentObject')
               setappdata(hc,'CurrentObject',h);
            end
            return
   else
            return
   end

   %***************************************

   tud.Axes           = axe;
   tud.StartPoint     = cp;
   tud.AxesPoint      = cp;
   tud.Selection      = select;

   if uic
      tud.EraseMode   = '';
   else
      tud.EraseMode   = get( h , 'EraseMode' );
   end



   ptr = tud.ActPointer;
   for ff = fieldnames(ptr)'
       ptr = setfield(ptr,ff{1},get(fig,ff{1}));
   end

   tud.ActPointer = ptr;

   tud.WindowButtonMotionFcn = get( fig , 'WindowButtonMotionFcn' );
   tud.WindowButtonUpFcn     = get( fig , 'WindowButtonUpFcn'     );

   form = sprintf( '%%.%.0fg',  ceil(abs(log(eps)/log(10)))+1 );
   hh   = sprintf(form,h);

   MoveFcn = sprintf('%s(''Move'',%s);%s',fcn,hh,tud.WindowButtonMotionFcn);
     UpFcn = sprintf('%s(''Up'',%s);%s',fcn,hh,tud.WindowButtonUpFcn);

   if ~uic

       set( hdl , 'EraseMode' , 'xor' , ...
                  'Clipping'  , 'off'       );

       ii = strcmp(get(hdl,'type'),'text');

       if any(ii)
          it = find(ii);
          set(hdl(ii),'units','data');
       end

   end

   setappdata(h,app,tud);

   set( fig , 'WindowButtonMotionFcn' , MoveFcn , ...
              'WindowButtonUpFcn'     ,   UpFcn        );
 
   set( fig , tud.MovPointer );

%******************************************
case 'MOVE'

   cp = get( tud.Axes , 'CurrentPoint' );

   cp = cp(1,[1 2]);

   dp = cp - tud.AxesPoint;
   
   tud.AxesPoint = cp;

   for hh = hdl(:)'
        typ = get(hh,'type');
       if any(strcmp(typ,{'text' 'light'}))
          set( hh , 'position' , get(hh,'position') + [ dp  0 ]  );
       elseif any(strcmp(typ,{'uicontrol'}))
          set( hh , 'position' , get(hh,'position') + [ dp  0  0 ]  );
       else           
          set( hh , 'xdata' , get(hh,'xdata') + dp(1) , ...
                    'ydata' , get(hh,'ydata') + dp(2)   );
       end
   end

   setappdata(h,app,tud);

   %-----------------------------------------------------
   % Try MotionFcn

   fcn = tud.MotionFcn;

   if isempty(fcn)
      return
   end

   fcn = cat( 2 , fcn , { h tud.AxesPoint tud.FigurePoint } );

   MsgF = evalfcn(h,fcn,tud.MotionEval,tud.MotionNarg,'Motion');

   if ~isempty(MsgF)
       fprintf(1,'%s %s\n',Msg0,MsgF);
   end
   
%******************************************
case 'UP'

   if uic 
      fig = tud.Axes;
   else
      fig = get( tud.Axes , 'parent' );
      set( hdl , 'EraseMode' , tud.EraseMode       )
   end

   set( fig , 'WindowButtonMotionFcn' , tud.WindowButtonMotionFcn , ...
              'WindowButtonUpFcn'     , tud.WindowButtonUpFcn        );

   set( fig , tud.ActPointer );


   %-----------------------------------------------------
   % Try UpFcn

   fcn = tud.UpFcn;

   if isempty(fcn)
      return
   end

   fcn = cat( 2 , fcn , { h tud.AxesPoint tud.FigurePoint } );

   MsgF = evalfcn(h,fcn,tud.UpEval,tud.UpNarg,'Up');

   if ~isempty(MsgF)
       fprintf(1,'%s %s\n',Msg0,MsgF);
   end

end 

%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = evalfcn(h,fcn,evl,narg,name)

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


%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [h,v] = hndlvarg(v,n,type,tag);

% HNDLVARG  Checks Input-varargin for Handle
%
%  [H,V] = HANDLVARG( V , n , Type , Tag )
%
%  Search in V{n(1)} .. V{n(2)} for Handle of Type and Tag
%
%  If no Handle found and Type =  'figure' | 'axes', 
%    GCF | GCA is returned 
%
%  see also: CHKHNDL
%

Nin = nargin;

h = [];

if Nin == 0
   v = cell(1,0);
   return
end


if ~iscell(v)
   error('Input must be a CellArray.');
end

nv = prod(size(v));

if Nin < 2
   n = [];
end


if Nin < 3
   type = 'figure';
end

chkin = { type };

if Nin == 4
   chkin = cat( 2 , chkin , {tag} );
end

%------------------------------------------------------------

if isempty(n)
   n = [ 1  nv ];
else
   n = n(:)';
   n = n( 1 : max(1,size(n,2)) );
   n = min( n , nv );
   if size(n,2) == 1
      n = [ 1  n ];
   end
end


%-------------------------------------------------------------

ok = 0;

for ii = n(1) : n(2)

   ok = chkhndl(v{ii},chkin{:});

   if ok
      break
   end

end

%-------------------------------------------------------------

if ok

  h = v{ii};

  v(ii) = [];   % Reduce v

  return

end

if ~any( strcmp( type , { 'figure' 'axes' } ) )
   return
end

h = get( 0 , 'currentfigure' );
if isempty(h)
   return
end

if strcmp( type , 'axes' )
   h = get( h , 'currentaxes' );
end

ok = chkhndl(h,chkin{:});

if ~ok
   h = [];
end



%***************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [ok,msg] = chkhndl(h,type,tag);

% CHKHNDL(H,Type,Tag)  Checks, if H is a Handle of specified Type and Tag
%
%  Tag   CharArray or CellStringArray which Tags, the Handle has to be.
%          The Wildcard '*' can be used at Begin and/or End.
%

Nin = nargin;

ok  = 0;
msg = [];


if Nin == 0
   return
end

ok = ( isnumeric(h) & ( prod(size(h)) == 1 ) );

if ~ok
   return
end

ok = ishandle(h);
if ~ok | ( Nin < 2 )
   return
end

%-------------------------------------------------------------------------
% Check with Type

is_type = ( ischar(type)  &  ...
            ( isempty(type)  |  ( prod(size(type)) == size(type,2) ) ) );

if is_type

  ok = isempty(type);
  if ~ok
      ok = strcmp( get(h,'type') , type  );
  end

  if ~ok | ( Nin < 3 )
     return
  end

elseif ( Nin == 3 )

  msg = 'Input Type must be a String.';
  ok  = 0;
  return

else
 
  tag = type;
  
end


%-------------------------------------------------------------------------
% Check Tag

if ischar(tag)
   tag = cellstr(tag);
end

is_tag = iscellstr(tag);
if is_tag
   try
     cat(2,tag{:});
   catch
     msg = 'Input Tag must be a CharArray or CellStringArray.';
     ok  = 0;
     return
   end
end   


%-------------------------------------------------------------------------
% Check with Tag

 t = get(h,'tag');

nt = size(t,2);

ok = 0;

for tt = tag(:)'

    tok = strcmp( tt{1} , t );

    if ~tok  &  ~isempty(tt{1})  &  ~isempty(tag)

       n2 = size(tt{1},2);

       wc = ( double(tt{1}) == double('*') );

       if wc(1) & ( n2 == 1 )

          tok = 1;

       elseif wc(1) & wc(n2)

          tok = ~isempty( findstr(t,tt{1}(2:n2-1)) );

       elseif wc(1)

          n   = min(nt,n2-1);
          tok = strcmp(tt{1}(2:n2),t(nt-n+1:nt));

       elseif wc(n2)

          n   = min(nt,n2-1);
          tok = strcmp(tt{1}(1:n2-1),t(1:n));
          
       end

    end

    ok = ( ok | tok );

end

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

switch mode

case 'def'

c = [ ...
'ooo.............'
'o##oo...........'
'o####oo.........'
'.o#####oo.......'
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
'..............oo'];

case 'mov'

c = [ ...
'ooo.......###...'
'o##oo.....##....'
'o####oo...#.#...'
'.o#####oo....#.#'
'.o#######oo...##'
'..o########oo###'
'..o##########o..'
'...o#####oooo...'
'...o#####o......'
'....o##oo#o.....'
'....o##o.o#o....'
'.....o#o..o#o...'
'.....o#o...o#o..'
'......o.....o#o.'
'.............o#o'
'..............oo'];

case 'rot'

c = [ ...
'ooo.........###.'
'o##oo......#...#'
'o####oo....#...#'
'.o#####oo..#...#'
'.o#######oo.###.'
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
'..............oo'];

end


c = double(c);

c(find(c==double('.'))) = NaN;
c(find(c==double('#'))) = 1;
c(find(c==double('o'))) = 2;


c = struct( 'Pointer'             , { 'custom' }  , ...
             'PointerShapeCdata'  , { c }         , ...
            'PointerShapeHotSpot' , { [ 1  1 ] }       );
