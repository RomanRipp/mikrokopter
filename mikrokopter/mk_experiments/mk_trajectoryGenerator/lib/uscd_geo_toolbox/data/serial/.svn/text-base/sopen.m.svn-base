function [dev,msg] = sopen(port,baud,dps,flow,varargin)

% SOPEN  Creates and Opens Serial Port Object
%
% [Obj,Msg] = SOPEN( PortName , BaudRate , DPS , FlowControl )
%
%------------------------------------------------------------------------
%
% PortName    COM1 .. COM*
%                /dev/tty*
%               S0 ..   S*   expanded with /dev/tty
%             USB0 .. USB*   expanded with /dev/tty
%
% Serial Port Object will not open if PortName is preceded by "-"
% Serial Port Object will     open if PortName is preceded by "+" (default)
%
%------------------------------------------------------------------------
%
% Use following Property-Value-Pairs to define additional Properties
%  of a Serial Port Object:
%
%  SOPEN( ... , Property , Value , ... )
%
% SOPEN will use an existing Serial Port Object with the same PortName.
%  (and a same Tag, if that property is defined by following Inputs)
%  The existing Port will closed, it's parameter reset and reopened.
%
%------------------------------------------------------------------------
% Basic Serial Port Parameters:
%
% DPS = '[DataBits][Parity][StopBits]',  3-Character-String,  default: '8n1'
%
% DataBits: {8} | 7
% Parity:   {n} | e | o | m | s ==  [ {none} | odd | even | mark | space ]
% StopBits: {1} | 0
%
% FlowControl: {n} | h | s      ==  [ {none} | hardware | software ]
%
%------------------------------------------------------------------------
%
% Default Parameter
%
%               BaudRate: 9600
%               DataBits: 8
%                 Parity:   none
%               StopBits: 1
%            FlowControl: none
%
%             Terminator: CR/LF
%                TimeOut:    3
% BytesAvailableFcnCount: 16
%
%------------------------------------------------------------------------
%
% See also: serial, serial/fopen, serial/fclose, serial/get, serial/set
%


Nin  = nargin;
Nout = nargout;

dev = [];
msg = '';

%---------------------------------------------------------
% Default Settings

if Nin < 2, baud = []; end
if Nin < 3, dps =  ''; end
if Nin < 4, flow = []; end

vin = {      'Terminator' , 'CR/LF' , ...
             'Timeout'    , 3       , ...
 'BytesAvailableFcnCount' , 16        };

nop = 0;

%******************************************************************

msg = {};

if Nin < 1
   if isunix
      port = 'S0'; 
   else
      port = 'COM1';
   end
end

%---------------------------------------------------------
% Check Port

s = size(port); p = prod(s);

ok = ( ~isempty(port) & ischar(port) & ( p == s(2) ) );
if ok
   if any( port(1) == '+-' )
      ok = ( s(2) > 1 );
     nop = ( port(1) == '-' ); %  No_Open
      if ok
         port = port(2:p);
      end
   end
end

if ~ok
    msg = cat(1,msg,{'Port must be a non-empty string.'});
end

%---------------------------------------------------------
% Check BaudRate

if ~isempty(baud)

    s = size(baud); p = prod(s);

    if ( ischar(baud) & ( p == s(2) ) );
       varargin = cat(1,{baud;dps;flow},varargin(:));
       baud = []; dps = ''; flow = '';
    else
       ok = ( isnumeric(baud) & ( prod(size(baud)) == 1 ) );
       if ok
          ok = ( isfinite(baud) & ( baud > 0 ) );
       end
       if ~ok
           msg = cat(1,msg,{'BaudRate must be a single positive Numeric.'});
       end
    end

end

%---------------------------------------------------------
% Check DataBits / Parity / StopBits

if ~isempty(dps)
    ok = ( ischar(dps) & isequal(size(dps),[1 3]) );
    if ok
       dps = lower(dps);
       ok = ( any( dps(1) == '78' ) & ...
              any( dps(2) == 'noems' ) & ...
              any( dps(3) == '01' )           );
    end
    if ~ok
        msg = cat(1,msg,{'DPS must be a valid 3-Character-String.'});
    end
end

%---------------------------------------------------------
% Check FlowControl

if ~isempty(flow)
    ok = ischar(flow);
    if ok
       flow = lower(flow(1));
       ok = any( flow == 'nhs' );
    end
    if ~ok
        msg = cat(1,msg,{'FlowControl must be a valid String.'});
    end
end

%---------------------------------------------------------
% Property-Value-Pairs

if ~isempty(varargin)

    nv = prod(size(varargin));
    ok = ( mod(nv,2) == 0 );
    if ok
       ok = chkcstr(varargin(1:2:end));
    end
    if ~ok
        msg =  cat( 1 , msg , ...
                       { 'Following Inputs must be Property-Value-Pairs.'} , ...
                       { 'Properties must be Strings.'} );
    end

end

%---------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    msg = sprintf('Invalid Inputs.\n%s',msg);
    if Nout < 2, error(msg), end
    return
end

msg = '';

%******************************************************************

if isempty(baud), baud = 9600; end
if isempty(dps),  dps  = '8n1'; end
if isempty(flow), flow = 'none'; end
    
%******************************************************************
% Check UNIX-PortName

if any( port(1) == 'SU' ) & ~any( port == filesep )
   port = [ '/dev/tty' port ];
end

%---------------------------------------------------------
% Check for Property "Tag" in Inputs

opt = {};

ok = strcmp( upper(varargin(1:2:end)) , 'TAG' );
if any(ok)
   ok = 2 * max(find(ok));
   if chkstr(varargin{ok})
      opt = { 'Tag'  varargin{ok} };
   end
end

%---------------------------------------------------------
% Check for existing Object (with Tag)
%  Try lower and upper-case PortNames

for prt = { port  lower(port) upper(port) }
    dev = instrfind('Port',prt{1},opt{:});
    if ~isempty(dev)
        port = prt{1};
        break
    end
end

%---------------------------------------------------------
% Create Object

if isempty(dev)
   try
      dev = serial(port);
   catch
      msg = sprintf('Error create Serial Object: SERIAL( %s ).\n%s', ...
                     port,lasterr);
      if Nout < 2, error(msg), end
      return
   end
end

%---------------------------------------------------------
% Close existing Object

if strcmp(get(dev,'status'),'open')
   try
      fclose(dev);
   catch
      msg = sprintf('Error close Serial Object.\n%s',lasterr);
      if Nout < 2, error(msg), end
      return
   end
end

%---------------------------------------------------------
% Set Parameter

try
   set( dev , 'BaudRate' , baud , ...
              'DataBits' , eval(dps(1),'8') , ...
              'Parity'   ,      dps(2)      , ...
              'StopBits' , eval(dps(3),'1') , ...
           'FlowControl' , flow , ...
                  vin{:} , varargin{:} );
catch
   msg = sprintf('Error set Properties of Serial Object.\n%s',lasterr);
   if Nout < 2, error(msg), end
   return
end

%---------------------------------------------------------
% Open Port

if nop
   return
end

try
   fopen(dev);
catch
   msg = sprintf('Error open Serial Object.\n%s',lasterr);
   if Nout < 2, error(msg), end
   return
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
   str = cellstr(str);
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


