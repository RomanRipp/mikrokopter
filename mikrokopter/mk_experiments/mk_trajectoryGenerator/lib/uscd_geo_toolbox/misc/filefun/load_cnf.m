function [Msg,MsgC,cnf,cok,cfg] = load_cnf(file,ini,varargin);

% LOAD_CNF  Loads and Checks ConfigurationFile, returns Structure
%
% [Msg,MsgCNF,CNF,INI_OK,V] = LOAD_CNF( ConfigFile , INI , V )
%
% INI = { KeyWord  {DefaultValue}  [Type]  [Size]  [Optional] }
%
% Loads the KeyWords, defined in 1. Column of INI, 
%  from ConfigFile using LOAD_ASC.
%
%  The Elements of the 2. Column (DefaultValue) must be a CellArray 
%    of Length 1 or a empty CellArray. The Element of the nonempty
%    {DefaultValue} is used if  a Keyword is not found in the ConfigFile, 
%     or the Value is incorrect,
%
% If INI contains more then 2 Columns for { Type  Size Optional }
%  or a StructureInput V is given, CHECK_CNF is used to check the KeyWords.
%
%  Type = { 'none'  'numeric'  'string'  
%           'cellstr'  'cellvec'  'object'   <ClassName> }
%
%  Optional =  CellArray of Values, which can be equal to the KeyWordValue
%
%  Specials for empty Type or 'none':
%
%        empty Optional: any KeyWordValue matching Size ok
%     nonempty Optional: only OptValues ok
%
%  Specials for empty Size or NaN's:
%
%            empty Size:  any Size Ok
%         NaN's in Size:  all Lengths in this Dimensions Ok
%
%  If an StructureInput V is given, the Value of CNF will replaced by the
%    Values of V before check.
%
%  V can be a Structure, a CellArray, defines a Structure using struct(V{:}),
%    or multiple FieldName-Value-Inputs.  
%
% INI_OK : -1  FieldName of INI doesn''t exist
%           0  Invalid Value
%           1  Value ok
%
%--------------------------------------------------------------------------
%
%  Use following Syntax to define the KeywordMarker <KM> and CommentMarker <CM>
%   which are used in the ConfigFile:
%
% [Msg,MsgCNF,CNF] = LOAD_CNF( ConfigFile , INI , '-k <KM>', '-c <CM>' , V )
%
% KeywordStart and KeywordEnd can be separated by Blanks (" ") in KM
%
%   default: <KM> ==  '# :'
%            <CM> ==  '%'
%
%--------------------------------------------------------------------------
%
%  see also:  LOAD_ASC, CHECK_CNF
%


Msg  = '';
MsgC = '';

Msg0 = 'LOAD_CNF: ';

cnf = [];
cfg = [];
cok = [];

Nin = nargin - 2;

if Nin < 0
   Msg = sprintf('%sNot enough Input Arguments.',Msg0);
   return
end

%-----------------------------------------------

Msg = cell(0,1);

if ~chkstr(file,1)
    Msg = cat( 1 , Msg , {'Input File must be a String.'} );
end

ok = ( iscell(ini)  &  ~isempty(ini)  &  ( ndims(ini) == 2 ) );

if ~ok

   Msg = cat( 1 ,  Msg , {'INI must be a CellArray: { KeyWord {DefaultValue} [Type] [Size] [Optional] }.'} );

else

   ok = chkcstr(ini(:,1));
   if ok
      ok = ~any(strcmp(ini(:,1),''));
   end
   if ~ok
       Msg = cat( 1 ,  Msg , {'Elements of 1. Column of INI must be nonempty Strings (KeyWord).'} );
   end

   if size(ini,2) > 2
      for ii = 1 : size(ini,1)
          if ~( ( iscell(ini{ii,2}) & ( prod(size(ini{ii,2})) == 1 ) ) | isempty(ini{ii,2}) )
              Msg = cat( 1 ,  Msg , {'Elements of 2. Column of INI must be empty or a single element CellArray.'} );
              break
          end
      end
   end

end

%--------------------------------------------------------------------

if ~isempty(Msg)
    Msg = sprintf('%s\n',Msg{:});
    Msg = sprintf('%sInvalid Inputs.\n%s',Msg0,Msg);
    cok = NaN * zeros(size(ini,1),1);
    return
end

Msg = '';

%--------------------------------------------------------------------

if size(ini,2) < 2
   ini      = ini(:,[1 1]);
   ini(:,2) = { {} };       % DefaultValue
end

%**************************************************
% Check varargin for Keyword & Comment - Marker

km = {'#' ':'};  % KeywordMarker
cm = '%';        % CommentMarker

ok = zeros(1,2);

for ii = 1 : min(2,Nin)

  val = varargin{ii};
   nv = prod(size(val));

  ok(ii) = ( ischar(val) & ( nv >= 3 ) & ( nv == size(val,2) ) );
 
  if ok(ii)
     ok(ii) = any( strcmp( val(1:2) , { '-k'  '-c' } ) );
     if ok(ii)
        key = val(1:2);
        val = separate(val(3:nv),' ');
        ok(ii) = ( ~isempty(val) & ( prod(size(val)) <= 1+strcmp(key,'-k') ) );
     end
     if ok(ii)
        switch  key 
           case '-k'
                  km = val;
           case '-c'
                  cm = val{1};
        end
     end
  end

end

if any(ok)
   varargin(find(ok)) = [];
end

%**************************************************
% Load KeyWords

try
  [Msg,key] = load_asc(file,ini(:,1),km,cm);
catch
   Msg = lasterr;
end


if ~isempty(Msg)
   Msg = sprintf('%sError call LOAD_ASC(%s).\n%s',Msg0,file,Msg);
   return
end


%**************************************************
% Check 

ii = strcmp( key(:,3) , 'none' );
if any(ii)
   ww = warnstat;
   if ~strcmp(ww,'off')
       warning('on');
       ii = find(ii);
       wrn = sprintf('%s ',key{ii,1});
       wrn = sprintf('Error evaluate Keywords in "%s":\n%s%s', ...
                      file,char(32*ones(1,9)),wrn);
       fprintf(1,'\n%s',char(7));
       warning(wrn);
       fprintf(1,'\n%s',char(7));
       warning(ww);
    end
end

nk = size(key,1);

MsgC    = cell(nk,1);
MsgC(:) = {''};

cok     = ones(nk,1);

for ii = 1 : nk

    m = '';

    if isempty( key{ii,3} );
       m = sprintf('Keyword "%s" not found.',key{ii,1});
       cok(ii) = -1;
    elseif strcmp( key{ii,3} , 'none' )
       m = sprintf('Invalid Syntax in Keyword "%s".',key{ii,1});
       cok(ii) = 0;
    end

    if isempty(m)
       ini{ii,2} = key(ii,2);
    elseif ~isempty(ini{ii,2});
       % Use DefaultValue
       m = '';
    end

      ok(ii) = isempty(m);
    MsgC{ii} = m;

end
% ii

ok = ( abs(cok) == 1 );

if any(~ok)
   MsgC = sprintf('%s\n',MsgC{find(~ok)});
   return
end

MsgC = '';

%--------------------------------------------------
% Build Structure

cnf = ini(:,[1 2])';
cnf = struct( cnf{:} );

if ( size(ini,2) < 3 ) & isempty(varargin)
   return
end

%--------------------------------------------------
% Check Structure

try
  [Msg,MsgC,cnf,cok,cfg] = check_cnf( cnf , ini(:,[1 (3:size(ini,2))]) , varargin{:} );
catch
   Msg = lasterr;
end

if ~isempty(Msg)
    Msg = sprintf('%sError call CHECK_CNF.\n%s',Msg0,Msg);
end


%****************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = separate(str,sep);

c = cell(1,0);

if isempty(str)
   return
end

n  = size(str,2);

is = ( double(str) == double(sep) );

if all(is)
   return
end

%---------------------------------------------

i0 = ~is(1);
i1 = ~is(n);

is = cat( 2 , ones(1,i0) , is , ones(1,i1) );

is = find( is );

is = is(:);

ni = size(is,1) - 1;

if ni == 0 
   return
end
   
%---------------------------------------------
% [ Start  End ]

ind = ( 1 : ni ) ;

is  = cat( 2 , is(ind)+1 , is(ind+1)-1 ) - i0;

%---------------------------------------------
% Take care for duplicate Seperators

if ( double(sep) == char(32) )
   is( find( is(:,1) > is(:,2) ) , : ) = [];
end

%---------------------------------------------

ni = size(is,1);

if ni == 0
   return
end

%---------------------------------------------

c = cell(1,ni);

for ii = 1 : ni
    c{ii} = str( is(ii,1) : is(ii,2) );
end


%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ww = warnstat

% WARNSTAT  Returns global WarningStatus
%
%  WARNSTAT returns the Status of WARNING
%
% Matlab R<13   WARNING
% Matlab R>12   WARNING for Identifier ALL
%

ww = warning;

if isstruct(ww)   % New Matlab R>12 Syntax
   try
      id = strcmp({ww.identifier},'all');
      if any(id)
         id = find(id);
         ww = ww(id(1)).state;
      else
         ww = '';
      end
   catch
      ww = '';
   end
elseif ~chkstr(ww)
   ww = '';
end

%*********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


%****************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
