function [msg,cnf,ini] = dbset(varargin)

% DBSET  Set Configuration for DataBase-Access
%
% DBSET( Parameter , Value , ... )
%
% DBSET( ParameterStruct )
%
% DBSET without any Input displays the actual Configuration
%
% [Msg,Config,Init] = DBSET( ... ) returns the ErrorMessage,
%   the actual Configuration and the Initilisation.
%
% Use DBROOT to set the Parameter "Root"
%
% The Parameters:
%
%          Root: Root Directory to search for DataFiles
%     Separator: Seperator for coded Identifer
%     Extension: default Extensions for DataFiles
%        Spacer: Spacer in DataFileName before FileNumber
%        Digits: Number of Digits for FileNumber in FileName
%
%         Dummy: DummyValues in DataFile
%
%      MaxLines: Maximum Number of HeaderLines
%     MaxLength: Maximum Length of HeaderLine
%        Marker: Marker at Begin of a HeaderLine
%    Assignment: AssignmentCharacter for HeaderValue
%      Comments: Markers for Comments in HeaderLine 
%    ColumnKeys: Keywords for Variable-Columns
%
%       Verbose: Display Messages in MatlabCommand if "on"
%       LogFile: Name of LogFile
%       LogMode: Log Messages in LogFile if "on"
%         LogID: FileIdentifer of LogFile 
%
% See in DBFILE how a DataFileName is build using 
%  Root, Spacer, Digits and Extension.
% 
% The DefaultValues:
%
%          Root: {''}
%     Separator: ':'
%     Extension: {'dat'  'edt'  'asc'  'raw'  'cal'  '*'}
%        Spacer: '_'
%        Digits: [3 2 1]
%         Dummy: [-9999 -9.9999e+03 1.0000e+32]
%      MaxLines: 256
%     MaxLength: 512
%        Marker: [1x0 char]
%    Assignment: '='
%      Comments: {'%'}
%    ColumnKeys: {'columns','columns cont.'}
%       Verbose: 'on'
%       LogFile: [1x0 char]
%       LogMode: 'off'
%         LogID: -1
%
% see also: DBGET, DBROOT
%

Nin  = nargin;
Nout = nargout;

msg = cell(0,1);

%---------------------------------------
% Initialisation

[ini,cnf] = defaults;

def = cnf;

ini = ini(:,[1 3 4 5]);

%---------------------------------------
% Check for existing ApplicationData

app = 'DBASE';

ok = isappdata(0,app);
if ok
   db = getappdata(0,app);
   ok = isequal(fieldnames(db),fieldnames(cnf));
end

if ok
   [m1,m2,cnf,ok] = check_cnf(cnf,ini,db);
   jj = ( ok == 0 );
   if any(jj)
      jj = find(jj);
      for ii = jj(:)'
          cnf = setfield( cnf , ini{ii,1} , getfield(def,ini{ii,1}) );
      end
   end
else
   [m1,m2,cnf,ok] = check_cnf(cnf,ini);
   if ~isempty(m1)
       msg = cat( 1 , msg , { rmblank(m1,2) } );
   end
   if ~isempty(m2)
       msg = cat( 1 , msg , { rmblank(m2,2) } );
   end
   if ~isempty(msg)
       msg = sprintf('%s\n',msg{:});
       msg = sprintf('Invalid DefaultConfiguration.\n%s',msg);
       if Nout == 0
          error(msg)
       end
       return
   end
end

%------------------------------------
% Check actual Logging

if strcmp(cnf.LogMode,'off')

   % Check if LogFile is open
   if ~isempty(fopen(cnf.LogID))   
       try, fclose(cnf.LogID); end
       cnf.LogID = -1;
   end

elseif isempty(fopen(cnf.LogID))

   % Check if LogFile is not open
   if ~isempty(cnf.LogFile)
       cnf.LogID   = fopen(cnf.LogFile,'a');
       cnf.LogFile = fopen(cnf.LogID);
   end

end

%------------------------------------

setappdata(0,app,cnf);
    
%---------------------------------------
% Check InputConfiguration

if Nin > 0

   if isequal(varargin,{0})
      varargin = {def};
   end

   %------------------------------------
   % Store actual Configuration

    org = cnf;

   %------------------------------------
   % Set / Check InputParameter

   for ii = 1 : Nin
       if iscell(varargin{ii}) & ~( prod(size(varargin{ii})) == 1 )
          varargin{ii} = varargin(ii);
       end
   end

   [m1,m2,cnf] = check_cnf(cnf,ini,varargin{:});

   if ~isempty(m1)
       msg = cat( 1 , msg , { sprintf('Invalid Inputs.\n%s',rmblank(m1,2)) } );
   end
   if ~isempty(m2)
       msg = cat( 1 , msg , { sprintf('Invalid Parameter.\n%s',rmblank(m2,2)) } );
   end

   %------------------------------------

   if isempty(msg)
      msg = '';
   else
      msg = sprintf('%s\n',msg{:});
      if Nout == 0
         error(msg)
      end
      return
   end

   %------------------------------------
   % Check LogFile

   cnf.LogID = org.LogID;

   log_on  = ~isempty(fopen(cnf.LogID));
      file = cnf.LogFile;
   newfile = ~strcmp( file , org.LogFile );

   if log_on & ( strcmp(cnf.LogMode,'off') | newfile )
      try, fclose(cnf.LogID); end
      cnf.LogID = -1;
      log_on = 0;
   end

   if strcmp(cnf.LogMode,'on') & isempty(fopen(cnf.LogID)) & ~isempty(file)
      cnf.LogID   = fopen(cnf.LogFile,'a');
      cnf.LogFile = fopen(cnf.LogID);
      if isempty(cnf.LogFile)
         dbmsg(sprintf('Can''t open LogFile "%s".',file),1,cnf);
      end
   end

   %------------------------------------

   setappdata(0,app,cnf);

end

msg = '';

if Nout == 0
   clear msg
   disp(' ')
   disp(cnf)
end

%***********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ini,cnf] = defaults(app)

sets = { 'off'  'on' };
sep  =  ':;,';
val  = [ -9999.0  -9999.9  1e32 ];
ext  = { 'dat'  'edt'  'asc'  'raw'  'cal'  '*' };
spc  = '_';
dgt  = [ 3  2  1 ];

cls ={'columns','columns cont.'};

cm   = '%';

ini = { ...

'Root'       { {''}  }  'cellstr' [ NaN 1  ] {}    {'PathRoot'}
'Separator'  {  ':'  }  'char'    [  1  1  ] {}    {'Separator in Identifers'}
'Extension'  {  ext  }  'cellstr' [  1 NaN ] {}    {'Extensions for DataFiles'}
'Spacer'     {  spc  }  'char'    [  1 NaN ] {}    {'Spacer for Number in FileName'}
'Digits'     {  dgt  }  'numeric' [  1 NaN ] {}    {'Digits for FileNumber'}
'Dummy'      {  val  }  'numeric' [  1 NaN ] {}    {'DummyValues in DataFile'}
'MaxLines'   {  256  }  'uint16'  [  1  1  ] {}    {'Max. Number of HeaderLines'}
'MaxLength'  {  512  }  'uint16'  [  1  1  ] {}    {'Max. Length of HeaderLine'}
'Marker'     { ''    }  'char'    [  1 NaN ] {}    {'Marker for HeaderLines'}
'Assignment' { '='   }  'char'    [  1  1  ] {}    {'Character to assign HeaderValues'}
'Comments'   { {'%'} }  'cellstr' [  1 NaN ] {}    {'Characters to mark a Comment'}
'ColumnKeys' { cls   }  'cellstr' [  1 NaN ] {}    {'Keywords for Variable-Columns'}
'Verbose'    { 'on'  }  'none'    [  1 NaN ] sets  {'VerboseMode, display in MatlabCommand'}
'LogFile'    { ''    }  'char'    [  1 NaN ] {}    {'Name of LogFile'}
'LogMode'    { 'off' }  'char'    [  1 NaN ] sets  {'Mode of Logging'}
'LogID'      {  -1   }  'numeric' [  1  1  ] {}    {'ID of LogFile (readonly)'}

};

% 'Delimiter'  {  ':'  }  'char'    [ 1  1  ] {}    {'?Delimiter in Columns?'}

if nargout > 1

  cnf = permute( ini(:,[1 2]) ,  [ 2  1 ] );

  cnf = struct(cnf{:});

end
