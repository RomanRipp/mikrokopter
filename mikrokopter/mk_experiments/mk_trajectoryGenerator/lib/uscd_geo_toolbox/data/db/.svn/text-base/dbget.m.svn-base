function varargout = dbget(varargin)

% DBGET  Get Configuration for DataBase-Access
%
%  V = DBGET( P )
%
% Returns the ParameterStructure V which corresponds
%  to the ParameterNames in the CellStringArray P.
%
%  V = DBGET  returns the full ParameterStructure.
%
% [ V1 , V2 , V3 , ... ] = DBGET( P1 , P2 , P3 , ... )
%
% Returns the Parameter V* corresponding to the 
%  ParameterNames P*.
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
% see also: DBSET, DBROOT
%

Nin  = nargin;
Nout = nargout;

Nout = max(1,Nout);

varargout = cell(1,Nout);

app = 'DBASE';

if isappdata(0,app)
   db = getappdata(0,app);
else    
   [msg,db] = dbset;
end


if ( Nin == 0 )
   varargout{1} = db;
   return
end

%--------------------------------------------------
% Check for single CellStringArray or Strings

one = ( Nin == 1 );
if one
   [one,v] = chkcstr(varargin{1},1);
end
if ~one
    [ok,v] = chkcstr(varargin);
    if ~ok
        error('Inputs must be Strings or a CellArray of Strings.')
    end
end

%--------------------------------------------------
% Multiple Outputs

fld = fieldnames(db);
lf  = lower(fld);

if ~one
    v = lower(v);
    for ii = 1 : min(Nin,Nout)
        jj = strcmp(lf,v{ii});
        if any(jj)
           jj = find(jj);
           varargout{ii} = getfield(db,fld{jj});
        end
    end
    return
end

%--------------------------------------------------
% Single Output

w = v(:)';

w = w([1 1],:);

w(2,:) = { {[]} };

n = size(w,2);

v = lower(v);

for ii = 1 : n
    jj = strcmp(lf,v{ii});
    if any(jj)
       jj = find(jj);
       w{2,ii} = { getfield(db,fld{jj}) };
    end
end

try
   w = struct(w{:});
catch
   error('Invalid Names for Parameter.');
end

varargout{1} = w;
