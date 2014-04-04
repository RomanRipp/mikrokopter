function f = save(s,f,varargin)

% SAVE  Save structure fields
%
% SAVE( Structure , FileName , [Options] )
%
% Save the fields of the scalar Structure ([1x1] struct)
%  as individual variables within the file FileName.
%
% Corresponding Syntax to GENERAL/SAVE;
%
% GENERAL: SAVE(FileName,'-STRUCT','S') % Input is Name of Variable
%
%  STRUCT: SAVE(S,FileName)             % Input is Variable
%      or  SAVE(FileName,S)
%
%---------------------------------------------------------
%
% FileName and Options must be Strings, 
%  Options start with the Character '-' 
%
% see under GENERAL/SAVE for valid Options
%
%---------------------------------------------------------
%
% A missing FileName:  SAVE( Structure , [Options] )
%  opens UIGETFILE to select a FileName
%
% FileName = SAVE( Structure , ... ) retuns the FileName
%
%---------------------------------------------------------
%
% See also: GENERAL/SAVE
%
%---------------------------------------------------------
%
% Following Fieldnames conflicts with this function:
% 
% savethestructure, structurematlabfilename, 
% structureinitialisation, structuresaveoptions,
% structurefieldnumbers, structureloopindex,
%
% function, eval, save, break, return,
% for, while, if, else, elseif, end,
% try, catch, switch, case, otherwise
%
% and other language constructs
%

Nin  = nargin;
Nout = nargout;

if Nin < 2
   f = '';
end

if ~isstruct(s)
    n = s;
    s = f;
    f = n;
end

if ~( isstruct(s) & ( prod(size(s)) == 1 ) )
    error('Skalar Structure Variable required.');
end

v = varargin;

if ~isempty(f)
    if ~( ischar(f) & ( prod(size(f)) == size(f,2) ) )
        error('FileName must be a String.');
    end
    if f(1) == '-'
       v = cat( 1 , {f} , v(:) );
       f = '';
    end
end

if ~isempty(v)
    if ~chkcstr(v)
        error('Options must be Strings.');
    end
end

if isempty(f)
   [f,p] = uigetfile('*.mat','Save Structure as ...');
   if ~ischar(f) | isempty(f)
       f = '';
       if Nout == 0, clear f, end
       return
   end
   f = fullfile(p,f);
end

s = cat( 2 , fieldnames(s) , struct2cell(s) );

try
   savethestructure(f,s,size(s,1),v(:));
catch
   error(sprintf('Error save Fields.\n%s',lasterr));
end

if Nout == 0, clear f, end

%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function savethestructure( structurematlabfilename , structureinitialisation , ...
                           structurefieldnumbers , structuresaveoptions )

for structureloopindex = 1 : structurefieldnumbers

    eval([  structureinitialisation{structureloopindex,1} ' = ' ...
           'structureinitialisation{structureloopindex,2};'         ]);

end

save( structurematlabfilename , structureinitialisation{:,1} , ...
      structuresaveoptions{:} );

%********************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

