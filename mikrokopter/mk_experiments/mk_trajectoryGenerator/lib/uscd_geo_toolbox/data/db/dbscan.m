function [file,ini,msg] = dbscan(ident)

% DBSCAN  Check for DataBaseFiles and their Variables
%
%
%  [File,INI,Msg] = DBSCAN(IDENT)
%
% input  :      ident           : string containing filenames
%                                 either 'clear' or 'coded' 
%
% Coded Ident: [Base]:[Name]:[Type]:[Numbers]
%
% Check for Files like:
%
%  Root/[Base]/Name/Type/Name_???.###
%
% see also: DBFILE, DBCHECK
%

Nin  = nargin;
Nout = nargout;

if Nin < 1
   ident = '*:*:*';
end

%*********************************************************************
% Get Files

[file,msg] = dbfile(ident);

if ~isempty(msg) & ( Nout < 3 )
    dbmsg(msg,1);
end

if isempty(file) | ( Nout < 2 )
   ini = dbident;
   return
end

%-----------------------------------------------------------------
% get keyword numbers of wanted informations

[ini,ok] = dbcheck(file);

if ~any(ok)
    file = {};
end

if ~all(ok)
    file = file(find(ok));
end

