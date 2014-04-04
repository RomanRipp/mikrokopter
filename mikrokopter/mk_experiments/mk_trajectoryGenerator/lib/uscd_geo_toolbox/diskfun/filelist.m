function [msg,c] = filelist(varargin)

% FILELIST  creates FileList
%
% [Msg , C ] = FILELIST(Pfad,'-r','.ext1','.ext2', ...)
%
%  C    = { Dir  Name  Bytes Date  DateNum };
%    
%  note:  give '.'  or no Extension to list all Extensions
%
% see also: DIRCONT, DIRINFO, DIRHIST
%

 
%--------------------------------------------------------

Nout = nargout;

msg = '';
  c = cell(0,5);

nl  = char(10);

fs  = filesep;

%****************************************************************************
% Get and Check Inputs

try
   [msg,Pfad,ext,recurse] = checkin(varargin{:});
catch
   msg = cat(2,'Error call FILELIST:CHECKIN',nl,lasterr);
end


%****************************************************************************
% Create FileList

if isempty(msg)

   try   
       c = getlist(Pfad,ext,recurse);
   catch
     msg = cat(2,'Error call FILELIST:GETLIST',nl,lasterr);
   end

end

%****************************************************************************
% Return

if ~isempty(msg) & ( Nout == 0 )

  fprintf(nl)
  fprintf(strrep(strrep(msg,'\','\\'),'%','%%'))
  fprintf([ nl  nl ]);

  msg = '';

  return

end

%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = getlist(Pfad,ext,recurse)

nl  = char(10);

fs  = filesep;

c = cell(0,5);  % { Pfad  Name  Bytes  Date  DateNum }

for e = ext(:)'

    d = dircont(Pfad,cat(2,'*',e{1}));  % { IsDir  Date  Bytes  Name  DateNum };

    if ~isempty(d)

       d = d( find( ( cat(1,d{:,1}) == 0 ) ) , : );

       if ~isempty(d)

          d(:,1) = { Pfad };

          c = cat( 1 , c , d(:,[1 4 3 2 5])  );

       end

     end

end

if ~recurse
   return
end

%----------------------------------------------------
% Recurse

d = dircont(Pfad,'/');

if ~isempty(d)

   bad = ( strcmp( d(:,4) , '..' )  |  strcmp( d(:,4) , '..' ) );

   d = d( find( ( cat(1,d{:,1}) == 1 ) & ( ~bad ) ) , : );

   if ~isempty(d)

      for ii = 1 : size(d,1)

          p = cat( 2 , Pfad , d{ii,4} , fs );

          c = cat( 1 , c , getlist(p,ext,recurse) );

       end

    end
   
end



%****************************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,Pfad,ext,recurse] = checkin(Pfad,varargin);


% CHECKIN  Check of Inputs 


msg = '';
nl  = char(10);

fs  = filesep;

Nin = nargin;

p0 = cd;

if isempty(p0)
   p0 = fs;
end

%*********************************************************
% Defaults

recurse = 0;

ext = cell(1,0);

ext_all = char( ones(1,0) );

%*********************************************************
% Get Pfad

if Nin < 1
  Pfad = '';
end

if isempty(Pfad)
   Pfad = p0;
end

if ~( ischar(Pfad) &  ( size(Pfad,2) == prod(size(Pfad)) ) )
  msg = 'Pfad must be a String.';
  return
end

%*********************************************************
% Check Pfad

Pfad = cat( 2 , Pfad , fs(1:(end*(~strcmp(Pfad(end),fs)))) );

if ~( exist(Pfad,'dir') == 7 )
   msg = (['Directory doesn''t exist. ' Pfad ])
   return
end

d = dir(Pfad);

if isempty(d)
   msg = (['Can''t read Directory: ' Pfad ])
   return
end

%*********************************************************
% Get Options from VarArg


if Nin <=  1

  ext = { ext_all };

  return

end

 
  VarIn = varargin;
  VarIn = VarIn(:);

  Nin = size(VarIn,1);

  msgK = '';

  for ii = 1 : Nin

    if ischar(VarIn{ii})  &  ~isempty(VarIn{ii})  &  ...
       ( size(VarIn{ii},2) == prod(size(VarIn{ii})) )

       %--------------------------------------------------------
       recurse = ( recurse | strcmp(VarIn{ii},'-r') );

       %--------------------------------------------------------
       if strcmp( VarIn{ii}(1) , '.' );
          ext = cat(2,ext,VarIn(ii));
       end

    end 
    % if char

  end
  % ii

  %---------------------------------------------------------------------------
  get_all = isempty(ext);
  if ~get_all
      get_all = any( strcmp( ext , '.' ) );
  end
 
  if get_all

     ext = { ext_all };

  end


