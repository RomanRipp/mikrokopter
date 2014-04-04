function [dat,ini,varargout] = dbload(ident,vars)

% DBLOAD  Load DataBase files.
%
% [DAT,INI] = DBLOAD(IDENT,VARS)
%
% IDENT      : string containing datafilenames or identifers
%               see: DBFILE
%
% VARS       : string or cellstring containing variables identifers
%               see: DBIDENT
%
%              empty: load all Variables
%              ZERO:  load HeaderVariables only
%              ONE:   load   DataVariables only
%
% DAT  VariableStructure
% INI  IdentStructure for Variables and DataFiles
%
% In case of multiple Datafiles the HeaderVariables are concatenated 
%  along the 1. Dimension (individualy rows), the DataVariables along
%  the 2. Dimension (individualy columns).
%
% The Values of HeaderVariables will be converted by their Format if
%  the Variables are defined in DBINIT.
%
% The requested variables are returned in individualy outputs
%  if more then 2 outputs are requested:
%
% [V1,V2,V3, ... ] = DBLOAD( IDENT , 'V1:V2:V3:...' )
% [V1,V2,V3, ... ] = DBLOAD( IDENT , {'V1' 'V2' 'V3' ... } )
%
%
% The Configurations for the DataFileNames, IdentiferSeparator etc.
%  can be modified using DBSET
%
% see also: DBFILE, DBIDENT, DBINIT, DBREAD, DBSCAN, DBSET
%
%


Nin  = nargin;
Nout = nargout;

if Nin < 1
   ident = '*:*:*';
end

if Nin < 2
   vars = {};
end

varargout = cell(1,Nout-2);

%*********************************************************************
% Check Inputs

msg = cell(0,1);

[ok,ident] = chkcstr(ident);

if ~ok
    msg = cat(1,msg,{'Identifer must be a CharacterArray or CellArray of Strings.'});
end

if ~( isempty(vars) | isequal(vars,0) | isequal(vars,1) | ...
       chkstr(vars) | chkcstr(vars)  )
    msg = cat(1,msg,{'Variables must be a CharacterArray or CellArray of Strings.'});
end

%-----------------------------------------------------------------

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%*********************************************************************
% Check Files and VariableIdentifer

if isempty(vars) | isnumeric(vars)

   [file,ini,msg] = dbscan(ident);

   if ~isempty(vars)
       ok = find( ini.Type == 1+isequal(vars,1) );
       
      fld = fieldnames(ini);
      for f = fld(:)'
          v   = getfield(ini,f{1});
          v   = v(ok,:);
          ini = setfield(ini,f{1},v);
      end

   end

else

   % Get VariableIdentifers

   ini = dbident(vars);

   % Check for Duplicate Names

   ini = chknam(ini,0);

   % Get Files

   [file,msg] = dbfile(ident);

end

if ~isempty(msg)
    dbmsg(msg,1);
end

if isempty(file) | isempty(ini.ID)
   dat = [];
   return
end

nvar = size(ini.ID,1);

%-----------------------------------------------------------------
% Prepare

nf = prod(size(file));

ok  = zeros(nvar,nf);

nd = size(ini.ID,1);

dat = cell(nd,1);

% FormatStrings which defines CellArrays
cfrm = { '%s' 'list' 'cell' };

for ii = 1 : nd
    if ini.Type(ii) == 1
       if    ( ini.Length(ii) < 0 ) | ...
          any( strcmp( ini.Format{ii} , cfrm ) )
          dat{ii} = {};
       else
          dat{ii} = zeros(0,max(1-isnan(ini.ID(ii)),ini.Length(ii)));
       end
    else
       dat{ii} = zeros(0,0);
    end
end

%-----------------------------------------------------------------
% Common Inputs for DBREAD


cnf = dbget;

cid = dbident(cnf.ColumnKeys);
cid = cid.ID;

%-----------------------------------------------------------------
% Loop over Files

typ = ini.Type;

v = sprintf(['%s' cnf.Separator],ini.Name{:});
v = v( 1 : end-1 );

dbmsg([v char(10)],1);

for ff = 1 : nf

    dbmsg(file{ff});
      
    [msg,ini] = dbread(file{ff},ini,cid,cnf);

    if all( ini.Type == 0 )
       msg = cat( 2 , msg , 'No Data.' );
    end

    if ~isempty(msg)
        dbmsg(sprintf(': %s',msg));
    end        

    dbmsg(char(10));

    ok(:,ff) = ~( ini.Type == 0 );

    if any(ok(:,ff))
 
       %-------------------------------------------------------------------------------
       % Header
       %-------------------------------------------------------------------------------

       tp0 = typ;

       jj = ( ini.Type == 1 );
 
       if any(jj)
          jj = find(jj);
          typ(jj) = typ(jj) + double( ( typ(jj) == 0 ) | ( typ(jj) == 2 ) );  % 0/2 --> 1/3
          for ii = jj(:)'
              v = dat{ii};
              w = ini.Value{ii};
              if tp0(ii) == 0
                 v = permute(v,[2 1]);
              end
              n = size(v,2);
              if ~iscell(v) 
                  if ( typ(ii) == 3 )  % Data before
                     vv = cell(n,1);
                     for iv = 1 : n
                         vv(iv) = {v(:,n)};
                     end
                     v = vv;
                  elseif ischar(w)
                     if isempty(v)
                        v = {};
                     elseif tp0(ii) == 0
                        v = cell(size(v,1),1);
                        v(:) = {''};
                     else
                        v = cellstr(v);
                     end
                  end
              end

              if iscell(v)
                 v = cat( 1 , v , {w} );
              else
                 m = size(w,2);
                 if     n < m
                        v = cat( 2 , v , NaN * zeros(size(v,1),m-n) );
                 elseif m < n
                        w = cat( 2 , w , NaN * zeros(size(w,1),n-m) );
                 end
                 v = cat( 1 , v , w );
              end
              dat{ii} = v;
          end
       end

       %-------------------------------------------------------------------------------
       % Data
       %-------------------------------------------------------------------------------

       jj = ( ini.Type == 2 );

       if any(jj)
          jj = find(jj);
          typ(jj) = typ(jj) + 2 * double( ( typ(jj) == 0 ) | ( typ(jj) == 1 ) );  % 0/1 --> 2/3
          for ii = jj(:)'
              v = dat{ii};
              w = ini.Value{ii};
              n = size(v,1);
              if ~iscell(v) & ( typ(ii) == 3 ) % Numerical Header before
                  vv = cell(n,1);
                  for iv = 1 : n
                      vv(iv) = {v(iv,:)};
                  end
                  v = vv;
              end
              if iscell(v)
                 v = cat( 1 , v , {w} );
              else
                 m = size(w,1);
                 if     n < m
                        v = cat( 1 , v , NaN * zeros(m-n,size(v,2)) );
                 elseif m < n
                        w = cat( 1 , w , NaN * zeros(n-m,size(w,2)) );
                 end
                 v = cat( 2 , v , w );
              end
              dat{ii} = v;
          end
       end

       %-------------------------------------------------------------------------------
       % Nothing
       %-------------------------------------------------------------------------------

       jj = ( ini.Type == 0 );

       if any(jj)
          jj = find(jj);
          for ii = jj(:)'
              v = dat{ii};
              if iscell(v)
                 v = cat( 1 , v , {[]} );
              elseif typ(ii) == 1
                 v = cat( 1 , v , NaN * ones(1,size(v,2)) ); % Add a Row
              else   % any( typ(ii) == [ 0  2 ] )
                 v = cat( 2 , v , NaN * ones(size(v,1),1) ); % Add a Column
              end
              dat{ii} = v;
          end
       end

    end

end

%----------------------------------------------------------------
% Check for individual VariableOutput
%  (old RODBLOAD-Statement)

if ( Nout > 2 )

   nd = size(dat,1);

   for ii = 3 : min(Nout,nd)
       varargout{ii-2} = getout(dat{ii});
   end

   if nd >= 2
      ini = getout(dat{2});
   else
      ini = [];
   end

   if nd >= 1
      dat = getout(dat{1});
   else
      dat = [];
   end

   return

end

%----------------------------------------------------------------
% Build OutputStructure

for ii = 1 : nd
    dat{ii} = dat(ii);
end

dat = cat( 2 , ini.Name , dat );

dat = permute(dat,[2 1]);

dat = struct(dat{:});

if Nout == 2

   ini.Type  = typ;
   ini.Value = {};

   ok = any(ok,1);

   if ~all(ok)
       if ~any(ok)
           file = {};
       else
           ok = find(ok);
           file = file(ok);
       end
   end

   ini.File = file;

end


%*******************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function v = getout(v)

if isnumeric(v) | ischar(v)
   return
end

if iscellstr(v)
   v = char(v);
end

