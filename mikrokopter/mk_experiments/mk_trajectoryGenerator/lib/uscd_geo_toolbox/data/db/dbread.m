function [msg,ini] = dbread(fid,ini,cid,cnf)

% DBREAD  Read Variables from DataBaseFile
%
% [Msg,INI] = DBREAD( File , INI )
%
% Empty INI.ID returns VariableInfo only, 
%  don't read Data!!!
%
% File can be a valid FileIdentifer by FOPEN.
%


msg  = '';

%-------------------------------------------------------

Nin  = nargin;

if Nin < 2
   ini = [];
end

if Nin < 3
   cid = [];
end

if Nin < 4
   cnf = [];
end

%-------------------------------------------------------

if isempty(ini)
   ini = dbident;
elseif ~isempty(ini.Type)
   ini.Type = zeros(size(ini.Type));
end


if isempty(cnf)
   cnf = dbget;
end

if isempty(cid)
   cid = dbident(cnf.ColumnKeys);
   cid = cid.ID;
end

%-------------------------------------------------------
% Check File

new = 0;

ok = ( isnumeric(fid) & ( prod(size(fid)) == 1 ) );
if ok
   file = fopen(fid);
   if isempty(file)
      msg = 'Invalid FID. ';
      return
   end
   if feof(fid)
      return
   end
elseif ~chkstr(fid,1)
   msg = 'Invalid FID or FileName. ';
   return
else
   file = fid;
   fid  = fopen(file,'r');
   if isempty(fopen(fid))
      msg = 'Can''t open File. ';
      return
   end
   new = 1;
end

%-------------------------------------------------------

if feof(fid)
   if new, fclose(fid); end
   return
end

nmax = double(cnf.MaxLines);   % Max. Number of HeaderLines
lmax = double(cnf.MaxLength);  % Max. Length of HeaderLine

vnan = cnf.Dummy;

km   = cnf.Marker;
sp   = cnf.Assignment;
cm   = cnf.Comments;

acc  = 1e-10;  % Accuracy

%----------------------------------------------------------
% Read HeaderLines
%

  head  = {};
  bb    = '';
  nn    = 0; 
  cc    = 0;

  fpos  = 0;

while ( cc == 0  ) & ( nn < nmax ) & ~feof(fid)

    % Goto end of previous HeaderLine
    %  ( Matlab7 jumps 3 Characters  at read Data like:
    %    [dat,cc] = fscanf(fid,'%f',1);
    %    at the Line: "InstrDepth   = ... " !!! )

    fseek(fid,fpos,'bof');  

    % read header line, max. LMAX Characters

    [bb,tt] = fgets(fid,lmax);

    fpos = ftell(fid);

    if ~isequal(bb,-1)

        nn = nn + 1;

        if isempty(tt) & ~feof(fid)
           msg = 'Long Lines. ';
           break
        end
 
        if ~isempty(tt)
            bb = bb( 1 : end-prod(size(tt)) );
        end

        if ~isempty(bb)

            %--------------------------------------------------
            % Check for valid Characters

            if ~all( ( bb ==  9 ) |  ...
                     ( bb == 10 ) |  ...
                     ( bb == 13 ) |  ...
                     (  28 <= bb  &   bb <= 126 ) | ...
                     ( 160 <= bb  &   bb <= 255 )          )
                msg = 'Invalid Characters. ';
                break
            end

            bb = rmblank(bb,2+i);

            %--------------------------------------------------
            % Check for KeyWordMarker

            if ~( isempty(bb) | isempty(km) )
                sk = size(km,2);
                if ~strcmp(bb(1:min(end,sk)),km)
                   bb = '';
                else
                   bb = bb( sk+1 : end );
                   bb = rmblank(bb,2+i);
                end
            end

            %--------------------------------------------------
            % Check for Letter and Comments

            if ~isempty(bb)
                if ~isletter(bb(1))
                    bb = '';
                elseif ~isempty(cm)
                    for c = cm(:)'
                        if ~isempty(c{1})
                            nc = size(c{1},2);
                            if size(bb,2) < nc
                               break
                            end
                            ii = findstr(bb,c{1});
                            if ~isempty(ii)
                                bb = bb( 1 : ii-1 );
                                if isempty(bb)
                                   break
                                end
                            end
                        end
                    end
                end
            end

            %--------------------------------------------------
            % Check for AssignmentCharacter

            if ~isempty(bb)
                if ~any( bb == sp )
                    bb = '';
                else
                    bb = rmblank(bb,2-i);
                end
            end

            if ~isempty(bb)
                head = cat(1,head,{bb});
            end

        end

        if ~feof(fid)
            % read One Data; if this fails continue with header lines
            [dat,cc] = fscanf(fid,'%f',1);
        end

    end

end

%----------------------------------------------------------

if isempty(msg)
   if ( nn == nmax ) & ( cc == 0 )
      msg = 'To much Lines. ';
   elseif isempty(head)
      msg = 'No Header. ';
   else
      head = dbident(head);
      if isempty(head.ID)
         msg = 'Invalid Header. ';
      end
   end
end

if ~isempty(msg)
    if new, fclose(fid); end
    return
end

%----------------------------------------------------------
% Get Number of Columns from 1. DataLine

ncol = 0;

if ~( ( cc == 0 ) | feof(fid) )

    fseek(fid,fpos,'bof');              % Goto End of Header

    % Read over empty Lines
    while ~feof(fid) & ( nn < nmax ) 
           bb = fgetl(fid);
           if ~isempty(bb)
              break
           end
           nn = nn + 1;
    end

    if ( nn == nmax ) & isempty(bb)
        msg = 'To much Lines. ';
        if new, fclose(fid); end
        return
    end

    if ~( feof(fid) | isequal(bb,-1) )
        ncol = prod(size(sscanf(bb,'%f')));  % 1. DataLine
    end

    fseek(fid,fpos,'bof');              % Goto End of Header

end

%----------------------------------------------------------
% Get ID's of DataColumns if Data exist

data = dbident;

fld  = fieldnames(data);
fld  = fld(:)';

if ~( ncol == 0 )

    isc = zeros(1,0);

    for ii = cid(:)'
        jj = ( head.ID == ii );
        if any(jj)
           jj = find(jj);
           isc = cat( 2 , isc , max(jj) );
        end
    end

    if isempty(isc)

       msg = 'No Columns. ';

    else

       for ic = isc
           cc = head.Value{ic};
           if ~isempty(cc)
               cc = rmblank(cc,2);
               dd = eval('dbident(cc)','[]');
               if ~isempty(dd)
                   for f = fld
                       data = setfield( data , f{1} , ...
                          cat(1,getfield(data,f{1}),getfield(dd,f{1})) );
                   end
               end
           end
       end

       if isempty(data.ID)
          msg = 'No valid Columns. ';
       else
          nd   = size(data.ID,1);
          if ~( nd == ncol )
              msg = 'Differing Columns. ';
          end
          if nd > ncol    
             for f = fld
                 v    = getfield(data,f{1});
                 data = setfield(data,f{1},v(1:ncol));
             end
          end
       end
    end

end

%----------------------------------------------------------
% Return HeaderID's, Names and Type

if  isempty(ini.ID)

    head.Type = 1 + 0 * head.Type;
    data.Type = 2 + 0 * data.Type;

    for f = fld
        ini = setfield( ini , f{1} , ...
                cat(1,getfield(head,f{1}),getfield(data,f{1})) );
    end
 
    if new, fclose(fid); end

    return

end

%----------------------------------------------------------
% Find requested Variables in Header or Data

nh = prod(size(head.ID));
nd = prod(size(data.ID));

hok = zeros( max(nh,1) , 1 );
dok = zeros( max(nd,1) , 1 );

ish = ( nh > 0 );
isd = ( nd > 0 );

for ii = 1 : prod(size(ini.ID))

    hid = ( ini.ID(ii) == head.ID );                 % HeaderID
    hnm = strcmp( ini.Name{ii} , head.Name );        % HeaderName
    jj = ( hid & hnm );
    if any(jj)
       jj = max(find(jj));                          % take last
       hok(jj) = ii;
    else
       if isd
          did = ( ini.ID(ii) == data.ID );           % DataID
          dnm = strcmp( ini.Name{ii} , data.Name );  % DataName
          jj = ( did & dnm );
       end
       if any(jj)
          jj = min(find(jj));                       % take first
          dok(jj) = ii;
       else
          if any(hnm)                      % HeaderName, take last
             jj = max(find(hnm));
             hok(jj) = ii;
          elseif any(hid)                  % HeaderID, take last
             jj = max(find(hid));
             hok(jj) = ii;
          elseif isd
             if any(dnm)                   % DataName, take first
                jj = min(find(dnm));
                dok(jj) = ii;
             elseif any(did)               % DataID, take first
                jj = min(find(did));
                dok(jj) = ii;
             end
          end
       end
    end
end

%----------------------------------------------------------
% Read and Reshape Data / Check Data for Dummys

if ~feof(fid) & ( ncol > 0 ) & any(dok)

    f0 = ftell(fid);

    dat = fscanf(fid,'%f');

    if ~feof(fid)

        f1 = ftell(fid); fseek(fid,0,'eof');
        fe = ftell(fid); fseek(fid,f1,'bof');

        pp = 100 * ( f1- f0 ) / ( fe - f0 ) ;

        msg = cat( 2 , msg , sprintf( 'Read %.0f%%. ', floor(pp) ) );

    end

    % Reshape data

    nrow = floor( prod(size(dat)) / ncol );

    dat = reshape(dat(1:(nrow*ncol)),ncol,nrow)';

    % replace dummy values with NaN's
      isn = zeros(size(dat));
      for v = vnan(:)'
          isn = ( isn | ( abs(dat - v) < acc ) );
      end
      isn = find(isn);
      dat(isn) = NaN;

end

%----------------------------------------------------------

ini.Value(:) = {[]};

%----------------------------------------------------------
% Get requested HeaderVariables

    jj = ~( hok == 0 );
    if any(jj)
       jj = find(jj);  % VarIndex
       for ii = jj(:)'
           kk = hok(ii);
           ini.Value{kk} = head2val( head.Value{ii} , ini.Format{kk} , ...
                                      ini.Length(kk) , vnan );
            ini.Type(kk) = 1;
       end
    end

%----------------------------------------------------------
% Get requested DataVariables

    jj = ~( dok == 0 );
    if any(jj)
       jj = find(jj);  % VarIndex
       for ii = jj(:)'
           kk = dok(ii);
           ini.Value{kk} = dat(:,ii);
            ini.Type(kk) = 2;
       end
    end

%----------------------------------------------------------

if new, fclose(fid); end

