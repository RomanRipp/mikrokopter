function [msg,head,cols,dat,time] = read_csa(file,comm)

% READ_CSA  Reads CSA-Ascii-DataFiles
%
% [Msg,Header,Columns,Data,Time] = READ_CSA( FileName )
%
% CSA format is an ASCII format in which the first 
%  five lines of the file contain the header and data information. 
%
% Lines 1 through 3 contain the static header variable:
%   names, units, and values; 
% lines 4 and 5 contain the dynamic variable names and units. 
%
% Msg contains ErrorMessages about invalid Data etc.
%
% Header  = { Name  Unit  Value }  [ K by 3 ] CellArray
% Columns = { Name Unit }          [ M by 2 ] CellArray
% Data    =  DataMatrice           [ N by M ] Numeric
%
% READ_CSA( FileName , CommentMarker ) defines the Marker
%  for Lines with Comments in header or data information.
%
% default: CommentMarker = '%%'
%
% See also: READ_ASC
%

lmax = 512;   % Max Length of Line
nmax = 64;    % MaxNumber of Lines for Header and Data

if nargin < 2
   comm = '%%';   % CommentMarker
end

msg  = '';
head = cell(0,3);
cols = cell(0,2);
dat  = [];
time = [];

%---------------------------------------------------------------
% Check, Open File

if ~chkstr(file)
    error('File must be a String.');
end

if ~( chkstr(comm) | isempty(comm) )
    error('CommentMarker must be a String.');
end


if ~( exist(file,'file') == 2 )
    msg = 'File does''nt exist.';
    return
end

fid  = fopen(file,'r');
if isempty(fopen(fid))
   msg = 'Can''t open File.';
   return
end

nrow = 0;

%---------------------------------------------------------------
% Read Header

[msg,head,nrow] = readhead(fid,3,lmax,nmax,nrow,comm);

if ~isempty(msg)
    msg = sprintf('Error read Parameter for Header.\n%s',msg);
    fclose(fid);
    return
end

[head,it] = head2val(head);

%---------------------------------------------------------------
% Read DataColumns

[msg,cols,nrow] = readhead(fid,2,lmax,nmax,nrow,comm);

if ~isempty(msg)
    msg = sprintf('Error read Parameter for Columns.\n%s',msg);
    fclose(fid);
    return
end

%---------------------------------------------------------------
% Read 1. DataLine

nc = size(cols,1);  % Number of Columns

[msg,fln,nrw] = readhead(fid,-nc,lmax,2*nmax,nrow,comm);

if ~isempty(msg)
    msg = sprintf('Error read 1. DataLine.\n%s',msg);
    fclose(fid);
    return
end

if nargout < 4
    fclose(fid);
    return
end

msg = cell(0,1);

if nrow+1 < nrw
   msg = cat(1,msg,{sprintf('Skip %.0f Rows.',nrw-(nrow+1))});
end

%---------------------------------------------------------------
% Read Data

f0 = ftell(fid);

dat = fscanf(fid,'%f');

if ~feof(fid)

    f1 = ftell(fid); fseek(fid,0,'eof');
    fe = ftell(fid); fseek(fid,f1,'bof');

    pp = 100 * ( f1 - f0 ) / ( fe - f0 ) ;

    msg = cat(1,msg,{sprintf('Read %.0f%% of Data.',floor(pp))});

end

fclose(fid);

%---------------------------------------------------------------
% Reshape data

ncol = size(fln,2);

nrow = floor( prod(size(dat)) / ncol );

if nrow > 0

   dat = reshape(dat(1:(nrow*ncol)),ncol,nrow)';

else

   dat = zeros(0,ncol);

end

dat = cat( 1 , fln , dat );

%---------------------------------------------------------------


if any(it) & ( nargout >= 5 )

   i0 = strcmp(head(:,1),'start_time');
   ii = strcmp(head(:,1),'samp_interv');
   if ~any(ii)
       ii = strwcmp(head(:,1),'*interv*');
   end
   
   if ~( sum(i0) == 1 ) & ( sum(ii) == 1 )

       warning('Can''t find "start_time" AND "samp_interv".');

   else

       i0 = find(i0);
       ii = find(ii);

       try
          t0 = datenum(head{i0,3});
       catch
          t0 = [];
       end

       if ~( isnumeric(t0) & ( prod(size(t0)) == 1 ) )
            warning('Invalid Value for StartTime.');
            t0 = [];
       end

       sc = [ 1 60 60 24 ];
       dt = head{ii,3};
       if ~( isnumeric(dt) & ( prod(size(dt)) == 1 ) )
           warning('Invalid Value for Intervall.');
           dt = [];
       else
           un = lower(head{ii,2});
           kk = 1 * strwcmp(un,'*sec*') + ...
                2 * strwcmp(un,'*min*') + ...
                3 * strwcmp(un,'*hour*') + ...
                4 * strwcmp(un,'*day*'); 
           if ( kk == 0 )
              warning('Can''t get Units for Intervall.');
              dt = [];
           else
              dt = dt * prod(sc(1:kk)) / prod(sc);
           end
       end

       if ~( isempty(t0) | isempty(dt) )

           time = t0 + dt * (0:(size(dat,1)-1))';

       end   

   end

end

%---------------------------------------------------------------

if isempty(msg)
   msg = '';
else
   msg = sprintf('%s\n',msg{:});
end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,dat,nrow] = readhead(fid,cnt,lmax,nmax,nrow,comm)

msg = '';

read_dat = ( cnt <= 0 );

if read_dat
   dat = zeros(0,cnt);
   if cnt == 0
      cnt = inf;
   else
      cnt = abs(cnt);
   end
else
   dat =  cell(0,cnt);
end

sc = size(comm,2);

col = 0;

while ~feof(fid) & ( nrow < nmax )  & ( col < cnt )

   [bb,tt] = fgets(fid,lmax);

   if ~isequal(bb,-1)

        nrow = nrow + 1;

        if isempty(tt) & ~feof(fid)
           msg = 'Long Lines.';
           return
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
                msg = 'Invalid Characters.';
                return
            end

            bb = rmblank(bb,2+i);

        end

        if ~isempty(bb) & ~( sc == 0 )
            if strcmp( bb(1:min(sc,end)) , comm )
               bb = '';
            end
        end

        if ~isempty(bb)

            if read_dat

               [v,cc] = sscanf(bb,'%f');

               if ~isempty(v)
                   if isfinite(cnt) & ~( cc == cnt )
                      msg = 'Invalid Number of Columns.';
                   else
                      dat = v(:)';
                   end
                   return
               end

            else

               bb(find(bb==9)) = 32;  % TAB !!!

               bb = sepname(bb,NaN,' ');

               col = col + 1;

               if col == 1

                  dat = bb(:);
                  dat = dat(:,ones(1,cnt));
                  dat(:,2:cnt) = {''};

               else

                  if ~( size(bb,2)  == size(dat,1) )
                       msg = 'Invalid Number of Parameter.';
                       return
                  else
                       dat(:,col) = bb(:);
                  end

               end

            end

        end

   end

end

%*************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function  [head,ok] = head2val(head)

% Transform String of HeaderValue to Value
%  valid if String contains single numeric Value
%
% Try to match DateTimeSyntax like: 
%
%   YYYYMMDDhhmmss --> [ YYYY MM DD hh  mm ss ]
%
 
n = size(head,1);

ok = zeros(n,1);

for ii = 1 : n

    str = head{ii,3};

    val = sscanf(str,'%f');

    if prod(size(val)) == 1

       prp = head{ii,1};
       uni = head{ii,2};

       prp = lower(prp);

       ok(ii) = ( isequal(size(uni),size(str)) & ( size(prp,2) >= 4 ) );
       if ok(ii)
          ok(ii) = ~( isempty(findstr(prp,'date')) & ...
                      isempty(findstr(prp,'time'))       );
       end

       if ok(ii)
          iv = ~( diff(uni) == 0 );
          if any(iv)
             nv = sum(iv) + 1;
             iv = find(iv);
            str = insert(str,iv,' ');
              v = sscanf(str,'%f',[1 nv]);
             if isequal(size(v),[1 nv])
                val = v;
             end
          end
       end

       head{ii,3} = val;

    end

end
