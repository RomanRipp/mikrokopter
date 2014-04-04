function [n,ok] = numfiles(pfad,ext,f)

% NUMFILES  Numerate Files by ascending Order of Date
%
% [NN,Ok] = NUMFILES( Directory , WildCard , [Number] )
%
% NN = { OldName  NewName }
%
% Ok =  1  succesfull
%       0  error rename Files
%      -1  File with NewName allready exists
%
% NewName = sprintf( '%<Number>.<Number>d_%s' , OldName )
%

n  = cell(0,2);
ok = zeros(0,1);

Nin = nargin;

if Nin < 2
   error('Not enough InputArguments');
end

if ~( chkstr(pfad,1) & chkstr(ext,1) )
    error('Inputs must be nonempty Strings.');
end

if ~( exist(pfad,'dir') == 7 )
    error('Directory doesn''t exist.');
end

if Nin < 3
   f = [];
else
   ok = ( isnumeric(f) & ( prod(size(f)) == 1 ) );
   if ok
      ok = ( ( f > 0 ) & ( mod(f,1) == 0 ) );
   end
   if ~ok
       error('Format must be a single positive Integer.');
   end
end

if ~isunix
    warning('UNIX-Enviroment required.');
    return
end

%***********************************************************
% Get Files

if ~any( ( ext == '*' ) | ( ext == '?' ) )
    ext = cat(2,'*',ext);
end

fs = filesep;

pfad = cat( 2 , pfad , fs(1:(end*(~strcmp(pfad(end),fs)))) );

fprintf(1,'%s: ',[pfad ext]);

d = dir([pfad ext]);

if ~isempty(d)
    d = d(find(~cat(1,d.isdir)));
end

if isempty(d)
   fprintf(1,'%s\n','No Files found.');
   return
end

m = prod(size(d));

fprintf(1,'%.0f Files\n',m);

if isempty(f)
   f = floor(log(m)/log(10)) + 1;
end

f = sprintf('%s%.0f.%.0fd_','%',f,f);

n = {d.name};
d = datenum(char(ger2eng({d.date})));

[d,si] = sort(d);

n = n(si);

d = round(datevec(d));

ds = diff(d(:,6),1,1);
if any( ds <= 0 )
   ds = max(1,ds);
   d(:,6) = cumsum(cat(1,d(1,6),ds)) + ( d(m,6) - (d(1,6)+sum(ds)) );
   d = datenum(d(:,1),d(:,2),d(:,3),d(:,4),d(:,5),d(:,6));
   d = round(datevec(d));
end

%***********************************************************
% Rename Files

n = n(:);
n = n(:,[1 1]);

ok = zeros(m,1);

mv = 'mv "%s" "%s"';
tc = 'touch -t %2.2d%2.2d%2.2d%2.2d %4.0f .%2.2d "%s"'; % -t MMDDhhmm [ [CC]YY] [.ss] ]

d = d(:,[2 3 4 5 1 6]); %  [ MM DD hh mm YYYY ss ];

d = round(d);

for ii = 1 : m

    nr = sprintf(f,ii);

    ok(ii) = strncmp(n{ii,1},nr,size(nr,2));

    if ok(ii)

        n{ii,2} = n{ii,1};

    else

        n{ii,2} = [ nr  n{ii,1} ];

    end

    fprintf(1,'%s ---> %s  ...  ',n{ii,:});

    n{ii,1} = [ pfad  n{ii,1} ];
    n{ii,2} = [ pfad  n{ii,2} ];

    if ~ok(ii)

        ok(ii) = -1 * ( exist(n{ii,2},'file') == 2 );

        if ~ok(ii)
             c    = sprintf(mv,n{ii,:})
%            [s,w] = unix(c);
 s = 0;
           ok(ii) = ( s == 0 );
           if ok(ii)
              ok(ii) = 1;( exist(n{ii,2},'file') == 2 );
              if ok(ii)
                  c    = sprintf(tc,d(ii,:),n{ii,2})
%                 [s,w] = unix(c);
              else
                  w = 'New File doesn''t exist.';
              end
           end
        end

    end

    if     ( ok(ii) == 1 )
       fprintf(1,'%s\n','ok');
    elseif ( ok(ii) == 0 )
       fprintf(1,'%s\n%s\n','error',w);
    elseif ( ok(ii) == -1 )
       fprintf(1,'%s\n','file exists');
    end

end

