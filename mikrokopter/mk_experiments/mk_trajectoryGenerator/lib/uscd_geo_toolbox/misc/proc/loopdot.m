function cl = loopdot(sc,n,z,force,c0)

% LOOPDOT  Shows the Progress of a Loop in the MatlabCommandWindow
%
%------------------------------------------------------
% 1. Initialising (Intro)
%
%  LOOPDOT( Scale , LoopLength , Name )
%
%    Scale   = [ DotsPerRow  MaxCountsPerRow ];
%   
%    Name    Name of Progress in Intro, optional
%
%------------------------------------------------------
% 2. Show Progress
%
%  LOOPDOT( Scale , LoopLength , Count )
%
%    Scale   = [ DotsPerRow  MaxCountsPerRow ];
%    Scale   = [ DotsPerRow  MaxRowNumber+i ];
%
%  if the Value for DotsPerRow is negative, the Number 
%   of proceeded Counts after each Row will displayed
%
%  NOTE: for faster proceeding of LOOPDOT (no Check of Inputs),
%          use an 4. Input, like: 
% 
%        LOOPDOT( Scale , LoopLength , Count , 0 )
%
%  Give the Number of Blanks before each Line as 4. Input:
%
%        LOOPDOT( Scale , LoopLength , Count , NBlank )
%        
%------------------------------------------------------
% 3. Show Progress with Elapsed Time
%
%  CL = LOOPDOT( ... ) returns the actual Clock: [ YY MM DD hh mm ss ]
%
%  Use this Output from Initialisation to diplay a Information
%   about Elapsed Time during Progress:
%
%  LOOPDOT( Scale , LoopLength , Count , 1 , CL )
%
%  LOOPDOT( Scale , LoopLength , Count , 1 , CL_NUM )
%
%
%  In that case, the Value of DotsPerRow can be positive.
%
%------------------------------------------------------
% 4. Examples
%
%------------------------------------------------------
%% Basic
%
%    n = 1000;
%   sc = [ 50  100 ];
%
%   loopdot(sc,n,'Calculate SIN')
%
%   for ii = 1 : n
%       sin(1:1000);
%       loopdot(sc,n,ii,1);
%   end
%
%------------------------------------------------------
%% + CountInfo
%
%    n = 1000;
%   sc = [ -50  100 ];
%
%   loopdot(sc,n,'Calculate SIN')
%
%   for ii = 1 : n
%       sin(1:1000);
%       loopdot(sc,n,ii,0);
%   end
%
%------------------------------------------------------
%% + TimeInfo
%
%    n = 1015;
%   sc = [ 50  100 ];
%
%   cl = loopdot(sc,n,'Calculate SIN')
%
%   for ii = 1 : n
%       sin(1:1000);
%       loopdot(sc,n,ii,6,cl);
%   end
%
%------------------------------------------------------
%% Max. 5 Rows
%
%    n = 1015;
%   sc = [ 50  5+i ];
%
%   cl = loopdot(sc,n,'Calculate SIN')
%
%   for ii = 1 : n
%       sin(1:1000);
%       loopdot(sc,n,ii,6,cl);
%   end
%

%*************************************************
% Check Inputs if NO force

Nin  = nargin;
Nout = nargout;

if Nin < 4

  if Nin < 2
     error('Not enough InputArguments.');
  end

  if Nin < 3
     z = '';   % Name
  end

  msg = check_in(sc,n,z);

  if ~isempty(msg)
     error(msg)
  end

end

if Nin < 5
   c0 = [];
end

nc = prod(size(c0));
if ~( nc == 0 )
   if ~( isnumeric(c0) & any( nc == [ 1  6 ] ) )
       c0 = [];
   end
end


cl = clock;


%*************************************************

show_count = ( ( sc(1) < 0 ) | ~isempty(c0) );

nr = ~( imag(sc(2)) == 0 );  % Max. Numberof Rows
sc = abs(real(sc));

if nr
   sc(2) = n / sc(2);
   sc(2) = sc(1) * ceil( sc(2) / sc(1) );
end


cpp = floor(min(sc(2),n)/sc(1));  % CountsPerPoint

cpp = cpp + ( cpp == 0 );

cpl = cpp * sc(1);                % CountsPerLine

%*************************************************
% Intro

if ischar(z)

  if isempty(z)
     z = 'Run loop';
  end

  fprintf(1,'%s %.0f times' , z , n );

  if cpp > 1
     fprintf(1,', %.0f per Point' , cpp );
  end
   
  if n > cpp*sc(1)
     fprintf(1,', %.0f per Row' , cpl );
  end

  fprintf(1,'\n')

  if Nout == 0
     clear cl
  end

  return

end

%*************************************************

show_end = ( z == n );

show_dot = ( mod(z,cpp) == 0 );
show_nl  = ( mod(z,cpl) == 0 );
show_one = ( mod(z,cpl) == 1 );

% Number of Blanks to add if show_end

nb = sc(1) - floor( ( z - cpl * floor(z/cpl) ) / cpp );

nb = nb * show_end * (~show_nl);

% if show_end,keyboard,end
% nb = ( nb - 0*show_dot ) * show_end * (~show_nl);

show_nl = ( show_nl | show_end );

show_count = ( show_count & show_nl );

%-------------------------------------------------

if show_one & ( Nin >= 4 )
   if isnumeric(force) & ~isempty(force)
      nbl = ceil(abs(real(force(1))));
      if nbl > 0
         fprintf(1,'%s',char(32*ones(1,nbl)));
      end
   end
end

if show_dot
   fprintf(1,'.');
end

if show_count
   
   pt10  = floor( ( log(n) / log(10) ) + 1e2*eps ) + 1;
   zform = sprintf('%%%.0f.0f of %%%.0f.0f',pt10([1 1]));
   str   = sprintf(zform,z,n);

   fprintf( 1 , '%s %s  %3.0f%%', char(32*ones(1,nb)) , str , floor(z/n*100) );

   if ~( nc == 0 )
       if ( nc == 6 )
          c0 = datenum(c0(1),c0(2),c0(3),c0(4),c0(5),c0(6));
       end
       dt = datenum(cl(1),cl(2),cl(3),cl(4),cl(5),cl(6)) - c0;
       dt = dt * 24 * 3600;          % Elapsed sec
       et = ( n - z ) * ( dt / z );  % Estimated sec
       dt = ceil(dt); dm = fix(dt/60); ds = dt - 60*dm;
       et = ceil(et); em = fix(et/60); es = et - 60*em;
       if z == n
          fprintf(1,' %4.2d:%2.2d',dm,ds);
       else
          fprintf(1,' %4.2d:%2.2d / %3.2d:%2.2d',dm,ds,em,es);
       end
   end

end

if show_nl
   fprintf(1,'\n');
end

if Nout == 0
   clear cl
end

%****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function msg = check_in(sc,n,z);

% InputCheck

msg = '';
nl  = char(10);

  %-------------------------------------------------
  % Check Scale

  ok = ( isnumeric(sc) & ( prod(size(sc)) == 2 ) );
  if ok
     ok = all( ( abs(real(sc)) > 0 ) & ( mod(real(sc),1) == 0 ) );
  end

  if ~ok
     msg = 'Values for Scale must have 2 positive Integers.';
  end

  %-------------------------------------------------
  % Check Size

  ok = ( isnumeric(n) & ( prod(size(n)) == 1 ) );
  if ok
     ok = ( ( n > 0 ) & ( mod(n,1) == 0 ) ); 
  end

  if ~ok
     msg = [ msg nl(1:(end*(~isempty(msg)))) ...
            'Value for Size must be a positive Integer.' ];
  end

  %-------------------------------------------------
  % Check Name / Counter

  ok = ischar(z);
  if ok
     ok = ( ( prod(size(z)) == size(z,2) ) | isempty(z) );
     if ~ok 
        msg = [ msg nl(1:(end*(~isempty(msg)))) ...
                'Value for Name must be a String.' ];
     end
  else
     ok = ( isnumeric(z) & ( prod(size(z)) == 1 ) );
     if ok
        ok = ( ( 1 <= z ) & ( z <= n )  &  ( mod(z,1) == 0 ) );
     end
     if ~ok 
        msg = [ msg nl(1:(end*(~isempty(msg)))) ...
                'Value for Counter must be a positive Integer between [ 1  LoopLength ].' ];
     end
  end


