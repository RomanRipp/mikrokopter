function str = splitstr(str,splt,ins);

% SPLITSTR Split String
%
% SPLITSTR( String , Split , Insert )
%
% Split = RowLength + i*MaxLength
%
% RowLength  Maximum Length of Row,
% MaxLength  Maximum Length of String
%
% RowLength < 0  ==> call of RMBLANK
% 
% RowLength + OverSize/RowLength with  OverSize <= RowLength
%
% Insert     Strings, after each a NewLine will insertet 
%

Nin = nargin;

nl  = char(10);  % NewLine  (SplitCharacter)

%-----------------------------------------
% Check Inputs

if isempty(str)
   return
end

if Nin < 2
   splt = [];
end

if Nin < 3
   ins = [];
end

if ~chkstr(str)
    error('Input must be a String.');
end

if isempty(splt)
   splt = 60;
end

nmax = abs(imag(splt(1)));
splt =     real(splt(1));

if ~isempty(ins)
    [ok,ins] = chkcstr(ins);
    if ~ok
        ins = [];
    end
end

%-----------------------------------------
% Check for Quota

n = size(str,2);

is_quot = all( str([1 n]) == '"' );

if nmax == 0
   nmax =  n - 2*is_quot;
end

%-----------------------------------------
% RMBLANK if negative SPLIT

if splt <= 0
   if is_quot
      str = cat(2,str(1),rmblank(str(2:n-1),2),str(n));
   else
      str = rmblank(str,2);
   end
end

%-----------------------------------------
% Check with Length of String

nmax = nmax + sum( str == nl );

if all( n <= [ splt nmax ]+2*is_quot )
   return
end

app = '';

if n > nmax+2*is_quot

   app = ' ...';

   ind = cat( 2 , ( 1 : nmax+is_quot ) ,  n*ones(1,is_quot) );

   str = str(ind);

   n   = nmax + 2*is_quot;

   if splt <= 0
      if is_quot
         str = cat(2,str(1),rmblank(str(2:n-1),2-i),str(n));
      else
         str = rmblank(str,2-i);
      end
   end
   
end

%-----------------------------------------
% Check with Insert

if ~isempty(ins)

    n  = size(str,2);
    i1 = zeros(1,0);

    for cc = ins(:)'
        jj = findstr( str , cc{1} );
        if ~isempty(jj)
            i1 = cat( 2 , i1 , jj+size(cc{1},2)-1 );
        end
    end

    if ~isempty(i1)
        i1 = sort(i1);
       str = insert(str,i1,nl);
    end

end

%-----------------------------------------
% No SPLIT ==> return

splt = abs(splt);

ovs  = splt;
splt = floor(splt);
ovs  = ovs - splt;
ovs  = 1e-10 * round(ovs/1e-10);
ovs  = splt + ceil( ovs * splt );

if splt == 0
   if ~isempty(app)
       na  = size(app,2);
       ind = ( 1 : na ) + n-na-is_quot;
       str(ind) = app;
   end
   return
end

%-----------------------------------------
% Check Length of Segments between NewLine

n  = size(str,2);

i0 = find( str == nl );
lg = diff(cat(2,1,i0+1,n+1));

if ~any( lg > splt )
    if ~isempty(app)
        na  = size(app,2);
        ind = ( 1 : na ) + n-na-is_quot;
        str(ind) = app;
    end
    return
end

%-----------------------------------------
% Check for Characters to Split

cc = { '- '  ', '  '; '  '. '   ' ' };

i1 = cat( 2 , 0 , i0 , n+1 );
ok = 0 * i1;

ok(end) = NaN;

for ii = cc

    jj = findstr( str , ii{1} );

    if ~isempty(jj)

        len = size(ii{1},2);

        j0 =  1 + ( jj( 1) == 1 );

        j1 = size(jj,2);
        j1 = j1 - ( jj(j1) >= n-(len-1) );

        jj = jj( j0 : j1 );

        if ~isempty(jj)
     
            % NewLine   before                behind
            kk = ~( ( str(jj-1) == nl ) | ( str(jj+len) == nl ) );
 
            if any(kk)
               if ~all(kk)
                   kk = find(kk);
                   jj = jj(kk);
               end
               jj = jj + len - 1;        % End of Character to Split
               i1 = cat( 2 , i1 , jj );
               ok = cat( 2 , ok , ones(size(jj))+strcmp(ii{1},' ') );
            end

        end

    end

end

if all( ok == 0 )
   if ~isempty(app)
       na  = size(app,2);
       ind = ( 1 : na ) + n-na-is_quot;
       str(ind) = app;
   end
   return
end


[i1,si] = sort(i1);
 ok     = ok(si);

jj = ( diff(i1,1,2) == 0 );   % Same EndIndex !!!
if any(jj)
   jj = find(jj) + 1;
   i1(jj) = [];
   ok(jj) = [];
end

%-----------------------------------------
% Check for Correct NewLine-Statement
%
%  [ TAB NL Blank NBSP ] after NL
%    NL before NL in I1 

i0 = find( ok == 0 );
i0 = i0( 1+(i0(1)==1) : end-(i1(i0)==n) );

jj = ( ok(i0-1) == 0 );     % NewLine before

for cc = [ 9 10 32 160 ]    % Valid Characters behind NewLine
    jj = ( jj | ( str(i1(i0)+1) == cc ) );
end

if ~all(jj)
        jj  = find(~jj);
        i0  = i0(jj);
        jj  = i1(i0);
     i1(i0) = [];
     ok(i0) = [];
end
 
%-----------------------------------------
% Check for Correct Split-Statement
%
% No NewLine before or behind

%-----------------------------------------

m  = size(i1,2);

l1 = diff(cat(2,i1(1),i1),1,2);    % Length between Characters to Split, incl.

l1 = cumsum(l1,2);

ind = ( 2 : m );

%%% [ splt  ovs ]

while 1

   % Remove Length-Offset

   n1 = zeros(1,m);
   i0 = find( ok <= 0 );
   n1(i0) = 1;
   n1 = cumsum(n1,2);
   n1 = l1(i0(n1));

   ls = ( l1 - n1 ) .* ~isnan(ok);

   ns = floor(ls/splt);

   jj = ( ( ns(ind-1) == 0 ) & ( ns(ind) > 0 ) );

   if ~any(jj)
       break
   end

   jj = find(jj) + 1;

   % NewLine in Range of OverSize or Split after Blank in Range of OverSize
   kk = ( ( ( ok(jj+1) == 0 ) & ( ls(jj+1) <= ovs ) ) | ...
          ( ( ok(jj+1) == 1 ) & ( ls(jj+1) <= ovs ) & ( ok(jj) == 2 ) ) );

   %%% lk = 0*ok; lk(jj) = 1-2*kk; [i1 ; l1 ; ok ; n1 ; ls ; ns ; lk] 

   ok(jj) = -1;

   if any(kk)
        m    = m - sum(kk);
        ind  = ind( 1 : (m-1) );
         kk  = jj(find(kk));
      i1(kk) = [];
      l1(kk) = [];
      ok(kk) = [];
   end

end

ok = ( ok == -1 );
if any(ok)
   ok = find(ok);
   i1 = i1(ok);
   str = insert(str,i1,nl);
end

if ~isempty(app)
    n   = size(str,2);
    na  = size(app,2);
    ind = ( 1 : na ) + n-na-is_quot;
    str(ind) = app;
end

