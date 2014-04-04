function [str,val] = splitunit(str,ref);

% SPLITUNIT  Splits a HeaderValue in Value / Unit
% 
% Check for Parenthesis-Stements in Values like for Units
%

Nout = nargout;

val = '';

if isempty(str)
   return
end

if nargin < 2
   ref = cat( 1 , '|/\([{' , ...
                  '|/\)]}'       );
end

if isempty(ref)
   return
end

ok = zeros(size(str));

for cc = ref(1,:)
    ok = ( ok | ( str == cc ) );
end

if any(ok)

    ok = min(find(ok));
    sp = str(ok);
   val = str( ok+1 : end ); % Unit
   str = str( 1 : ok-1 );   % Value only

    if ~isempty(val) & ( Nout == 2 )
        val = rmblank(val,2);
        if ~isempty(val)
            jj = find( sp == ref(1,:) );
           val = val( 1 : ( end - double(val(end)==ref(2,jj)) ) );
           val = rmsurrd(val,ref(:,2:end));
        end
    end

end

str = rmblank(str,2);
