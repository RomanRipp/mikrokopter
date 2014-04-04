function cp(varargin)

% CP  Executes CP
%
% CP( Quelle , Ziel )
%
    

if nargin < 2
  error('CP Usage: cp Quelle Ziel ');
end

VarIn = varargin;

ok = iscellstr(VarIn);
if ok
   try
     str = cat(2,VarIn{:});
      ok = ( ischar(str)  &  ( prod(size(str)) == size(str,2) ) );
   catch
     ok = 0;
   end
end

if ~ok
  error('CP: Inputs must be Strings.');
end

comm = 'cp';

%----------------------------------------
% Check for MCOPY

chk = 'a:';
nc  = size(chk,2);

for vv = VarIn(:)'
    n = min( nc , size(vv{1},2) );
    if strcmp( vv{1}(1:n) , chk );
       comm = 'mcopy';
    end
end

%----------------------------------------

command = [ comm  ' '  strhcat(VarIn,' ')   ' &' ];

fprintf([ char(10) 'UNIX: '  command char([10 10]) ]);

[s,w] = unix(command);

if ~isempty(w)
  fprintf([ char(10) w char([10 10]) ]);
end
