function rcp(varargin)

% RCP  Executes RCP
%
% RCP( Quelle , Ziel )
%
    

if nargin < 2

  fprintf([ char(10) 'Usage: rcp Quelle Ziel '  char([10 10]) ]);

  return

end

VarIn = varargin;

command = ['scp ' strhcat(VarIn,' ')    ' &'];

fprintf([ char(10) 'UNIX: '  command char([10 10]) ]);

[s,w] = unix(command);

if ~isempty(w)
  fprintf([ char(10) w char([10 10]) ]);
end
