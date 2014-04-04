function man(str)

% MAN  Executes MANUAL
%
% MAN( Command )
%
    

if nargin == 0

  str = '';

end

if isempty(str)

  return

end

command = ['! man ' str  ];

fprintf([ char(10)  command char([10 10]) ]);

eval(command);

