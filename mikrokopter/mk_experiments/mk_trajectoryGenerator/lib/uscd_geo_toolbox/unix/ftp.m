function ftp(opt)

% FTP  Open's FTP-Connection
%
% FTP( Server )
%
    

if nargin == 0

  opt = '';

end

if ~( ischar(opt)  &  ( prod(size(opt)) == size(opt,2) ) )
  error('FTP: Input Server must be a String.')
end


command = ['ftp ' opt  ];

fprintf([ char(10) 'UNIX: '  command char([10 10]) ]);

unix(command);

