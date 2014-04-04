function txt = enc_dos(file,outfile)

% ENC_DOS  Display's bad Characters from DOS-ASCII-Files
%
%  EncodeText = ENC_DOS( DOS_ASCII_File )
%
%  further Options:
%
%  EncodeText = ENC_DOS( DOS_ASCII_File , OutFile)
%  
%    writes the EncodeText to OutFile
%

Nin = nargin;

txt = '';

fid = fopen(file,'r');

if fid == -1
  error(['Can''t open File: '  file ]);
end

bb = fread(fid,'char');

fclose(fid);

% Remove CR  !!!!!!!!!!!!!!!!!!!!!!!!!!!
bb( find(bb==13) ) = []; 


% Indize of Bad Characters

 ii = find( ~( (  bb ==  9 ) |  ...
               (  bb == 10 ) |  ...
               (  bb == 13 ) |  ...
               (  28 <= bb  &   bb <= 126 ) | ...
               ( 161 <= bb  &   bb <= 255 )       ) );

if isempty(ii)
   return
end

bad = bb(ii);

bb(ii) = double('~');

[bad,s_i] = sort(bad);


% Extract 20 Characters arround the Bad Characters
 
n2 = 10;

ii          = ii(s_i) - n2;

ii          = ii(:,ones(1,2*n2+1));
ii(:,2:end) = 1;

ii          = cumsum(ii,2);

bb = cat(1,bb,32*ones(n2,1));

txt         = bb(ii);

txt( find( txt ==  9 ) ) = 32;  % TAB --> Space
txt( find( txt == 10 ) ) = 32;  % NL  --> Space

txt        = char(txt);


btxt = sprintf('''%4.0f: '';',bad);
btxt = eval([ '['  btxt(1:end-1)  ']' ]); 

txt = cat(2,btxt,txt);


if Nin == 2
 
  fid = fopen(outfile,'wt');

  if fid == -1
    fprintf([ char([10 10]) 'Can''t open OutFile: ' outfile char([10 10]) ]);
  end

  ftxt = cat( 2 , txt , char(10)*ones(size(txt,1),1) );
  ftxt = ftxt';
  ftxt = ftxt(:)';

  ftxt = strrep(ftxt,'%','%%');
  ftxt = strrep(ftxt,'\','\\');
 
  fprintf( fid , ftxt );
 
  fclose(fid);

end
