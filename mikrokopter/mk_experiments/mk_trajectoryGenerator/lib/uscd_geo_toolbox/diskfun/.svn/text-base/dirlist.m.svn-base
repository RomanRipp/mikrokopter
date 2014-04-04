function liste = dirlist( liste , file );

% DIRLIST  returns String-List from Output of DIRINFO
%
% LISTE = DIRLIST( LIST )
%
% LIST is a 5-Column-CellStringArray from the 8. Column of 
%        Output from DIRINFO with Option '-l'.
%
% LIST = { Identifer  Date  Bytes  Name  'NDirectory  NFiles' } 
%
% LISTE is a String, creates from LIST, 
%        Rows separated by NewLine.
%
% DIRLIST( LIST , filename )  writes LISTE into File filename
%  
%


Nin  = nargin;
Nout = nargout;

nl = char(10);

if ~( iscellstr(liste)  &  ( size(liste,2) == 5 ) )
  fprintf([nl 'DIRLIST: ' ...
              'LIST must be a 5-Column-CellstringArray,' nl ...
              ' from the 8. Column of Output from ' ...
              'DIRINFO with Option ''-l''.'  nl ]);
  liste = '';
  return
end

if Nin == 2
  if ~ischar(file)
    fprintf([nl 'DIRLIST: filename must be a String.' nl ]);
    file = '';
  else
    file(find(abs(file)==32)) = [];
    if isempty(file)
     fprintf([nl 'DIRLIST: filename is empty.' nl ]);
    end
  end
else
    file = '';
end


fprintf([ nl 'DIRLIST: Build FileList, Please wait ... ' ]);

     s2 = size(liste,2);

  liste = [ liste  cell(size(liste,1),5) ];
  liste(:,s2+1) = {'  '};
  liste(:,s2+2) = {' '};
  liste(:,s2+3) = {'  '};
  liste(:,s2+4) = {'    '};
  liste(:,s2+5) = { nl };

  liste = liste(:,[1 s2+1 2 s2+2 3 s2+3 4 s2+4 5 s2+5 ]);
 
  liste = liste';
  liste = cat(2,liste{:});

  fprintf([ nl nl ]);

if ~isempty(file)
  fprintf([ 'DIRLIST: Write LISTE into File: ' ...
            strrep(file,'\','\\')  ' ... '  ]);

  fid = fopen(file,'wt');
  if fid == -1
    fprintf([ 'Error, can''t open File.' nl ]);
  else

    l = strrep(liste,'\','\\');
    l = strrep(l,'%','%%');

    fprintf(fid,l);

    fclose(fid);

    fprintf(nl);
  end

  fprintf(nl);

end


if Nout == 0
  liste = '';
end