function [msg,c,ok] = lowername(varargin)

% LOWERNAME  Renames Files to LowerCase-Names
%
% [Msg , C , OK ] = LOWERNAME(Pfad,'-r','.ext1','.ext2', ...)
%
%  C    = { Pfad OldName  NewName  OldFullFile  NewFullFile };
%
%  OK   = [  0    not converted, allready LowerCase
%            1        converted
%           -1    not converted, NewFullFile allready exist
%    
% [Msg , C ] = LOWERNAME( C )  reverts the Rename
%
%  note:  give '.'  or no Extension to list all Extensions
%
% required m-file:  FILELIST
%
% see also: DIRCONT, DIRINFO, DIRHIST, FILELIST
%

 
%--------------------------------------------------------

Nout = nargout;

msg = '';
  c =  cell(0,5);
 ok = zeros(0,1);

 c0 = c;

nl  = char(10);

%*****************************************
% Check for revert

revert = iscellstr(varargin{1});

if revert

   c = varargin{1};

   revert = ( ( size(c,2) == 5 ) & ...
              ( prod(size(c)) == size(c,1)*size(c,2) ) );

  if revert

      for ii = 1 : size(c,1)

        revert = ( revert & ...
                      strcmp( lower(c{ii,3}) , lower(c{ii,2}) ) & ...
                      strcmp( c{ii,4} , cat( 2 , c{ii,[1 2]} ) ) & ...
                      strcmp( c{ii,5} , cat( 2 , c{ii,[1 3]} ) )       );

      end

  end

  if revert

      c = c(:,[1 3 2 5 4]);

   else

      msg = 'First Input to Revert must be a 5-Column-CellStringArray by LOWERNAME.';

   end

else

  %*****************************************
  % Get FileList

  try

     [msg,c] = filelist(varargin{:});  % { Pfad  Name  Bytes  Date  DateNum }

  catch
   
      msg = [ 'Error call FILELIST.' nl lasterr ];

  end

  if ~isempty(msg)

    c = c0;

  end

end


%*****************************************
if ~isempty(msg)

 if ( Nout == 0 )

    fprintf(nl)
    fprintf(strrep(strrep(msg,'\','\\'),'%','%%'))
    fprintf([ nl  nl ]);

    msg = '';

  end

  return

end

if isempty(c)
   return
end
 
n = size(c,1);

%*****************************************
% Convert Names

n = size(c,1);

if ~revert

  for ii = 1 : n

    c{ii,3} = lower(c{ii,2});

    c{ii,4} = cat( 2 , c{ii,[1 2]} );   % Full OrgName
    c{ii,5} = cat( 2 , c{ii,[1 3]} );   % Full NewName

  end

end

ok = zeros(n,1);

for ii = 1 : n

     ok(ii) = ~strcmp( c{ii,2} , c{ii,3} );

     if ok(ii)

        ok(ii) = 1 - 2 * any( strcmp( c{ii,5} , c(:,4) ) );

     end

     if ok(ii) == 1

        fprintf('%s   %s ---> %s   ... ', c{ii,[1 2 3]} );
 
        command = [ 'mv '  c{ii,4}  '  ' c{ii,5} ];

        [s,w] = unix(command);

        if ~isempty(w) | ~( s == 0 )

           msg = cat( 2 , msg , nl(1:(end*(~isempty(msg)))) , ...
                     'Error call UNIX: ' , command ,  ...
                      nl(1:(end*(~isempty(w)))) , w , nl );

           fprintf('error')

        else

           fprintf('ok')

        end

        fprintf(nl)

     end

end

if ~isempty(msg) & ( Nout == 0 )

    fprintf(nl)
    fprintf(strrep(strrep(msg,'\','\\'),'%','%%'))
    fprintf([ nl  nl ]);

    msg = '';

end
