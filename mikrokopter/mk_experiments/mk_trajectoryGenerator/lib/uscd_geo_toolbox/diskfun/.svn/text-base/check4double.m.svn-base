function [dd,txt] = check4double(pfad);

% CHECK4DOUBLE   Checks Directory recursively for duplicate Files
%
%  CHECK4DOUBLE( DIRECTORY )
%    
%   get recursively the DirectoryStructure from DIRECTORY, 
%    displays this Structure and checks for duplicate Files.
%
%   CHECK4DOUBLE without an Output will diplay the duplicate
%    Files in Matlab's CommandWindow.
%
%   [ DoubleList , Text ] = CHECK4DOUBLE( DIRECTORY )
% 
%  returns a CellArray DoubleList, containing SubCells with 
%    { Date  Bytes  FileName  Name }  of the duplicate Files  
%   and  the Text about that.
%
%   to get a String from DoubleList, use:
%  
%   >> List = strhcat( permute( cat(1,DoubleList{:}) , [ 2  1 ] ) , ...
%                      '   ' , size(DoubleList,2) );
%
%
% CHECK4DOUBLE  requires DIRINFO
%
%  see also:  DIRINFO  DIRLIST
%

dd  = {};
txt = '';

nl = char(10);

fs = filesep;

Nout = nargout;

is_txt = any( Nout == [ 0  2 ] );
fprintf([ nl ...
         'Get DirectoryStructure: ' ... 
         char([10 10])  ]);

 d = dirinfo(pfad,'-r','-l');

fprintf([ nl ...
         'Check for duplicate Files ... ' ]);

 % Get List

 d = d{8};


 % Get { Date  Bytes  FileNames }

 d = d( find(strcmp('F',d(:,1))) , [ 2  3  4 ] );

 N = size(d,1);

 d(:,4) = {''};  % FileName without Path

 for ii = 1 : N

    jj = findstr(d{ii,3},fs);
    if ~isempty(jj)
      d{ii,4} = d{ii,3}( max(jj)+size(fs,2) : end );
    else
      d{ii,4} = d{ii,3}; 
    end

 end


 ok = zeros(N,3);

 for ii = 1 : N

   if ~ok(ii,1)
         jj    = find( strcmp(d{ii,4},d(:,4)) );
      ok(jj,1) = ii;
      ok(jj,2) = length(jj);
      ok(ii,3) = length(jj);
   end

 end

 is_2 = find( ok(:,2) > 1 );
 
 ok = ok(is_2,:);
  d =  d(is_2,:);

 [ok(:,1),s_i] = sort(ok(:,1));

   ok(:,[2 3]) = ok(s_i,[2 3]);

             d = d(s_i,:);


 is_2 = find( ok(:,3) > 0 );
   n2 = size(is_2,1);

 if Nout >= 1
   dd = cell(n2,1);
 end


 for ii = 1 : n2

    jj = find( ok(:,1) == ok(is_2(ii),1) );
   
    if Nout >= 1
        dd{ii} = d(jj,:);
    end

    if is_txt

      txt1 = cell(size(jj,1),6);
      txt1(:,[1 3 5]) = d(jj,[1 2 3]);
      txt1(:,2)       = { char(32*ones(1,2)) };
      txt1(:,4)       = { char(32*ones(1,3)) };
      txt1(:,6)       = { nl           };

      txt1            = txt1';
      txt             = cat(2,txt,d{is_2(ii),4},nl,txt1{:},nl);
    end

 end

fprintf(['done' ... 
         char([10 10])  ]);

 if Nout == 0 
   fprintf([ strrep(txt,'\','\\') nl ]);
 end

