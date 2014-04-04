function str = remsurrd(str,cc,rm)

% RMSURRD  Remove surrounding Characters from String
%
%  STR = RMSURRD( String , CC )
%
%   CC is an 2-Row-CharacterArray for the characters to remove
%    from String if they match the Begin  and End of String.
%
%    default: CC = [ '|/\([{'     % Characters at Begin
%                    '|/\)]}' ];  % Characters at End
%
%  STR = RMSURRD( String , CC , BL )
%  
%   Removes from String following BlankCharacters from BL after
%    the removal of CC-Characters, like RMBLANK.
% 
%    default: BL = char([32 13 10 9 0]); % [ Space CR LF TAB ZERO ]
%
%  Use "-1" for CC or BL to use the DefaultSettings.
%
% see also: RMBLANK
%
% Example:
%
%  rmsurrd('{[1 2 3]}')          % returns '1 2 3'
%  rmsurrd('{[1 2 3]}','{}')     % returns '[1 2 3]'   % Only '{}' removed
%
%  rmsurrd('{[1 2 3] }')         % returns '1 2 3'     %    Blanks removed
%  rmsurrd('{[1 2 3] }',-1,'')   % returns '[1 2 3] '  % No Blanks removed
%  rmsurrd('{[1 2 3]*}',-1,'*')  % returns '1 2 3'     % '*' removed
% 
%

Nin = nargin;

%**********************************************************

if Nin < 1
   error('Input String is missing.');
end

if Nin < 2
   cc = -1;
end

if Nin < 3
   rm = -1;
end

%**********************************************************

msg = cell(0,1);

if ~( ischar(str) & ( prod(size(str)) == size(str,2) ) )
    msg = cat(1,msg,{'1. Input must be a String'});
end

if isequal(cc,-1)
   cc = cat( 1 , '|/\([{' , ...
                 '|/\)]}'       );
elseif ~( ischar(cc) & ( ndims(cc) == 2 ) )
    msg = cat(1,msg,{'2. Input must be a CharacterArray.'});
end

if isequal(rm,-1)
     rm = char([32 13 10 9 0]);
elseif ~( ischar(rm) | isempty(rm) )
    msg = cat(1,msg,{'3. Input must be a Characters'});
elseif ~isempty(rm)
     rm = rm(:)';
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%**********************************************************

if isempty(cc) | isempty(str)
   return
end

if size(cc,1) == 1
   if size(cc,2) == 2
      cc = cc';
   else
      cc = cc([1 1],:);
   end
end

%**********************************************************

n = size(str,2);

while 1

   if n <= 0
      break
   end

   if ~any( ( str(1) == cc(1,:) ) & ( str(n) == cc(2,:) ) )
       break
   end

   str = str( 2 : n-1 );

     n = n - 2;
   
   if ~isempty(rm)
       bl = zeros(1,n);
       for ii = rm
          bl = ( bl | ( str == ii ) );
       end
       if any(bl)
          bl = double(bl);
          i0 = 1 + sum(cumprod(bl,2),2);
          bl = bl(n:-1:1);
          i1 = n - sum(cumprod(bl,2),2);
         str = str(i0:i1);
           n = i1 - i0 + 1;
       end
   end

end