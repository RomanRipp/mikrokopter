function [msg,ind,n] = cyclind(si,ind,mode)

% CYCLIND  Returns cyclic Index
%
% [Msg,Index,N] = CYCLIND( Size , Index , Mode )
%
%  Size    Size of Object, Index refers too
%  Index   CellArray of IndexVectors to convert
%  N       Vector with Number of Elements per IndexVector
%

msg = '';

if nargin < 3
   mode = 0;
end

mode = isequal(mode,1);

ind = ind(:);

ni = size(ind,1);

if ni == 1

   si = prod(si(:));

else

   si = si(:);

   si = cat( 1 , si , ones(ni-size(si,1),1) );

end

n = si;

%-------------------------------------------------------------
% Check IndexVectors


for ii = 1 : ni

    if si(ii) == 0
       if mode & ~( isequal(ind{ii},':') & ischar(ind{ii}) )
           si(ii) = 1;
       else
          ind{ii} = [];
       end
    end

    if ~isempty(ind{ii})

%       if isequal(ind{ii},':') & ischar(ind{ii})
%          ind{ii} = ( 1 : si(ii) )';
%       else

        if ~( isequal(ind{ii},':') & ischar(ind{ii}) )

          ok = isnumeric(ind{ii});
          if ok & ~isempty(ind{ii})
             ok = all( isfinite(ind{ii}) );
          end

          if ~ok
             msg = 'Subscript indices must be finite numeric values.';
             return
          end

          ind{ii} = ind{ii} - si(ii)*floor(ind{ii}/si(ii));  % [ 0 .. si(ii)-1 ]

          ind{ii} = fix(ind{ii}) + si(ii) * ( ind{ii} == 0 );

            n(ii) = prod(size(ind{ii}));

        end

    end

end

%-------------------------------------------------------------
% Fill missing Dimensions

for ii = ni+1 : size(si,1)

    ind = cat( 2 , ind , {(1:si(ii))} );

end

