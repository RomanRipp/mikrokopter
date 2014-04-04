function c = cat(dim,varargin)

% CAT  Concatenate CYCL-arrays
%
% B = CAT( DIM , A1 , A2 , A3 , ... );
%

msg = '';
nl  = char(10);

Nin = nargin-1;

%*****************************************************************
% Check Inputs

%-----------------------------------------------------------------
% Check DIM

ok = ( strcmp(class(dim),'double') & ( prod(size(dim)) == 1 ) );
if ok
   ok = ( ( dim > 0 ) & ( mod(dim,1) == 0 ) );
end

if ~ok
    error('Value for DIM must be a single positive Integer.');
end

if Nin == 1
   c = varargin{1};
   return
end

%-----------------------------------------------------------------
% Check Class of Objects, get Size

si =  cell(Nin,1);
ns = zeros(Nin,1);

cl =  cell(Nin,1);   % Class of Values

for ii = 1 : Nin
    cl{ii} = class(varargin{ii});
    if strcmp(cl{ii},'cycl');
       varargin{ii} = varargin{ii}.VALUE;
             cl{ii} = class(varargin{ii});
    end
    si{ii} = size(varargin{ii});
    ns(ii) = prod(size(si{ii}));
end

if ~all(strcmp(cl{1},cl))
    error('Inputs must be of same Class.');
end

if ~isempty(msg)
    error(msg);
end

%-----------------------------------------------------------------
% Check Size

nn = max(max(ns),dim);

chk = ( ns < nn );
if any(chk)
   chk = find(chk);
   for ii = chk(:)'
       si{ii} = cat( 2 , si{ii} , ones(1,nn-ns(ii)) );
   end
end

si = cat(1,si{:});

chk        = si;
chk(:,dim) = 0;
chk        = diff(chk,1,1);

if ~all( chk(:) == 0 )
    error('CAT arguments dimensions are not consistent.');
end

%*****************************************************************
% Concatenate

ind = cell(1,nn);

ind(:)   = {':'};

cs = cumsum(cat(1,0,si(:,dim)));

ind(dim) = { cat( 2 , (1:si(1,dim)) , ones(1,cs(end)-si(1,dim))) };

c = varargin{1}(ind{:});

for ii = 2 : Nin

    ind(dim) = { ( 1 : si(ii,dim) ) + cs(ii) };

    c(ind{:}) = varargin{ii};

end

c = cycl(c);