function cdp(sub)

% CDP  Change to Subdirectory of Matlab/Project
%
% CDP( SubDirectory )
%

pfad = '/d1/project/matlab/';

if nargin < 1
   sub = '';
elseif isempty(sub)
   sub = '';
elseif ~( ischar(sub) & ( prod(size(sub)) == size(sub,2) ) )
   error('SubDirectory must be a String.')
end

p = fullfile(pfad,sub);

if ~( exist(p,'dir') == 7 )
    error(['Directory doesn''t exist: ' p ]);
end

cd(p)
