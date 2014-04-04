function c = rgb2grey(c)

% RGB2GREY   Converts between RGB- and GREY- Colors
%
% GREY = RGB2GREY(  RGB )
% 
%  RGB = RGB2GREY( GREY )
%
% RGB :  N by 3   , N Colors
% GREY:  N by 4   , [ Grey dR dG dB ]; 
%                     Grey = 0.299*R + 0.587*G + 0.114*B
% 
%

s2 = size(c,2);

if ~( isnumeric(c) & ( ndims(c) == 2 ) & ...
      any( s2 == [ 3  4 ] ) )
    error('Input must be a numeric Matrice with 3 or 4 Columns.');
end

if isempty(c)
   c = ones(0,4-(s2==4));
   return
end

if strcmp(class(c),'uint8')
   c = double(c) / 255;
end

sc = [ 0.299  0.587  0.114 ];

sk = ( 1 - sc );

s1 = size(c,1);

if s2 == 3
%  RGB --> GREY

   y = c * sc(:);
   c = c - y(:,[1 1 1]);
   c = 0.5 * ( c ./ sk(ones(1,s1),:)  + 1 );

   c = cat( 2 , y , c );

else
% GREY --> RGB

   c = ( 2*c(:,2:4) - 1 ) .* sk(ones(1,s1),:) + c(:,[1 1 1]);

end

 
   
