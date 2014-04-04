function [x,y,z] = read_stl(file);

% READ_STL  Reads Data from STL-Format (triangulated)
%
%     XYZ = READ_STL( File );   XYZ = [ 3 by NFaces by XYZ ]
%
% [X,Y,Z] = READ_STL( File ); X,Y,Z = [ 3 by NFaces ]
%
%
% Example for Visualisation:
%
%  patch( 'xdata' , X , 'ydata' , Y , 'zdata' , Z , , ...
%         'facecolor' , 'none' , 'edgecolor' , 'k' );
%
% See also: WRITE_STL, READRODB (required)
%


[m,h,x] = read_asc(file);

if ~isempty(m)
    error(m);
end

x = permute(x,[2 1]);

x = reshape(x,size(x,1),4,size(x,2)/4);

x(:,1,:) = [];           % [ XYZ by 3 by NFaces ]

x = permute(x,[2 3 1]);  % [ 3 by NFaces by XYZ ]

if nargout == 1
   return
end

z = x(:,:,3);
y = x(:,:,2);
x = x(:,:,1);
