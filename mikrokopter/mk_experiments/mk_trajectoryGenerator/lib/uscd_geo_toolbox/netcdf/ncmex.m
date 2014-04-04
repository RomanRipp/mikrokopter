function varargout = ncmex(varargin)

% NCMEX    Driver for NetCDF C-Language interface.
%
%  NCMEX('action', ...) performs the specified NetCDF action.
%   Variables are returned as multi-dimensional arrays whose
%   dimensions are arranged in the left-to-right order defined
%   in 'vardef' or retrieved by 'varinq'.  No pre-put or post-get
%   permutation of dimensions is required.  The base-index for
%   slabs is zero (0), and -1 can be used to specify the remaining
%   count along any variable direction from the starting point.
%
% NCMEX('USAGE')
% 
%*******************************************************************
%HANDLE FILE
%
%  [cdfid, rcode ] = ncmex('CREATE', 'path',CMode )
%  [cdfid, status] = ncmex('OPEN'  , 'path', Mode );  
%
%   CMode:  'clobber' | 'noclobber' 
%    Mode:  'write'   | 'nowrite'
%
%
%   status        = ncmex('REDEF' , cdfid)
%   status        = ncmex('ENDEF' , cdfid)
%   status        = ncmex('SYNC'  , cdfid)
%   status        = ncmex('ABORT' , cdfid)
%   status        = ncmex('CLOSE' , cdfid)
%
%  [ndims, nvars, natts, recdim, status] = ncmex('INQUIRE', cdfid)
%
%
%
%*******************************************************************
%STORE DIMENSION
% 
%  status = ncmex('DIMDEF'   , cdfid, 'name', length)
%  status = ncmex('DIMRENAME', cdfid, 'name')
%
%NOTE:  Use   length = 'UNLIMITED'   for Unlimited Dimensions
%
%
%RETRIEVE DIMENSION
% 
%  [ dimid , rcode         ] = ncmex('DIMID' , cdfid, 'dname')
%  ['dname', length, status] = ncmex('DIMINQ', cdfid,  dimid )
%
%
%NOTE:  Use   'dname'  instead of  "dimid" to specify the Dimension
%
%*******************************************************************
%STORE VARIABLE
%
%  status = ncmex('VARDEF'   , cdfid, 'name', datatype, ndims, [dim])
%  status = ncmex('VARPUT'   , cdfid, varid , start, count, value, autoscale)
%  status = ncmex('VARPUT1'  , cdfid, varid , coords, value, autoscale)
%  status = ncmex('VARPUTG'  , cdfid, varid , start, count, stride, [], value, autoscale)
%  status = ncmex('VARRENAME', cdfid, varid, 'name')
%
%
%RETRIEVE VARIABLE
%
%  [ varid , rcode  ] = ncmex('VARID'  , cdfid, 'vname')
%  ['vname', datatype, ndims, [dim] , natts, ...
%    natts , status ] = ncmex('VARINQ' , cdfid, varid )
%
%  [ value , status ] = ncmex('VARGET' , cdfid, varid, start, count, autoscale)
%  [ value , status ] = ncmex('VARGET1', cdfid, varid, coords, autoscale)
%  [ value , status ] = ncmex('VARGETG', cdfid, varid, start, count, stride, [], autoscale)
%
%
%  DataTypes:  1. 'byte'      [        -128        127 ]      'int8'
%              2. 'char'      [           0        255 ]     'uint8'
%              3. 'short'     [      -32768      32767 ]     'int16'
%              4. 'long'      [ -2147483648 2147483647 ]     'int32'
%              5. 'float'       32 bit floating point
%              6. 'double'      64 bit floating point
% 
%
%
%NOTE:  Use   'vname'  instead of  "varid" to specify the Variable
%        
%
%*******************************************************************
%STORE ATTRIBUTE
%
%  status = ncmex('ATTPUT'   , cdfid, varid, 'aname', datatype, len, value) 
%  status = ncmex('ATTCOPY'  , incdf, invar, 'aname', outcdf, outvar)
%  status = ncmex('ATTRENAME', cdfid, varid, 'aname', 'newname')
%  status = ncmex('ATTDEL'   , cdfid, varid, 'aname')
%
%
%RETRIEVE ATTRIBUTE
%
%  [datatype, len, status] = ncmex('ATTINQ' , cdfid, varid, attid )
%  ['aname' ,      status] = ncmex('ATTNAME', cdfid, varid, attid )
%  [ value  ,      status] = ncmex('ATTGET' , cdfid, varid, attid )
%
% 
%NOTE:  Use   varid = 'global'  (-1)   for Global Attributes
%
%
%
%*******************************************************************
%STORE RECORD
%
%  status = ncmex('RECPUT', cdfid, recnum, [data], autoscale, recdim)
%
%
%RETRIEVE RECORD
%
%  [ [recvarids], [recsizes], status] = ncmex('RECINQ', cdfid, recdim)
%  [ [data]     ,             status] = ncmex('RECGET', cdfid, recnum, autoscale, recdim)
% 
%
%
%*******************************************************************
%OPTIONS
%
%  len          = ncmex('TYPELEN', datatype)
%  old_fillmode = ncmex('SETFILL', cdfid, fillmode)
% 
%  old_ncopts   = ncmex('SETOPTS', ncopts)
%  ncerr        = ncmex('ERR')
%  code         = ncmex('PARAMETER', 'NC_...')
%
%
%
%
%*******************************************************************
%NOTE
%
%   1. The rcode is always zero.
%   2. The dimid can be number or name.
%   3. The varid can be number or name.
%   4. The attid can be name or number.
%   5. The operation and parameter names are not case-sensitive.
%   6. The cmode defaults to 'NC_NOCLOBBER'.
%   7. The  mode defaults to 'NC_NOWRITE'.
%   8. The value -1 determines length automatically.
%   9. The operation names can prepend 'nc'.
%  10. The parameter names can drop 'NC_' prefix.
%  11. Dimensions: Matlab (i, j, ...) <==> [..., j, i] NetCDF.
%  12. Indices and Identifiers are ZERO-Based.
%  13. One-Dimensional slabs are returned as COLUMN-Vectors.
%  14. Scaling can be automated via 'scale_factor' and 'add_offset'.
%

% Copyright (C) 1992-1997 Dr. Charles R. Denham, ZYDECO.
% All Rights Reserved.

if nargin < 1, help ncmex, return, end

%------------------------------------------------------
% Looking for Mex-file Interface

fcn = '';

ff = { 'mexcdf'  'mexnc' };  % Basic name of MEX-File

ext = mexext;

frm = '%s.%s';

for f = ff
    if exist(sprintf(frm,f{1},ext),'file') == 3     % True for MexFile
       fcn = f{1};
       break
    end
end

%-----------------------------------------------------
if isempty(fcn)
%-----------------------------------------------------
% Check by Version- / Release- Number
%  Start with actual Setting, descending

   id = [ 00  71 70 65 61 60 53 7  6  5
          00  14 14 13 12 12 11 14 13 11 ];

   v = version;

   id(1,1) = eval(v([1 3]),'0');

   r = ( v == 'R' );
   if any(r)
      r       = min(find(r)) + [ 1  2 ];
      id(2,1) = eval(v(r),'0');
   end

   fm = '%s%.0f';

   for f = ff
       for ii = find( ~( id(:,1)' == 0 ) )
           for jj = find( id(ii,:) <= id(ii,1) )
                n = sprintf(fm,f{1},id(ii,jj));
                if exist(sprintf(frm,n,ext),'file') == 3   % True for MexFile
                   fcn = n;
                   break
                end
           end
           if ~isempty(fcn), break, end
       end
       if ~isempty(fcn), break, end
   end

   if isempty(fcn)
      ff = upper(ff);
      msg = sprintf('%s / %s  for Version %.0f / Release %.0f.', ...
                     ff{:},id(1,1),id(2,1));
      error(sprintf('Cann''t find MexCDF-Interface: .%s\n%s',upper(ext),msg));
   end

%-----------------------------------------------------
end
%------------------------------------------------------

% The "record" routines are emulated.

op = lower(varargin{1});
if any(findstr(op, 'rec'))
   fcn = op;
   if ~strcmp(fcn(1:2), 'nc')
      fcn = ['nc' fcn];
   end
   varargin{1} = [];
end

% Matlab-5 comma-list syntax.

if nargout > 0
   varargout = cell(1, nargout);
   [varargout{:}] = feval(fcn, varargin{:});
  else
   feval(fcn, varargin{:});
end
