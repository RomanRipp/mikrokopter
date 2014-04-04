function x = meannan(x,dim)

% MEANNAN    Column mean with missing data.
%
%  Y = MEANNAN(X) returns the mean of each column of X as a row vector
%  where missing data values are encoded as NaNs. For vectors, MEANNAN(X)
%  returns the mean value of the elements in X.
%  Y = MEANNAN(X,DIM) averages over dimension DIM (only version 5 or higher)

%  This code has been suggested by Douglas M. Schwarz (schwarz@kodak.com) 
%  in the news group comp.soft-sys.matlab.
%
%  C. Mertens, IfM Kiel
%  $Revision: 1.1 $ $Date: 1995/03/08 14:27:31 $
%
% added ~isempty	G.Krahmann, IfM Kiel, Oct 1995
% added compatibility for version 5.	G.Krahmann, LODYC Paris, Jul 1997

% check for version

 
v=version;
if any(strcmp(v(1),{ '3' '4'}))

    dim = 0;

else

  if nargin==1
    dim = find( size(x) > 1 );
    if isempty(dim)
      dim = 1;
    else
      dim = dim(1);
    end
   end

end



      ok   = ~isnan(x);

  if dim == 0

      bad  = find( ~ok );
    x(bad) = 0*bad; 

    ok     = sum(ok);
    x      = sum(x) ./ ( ok + ( ok == 0 ) );

      bad  = find(ok==0);
    x(bad) = NaN*bad;

  else

    x(find(ok==0)) = 0;
 
    ok     = sum(ok,dim);
    x      = sum(x,dim) ./ ( ok + ( ok == 0 ) );

    x(find(ok==0)) = NaN;

  end


