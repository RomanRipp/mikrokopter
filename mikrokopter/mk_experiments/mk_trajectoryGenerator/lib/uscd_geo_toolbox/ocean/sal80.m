function s = sal80(r,t,p)

% SAL80   Converts Conductivity to Salinity
%
% S = SAL80(C,T,P);
%
%-----------------------------------------------------------------------------
%     Calculates Salinity from Conductivity, Temperature and Pressure
%
%     Input parameters:
%            C  -  Conductivity [mS/cm]
%            T  -  Temperature  [IPTS-68]
%            P  -  Pressure     [dbar]
%
%     Output parameters:
%            S  -  Salinity  [PSU-78]
%
%
%     References:  UNESCO Report No. 37, 1981
%                  Practical Salinity Scale 1978: E.L. Lewis,IEEE Ocean
%                  Engineering, Jan., 1980
%
%---------------------------------------------------------------------

%
%  S. Chiswell 1991
%

      r  = r / 10;  % [S/m]

      r0 = 4.2914;
      tk = 0.0162;

      a  = [ 2.070e-05 -6.370e-10  3.989e-15 ];
      b  = [ 3.426e-02  4.464e-04  4.215e-01 -3.107e-3];

      aa = [ 0.0080 -0.1692 25.3851 14.0941 -7.0261  2.7081];
      bb = [ 0.0005 -0.0056 -0.0066 -0.0375  0.0636 -0.0144];

      cc = [ 6.766097e-01  2.00564e-02  1.104259e-04  -6.9698e-07  1.0031e-09];

      rt = cc(1) + cc(2)*t + cc(3)*t.*t + cc(4)*t.*t.*t + cc(5)*t.*t.*t.*t;

      rp = p .* (a(1) + a(2)*p + a(3)*p.*p);
      rp = 1 + rp./ (1 + b(1)*t + b(2)*t.*t + b(3)*r/r0 + b(4)*r/r0.*t);

      rt = r ./ (r0*rp.*rt);

      art = aa(1);
      brt = bb(1);
      for ii = 2 : 6
         rp  = rt.^((ii-1.0)/2.0);
         art = art + aa(ii).*rp;
         brt = brt + bb(ii).*rp;
      end

      rt = t - 15.0;

      s = art + (rt./(1+tk.*rt)).*brt;