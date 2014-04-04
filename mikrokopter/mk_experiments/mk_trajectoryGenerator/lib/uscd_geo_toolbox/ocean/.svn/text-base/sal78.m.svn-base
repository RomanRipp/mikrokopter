function X = sal78(P,T,CND,C1535,M)

% SAL78  Converts between Salinity and Conductivity
%
% SAL = SAL78( P , T68 , CND , C1535 , 0 )
% CND = SAL78( P , T68 , SAL , C1535 , 1 )
%
% FUNCTION SERVES FOR TWO PURPOSES (last Input):
%	1:	CONVERT CONDUCTIVITY TO SALINITY (0)
%	2:	CONVERT SALINITY TO CONDUCTIVITY (1)
%
%	ALGORITHMS RECOMMENDED BY JPOTS USING THE 1978
%	PRACTICAL SALINITY SCALE(IPSS-78) AND IPTS-68
%	FOR TEMPERATURE.
%
%       NOTE: The Conversion from IPTS-68 to ITS90 is:
%              T90 = 0.99976 * T68
%              T68 = 1.00024 * T90
%
% SAL78 COMPUTES EITHER 
%     THE CONDUCTIVITY RATIO    if C1535 = 1.0
%  or THE ABSOLUTE CONDUCTIVITY if C1535 = 42.914 = C(T=15,S=35)
%
%	UNITS:
%		PRESSURE        P         DBARS
%		TEMPERATURE     T         DEG.C.
%		SALINITY        S         NSU
%
%	RETURNS ZERO FOR CND < 0.0005  AND  last 0
%	RETURNS ZERO FOR SAL < 0.02    AND  last 1
%
% CHECKVALUES:
%
%           SAL78 = 1.888091
%       FOR   SAL =    40 NSU
%               T =    40 DEG
%               P =   10000 DBARS
%
%             CND = SAL78(10000,40,40,1,1) = 1.888091
%
%           SAL78 = 39.99999
%       FOR   CND = 1.888091
%               T =      40 DEG C.
%               P =   10000 DBARS
%
%             SAL = SAL78(10000,40,1.888091,1,0) = 39.999996
%


P = P/10 ;

if nargin == 3
	C1535 = 42.914 ;
	M = 0 ;
end

%ZERO SALINITY TRAP
if M == 0 
     zerocnd = find(CND < 5e-4) ; 
else
     zerocnd = find(CND < 0.02) ;
end

%SELECT BRANCH FOR SALINITY (M=0) OR CONDUCT.(M=1)
DT = T - 15.0 ;

if M == 0

     %CONVERT CONDUCTIVITY TO SALINITY
     R = CND/C1535 ;
     RT = R./(rt35(T).*(1.0 + c(P)./(b(T) + a(T).*R))) ;
     RT = sqrt(abs(RT)) ;
     %SALINITY RETURN
     X = sal(RT,DT) ;
     X(zerocnd) = zeros(size(zerocnd)) ;
else
     %CONVERT SALINITY TO CONDUCTIVITY
     %FIRST APPROXIMATION
     RT = sqrt(CND/35.0) ;
     SI = sal(RT,DT) ;
     [m,n] = size(CND) ;
     for i=1:m
          for j=1:n
               N = 0 ;
               DELS = 1 ;
               %ITERATE (MAX 10 ITERAT.) TO INVERT SAL POLYNOMIAL
               %FOR sqrt(RT)
               while (DELS > 1.e-4 & N < 10)
                    RT(i,j) = RT(i,j) + (CND(i,j) - SI(i,j))/dsal(RT(i,j),DT(i,j)) ;
                    SI(i,j) = sal(RT(i,j),DT(i,j)) ;
                    N = N + 1 ;
                    DELS = abs(SI(i,j) - CND(i,j)) ;
               end
          end
     end
     %COMPUTE CONDUCTIVITY RATIO 
     RTT = rt35(T) .* RT .* RT ;
     AT = a(T) ;
     BT = b(T) ;
     CP = c(P) ;
     CP = RTT.*(CP + BT) ;
     BT = BT - RTT.*AT ; 
     R = sqrt(abs(BT.*BT+4.0*AT.*CP)) - BT ;
     %CONDUCTIVITY RETURN
     X = 0.5*C1535*R./AT ;
     X(zerocnd) = zeros(size(zerocnd)) ;
end


%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function D = dsal(XR,XT)
%DSAL	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM, Kiel

	D=((((13.5405*XR-28.1044).*XR+42.2823).*XR+50.7702).*XR ...
                  -0.1692)+(XT./(1.0+0.0162*XT)).*((((-0.0720*XR+0.2544).*XR ...
                  -0.1125).*XR-0.0132).*XR-0.0056) ;

%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function rt = rt35(XT)
%RT35	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM Kiel

	rt=(((1.0031E-9*XT-6.9698E-7).*XT+1.104259E-4).*XT ...
             +2.00564E-2).*XT+0.6766097 ;


%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function S = sal(XR,XT)
%SAL	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM Kiel

	S = ((((2.7081*XR-7.0261).*XR+14.0941).*XR+25.3851).*XR ...
     	    -0.1692).*XR+0.0080+(XT./(1.0+0.0162*XT)).*(((((-0.0144*XR+ ...
     	     0.0636).*XR-0.0375).*XR-0.0066).*XR-0.0056).*XR+0.0005) ;



%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function x = a(XT)
%ACOEF	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM Kiel

x=-3.107E-3*XT+0.4215 ;

%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function x = b(XT)
%BCOEF	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM Kiel

x = (4.464E-4*XT+3.426E-2).*XT+1.0 ;

%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function x = c(XP)
%CCOEF	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM Kiel

x = ((3.989E-12*XP-6.370E-8).*XP+2.070E-4).*XP ;
