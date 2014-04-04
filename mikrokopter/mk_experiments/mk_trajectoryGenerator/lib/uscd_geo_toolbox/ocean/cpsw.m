function cp = cpsw(p,t,s);

% CPSW specific heat of sea-water at constant pressure
%
%      as described in Landolt-Boernstein, New Series V/3 a, p. 247-250.
%
%        cp = CPSW(P,T,S);
%
%        with matrices  S  : salinity in psu
%                       T  : temperature in deg C
%                       P  : pressure in dbar
%                       cp : specific heat at constant pressure in J/(K*kg)
%
%        check value :  cpsw(    0,40.0,40) = 3980.1 J/(K*kg)
%                       cpsw(10000,40.0,40) = 3849.5 J/(K*kg)

%        by Ulf Garternicht in August 1994, IfM Kiel
% changed input order 	G.Krahmann, Sep 1995

% init some variables

s15=s.^(3/2);

t2=t.^2;
t3=t.^3;
t4=t.^4;

p2=p.^2;
p3=p.^3;

% compute cp(s,t,0)

a = [4.2174e3,-3.720283,0.1412855,-2.654387e-3,2.093236e-5];
b = [-7.643575,0.1072763,-1.38385e-3];
c = [0.1770383,-4.07718e-3,5.148e-5];

cpst0 =           a(1)+a(2)*t+a(3)*t2+a(4)*t3+a(5)*t4 ...
        + s .*   (b(1)+b(2)*t+b(3)*t2) ...
        + s15 .* (c(1)+c(2)*t+c(3)*t2);

% compute cp(s,t,p)

a = [-4.9592e-2,1.45747e-3,-3.13885e-5,2.0357e-7,1.7168e-9];
b = [2.4931e-6,-1.08645e-7,2.87533e-9,-4.0027e-11,2.2956e-13];
c = [-5.422e-11,2.6380e-12,-6.5637e-14,6.136e-16];
d = [4.9247e-4,-1.28315e-5,9.802e-8,2.5941e-9,-2.9179e-11];
e = [-1.2331e-5,-1.517e-7,3.122e-9];
f = [-2.9558e-8,1.17054e-9,-2.3905e-11,1.8448e-13];
g = [9.971e-10];
h = [5.540e-13,-1.7682e-14,3.513e-16];
k = [-1.4300e-15];

cp = cpst0 + p         .* (a(1)+a(2)*t+a(3)*t2+a(4)*t3+a(5)*t4) ...
           + p2        .* (b(1)+b(2)*t+b(3)*t2+b(4)*t3+b(5)*t4) ...
           + p3        .* (c(1)+c(2)*t+c(3)*t2+c(4)*t3) ...
           + s   .* p  .* (d(1)+d(2)*t+d(3)*t2+d(4)*t3+d(5)*t4) ...
           + s15 .* p  .* (e(1)+e(2)*t+e(3)*t2) ...
           + s   .* p2 .* (f(1)+f(2)*t+f(3)*t2+f(4)*t3) ...
           + s15 .* p2  * g ...
           + s   .* p3 .* (h(1)+h(2)*t+h(3)*t2) ...
           + k * s15 .* p3 .* t;
