%-------------------------------------------------------------------------%
% The authors will be thankful if the users of this code reference the work
% where this code was presented:
% V.A.Lacerda, R.M.Monaro, D.Campos-Gaona, R.Pena-Alzola, D.V.Coury 
% "An Approximated Analytical Model of Pole-to-ground Faults in Symmetrical
% Monopole MMC-HVDC Systems" IEEE Journal on Emerging and Selected Topics
% in Power Electronics, 2020.

% This version was run on Matlab 9.8.0 (R2020a)
% Author: Vinícius A. Lacerda, University of São Paulo
% contact: vinicius.albernazlacerda@gmail.com
% Last version: 26 may 2020
%-------------------------------------------------------------------------%

% ---------- Deriving the expressions for resonant frequencies -----------%


% Calculating the equivalent reactance "seen" by both Cb1s
% Neglecting system losses (resistances)
% Defining reactances:
syms Cb1 Cb2 La1 La2 Ldc1 Ldc2 w positive
Xa1 = w*La1;
Xa2 = w*La2;
Xb1 = -1/(w*Cb1);
Xb2 = -1/(w*Cb2);
Xdc1 = w*Ldc1;
Xdc2 = w*Ldc2;
Xya1 = (Xb1*Xa1)/(2*Xb1 + Xa1); % For Delta-Y transformation
Xyb1 = (Xb1*Xb1)/(2*Xb1 + Xa1); % For Delta-Y transformation
Xya2 = (Xb2*Xa2)/(2*Xb2 + Xa2); % For Delta-Y transformation
Xyb2 = (Xb2*Xb2)/(2*Xb2 + Xa2); % For Delta-Y transformation
Xp2 = (Xdc2+Xya2)*(Xyb2)/(Xdc2+Xya2+Xyb2);

% The Thevenin equivalent reactance:
Xth = Xa1*(Xp2 + 2*Xdc1 + Xdc2 + Xya2)/(Xa1 + Xp2 + 2*Xdc1 + Xdc2 + Xya2);
Xth = ((La1*La2 + 2*La1*Ldc1 + 2*La1*Ldc2)*w + (- 2*Cb2*La1*Ldc2^2 - 2*Cb2*La1*La2*Ldc1 - 2*Cb2*La1*La2*Ldc2 - 4*Cb2*La1*Ldc1*Ldc2)*w^3 + (Cb2^2*La1*La2*Ldc2^2 + 2*Cb2^2*La1*La2*Ldc1*Ldc2)*w^5)/((Cb2^2*La2*Ldc2^2 + 2*Cb2^2*La2*Ldc1*Ldc2 + Cb2^2*La1*La2*Ldc2)*w^4 + (- 2*Cb2*Ldc2^2 - 4*Cb2*Ldc1*Ldc2 - 2*Cb2*La2*Ldc2 - 2*Cb2*La2*Ldc1 - 2*Cb2*La1*Ldc2 - Cb2*La1*La2)*w^2 + 2*Ldc2 + 2*Ldc1 + La2 + La1);

% In resonance, the equivalent reactance is zero
sol = solve(Xth==0,w);

% This yields two solutions:
w' = ((Ldc2^2 - (Ldc2^4 + 4*Ldc1*Ldc2^3 + 4*Ldc1^2*Ldc2^2 + La2^2*Ldc1^2)^(1/2) + La2*Ldc1 + La2*Ldc2 + 2*Ldc1*Ldc2)/(Cb2*La2*Ldc2^2 + 2*Cb2*La2*Ldc1*Ldc2))^(1/2)
w'' = (((Ldc2^4 + 4*Ldc1*Ldc2^3 + 4*Ldc1^2*Ldc2^2 + La2^2*Ldc1^2)^(1/2) + Ldc2^2 + La2*Ldc1 + La2*Ldc2 + 2*Ldc1*Ldc2)/(Cb2*La2*Ldc2^2 + 2*Cb2*La2*Ldc1*Ldc2))^(1/2)
% where w' is an rough approximation of the low-frequencies
% and w'' is a good approximation of the high-frequency related to Cb1
% As all parameters of w'' depends on the circuit seen by Cb1, we will call:
% w2 = w''.
% Using the equivalent frequencies (please refer to ODEs_to_Laplace.m)
% and approximating (Ldc + Ldc1 to 2*Ldc), we can write w1 as:
w1 = sqrt(wa2 + 1/2*wdc2 + 1/2*w2dc2 + (1/4*(Ldc1^2/Ldc2^2)*wdc2^2 + wa2^2)^(1/2));

% The same procedure can be done for Cb2, obtaining other resonant frequency:
w2 = sqrt(wa1 + 1/2*wdc1 + 1/2*w1dc1 + (1/4*(Ldc2^2/Ldc1^2)*wdc1^2 + wa1^2)^(1/2));

% Now, using Vieta's formula, w3 and w4 can be obtained.
% By Vieta's formula we know that:
r1 + r2 + r3 + r4 = p8/p9;
r1 * r2 * r3 * r4 = p0/p9;
% where r1,r2,r3,r4 are the roots (w1^2, w2^2, w3^2, w4^2)
% and p0, p8, p9 are terms of P(s)
% As we already know r1 and r2:
r3 + r4 = p8/p9 - r1 - r2
r3 * r4 = p0/p9 * 1/r1 * 1/r2 
% Thus, r3 and r4 are found by joining the two equations and solving
% the second-order polynomial.
w3 = sqrt(-a1/2 -a2/2 + a3/2 + sqrt((a1 + a2 - a3)^2 - a4)/2);
w4 = sqrt(-a1/2 -a2/2 + a3/2 - sqrt((a1 + a2 - a3)^2 - a4)/2);
%where
% a1 = (((Ldc2/Ldc1)^2)/4 *wdc1^2 + wa1^2)^(1/2) + w1dc1/2 + wdc1/2  + wa1;
% a2 = (((Ldc1/Ldc2)^2)/4 *wdc2^2 + wa2^2)^(1/2) + w2dc2/2 + wdc2/2  + wa2;
% a3 = 2*wdc1 + 2*wdc2 + 2*wa1 + 2*wa2 + Ldc2/Ldc1*wdc1 + Ldc1/Ldc2*wdc2;
% a4 = 4*(2*wa1*wa2*(wdc1*w2dc2 + wdc2*w1dc1) + w1dc1*w2dc2*(wa1*wdc2 + wa2*wdc1))/(a1*a2);
