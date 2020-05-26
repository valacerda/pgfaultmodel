%-------------------------------------------------------------------------%
% The authors will be thankful if the users of this code reference the work
% where this code was presented:
% V.A.Lacerda, R.M.Monaro, D.Campos-Gaona, R.Pena-Alzola, D.V.Coury 
% "Approximated Analytical Model of Pole-to-ground Faults in Symmetrical
% Monopole MMC-HVDC Systems" IEEE Journal on Emerging and Selected Topics
% in Power Electronics, 2020.

% This version was run on Matlab 9.8.0 (R2020a)
% Author: Vinícius A. Lacerda, University of São Paulo
% contact: vinicius.albernazlacerda@gmail.com
% Last version: 26 may 2020
%-------------------------------------------------------------------------%

%----- Build circuit ODEs using Kirchhoff's voltage and current laws -----%
% Hypotheses
 % 1 - The converter is considered a constant DC Voltage source
 % 2 - The converter internal voltage and DC voltage is considered balanced
 % 3 - Vdc is equal in both converters (Idc is added later)
% We define Ra = 2*Rarm/3 and La = 2*Larm/3 to reduce the expressions

% Variables and parameters
syms Vdc Rf Ra1 La1 Cb1 Cb2 Ra2 La2 Ldc1 Ldc2 Rdc1 Rdc2 positive
syms ic1 ibp1 ibn1 ifp1 ifn1 vbp1 vbn1 ic2 ibp2 ibn2 ifp2 ifn2 vbp2 vbn2
syms dic1 dibp1 dibn1 difp1 difn1 dic2 dibp2 dibn2 difp2 difn2

% KCLs
ifp1 = ic1 + ibp1;  difp1 = dic1 + dibp1;
ifn1 = ic1 + ibn1;  difn1 = dic1 + dibn1;
ifp2 = ic2 + ibp2;  difp2 = dic2 + dibp2;
ifn2 = ic2 + ibn2;  difn2 = dic2 + dibn2;
ift = ifp1 + ifp2;

% Capacitors' voltages
dvbp1 = -ibp1/Cb1;
dvbn1 = -ibn1/Cb1;
dvbp2 = -ibp2/Cb2;
dvbn2 = -ibn2/Cb2;

% KVLs
eqn1 = -Vdc + Ra1*ic1 + La1*dic1 + vbp1 + vbn1 == 0;
eqn2 = -vbp1 + Rdc1*ifp1 + Ldc1*difp1 + Rf*ift == 0;
eqn3 = -vbn1 + vbn2 + (Rdc1+Rdc2)*ifn1 + (Ldc1+Ldc2)*difn1 == 0;
eqn4 = vbp2 - Rf*ift - Rdc2*ifp2 - Ldc2*difp2 == 0;
eqn5 = Vdc - vbn2 - vbp2 - Ra2*ic2 - La2*dic2 == 0;
eqn6 = Ra1*ic1 + La1*dic1 + Rdc1*ifp1 + Ldc1*difp1 - Rdc2*ifp2 - Ldc2*difp2 - Ra2*ic2 - La2*dic2 + ifn1*(Rdc1+Rdc2) + difn1*(Ldc1+Ldc2) == 0;

% Eliminate ifs using KCLs
eqn1 = vbn1 - Vdc + vbp1 + La1*dic1 + Ra1*ic1 == 0;
eqn2 = Rf*(ibp1 + ibp2 + ic1 + ic2) - vbp1 + Ldc1*(dibp1 + dic1) + Rdc1*(ibp1 + ic1) == 0;
eqn3 = vbn2 - vbn1 + (dibn1 + dic1)*(Ldc1 + Ldc2) + (Rdc1 + Rdc2)*(ibn1 + ic1) == 0;
eqn4 = vbp2 - Rf*(ibp1 + ibp2 + ic1 + ic2) - Ldc2*(dibp2 + dic2) - Rdc2*(ibp2 + ic2) == 0;
eqn5 = Vdc - vbn2 - vbp2 - La2*dic2 - Ra2*ic2 == 0;
eqn6 = La1*dic1 - La2*dic2 + Ra1*ic1 - Ra2*ic2 + (dibn1 + dic1)*(Ldc1 + Ldc2) + (Rdc1 + Rdc2)*(ibn1 + ic1) + Ldc1*(dibp1 + dic1) - Ldc2*(dibp2 + dic2) + Rdc1*(ibp1 + ic1) - Rdc2*(ibp2 + ic2) == 0;

% Solve the system of KVLs
eqns = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6];
S = solve(eqns, [dic1 dic2 dibp1 dibn1 dibp2 dibn2]);

dic1 = -(vbn1 - Vdc + vbp1 + Ra1*ic1)/La1;
dibp1 = -(Ldc1*Vdc - La1*vbp1 - Ldc1*vbn1 - Ldc1*vbp1 + La1*Rdc1*ibp1 + La1*Rdc1*ic1 - Ldc1*Ra1*ic1 + La1*Rf*ibp1 + La1*Rf*ibp2 + La1*Rf*ic1 + La1*Rf*ic2)/(La1*Ldc1);
dibn1 = -(Ldc1*Vdc + Ldc2*Vdc - La1*vbn1 + La1*vbn2 - Ldc1*vbn1 - Ldc2*vbn1 - Ldc1*vbp1 - Ldc2*vbp1 + La1*Rdc1*ibn1 + La1*Rdc2*ibn1 + La1*Rdc1*ic1 - Ldc1*Ra1*ic1 + La1*Rdc2*ic1 - Ldc2*Ra1*ic1)/(La1*(Ldc1 + Ldc2));
dic2 = -(vbn2 - Vdc + vbp2 + Ra2*ic2)/La2;
dibp2 = -(Ldc2*Vdc - La2*vbp2 - Ldc2*vbn2 - Ldc2*vbp2 + La2*Rdc2*ibp2 + La2*Rdc2*ic2 - Ldc2*Ra2*ic2 + La2*Rf*ibp1 + La2*Rf*ibp2 + La2*Rf*ic1 + La2*Rf*ic2)/(La2*Ldc2);
dibn2 = -dic1 - dibn1 - dic2;
