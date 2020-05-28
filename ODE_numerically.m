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

%--------------- Numeric solution of the system of ODEs ------------------%

%% Parameters and constant values
Vdc = 640; Idc = 0; Rf = 5;

% --- Converter 1 ---
La1 = 3*21.2e-3; Ra1 = 3*0.2233; Cb1 = 4e-6; Ldc1 = 115e-3; Rdc1 = 10;

% --- Converter 2 ---
La2 = 3*28.266e-3; Ra2 = 3*0.295; Cb2 = 3e-6; Ldc2 = 180e-3; Rdc2 = 20;

% --- ODE parameters ---
% Time span to solve the system numerically [in seconds]
tspan = [0 10e-3];

% Initial values for variables (ic1, ibp1, ibn1, ic2, ibp2, ibn2, vbp1, vbn1, vbp2, vbn2)
y0 = [0 0 0 0 0 0 Vdc/2 Vdc/2 Vdc/2 Vdc/2]; 

%% Solving numerically the system of ODEs
[tode,yode] = ode45(@(t,y) odePGfault(t,y,La1,La2,Ra1,Ra2,Rdc1,Rdc2,Ldc1,Ldc2,Cb1,Cb2,Rf,Vdc), tspan, y0);

% Interpolate 
tint = linspace(0,tspan(2),501);
ic1 = interp1(tode, yode(:,1), tint);
ibp1 = interp1(tode, yode(:,2), tint);
ibn1 = interp1(tode, yode(:,3), tint);
ic2 = interp1(tode, yode(:,4), tint);
ibp2 = interp1(tode, yode(:,5), tint);
ibn2 = interp1(tode, yode(:,6), tint);
vbp1 = interp1(tode, yode(:,7), tint);
vbn1 = interp1(tode, yode(:,8), tint);
vbp2 = interp1(tode, yode(:,9), tint);
vbn2 = interp1(tode, yode(:,10), tint);

ifp1 = ic1+ibp1;
ifp2 = ic2+ibp2;
ifn1 = ic1+ibn1;
ifn2 = ic2+ibn2;
