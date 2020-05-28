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

%---------- Description of system's ODEs for numeric solution ------------%

function dydt = odePGfault(~,y,La1,La2,Ra1,Ra2,Rdc1,Rdc2,Ldc1,Ldc2,Cb1,Cb2,Rf,Vdc)
% Variables
ic1 = y(1); ibp1 = y(2); ibn1 = y(3); ic2 = y(4); ibp2 = y(5);
ibn2 = y(6); vbp1 = y(7); vbn1 = y(8); vbp2 = y(9); vbn2 = y(10);

% System of ODEs
dydt = zeros(10,1);
dydt(1) = -(vbn1 - Vdc + vbp1 + Ra1*ic1)/La1;
dydt(2) = -(Ldc1*Vdc - La1*vbp1 - Ldc1*vbn1 - Ldc1*vbp1 + La1*Rdc1*ibp1 + La1*Rdc1*ic1 - Ldc1*Ra1*ic1 + La1*Rf*ibp1 + La1*Rf*ibp2 + La1*Rf*ic1 + La1*Rf*ic2)/(La1*Ldc1);
dydt(3) = -(Ldc1*Vdc + Ldc2*Vdc - La1*vbn1 + La1*vbn2 - Ldc1*vbn1 - Ldc2*vbn1 - Ldc1*vbp1 - Ldc2*vbp1 + La1*Rdc1*ibn1 + La1*Rdc2*ibn1 + La1*Rdc1*ic1 - Ldc1*Ra1*ic1 + La1*Rdc2*ic1 - Ldc2*Ra1*ic1)/(La1*(Ldc1 + Ldc2));
dydt(4) = -(vbn2 - Vdc + vbp2 + Ra2*ic2)/La2;
dydt(5) = -(Ldc2*Vdc - La2*vbp2 - Ldc2*vbn2 - Ldc2*vbp2 + La2*Rdc2*ibp2 + La2*Rdc2*ic2 - Ldc2*Ra2*ic2 + La2*Rf*ibp1 + La2*Rf*ibp2 + La2*Rf*ic1 + La2*Rf*ic2)/(La2*Ldc2);
dydt(6) = (vbn1 - Vdc + vbp1 + Ra1*ic1)/La1 + (Ldc1*Vdc + Ldc2*Vdc - La1*vbn1 + La1*vbn2 - Ldc1*vbn1 - Ldc2*vbn1 - Ldc1*vbp1 - Ldc2*vbp1 + La1*Rdc1*ibn1 + La1*Rdc2*ibn1 + La1*Rdc1*ic1 - Ldc1*Ra1*ic1 + La1*Rdc2*ic1 - Ldc2*Ra1*ic1)/(La1*(Ldc1 + Ldc2)) + (vbn2 - Vdc + vbp2 + Ra2*ic2)/La2;
dydt(7) = -ibp1/Cb1;
dydt(8) = -ibn1/Cb1;
dydt(9) = -ibp2/Cb2;
dydt(10) = -ibn2/Cb2;
end
