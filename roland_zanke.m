%% This function calculates first and second Einstein integrals (Einstein, 1950)
% Based on the method suggested by Roland and Zanke (Abad et al., 2006)
% Discussion in Journal of Hydraulic Engineering-ASCE.
% For details please see Appendix 1-6 of
% Kaveh Zamani, Fabian Bombardelli and Babak Kamrani-Moghaddam (2016)
% "A comparison of current methods for the evaluation of Einstein’s integrals"
% Technical note in ASCE Journal of Hydraulic Engineering
%% Implemented by Kaveh Zamani at UC Davis, Department of Civil and Environmental Engg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of variables:
% Input
% z: Rouse dimensionless number
% E: Relative-bedload layer thickness
% Output
% J1: First Einstein integral
% J2: Second Einstein integral

function [J1,J2]=roland_zanke(z,E)

J1= (1/(z-1))*(((1-E)^z)/E^(z-1))-...
    (z/(z-1))*((1/(z-2))*(((1-E)^(z-1))/E^(z-2))-...
    ((z-1)/(z-2))*((1/(z-3))*(((1-E)^(z-2))/E^(z-3))-...
    ((z-2)/(z-3))*((z-3)*pi/(sin((z-3)*pi))- E^(4-z)/(4-z))...
    ));
GMA = 0.5772156649015328606065121;
FZ = (1-GMA)- log(abs(4-z))+ 1/(3-z)+ 0.5/(4-z) + 1/(24*(4-z)^2);
J2z3 = (2-z)*pi*FZ/sin((z-2)*pi) - E^(3-z)*log(E)/(3-z)+ (E^(3-z))/((3-z)^2);
J1z2 = (1/(z-2))*(((1-E)^(z-1))/(E^(z-2)))...
    - ((z-1)/(z-2))*((z-2)*pi/sin((z-2)*pi) - (E^(3-z))/(3-z));
Croshe = (1/(z-2))*(log(E)*((1-E)^(z-1))/(E^(z-2))- (z-1)*J2z3*J1z2);
J2= (1/(z-1))*(log(E)*((1-E)^z)/(E^(z-1))- z*(Croshe)+J1);

return
