%% This function calculates first and second Einstein integrals (Einstein, 1950)
% Based on the method suggested by Tatsuaki Nakato (Nakato, 1984)
% Journal of Hydraulic Engineering-ASCE.
% For details please see Appendix 1-1 of
% Kaveh Zamani, Fabian Bombardelli and Babak Kamrani-Moghaddam (2016)
% "A comparison of current methods for the evaluation of Einstein’s integrals"
% Technical note in ASCE Journal of Hydraulic Engineering
%% Implemented by Kaveh Zamani at UC Davis, Department of Civil and Environmental Engg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of variables:
% Input
% Z: Rouse dimensionless number
% E: Relative-bedload layer thickness
% Output
% J1: First Einstein integral
% J2: Second Einstein integral
function [J1,J2]=nakato(Z,E)

A =  E;
E = 0.5;

if (abs(Z-1)<0.01)
    F1 = log(E)-log(A);
    G1 = (log(E)^2-log(A)^2)/2;
    F2 = Z/(Z-2)*(E^(2-Z)-A^(2-Z));
    G2 = Z*E^(2-Z)/(Z-2)*(log(E)-1/(2-Z))-Z*A^(2-Z)/(Z-2)*(log(A)-1/(2-Z));
    F3 = Z*(Z-1)/2/(3-Z)*(E^(3-Z)-A^(3-Z));
    G3 = Z*(Z-1)*E^(3-Z)/2/(3-Z)*(log(E)-1/(3-Z))-Z*(Z-1)*A^(3-Z)/2/(3-Z)*(log(A)-1/(3-Z));
elseif (abs(Z-2)<0.01)
    F1 = 1/(1-Z)*(E^(1-Z)-A^(1-Z));
    G1 = E^(1-Z)/(1-Z)*(log(E)-1/(1-Z))-A^(1-Z)/(1-Z)*(log(A)-1/(1-Z));
    F2 = -2*(log(E)-log(A));
    G2 = -(log(E)^2)+log(A)^2;
    F3 = Z*(Z-1)/2/(3-Z)*(E^(3-Z)-A^(3-Z));
    G3 = Z*(Z-1)*E^(3-Z)/2/(3-Z)*(log(E)-1/(3-Z))-Z*(Z-1)*A^(3-Z)/2/(3-Z)*(log(A)-1/(3-Z));
elseif (abs(Z-3)<0.01)
    F1 = 1/(1-Z)*(E^(1-Z)-A^(1-Z));
    G1 = E^(1-Z)/(1-Z)*(log(E)-1/(1-Z))-A^(1-Z)/(1-Z)*(log(A)-1/(1-Z));
    F2 = Z/(Z-2)*(E^(2-Z)-A^(2-Z));
    G2 = Z*E^(2-Z)/(Z-2)*(log(E)-1/(2-Z))-Z*A^(2-Z)/(Z-2)*(log(A)-1/(2-Z));
    F3 = 3*(log(E)-log(A));
    G3 = 3*(log(E)^2-log(A)^2)/2;
else
    F1 = 1/(1-Z)*(E^(1-Z)-A^(1-Z));
    G1 = E^(1-Z)/(1-Z)*(log(E)-1/(1-Z))-A^(1-Z)/(1-Z)*(log(A)-1/(1-Z));
    F2 = Z/(Z-2)*(E^(2-Z)-A^(2-Z));
    G2 = Z*E^(2-Z)/(Z-2)*(log(E)-1/(2-Z))-Z*A^(2-Z)/(Z-2)*(log(A)-1/(2-Z));
    F3 = Z*(Z-1)/2/(3-Z)*(E^(3-Z)-A^(3-Z));
    G3 = Z*(Z-1)*E^(3-Z)/2/(3-Z)*(log(E)-1/(3-Z))-Z*(Z-1)*A^(3-Z)/2/(3-Z)*(log(A)-1/(3-Z));
end

fun1 = @(y,rouse)((1-y)./y).^Z ;
F4 = quadgk(@(y)fun1(y,Z),E,1);
fun2 = @(y,rouse)(((1-y)./y).^Z).*log(y) ;
G4 = quadgk(@(y)fun2(y,Z),E,1);
J1 = F1+F2+F3+F4;
J2 = G1+G2+G3+G4;

return
