%% This function calculates first Einstein integral
% Based on the method suggestion by Vantankhah and Velayati discussion on
% Journal of Hydraulic Engineering-ASCE.
% Kaveh Zamani, Fabian Bombardelli and Babak Kamrani-Moghadam (2016)
% "A comparison of current methods for the evaluation of Einstein’s integrals"
% Technical note in ASCE Journal of Hydraulic Engineering
%% Implemented by Kaveh Zamani at WRL UNSW Sydney, Department of Civil and Environmental Engg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of variables:
% Input
% Z: Rouse dimensionless number
% E: Relative-bedload layer thickness
% Output
% J1_8: First Einstein integral via equation 8
% J2_9: First Einstein integral via equation 9

function [J1_8,J1_9]=vatan(z,E)

Es = E/(1-E);

%%
term0 = 1 - Es^(1-z)/(1+Es);
term1 = - z*(Es^(1-z)/(1-z) + 1/(1+z));
term2 =   z*(Es^(2-z)/(2-z) + 1/(2+z));
term3 = - z*(Es^(3-z)/(3-z) + 1/(3+z));
term4 =   z*(Es^(4-z)/(4-z) + 1/(4+z));
term5 = - z*(Es^(5-z)/(5-z) + 1/(5+z));

J1_8 = term0+term1+term2+term3+term4+term5+0.591*z*(Es^(6-z)/(6-z)+1/(6+z));
J1_9 = term0 -z*Es^(1-z)/(1-z)+1.017*z*Es^(2-z)/(2-z)-0.595*z*(Es^(2.74-z)-1)/(2.74-z)-0.6*z/(0.84+z);

return