%% This function calculates first and second Einstein integrals
% Based on the method suggestion by Rajesh Srivastava (Abad et al., 2006)
% Discussionin Journal of Hydraulic Engineering-ASCE.
% For details please see Appendix 1-5 of
% Kaveh Zamani, Fabian Bombardelli and Babak Kamrani-Moghaddam (2016)
% "A comparison of current methods for the evaluation of Einstein’s integrals"
% Technical note in ASCE Journal of Hydraulic Engineering Vol. 143, Issue 4 
%% Implemented by Kaveh Zamani at UC Davis, Department of Civil and Environmental Engg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of variables:
% Input
% Z: Rouse dimensionless number
% E: Relative-bedload layer thickness
% Output
% J1: First Einstein integral
% J2: Second Einstein integral

function [J1,J2]=srivastava(Z,E)

Es = E/(1-E);
if (abs(Z-1)<0.01)
    
    j1term1 = -log(Es);
    j1term2 = 2.061*(Es^(2-Z)-1)/(2-Z);
    j1term3 = 1.385*(Es^(2.6-Z)-1)/(2.6-Z);
    j1term4 = 0.3327/(0.6703+Z);
    
    j2term1 = -(log(Es)^2)/2;
    j2term2 = -1.903*(Es^(2-Z)*(1-(2-Z)*log(Es))-1)/(2-Z)^2;
    j2term3 = 2.022*(E^(2.6-Z)*(1-(2.6-Z)*log(Es))-1)/(2.6-Z)^2;
    j2term4 = -0.2914/(1.652+Z);
    
elseif (abs(Z-2)<0.01)
    
    j1term1 = -(Es^(1-Z)-1)/(1-Z);
    j1term2 = 2.061*log(Es);
    j1term3 = 1.385*(Es^(2.6-Z)-1)/(2.6-Z);
    j1term4 = 0.3327/(0.6703+Z);
    
    j2term1 = (Es^(1-Z)*(1-(1-Z)*log(Es))-1)/(1-Z)^2;
    j2term2 = 1.903*(log(Es)^2)/2;
    j2term3 = 2.022*(E^(2.6-Z)*(1-(2.6-Z)*log(Es))-1)/(2.6-Z)^2;
    j2term4 = -0.2914/(1.652+Z);
    
elseif (abs(Z-2.6)<0.01)
    
    j1term1 = -(Es^(1-Z)-1)/(1-Z);
    j1term2 = 2.061*(Es^(2-Z)-1)/(2-Z);
    j1term3 = 1.385*log(Es);
    j1term4 = 0.3327/(0.6703+Z);
    
    j2term1 = (Es^(1-Z)*(1-(1-Z)*log(Es))-1)/(1-Z)^2;
    j2term2 = -1.903*(Es^(2-Z)*(1-(2-Z)*log(Es))-1)/(2-Z)^2;
    j2term3 = -2.022*(log(Es)^2)/2;
    j2term4 = -0.2914/(1.652+Z);
    
else
    
    j1term1 = -(Es^(1-Z)-1)/(1-Z);
    j1term2 = 2.061*(Es^(2-Z)-1)/(2-Z);
    j1term3 = -1.385*(Es^(2.6-Z)-1)/(2.6-Z);
    j1term4 = 0.3327/(0.6703+Z);
    
    j2term1 = ((Es^(1-Z))*(1-(1-Z)*log(Es))-1)/(1-Z)^2;
    j2term2 = -1.903*((Es^(2-Z))*(1-(2-Z)*log(Es))-1)/(2-Z)^2;
    j2term3 = 2.022*((E^(2.6-Z))*(1-(2.6-Z)*log(Es))-1)/(2.6-Z)^2;
    j2term4 = -0.2914/(1.652+Z);
    
end

J1 = j1term1 + j1term2 + j1term3 + j1term4;
J2 = j2term1 + j2term2 + j2term3 + j2term4;

return
