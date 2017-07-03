%% This function calculates first and second Einstein integrals (Einstein, 1950)
% Based on the method developed by
% Junke Guo and Pierre Y. Julien (2004); Journal of Hydraulic Engineering ASCE); for details look at the
% Appendix 1-2 of:
% Kaveh Zamani, Fabian Bombardelli and Babak Kamrani-Moghaddam (2016)
% "A comparison of current methods for the evaluation of Einstein’s integrals"
% Technical note in ASCE Journal of Hydraulic Engineering
%% Implemented by Kaveh Zamani at UC Davis, Department of Civil and Environmental Engg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of variables:
% Input
% Z: Rouse dimensionless number
% E: Relative-bedload layer thickness
% num_term: Number of terms in the partial sum of second infinite sum of
%           Equation (8) of Appendix 1.2 of Zamani et al. (2016)
% Output
% J1: First Einstein integral
% J2: Second Einstein integral

function [J1,J2]=guo(Z,E,num_term)

if(abs(Z-round(Z))>0.005)   
    F1_term2= 0;
    for k=1:10
        F1_term2= ((-1)^(k))/(k-Z)*((E/(1-E))^(k-Z)) + F1_term2;
        % Equation 9 in Zamani et al. (2016)
    end
    F1 = ((1-E)^Z)/(E^(Z-1)) - Z*F1_term2;
    J1 = Z*pi/sin(Z*pi) - F1 ;
    F2_term2 =0;
    for k=1:num_term
        ZZ=Z-k;
        FZmK =0;
        for i=1:num_term
            FZmK= ((-1)^(i))/(i-ZZ)*((E/(1-E))^(i-ZZ)) + FZmK;
        end
        phi_ZmK = ((1-E)^ZZ)/(E^(ZZ-1)) - ZZ*FZmK;
        F2_term2 = ((-1)^k)*phi_ZmK/(ZZ*(ZZ-1)) + F2_term2;
    end
    F2 = F1*(log(E)+1/(Z-1)) +Z*F2_term2;
    third_term = log(1+1.781*Z)-0.1361*Z/(1+1.284*Z)^2.15;
    J2 =Z*pi/sin(Z*pi)*...
        (pi*cot(Z*pi)- 1 - 1/Z + third_term)- F2;   
elseif (round(Z)==1)
    J1 = -1 - (log(E)-E);
    J2 = 1 - (E + log(E)^2/2 -E*log(E));
else
    n = round(Z);
    up_n = max(n-2,0);
    first_term = 0;
    for k=0:up_n
        first_term = first_term + ((-1)^k)*factorial(n)*(E^(k-n+1)-1)/(factorial(n-k)*factorial(k)*(n-k-1));
    end
    J1 = first_term + ((-1)^k)*(n*log(E)-E+1);
    first_term = 0;
    for k=0:up_n
        first_term = first_term + ((-1)^k)*factorial(n)/(factorial(n-k)*factorial(k))*...
            (E^(1+k-n)*log(E)/(n-k-1)+(E^(1+k-n)-1)/((n-k-1)^2));
    end
    J2 = first_term + ((-1)^n)*((n/2)*(log(E))^2 - E*log(E) + E - 1);
end % if

return
