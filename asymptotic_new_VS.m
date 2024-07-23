function [J1,J2,time_J1,time_J2]=asymptotic_new_VS(Z,E,n_1,n_2)

% 0.6IN FRONT OF TERM2 of J1 is hardwired
tic
term1 = 0;
term2 = 0;
term3 = 0;
term4 = 0;
term5 = 0;
term6 = 0;

switch n_1
    case 1
        term1 = E/(Z-1);
    case 2
        term1 = E/(Z-1);
        term2 = Z*E*E/((Z-1)*(Z-2));
    case 3
        term1 = E/(Z-1);
        term2 = Z*E*E/((Z-1)*(Z-2));
        term3 = 2*Z*E^3/((Z-3)*(Z-2)*(Z-1)); 
    case 4
        term1 = E/(Z-1);
        term2 = Z*E*E/((Z-1)*(Z-2));
        term3 = 2*Z*E^3/((Z-3)*(Z-2)*(Z-1)); 
        term4 = 6*Z*E^4/((Z-4)*(Z-3)*(Z-2)*(Z-1));
    case 5
        term1 = E/(Z-1);
        term2 = Z*E*E/((Z-1)*(Z-2));
        term3 = 2*Z*E^3/((Z-3)*(Z-2)*(Z-1)); 
        term4 = 6*Z*E^4/((Z-4)*(Z-3)*(Z-2)*(Z-1));
        term5 = 24*Z*E^5/((Z-5)*(Z-4)*(Z-3)*(Z-2)*(Z-1));
    otherwise
        term1 = E/(Z-1);
        term2 = Z*E*E/((Z-1)*(Z-2));
        term3 = 2*Z*E^3/((Z-3)*(Z-2)*(Z-1)); 
        term4 = 6*Z*E^4/((Z-4)*(Z-3)*(Z-2)*(Z-1));
        term5 = 24*Z*E^5/((Z-5)*(Z-4)*(Z-3)*(Z-2)*(Z-1));
        term6 = 120*Z*E^6/((Z-6)*(Z-5)*(Z-4)*(Z-3)*(Z-2)*(Z-1));
end

J1 = -(1/E-1)^Z*(-term1+term2-term3+term4-term5+term6); 
time_J1 = toc;
%% todo: 0.7
tic
term1 = 0;
term2 = 0;
term3 = 0;
term4 = 0;
term5 = 0;
term6 = 0;
 
switch n_2
    case 1
        term1 = (-Z*log(E)+log(E)-1)*E/((Z-1)^2); 
    case 2
        term1 = (-Z*log(E)+log(E)-1)*E/((Z-1)^2); 
        term2 = (-3*Z+2*Z*Z+2*Z*log(E)-3*Z*Z*log(E)+Z^3*log(E))*E*E/((Z-1)^2*(Z-2)^2);
    case 3
        term1 = (-Z*log(E)+log(E)-1)*E/((Z-1)^2); 
        term2 = (-3*Z+2*Z*Z+2*Z*log(E)-3*Z*Z*log(E)+Z^3*log(E))*E*E/((Z-1)^2*(Z-2)^2);
        term3 = ((-16*Z+13*Z^2-Z^4+12*Z*log(E)-22*Z^2*log(E)+12*Z^3*log(E)-2*Z^4*log(E))* E^3)/((-3+Z)^2 *(-2+Z)^2 *(-1+Z)^2) ;
    case 4
        term1 = (-Z*log(E)+log(E)-1)*E/((Z-1)^2); 
        term2 = (-3*Z+2*Z*Z+2*Z*log(E)-3*Z*Z*log(E)+Z^3*log(E))*E*E/((Z-1)^2*(Z-2)^2);
        term3 = ((-16*Z+13*Z^2-Z^4+12*Z*log(E)-22*Z^2*log(E)+12*Z^3*log(E)-2*Z^4*log(E))* E^3)/((-3+Z)^2 *(-2+Z)^2 *(-1+Z)^2) ;
        term4 = ((-180*Z+170*Z^2-5*Z^3-26*Z^4+5*Z^5+144*Z*log(E)-300*Z^2*log(E)+210*Z^3*log(E)-60*Z^4*log(E)+6*Z^5*log(E))*E^4)/ ...
                ((-4+Z)^2*(-3+Z)^2*(-2+Z)^2*(-1+Z)^2);
    case 5
        term1 = (-Z*log(E)+log(E)-1)*E/((Z-1)^2); 
        term2 = (-3*Z+2*Z*Z+2*Z*log(E)-3*Z*Z*log(E)+Z^3*log(E))*E*E/((Z-1)^2*(Z-2)^2);
        term3 = ((-16*Z+13*Z^2-Z^4+12*Z*log(E)-22*Z^2*log(E)+12*Z^3*log(E)-2*Z^4*log(E))* E^3)/((-3+Z)^2 *(-2+Z)^2 *(-1+Z)^2) ;
        term4 = ((-180*Z+170*Z^2-5*Z^3-26*Z^4+5*Z^5+144*Z*log(E)-300*Z^2*log(E)+210*Z^3*log(E)-60*Z^4*log(E)+6*Z^5*log(E))*E^4)/ ...
                ((-4+Z)^2*(-3+Z)^2*(-2+Z)^2*(-1+Z)^2);
        term5 = (2 *(1728*Z - 1838*Z^2 + 135*Z^3 + 385*Z^4 - 135*Z^5 + 13*Z^6 - 1440*Z*log(E) ...
                 + 3288* Z^2 *log(E) - 2700*Z^3*log(E) + 1020*Z^4*log(E) - 180*Z^5*log(E) + 12*Z^6*log(E))*E^5)/ ...
                ((-5 + Z)^2*(-4 + Z)^2*(-3 + Z)^2*(-2 + Z)^2*(-1 + Z)^2);
%     otherwise

end

J2 = -((-1+1/E)^Z)*(term1+term2+term3+term4-term5-term6);
time_J2 = toc;
end