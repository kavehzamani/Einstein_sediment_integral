%% This function calculates first and second Einstein integrals
% Based on the regression method developed by Jorge D. Abad and Marcelo H. Garcia (Abad et
% al. 2006; Journal of Hydraulic Engineering ASCE); for details look at the
% Appendix 1-3 of:
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

function [J1,J2]=abad(Z,E)

if (E <0.00001)
    error('Houston we have a problem! "E" is out of range!')
elseif(E < 0.015)
    c0=1.4852; c1=0.2025; c2=14.087; c3=20.918;  c4=-10.91; c5=2.034;  c6=-0.1345;
    b0=1.151;  b1=2.1787; b2=7.6572; b3=-0.2777; b4=-0.57;  b5=0.1424; b6=-0.0105;
elseif(E < 0.025)
    c0=1.2134; c1=1.9542; c2=10.613; c3=6.0002;  c4=-3.6259; c5=0.6938; c6=-0.0462;
    b0=1.1428; b1=2.4442; b2=4.2581; b3=-0.4713; b4=-0.1505; b5=0.0467; b6=-0.0036;
elseif(E < 0.035)
    c0=1.1409; c1=2.4266; c2=8.2541; c3=2.4058;  c4=-1.7617; c5=0.3474; c6=-0.0234;
    b0=1.1744; b1=2.4172; b2=3.0015; b3=-0.4405; b4=-0.049;  b5=0.0218; b6=-0.0018;
elseif(E < 0.045)
    c0=1.1138; c1=2.5982; c2=6.7187; c3=1.029;   c4=-1.001;  c5=0.2045; c6=0.0139;
    b0=1.2143; b1=2.364;  b2=2.3373; b3=-0.3955; b4=-0.0104; b5=0.0116; b6=-0.001;
elseif(E < 0.055)
    c0=1.1038; c1=0.6626; c2=5.6497; c3=0.3822;  c4=-0.6174; c5=0.1315; c6=-0.0091;
    b0=1.2574; b1=2.3159; b2=1.9239; b3=-0.3558; b4=0.0075;  b5=0.0064; b6=-0.0006;
elseif(E < 0.065)
    c0=1.102;  c1=2.6809; c2=4.864;  c3=0.0422;  c4=-0.3989; c5=0.0894; c6=-0.0063;
    b0=1.3023; b1=2.2773; b2=1.6411; b3=-0.3228; b4=0.0167;  b5=0.0035; b6=-0.0004;
elseif(E < 0.075)
    c0=1.1048; c1=2.6775; c2=4.2624; c3=-0.1487; c4=-0.2639; c5=0.0629; c6=-0.0045;
    b0=1.3486; b1=2.2481; b2=1.4351; b3=-0.2955; b4=0.0216;  b5=0.0017; b6=-0.0003;
elseif(E < 0.085)
    c0=1.1104; c1=2.6636; c2=3.787; c3=-0.2598; c4=-0.1757; c5=0.0454; c6=-0.0033;
    b0=1.3961; b1=2.2269; b2=1.2782; b3=-0.2728; b4=0.0243; b5=0.0005; b6=-0.0002;
elseif(E < 0.095)
    c0=1.1178; c1=2.6448; c2=3.4019; c3=-0.3254; c4=-0.1156; c5=0.0333;  c6=-0.0025;
    b0=1.445;  b1=2.2125; b2=1.1548; b3=-0.2536; b4=0.0258;  b5=-0.0002; b6=-0.0001;
elseif(E < 0.105)
    c0=1.1266; c1=2.6239; c2=3.0838; c3=-0.3636; c4=-0.0734; c5=0.0246; c6=-0.0019;
    b0=1.4952; b1=2.2041; b2=1.0552; b3=-0.2372; b4=0.0265; b5=-0.0008; b6=-0.00005;
else
    error('Houston we have a problem! E is out of range!')
end

J1 = 1/(c0+c1*Z+c2*Z*Z+c3*Z^3+c4*Z^4+c5*Z^5+c6*Z^6);
J2 = 1/(b0+b1*Z+b2*Z*Z+b3*Z^3+b4*Z^4+b5*Z^5+b6*Z^6);

J1 =  J1/(E/(1-E))^Z;
J2 = -J2/(E/(1-E))^Z;

return
