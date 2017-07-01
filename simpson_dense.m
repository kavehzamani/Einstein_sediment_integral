%% This function calculates first and second Einstein integrals
% Based on composite Simpson method (as known as cubic Simpson in some
% references)
% Appendix 2 of:
% Kaveh Zamani, Fabian Bombardelli and Babak Kamrani-Moghaddam (2016)
% "A comparison of current methods for the evaluation of Einstein’s integrals"
% Technical note in ASCE Journal of Hydraulic Engineering
%% Implemented by Kaveh Zamani at UC Davis, Department of Civil and Environmental Engg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of variables:
% Input
% Z: Rouse dimensionless number
% E: Relative-bedload layer thickness
% n_step: number of nodes for numerical computation of the integral 
% Output
% J1: First Einstein integral
% J2: Second Einstein integral

function [J1,J2]=simpson_dense(rouse,E,n_step)

if (n_step<4000)
    n_step = 4000;
    disp('Houston we have a problem! number of steps was not enough and set to 4000')
end

d_big = (1-E)/n_step;
y_big = E:d_big:1;
d_smal = d_big/3;
y_smal = E:d_smal:1;

fun1 = @(y,rouse)((1-y)./y).^rouse;
fun2 = @(y,rouse)(((1-y)./y).^rouse).*log(y);

y1 = fun1(y_smal,rouse);
y_n_1 = fun1(y_big,rouse);

y2 = fun2(y_smal,rouse);
y_n_2 = fun2(y_big,rouse);

J1= (2*sum(y1)+sum(y1(2:end-1))-sum(y_n_1))*d_big/8;
J2= (2*sum(y2)+sum(y2(2:end-1))-sum(y_n_2))*d_big/8;

return

