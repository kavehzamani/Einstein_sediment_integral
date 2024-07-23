function [J1,J2,time1,time2]=quadrature_machine(rouse,E)

tic
fun1 = @(y,rouse)((1-y)./y).^rouse ;
J1 = integral(@(y)fun1(y,rouse),E,1,'AbsTol',1e-15,'RelTol',1e-14);
time1 = toc;

tic
fun2 = @(y,rouse)(((1-y)./y).^rouse).*log(y) ;
J2 = integral(@(y)fun2(y,rouse),E,1,'AbsTol',1e-15,'RelTol',1e-14);
time2 = toc;

return