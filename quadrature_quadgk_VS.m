function [J1,J2,time_J1,time_J2]=quadrature_quadgk_VS(rouse,E,AbsTol,RelTol)
 tic
 fun1 = @(y,rouse)((1-y)./y).^rouse ;
 J1 = quadgk(@(y)fun1(y,rouse),E,1,'AbsTol',AbsTol,'RelTol',RelTol);
 time_J1 = toc;

 tic
 fun2 = @(y,rouse)(((1-y)./y).^rouse).*log(y) ;
 J2 = quadgk(@(y)fun2(y,rouse),E,1,'AbsTol',AbsTol,'RelTol',RelTol);
 time_J2 = toc;
return

