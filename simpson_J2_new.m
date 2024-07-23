function [J1,J2,time1,time2]=simpson_J2_new(rouse,E,n_step)

d_big = (1-E)/n_step;
d_smal = d_big/3;
y_smal = E:d_smal:1;
y_big = E:d_big:1;

tic
fun2 = @(y,rouse)(((1-y)./y).^rouse).*log(y);

y2 = fun2(y_smal,rouse);
y_n_2 = fun2(y_big,rouse);

A = sum(y2);
B = sum(y2(2:end-1));
C = sum(y_n_2);
J2 =(2*A+B-C)*d_big/8;
time2=toc;
tic
J1 = (2*sum(y2(1:end-1)./log(y_smal(1:end-1)))+sum(y2(2:end-1)./log(y_smal(2:end-1)))...
    - sum(y_n_2(1:end-1)./log(y_big(1:end-1))))*d_big/8; 
time1=toc;
return