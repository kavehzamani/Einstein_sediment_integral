function [J1,J2,time1,time2]=simpson_J1_new(rouse,E,n_step)

d_big = (1-E)/n_step;
d_smal = d_big/3;
y_smal = E:d_smal:1;
y_big = E:d_big:1;

tic
fun1 = @(y,rouse)((1-y)./y).^rouse;
y1 = fun1(y_smal,rouse);
y_n_1 = fun1(y_big,rouse);
A = 2*sum(y1);
B = sum(y1(2:end-1));
C = sum(y_n_1);
J1= (A+B-C)*d_big/8;
time1=toc;
tic
J2= (2*sum(y1(1:end-1).*log(y_smal(1:end-1)))+sum(y1(2:end-1).*log(y_smal(2:end-1)))...
    - sum(y_n_1(1:end-1).*log(y_big(1:end-1))))*d_big/8;
time2=toc;
return