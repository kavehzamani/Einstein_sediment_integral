%NOTE: THE n_step HAS TO BE MULTIPLE OF 3 TO MAKE IT WORK. 

function [J1,J2,time_J1,time_J2]=simpson_J1_J2_new_SV(rouse,E,n_step)
% J1 code Simpson 3/8
tic
a_b = (1-E)/n_step;
y_SV = E:a_b:1;

func_SV = @(y,rouse) ((1-y)./y).^rouse;
f = func_SV(y_SV,rouse);

Simpson3 = 0;
Simpson2 = 0;
for ii = 2:3:length(f)-1
    Simpson3 = Simpson3 + 3*(f(ii)+f(ii+1));
end
for jj = 4:3:length(f)-2
    Simpson2 = Simpson2 + 2*f(jj);
end

A = f(1);
B = Simpson3;
C = Simpson2;
D = f(end);

J1 = (3*a_b/8)*(A+B+C+D);

time_J1 = toc;
% J2 code Simpson 3/8
tic
func2_SV = @(y,rouse) ((1-y)./y).^rouse.*log(y);
f2 = func2_SV(y_SV,rouse);

Simpson2_3 = 0;
Simpson2_2 = 0;
for ii = 2:3:length(f2)-1
    Simpson2_3 = Simpson2_3 + 3*(f2(ii)+f2(ii+1));
end
for jj = 4:3:length(f2)-2
    Simpson2_2 = Simpson2_2 + 2*f2(jj);
end

A2 = f2(1);
B2 = Simpson2_3;
C2 = Simpson2_2;
D2 = f2(end);

J2 = -(3*a_b/8)*(A2+B2+C2+D2);
time_J2 = toc;
return