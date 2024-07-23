clc
clear variables

E = [0.1; 0.01; 0.001; 0.0001; 0.00001];
z = [0 0.00001 0.19101	0.38201	0.57301	0.76401	0.95501	1.14601	1.33701...
    1.52801	1.71901	1.91001	2.10101	2.29201	2.48301	2.67401	2.86501	...
    3.05601	3.24701	3.43801	3.62901	3.82001	4.01101	4.20201	4.39301	...
    4.58401	4.77501	4.96601];


%% Time for Simpson38
n_step = [500, 1000, 2000];
time_total_simpson38_J1 = zeros(1,length(n_step));
time_total_simpson38_J2 = zeros(1,length(n_step));
t_avg_simpson38_J1 = zeros(1,length(n_step));
t_avg_simpson38_J2 = zeros(1,length(n_step));

for nn = 1:length(n_step)
    for jj = 1:length(E)
        for kk = 1:length(z)
            [J1,J2,t_simp_SV_J1,t_simp_SV_J2] = simpson_J1_J2_new_SV(z(kk),E(jj),n_step(nn));
            time_total_simpson38_J1(nn) = time_total_simpson38_J1(nn) + t_simp_SV_J1;
            time_total_simpson38_J2(nn) = time_total_simpson38_J2(nn) + t_simp_SV_J2;
        end
    end
    t_avg_simpson38_J1(nn) = time_total_simpson38_J1(nn) / (length(E) * length(z));
    t_avg_simpson38_J2(nn) = time_total_simpson38_J2(nn) / (length(E) * length(z));
end



%% Time for Simpsondense
n_step = [500, 1000, 2000];
time_total_simpsondense_J1 = zeros(1,length(n_step));
time_total_simpsondense_J2 = zeros(1,length(n_step));

t_avg_simpsondense_J1 = zeros(1,length(n_step));
t_avg_simpsondense_J2 = zeros(1,length(n_step));

for nn = 1:length(n_step)
    for jj = 1:length(E)
        for kk = 1:length(z)
            [J1,J2,t_simp_SV_J1,t_simp_SV_J2] = simpson_dense(z(kk),E(jj),n_step(nn));
            time_total_simpsondense_J1(nn) = time_total_simpsondense_J1(nn) + t_simp_SV_J1;
            time_total_simpsondense_J2(nn) = time_total_simpsondense_J2(nn) + t_simp_SV_J2;
        end
    end
    t_avg_simpsondense_J1(nn) = time_total_simpsondense_J1(nn) / (length(E) * length(z));
    t_avg_simpsondense_J2(nn) = time_total_simpsondense_J1(nn) / (length(E) * length(z));
end

%% Time for Simpson J1
n_step = [500, 1000, 2000];

time_total_simpsonJ1_J1 = zeros(1,length(n_step));
time_total_simpsonJ1_J2 = zeros(1,length(n_step));

t_avg_simpsonJ1_J1 = zeros(1,length(n_step));
t_avg_simpsonJ1_J2 = zeros(1,length(n_step));

for nn = 1:length(n_step)
    for jj = 1:length(E)
        for kk = 1:length(z)
            [J1,J2,t_simp_SV_J1,t_simp_SV_J2] = simpson_J1_new(z(kk),E(jj),n_step(nn));
            time_total_simpsonJ1_J1(nn) = time_total_simpsonJ1_J1(nn) + t_simp_SV_J1;
            time_total_simpsonJ1_J2(nn) = time_total_simpsonJ1_J2(nn) + t_simp_SV_J2;
        end
    end
    t_avg_simpsonJ1_J1(nn) = time_total_simpsonJ1_J1(nn) / (length(E) * length(z));
    t_avg_simpsonJ1_J2(nn) = time_total_simpsonJ1_J2(nn) / (length(E) * length(z));
end

%% Time for Simpson J2
n_step = [500, 1000, 2000];

time_total_simpsonJ2_J1 = zeros(1,length(n_step));
time_total_simpsonJ2_J2 = zeros(1,length(n_step));

t_avg_simpsonJ2_J1 = zeros(1,length(n_step));
t_avg_simpsonJ2_J2 = zeros(1,length(n_step));

for nn = 1:length(n_step)
    for jj = 1:length(E)
        for kk = 1:length(z)
            [J1,J2,t_simp_SV_J1,t_simp_SV_J2] = simpson_J2_new(z(kk),E(jj),n_step(nn));
            time_total_simpsonJ2_J1(nn) = time_total_simpsonJ2_J1(nn) + t_simp_SV_J1;
            time_total_simpsonJ2_J2(nn) = time_total_simpsonJ2_J2(nn) + t_simp_SV_J2;
        end
    end
    t_avg_simpsonJ2_J1(nn) = time_total_simpsonJ2_J1(nn) / (length(E) * length(z));
    t_avg_simpsonJ2_J2(nn) = time_total_simpsonJ2_J2(nn) / (length(E) * length(z));
end

%% Time for Asymptotic
time_total_asymptotic_J1 = zeros(1,5);
time_total_asymptotic_J2 = zeros(1,5);

t_avg_asymptotic_J1 = zeros(1,5);
t_avg_asymptotic_J2 = zeros(1,5);

for nn = 1:5
    n_1 = nn;
    n_2 = nn;
    for jj = 1:length(E)
        for kk = 1:length(z)
            [J1,J2,time_J1,time_J2] = asymptotic_new_VS(z(kk),E(jj),n_1,n_2);
            time_total_asymptotic_J1(nn) = time_total_asymptotic_J1(nn) + time_J1;
            time_total_asymptotic_J2(nn) = time_total_asymptotic_J2(nn) + time_J2;
        end
    end
    t_avg_asymptotic_J1(nn) = time_total_asymptotic_J1(nn) / (length(E) * length(z));
    t_avg_asymptotic_J2(nn) = time_total_asymptotic_J2(nn) / (length(E) * length(z));
end

%% Time for quadrature_quadgk
AbsTol = [10.0e-8,10.0e-7,10.0e-6];
RelTol = [10.0e-7,10.0e-6,10.0e-5];

time_total_quadGK_J1 = zeros(length(AbsTol));
time_total_quadGK_J2 = zeros(length(AbsTol));

t_avg_quadGK_J1 = zeros(1,length(AbsTol));
t_avg_quadGK_J2 = zeros(1,length(AbsTol));

for nn = 1:length(AbsTol)
    abstol = AbsTol(nn);
    reltol = RelTol(nn);
    for jj = 1:length(E)
        for kk = 1:length(z)
            [J1,J2,time_J1,time_J2] = quadrature_quadgk_VS(z(kk),E(jj),abstol,reltol);
            time_total_quadGK_J1(nn) = time_total_quadGK_J1(nn) + time_J1;
            time_total_quadGK_J2(nn) = time_total_quadGK_J2(nn) + time_J2;
        end
    end
    t_avg_quadGK_J1(nn) = time_total_quadGK_J1(nn) / (length(E) * length(z));
    t_avg_quadGK_J2(nn) = time_total_quadGK_J2(nn) / (length(E) * length(z));

end

%% Time for quadrature_machine (exact)
%{
% time_total_quadm_J1 = 0;
time_total_quadm_J2 = 0;

AbsTol = [10.0e-15];
RelTol = [10.0e-14];

for nn = 1:length(AbsTol)
    abstol = AbsTol(nn);
    reltol = RelTol(nn);
    for jj = 1:length(E)
        for kk = 1:length(z)
            [J1,J2,time_J1,time_J2] = quadrature_quadgk_VS(z(kk),E(jj),abstol,reltol);
            time_total_quadm_J1 = time_total_quadm_J1 + time_J1;
            time_total_quadm_J2 = time_total_quadm_J2 + time_J2;
        end
    end
    t_avg_quadm_J1 = time_total_quadm_J1 / (length(E) * length(z));
    t_avg_quadm_J2 = time_total_quadm_J2 / (length(E) * length(z));
end
%}

%% Time for quadrature_machine (exact_1)
time_total_quadm1_J1 = 0;
time_total_quadm1_J2 = 0;

AbsTol = [10.0e-15];
RelTol = [10.0e-14];

for nn = 1:length(AbsTol)
    abstol = AbsTol(nn);
    reltol = RelTol(nn);
    for jj = 1:length(E)
        for kk = 1:length(z)
            [J1,J2,time_J1,time_J2] = quadrature_machine(z(kk),E(jj));
            time_total_quadm1_J1 = time_total_quadm1_J1 + time_J1;
            time_total_quadm1_J2 = time_total_quadm1_J2 + time_J2;
        end
    end
    t_avg_quadm1_J1 = time_total_quadm1_J1 / (length(E) * length(z));
    t_avg_quadm1_J2 = time_total_quadm1_J2 / (length(E) * length(z));
end


%% Plot of results

y = [t_avg_quadm1_J1 t_avg_quadm1_J2;
    t_avg_quadGK_J1(1) t_avg_quadGK_J2(1);
    t_avg_quadGK_J1(2) t_avg_quadGK_J2(2);
    t_avg_quadGK_J1(3) t_avg_quadGK_J2(3);
    t_avg_simpsonJ1_J1(1) t_avg_simpsonJ1_J2(1);
    t_avg_simpsonJ1_J1(2) t_avg_simpsonJ1_J2(2);
    t_avg_simpsonJ1_J1(3) t_avg_simpsonJ1_J2(3);
    t_avg_simpsonJ2_J1(1) t_avg_simpsonJ2_J2(1);
    t_avg_simpsonJ2_J1(2) t_avg_simpsonJ2_J2(2);
    t_avg_simpsonJ2_J1(3) t_avg_simpsonJ2_J2(3);
    t_avg_asymptotic_J1(1) t_avg_asymptotic_J2(1);
    t_avg_asymptotic_J1(2) t_avg_asymptotic_J2(2);
    t_avg_asymptotic_J1(3) t_avg_asymptotic_J2(3);
    t_avg_asymptotic_J1(4) t_avg_asymptotic_J2(4);
    t_avg_asymptotic_J1(5) t_avg_asymptotic_J2(5);];

y1 = [t_avg_asymptotic_J1(1) t_avg_asymptotic_J2(1);
    t_avg_asymptotic_J1(2) t_avg_asymptotic_J2(2);
    t_avg_asymptotic_J1(3) t_avg_asymptotic_J2(3);
    t_avg_asymptotic_J1(4) t_avg_asymptotic_J2(4);
    t_avg_asymptotic_J1(5) t_avg_asymptotic_J2(5);];

ylog = log10(y);

x = categorical({'GK Machine Precision','GK, tol=10e-8','GK, tol=10e-7','GK, tol=10e-6'...
    ,'CSimpson J1, 500','CSimpson J1, 1000','CSimpson J1, 2000','CSimpson J2, 500'...
    ,'CSimpson J2, 1000','CSimpson J2, 2000','Asymptotic 1-1','Asymptotic 2-2'...
    ,'Asymptotic 3-3','Asymptotic 4-4','Asymptotic 5-5'});
x = reordercats(x,{'GK Machine Precision','GK, tol=10e-8','GK, tol=10e-7','GK, tol=10e-6'...
    ,'CSimpson J1, 500','CSimpson J1, 1000','CSimpson J1, 2000','CSimpson J2, 500'...
    ,'CSimpson J2, 1000','CSimpson J2, 2000','Asymptotic 1-1','Asymptotic 2-2'...
    ,'Asymptotic 3-3','Asymptotic 4-4','Asymptotic 5-5'});

x1 = categorical({'Asymptotic 1-1','Asymptotic 2-2'...
     ,'Asymptotic 3-3','Asymptotic 4-4','Asymptotic 5-5'});
% x1 = reordercats(x,{'Asymptotic 1-1','Asymptotic 2-2'...
    % ,'Asymptotic 3-3','Asymptotic 4-4','Asymptotic 5-5'});

fig = figure('Units','inches');
fig.Position = [2, 3, 6.5, 4.2];

tiledlayout(1,4)
nexttile([1,3])

b = bar(x(1:end-5),y(1:end-5,:)*1000*1000, 'FaceColor','flat','BarWidth',1.3);
b(1).FaceColor = [40/255 40/255 40/255];
b(2).FaceColor = [180/255 180/255 180/255];

% xtips1 = b(1).XEndPoints - 0.3;
% ytips1 = b(1).YEndPoints + 0.2;
% labels1 = string(b(1).YData);
% h = text(xtips1, ytips1, labels1,'Rotation',90,'FontSize',8,...
%     'HorizontalAlignment','center','VerticalAlignment','top');

% xtips2 = b(2).XEndPoints - 0.3;
% ytips2 = b(2).YEndPoints + 0.2;
% labels2 = string(b(2).YData);
% text(xtips2, ytips2, labels2,'Rotation',90,'FontSize',8,'HorizontalAlignment','right')

% ylim([0,1e-3])
ylabel('Time $\mu s$','interpreter','latex','FontSize',11)
lgd = legend({'J1','J2'});
lgd.Box = 'off';
text(0,1.05,'(a)','FontSize',11,'Units','normalized')

nexttile(4)
b = bar(x1,y1*1000*1000, 'FaceColor','flat','BarWidth',1.3);
b(1).FaceColor = [40/255 40/255 40/255];
b(2).FaceColor = [180/255 180/255 180/255];

set(gca,'YAxisLocation','right')
ylabel('Time $\mu s$','interpreter','latex','FontSize',11)

text(0,1.05,'(b)','FontSize',11,'Units','normalized')

