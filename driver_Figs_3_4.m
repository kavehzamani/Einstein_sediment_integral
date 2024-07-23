clc; clear variables
close all;
format long

r1 = 0.1:0.1:0.9;
r2 = r1+1;
r3 = r2+1;
r4 = r3+1;
r5 = r4+1;
r6 = r5+1;
rouse = union(r1,union(r2,union(r3,union(r4,union(r5,r6)))));
E_rep = [0.0004,0.006,0.01,0.06];

n = length(rouse);
m = length(E_rep);
J1_exact = zeros(n,m);
J2_exact = zeros(n,m);

for i=1:n
   for j=1:m
   [J1_exact(i,j),J2_exact(i,j)] =quadrature_machine(rouse(i),E_rep(j));
   end
end

for i=1:n
   for j=1:m
   [J11(i,j),J21(i,j)] = asymptotic_new_VS(rouse(i),E_rep(j),1,1);
   [J12(i,j),J22(i,j)] = asymptotic_new_VS(rouse(i),E_rep(j),2,2);
   [J13(i,j),J23(i,j)] = asymptotic_new_VS(rouse(i),E_rep(j),3,3);
   [J14(i,j),J24(i,j)] = asymptotic_new_VS(rouse(i),E_rep(j),4,4);
   [J15(i,j),J25(i,j)] = asymptotic_new_VS(rouse(i),E_rep(j),5,5);
   end
end

err_J1_1 = 100*abs((J11-J1_exact)./J1_exact);
err_J2_1 = 100*abs((J21-J2_exact)./J2_exact);

err_J1_2 = 100*abs((J12-J1_exact)./J1_exact);
err_J2_2 = 100*abs((J22-J2_exact)./J2_exact);

err_J1_3 = 100*abs((J13-J1_exact)./J1_exact);
err_J2_3 = 100*abs((J23-J2_exact)./J2_exact);

err_J1_4 = 100*abs((J14-J1_exact)./J1_exact);
err_J2_4 = 100*abs((J24-J2_exact)./J2_exact);

err_J1_5 = 100*abs((J15-J1_exact)./J1_exact);
err_J2_5 = 100*abs((J25-J2_exact)./J2_exact);

%%
fontsize = 10;
numticks = 4;
fig1 = figure('Units','inches','Position',[0, 0.2, 6.5, 4.5]);
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','tight');
for i=1:m
    E= E_rep(i);
    
    ax1(i) = nexttile([1,1]); 
    semilogy(rouse,err_J1_1(:,i),'-','lineWidth',1.2,'markersize',2,'color',[0/256 0/256 0/256],'DisplayName','1 Term')
    hold on
    semilogy(rouse,err_J1_2(:,i),':','lineWidth',1.5,'color',[0/256 0/256 0/256],'DisplayName','2 Terms') 
    semilogy(rouse,err_J1_3(:,i),'--','lineWidth',1.2,'color',[90/256 90/256 90/256],'DisplayName','3 Terms')
    semilogy(rouse,err_J1_4(:,i),'-.','lineWidth',1.2,'markersize',2,'color',[90/256 90/256 90/256],'DisplayName','4 Terms')
    semilogy(rouse,err_J1_5(:,i),'-','lineWidth',1.2,'color',[180/256 180/256 180/256],'DisplayName','5 Terms')
    hold off
    labela = text(0.07,0.07,['E = ',num2str(E_rep(i))], 'Units','normalized','FontSize', fontsize);

    set(ax1(i),'FontSize',fontsize)

    if i == 1
        set(ax1(i),'xticklabel',[])
    elseif i == 2
        set(ax1(i),'xticklabel',[])      
    end

    % Ly = get(ax1(i),'YLim');
    % set(ax1(i),'XTick',linspace(Ly(1),Ly(2),numticks))


    % grid on
    % ax1(i).MinorGridLineWidth = 0.1;
    % ax1(i).MinorGridColor = [220/255, 220/255, 220/255];

end
ylabel(layout,'%','fontsize',fontsize)
xlabel(layout,'Rouse','fontsize',fontsize)
lgd = legend(ax1(i));
lgd.Box = 'off';
lgd.FontSize = 9;
lgd.ItemTokenSize(1) = 11;
lgd.ItemTokenSize(2) = 8;


%%
fontsize = 10;
fig1 = figure('Units','inches','Position',[0, 0.2, 6.5, 4.5]);
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','tight');
for i=1:m
    E= E_rep(i);
    
    ax1(i) = nexttile([1,1]); 
    semilogy(rouse,err_J2_1(:,i),'-','lineWidth',1.2,'markersize',2,'color',[0/256 0/256 0/256],'DisplayName','1 Term')
    hold on
    semilogy(rouse,err_J2_2(:,i),':','lineWidth',1.5,'color',[0/256 0/256 0/256],'DisplayName','2 Terms') 
    semilogy(rouse,err_J2_3(:,i),'--','lineWidth',1.2,'color',[90/256 90/256 90/256],'DisplayName','3 Terms')
    semilogy(rouse,err_J2_4(:,i),'-.','lineWidth',1.2,'markersize',2,'color',[90/256 90/256 90/256],'DisplayName','4 Terms')
    semilogy(rouse,err_J2_5(:,i),'-','lineWidth',1.2,'color',[180/256 180/256 180/256],'DisplayName','5 Terms')
    hold off
    labela = text(0.07,0.07,['E = ',num2str(E_rep(i))], 'Units','normalized','FontSize', fontsize);

    set(ax1(i),'FontSize',fontsize)

    if i == 1
        set(ax1(i),'xticklabel',[])
    elseif i == 2
        set(ax1(i),'xticklabel',[])    
    end  
end
ylabel(layout,'%','fontsize',fontsize)
xlabel(layout,'Rouse','fontsize',fontsize)
lgd = legend(ax1(i));
lgd.Box = 'off';
lgd.FontSize = 9;
lgd.ItemTokenSize(1) = 11;
lgd.ItemTokenSize(2) = 8;