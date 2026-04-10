function LSRM_RPP(CumInt,time,reg)
figure;
R=max(reg);
tiledlayout(1,R,"TileSpacing",'compact','Padding','compact');
for i=1:R
    nexttile;hold on;Fun_defaultAxes;
    id=reg==i;
    plot(time,CumInt(i,:),'LineWidth',1.9,'Color',[1 0.0745 0.65]); % Integral
    plot(time(id),1:sum(id),'LineWidth',1.5,'Color',[0 0 0],'LineStyle','-'); % Actual
    xlabel('Time'); ylabel('Num.Events');grid;
end

end