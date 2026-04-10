function LSRM_RPP_Transformed(CumInt,reg)
figure;
R=max(reg);
tiledlayout(1,R,"TileSpacing",'compact','Padding','compact');
for i=1:R
    nexttile;hold on;Fun_defaultAxes;
    id=reg==i;
    plot(1:sum(id),CumInt{i},'LineWidth',1.5,'Color',[1 0.0745 0.65]);
    plot(1:sum(id),1:sum(id),'LineWidth',1.5,'Color',[0 0 0],'LineStyle','--');
    xlabel('Transformed Time'); ylabel('Num.Events');grid;
    axis('tight');
end

end