function SRM_RPP_Transformed(CumInt,time)
figure;
Seq=1:length(time);
tiledlayout(1,1,"TileSpacing",'compact','Padding','compact');
nexttile;hold on;Fun_Decorat;
plot(Seq,CumInt,'LineWidth',1.9,'Color',[1 0.0745 0.65]); % Integral
plot(CumInt,CumInt,'LineWidth',1.5,'Color',[0 0 0],'LineStyle','-'); % Actual
xlabel('Transformed Time'); ylabel('Num.Events');grid;
axis('tight');

end

