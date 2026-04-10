function SRM_RPP(CumInt,time)
figure;
Seq=1:length(time);
tiledlayout(1,1,"TileSpacing",'compact','Padding','compact');
nexttile;hold on;Fun_Decorat;
plot(time,CumInt,'LineWidth',1.9,'Color',[1 0.0745 0.65]); % Integral
plot(time,Seq,'LineWidth',1.5,'Color',[0 0 0],'LineStyle','-'); % Actual
xlabel('Time'); ylabel('Num.Events');grid;

end

