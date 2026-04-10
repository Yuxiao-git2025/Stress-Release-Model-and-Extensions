function SRM_PlotResults(lambda,Otime,Omag,Mmin)
global evalNums
figure;
nexttile;hold on; Fun_defaultAxes;xlim('tight');
evalpts=linspace(min(Otime),max(Otime),evalNums)';
plot(evalpts, lambda, 'LineWidth', 1.6,'Color','b','LineStyle','-');
% A plot for Poisson rate (Mean occurrence)
yline(length(Otime)/(max(Otime)-min(Otime)),'Color',[.7 .7 .7],'LineWidth',1.4);  
xlabel('Time'); ylabel('\lambda(t)');

yyaxis right;ax=gca;
stem(Otime,Omag,'LineWidth',1,'LineStyle','-','Marker','o','Color','k', ...
    'MarkerSize',10,'BaseValue',Mmin);
ylim([Mmin,round(max(Omag)+2)])
ax.YAxis(1).Color='b';
ax.YAxis(2).Color='k';
ylabel('Magnitude');
end