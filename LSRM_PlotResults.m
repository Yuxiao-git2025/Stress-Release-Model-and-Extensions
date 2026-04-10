function LSRM_PlotResults(lambdai,Region,Otime,Omag,Mmin)
global evalNums
if isempty(evalNums)
    evalNums=200;
end
evalTime=linspace(min(Otime),max(Otime),evalNums)';
Nplt=size(lambdai,1);
if Nplt>=5
    error('Please try small numbers of region or checking input');
end


figure;
col=[.1 .1 .1;0.0745    0.6235    1;1  0.4118  0.1608;
    1 0.0745 0.651;.6 .6 .6;];
tiledlayout(1,Nplt,"TileSpacing","compact","Padding","compact");
for i=1:Nplt
    nexttile;hold on; Fun_defaultAxes;
    xlim([min(Otime) max(Otime)]); ylim([min(lambdai(:)) max(lambdai(:))]);
    yyaxis left;
    lami=lambdai(i,:);
    idzero=lami~=0;
    idreg=Region==i;
    plot(evalTime(idzero),lami(idzero), 'LineWidth', 1.6,'Color',col(i,:),'LineStyle','-');
    % A plot of Poisson rate (Mean Occurrence)
    yline(sum(idreg)/(max(Otime(idreg))-min(Otime(idreg))),'Color',[.7 .7 .7],'LineWidth',1.4);
    xlabel('Time'); ylabel('\lambda(t)');


    yyaxis right;
    stem(Otime(idreg),Omag(idreg),'LineWidth',1,'LineStyle','-','Marker','o','Color',col(i,:), ...
        'MarkerSize',10,'BaseValue',Mmin);
    ax=gca;
    ax.YAxis(1).Color='k';
    ax.YAxis(2).Color='k';
    ylabel('Magnitude');
    ylim([Mmin,round(max(Omag)+2)]);
end
set(gcf,'position',[150,50,1400,600]);
end


