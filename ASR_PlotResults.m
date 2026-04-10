function ASR_PlotResults(EventTime,EventMag)
global Mmin Mmax
figure;
if length(EventMag)>=50
    fprintf('# bvalue=%.2f\n',cal_blvalue(EventMag,Mmin,0.1,Mmax,0,0.75));
end
tiledlayout(1,1,"TileSpacing","compact","Padding","compact");
nexttile;Fun_defaultAxes;hold on;
plot( EventTime/1e3 ,EventMag , '.k' , 'MarkerSize' , 8)
xlim([0,max(EventTime/1e3)]);
xlabel( 'Years (\times10^3)' );
ylabel( 'Magnitude' );
yl=ylim;
ylim([min(EventMag) yl(2)]);
set(gcf,'position',[300,50,900,580]);
end