%% plot poi_sfr for each cell type for a given experiment
r=2; %cell type: 1=all; 2=grids; 3=conj; 4=noncj; 5=hd; 6=spatial; 7=theta; 8=int
s=1; %dataset: 1=8-oh inf; 2=8-oh sys; 3=diaz inf; 4=diaz sys; 5=dmso sys; 6=pbs inf; 7=pbs sys

figure; hold on
cle = [0 0 0; 1 0 0; .5 .5 .5; .75 .75 .75];
for c= 1:4
    inds = find(dat(s).crunched.cell(r,c).poi_sfr_p<0.05); %filter
    int = dat(s).crunched.cell(r,c).poi_sfr_inter(inds);
    slo = dat(s).crunched.cell(r,c).poi_sfr_slope(inds);

    sp = 0:2.5:30;
    clear lne
    for i = 1:length(slo)
        lne(i,:) = int(i) + slo(i).*sp;
    end
    plot(sp, nanmean(lne),'square','MarkerSize',8,'MarkerFaceColor',cle(c,:),'MarkerEdgeColor',cle(c,:));
    %plot(sp, nanmean(lne),'square','MarkerSize',13,'MarkerFaceColor','Color',cle(c,:),'LineWidth',20);
    errorbar(sp, nanmean(lne), nanstd(lne,[],1)/sqrt(size(lne,2)), 'Color',cle(c,:),'LineWidth',2)
    ylim([0 8]); xlim([-1 31])
end

title('Poisson speed versus firing rate','FontSize',18, 'Fontname', 'Times','FontWeight','Bold')
xlabel('Running Speed (cm/s)','FontSize',18, 'Fontname', 'Times','FontWeight','Bold')
ylabel('Firing rate (Hz)','FontSize',18, 'Fontname', 'Times','FontWeight','Bold')
set(gca,'fontsize',18, 'Fontname', 'Times','FontWeight','Bold')
axis square

%%
% r=2;
% close all
% figure; hold on
% for c = 1:4
%     inds = find(dat(1).crunched.cell(r,c).poi_sfr_p<0.05); %filter
%     int = dat(1).crunched.cell(r,c).poi_sfr_inter(inds);
%     slo = dat(1).crunched.cell(r,c).poi_sfr_slope(inds);
% 
%     sp = 0:5:40;
%     clear lne
%     
%     lam_int = nanmean(int);
%     lam_slo = nanmean(slo);
%     plot(sp, lam_int+sp*lam_slo,'Color',cle(c,:))
%     
% end