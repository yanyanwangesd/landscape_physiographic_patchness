
% this script uses the outputs of the landscape_physiographic_patchness.m
% to plot the statistics of multiple types of patch changes

% Author: Dr. Yanyan Wang (yanyan.wang@erdw.ethz.ch, or wangyanyan0607@hotmail.com)
% Last update: Aug. 30th, 2023


%%
types = [111   112   113   121   122   123   211   212   213   221   222   223   311   312   313   321   322   323  411 412 413 421 422 423];
type_merg = [1 2 3 1 5 6 7 8 9 7 11 12 13 14 15 13 17 18 19 20 21 19 23 24]; % not distinguish aspect for shallow slope

type_plot = unique(type_merg); % not distinguish aspect for shallow slope

mstart = 30400;
printstep = 5;
dt = 2e4;
stepplot = 1:1:520;

% plot patche change events, every 1 Myr
time = stepplot*dt*printstep/1e6; % convert model time to Myr
dtplot = 1; % intervals of model time to plot, Myr
bint = min(time):dtplot:max(time);

% group events into four elevation bands
bandtypes = 5; % the number of habitat types in each elevation band

born = zeros(length(bint)-1,4);
split = zeros(length(bint)-1,4);
merge = zeros(length(bint)-1,4);
migrate = zeros(length(bint)-1,4);
die = zeros(length(bint)-1,4);
npatch = zeros(length(bint)-1,4);
ninbin = zeros(length(bint)-1,4);

for k = 0:3
        
    for j = k*6+1:(k+1)*6
        if ismember(j, type_plot)
            
            dir = strcat('/landscape_physiographic_patchness/habitat_',num2str(targetH));            
            filename = strcat(dir,'/patch_evolution_statistics_habitat_',num2str(targetH),'.mat');
            load(filename)
                        
            for it = 1:length(bint)-1
                idt = time>bint(it)&time<=bint(it+1);
                born(it,k+1) = born(it,k+1)+sum(Nborn(idt));
                split(it,k+1) = split(it,k+1)+ sum(Nsplit(idt));
                merge(it,k+1) = merge(it,k+1)+sum(Nmerge(idt));
                migrate(it,k+1) = migrate(it,k+1)+ sum(Nmigrate(idt));
                die(it,k+1) = die(it,k+1)+sum(Ndie2(idt));
                npatch(it,k+1) = npatch(it,k+1)+sum(Npatch(idt));
                ninbin(it,k+1) = ninbin(it,k+1)+sum(idt);
            end
        end
        
    end

end

Npatchmean = npatch./ninbin*bandtypes;

%% plot figures

% 1) plot histograms of events
figure
for k=0:3
    subplot(4,1,k+1)
    by = horzcat(born(:,k+1),merge(:,k+1),split(:,k+1),die(:,k+1));
    %     by = by./Npatchmean(:,k+1)/10;
    hb1 = bar(bint(1:end-1)+0.5, by,'stacked');
    hb1(1).BarWidth = 1.;
    hb1(1).LineWidth = 0.1;
    hb1(2).LineWidth = 0.1;
    hb1(3).LineWidth = 0.1;
    hb1(4).LineWidth = 0.1;
    hb1(1).FaceColor = [178,171,210]/255;
    hb1(2).FaceColor = [230,97,1]/255;
    hb1(3).FaceColor = [253,184,99]/255;
    hb1(4).FaceColor = [94,60,153]/255;
    hold on
    hb3 = stairs(bint(1:end), [Npatchmean(:,k+1);Npatchmean(end,k+1)],'k-','LineWidth',3);
    legend('birth','merge','split','death','Mean patch number','Location','best','Orientation','Horizontal')
    xlabel('Model time (Myr)')
    ylabel('Number of event per Myr')
    if k==0
        title('elevation band 0-500 m')
    end
    if k==1
        title('elevation band 500-1500 m')
    end
    if k==2
        title('elevation band 1500-2500 m')
    end
    if k==3
        title('elevation band 2500-3000 m')
    end    
    set(gca,'FontSize',18)
    
end


% 2) plot recurrence time of events
figure
for k=0:3
    
    Npatchmean = npatch(:,k+1)./ninbin(:,k+1)*5;
    rct_born = 1./(born(:,k+1)./Npatchmean);
    rct_split = 1./(split(:,k+1)./Npatchmean);
    rct_merge = 1./(merge(:,k+1)./Npatchmean);
    rct_death = 1./(die(:,k+1)./Npatchmean);
    
    rct_all = 1./((die(:,k+1)+merge(:,k+1)+split(:,k+1)+born(:,k+1))./Npatchmean);
    
    subplot(4,1,k+1)
    hold on
    plot(bint(1:end-1)+0.5, rct_born,'k-s','MarkerSize',8,'MarkerFaceColor',[178,171,210]/255)
    hold on
    plot(bint(1:end-1)+0.5, rct_split,'k-s','MarkerSize',8,'MarkerFaceColor',[253,184,99]/255)
    hold on
    plot(bint(1:end-1)+0.5,rct_merge,'k-s','MarkerSize',8,'MarkerFaceColor',[230,97,1]/255)
    hold on
    plot(bint(1:end-1)+0.5,rct_death,'k-s','MarkerSize',8,'MarkerFaceColor',[94,60,153]/255)
    legend('birth','merge','split','death','Location','best','Orientation','Horizontal')
    xlabel('Model time (Myr)')
    ylabel('Recurrence time (Myr)')
    if k==0
        title('elevation band 0-500 m')
    end
    if k==1
        title('elevation band 500-1500 m')
    end
    if k==2
        title('elevation band 1500-2500 m')
    end
    if k==3
        title('elevation band 2500-3000 m')
    end    

    set(gca,'FontSize',18)
    box on
    
end






