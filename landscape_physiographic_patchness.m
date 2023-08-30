
%
%
% 1) The general algorithm
% This script reads topography (x, y, z) from LEMs and interpolates to a user-defined regular grid.
% The land surface is then represented by a series of elemental surfaces.
% Elevation, hillslope gradient, and hillslope aspect is calculated for
% each elemental surface. Elemental surfaces with parametric values falling
% in this range may be connected and form a contiguous patch of surface.
% Individual patches are indentified with image functions as objects.
%
% 2) Inputs: text files with three columns, column 1-2 stores x, y (meter),
%            column 3 stores z (meter). This script used outputs of the
%            landscape evolution model DAC(Divide And Capture) where the landscape
%            is characterized by rivers and divides.
%
% 3) Outputs: MATLAB data files (.mat) stores the number of different types
%             of events, and related patch identities of split and merge event.
%
% 4) Special attention: timestepping should be small so that the same patch
%             still overlap in space for consecutive steps.
%
%
% Author: Dr. Yanyan Wang (yanyan.wang@erdw.ethz.ch, or wangyanyan0607@hotmail.com)
% Last update: Aug. 30th, 2023



%% Basic model parameters of the landscape evolution model outputs

clear;

% text file directory
pathname  = '/landscape_physiographic_patchness/ASCII/';

% file name of the first step text file
mstart = 30400;

% file name incremental value
printstep = 5;

% timestep
dt = 2e4;

% steps to calculate
stepplot = 1:1:520;

% model domain size
XL = 250e3; % meter
YL = 100e3; % meter
dx = 100; % meter

% remeshing grid for interpolation
[xq,yq] = ndgrid(0:dx:XL, 0:dx:YL);


%% Physiographic parameters and range to define a habitat
z_cutoffs = [0 500 1500 2500 3000]; % elevation cutoffs
g_cutoffs = [0 3 20]; % hillslope gradient cutoffs
asp_cutoffs = [45 325]; % aspect cutoffs, aspect is defined as clockwise from 325 to 45 degree.

crita = 100; % critical area to be counted, [number of grid cells]

%% habitat lifespan
% 24 types of combinations of elevation, gradient and aspect of hillslope
% cutoffs.
types = [111   112   113   121   122   123   211   212   213   221   222   223   311   312   313   321   322   323  411 412 413 421 422 423];
type_merg = [1 2 3 1 5 6 7 8 9 7 11 12 13 14 15 13 17 18 19 20 21 19 23 24]; % not distinguish aspect for shallow slope

type_plot = unique(type_merg); % not distinguish aspect for shallow slope


% loop through all types of habitat types
for j = 1:length(type_plot)
    % examine the evolutionary of a target habitat
    targetH = type_plot(j); % the a target habitat
    
    dir = strcat('habitat_',num2str(targetH));
    
    % create a folder for data storage
    mkdir(dir);
    
    Nmerge = zeros(length(stepplot),1);
    Nsplit = zeros(length(stepplot),1);
    Ndie = zeros(length(stepplot),1);
    Nborn = zeros(length(stepplot),1);
    Nmigrate = zeros(length(stepplot),1);
    Npatch = zeros(length(stepplot),1);
    
    firststep = 1;
    for ip = firststep:length(stepplot)
        
        i = stepplot(ip);
        % 1) construct string of file names
        filestr = '00000000';
        istep = num2str(mstart+i*printstep);
        filestr(end-length(istep)+1:end) = istep;
        
        % 2-0) coordinates of river nodes
        fncoor = strcat(pathname, 'coordinates',filestr);
        coor = importdata(fncoor);
        xr = coor(:,1);
        yr = coor(:,2);
        zr = coor(:,3);
        
        % 2-1) divide coordinates and elvation
        fndivide = strcat(pathname,'no_channel_connection',filestr);
        divide = importdata(fndivide);
        divi = divide(:,1);% bounding node i
        divj = divide(:,2);% bounding node j
        
        fndivxy = strcat(pathname,'divides',filestr);
        divxy = importdata(fndivxy);
        
        goodiv = divxy(:,3)>zr(divi)&divxy(:,3)>zr(divj);
        
        xd = divxy(goodiv,1);
        yd = divxy(goodiv,2);
        zd = divxy(goodiv,3);
        
        clear divi
        clear divj
        clear divxy
        clear divide
        
        % interpolate to gridded data
        % tri = delaunay(X,Y) creates a 2-D Delaunay triangulation. 'tri' is a matrix representing the set of triangles that make up the triangulation.
        x = [xr; xd];
        y = [yr; yd];
        z = [zr; zd];
        % use the unique points
        [~,ia,~] = unique(x.^2+y.^2);
        x = x(ia);
        y = y(ia);
        z = z(ia);
        % interp to gridded data
        F = scatteredInterpolant(x,y,z,'natural');
        zq = F(xq,yq);
        
        % calculate dip angle and aspect
        [Nx,Ny,Nz] = surfnorm(xq,yq,zq);
        id = Nz<0;
        Nx(id) = -Nx(id);
        Ny(id) = -Ny(id);
        Nz(id) = -Nz(id);
        ASP = cart2pol(Nx,Ny);
        
        % this defines aspect to be North=0, East=90, South=180, West=270 in degree
        ASP = mod(90-ASP/pi*180,360);
        
        % calculate surface dipping angle
        DIP = acosd(abs(Nz)); % in degree
        
        % define habitat types
        id1 =  ASP<=asp_cutoffs(1)|ASP>asp_cutoffs(2); % facing to the north
        
        idz1 = zq>z_cutoffs(2)&zq<=z_cutoffs(3);
        idz2 = zq>z_cutoffs(3)&zq<=z_cutoffs(4);
        idz3 = zq>z_cutoffs(5);
        
        ids1 = DIP>g_cutoffs(2)&DIP<=g_cutoffs(3);
        ids2 = DIP>g_cutoffs(3);
        
        habitat_aspect = ones(size(xq));
        habitat_aspect(id1) = 2;
        
        habitat_elevation = ones(size(xq));
        habitat_elevation(idz1) = 2;
        habitat_elevation(idz2) = 3;
        habitat_elevation(idz3) = 4;
        
        habitat_slope = ones(size(xq));
        habitat_slope(ids1) = 2;
        habitat_slope(ids2) = 3;
        
        clear id1
        clear idz1
        clear idz2
        clear ids1
        clear ids2
        
        habitat = habitat_elevation*100+habitat_aspect*10+habitat_slope;
        
        % find habitat patches from image functions
        habitat_label = zeros(size(habitat));  % habitat label from merged
        habitat_patch_polygon = zeros(size(habitat));  % habitat patch polygon
        for ki = 1:numel(types)
            idk = habitat==types(ki);
            habitat_label(idk) = type_merg(ki);
        end
        
        % create folder under the habitat folder to store split and merge
        dir2 = strcat(dir,'/split_merge_event_list/');
        mkdir(dir2)
        
        %initialize merge event array
        mergename = strcat(dir2,'Step',num2str(ip),'_merge.mat');
        % the merge array first element is children
        mergeEvent = nan(50,6);
        
        %initialize split event array
        splitname = strcat(dir2,'Step',num2str(ip),'_split.mat');
        % the split array first element is parent
        splitEvent = nan(50,6);
        
        if ip==firststep
            habitat_label0 = habitat_label; % habitat label in the last step
            BW0 = habitat_label0==targetH;
            [B0, L0, nobj0, Aobj0] = bwboundaries(BW0,'holes');
            % deal with donut patches
            L02 = zeros(size(L0));
            L02(L0<=nobj0) = L0(L0<=nobj0);
            % filter small patches
            stats = regionprops('table',BW0,'Area');
            AREAobj0 = stats.Area;
            id = AREAobj0<crita;
            hb = 1:1:nobj0;
            id2 = ismember(L02(:),hb(id));
            L02(id2) = 0;
            
            Npatch(ip) = numel(unique(L02))-1;
            
            nobj0belowcrita = nobj0-Npatch(ip); % thoese should not be counted as death in the second step
            
            % initialize patch id, age, and patch center, allocate big space for
            % each array
            patchMaxN = nobj0*200;
            patchID0 = nan(length(stepplot),patchMaxN);
            patchID0(ip,1:nobj0) = 1:nobj0;
            patchAge0 = zeros(length(stepplot),patchMaxN);
            patchAge0(ip,1:nobj0) = 1; % age increment is one step
            patchCenterX0 = nan(length(stepplot),patchMaxN);
            patchCenterY0 = nan(length(stepplot),patchMaxN);
            patchArea0 = nan(length(stepplot),patchMaxN);   % number of pixels of each patch
            patchhistory0 = zeros(length(stepplot),patchMaxN); % the history of the patches, 0-initial step
            
            for ipk = 1:nobj0
                id = L02==ipk;
                patchCenterX0(ip,ipk) = mean(xq(id));
                patchCenterY0(ip,ipk) = mean(yq(id));
                patchArea0(ip,ipk) = sum(id(:));
            end
            
            patchID = patchID0;
            patchAge = patchAge0;
            patchCenterX = patchCenterX0;
            patchCenterY = patchCenterY0;
            patchArea = patchArea0;
            originalN = max(patchID0(:)); % habitat number of first step
            incpatch = 0; % array expansion number
            patchhistory = patchhistory0;
        else
            
            % patches current step
            BW = habitat_label==targetH;
            [~, L, nobj, ~] = bwboundaries(BW,'holes');
            L2 = zeros(size(L));
            L2(L<=nobj) = L(L<=nobj);
            % filter small patches
            crita = 100; % critical area to be counted as a new born, [n. of pixels]
            stats = regionprops('table',BW,'Area');
            AREAobj = stats.Area;
            id = AREAobj<crita;
            hb = 1:1:nobj;
            id2 = ismember(L2(:),hb(id));
            L2(id2) = 0;
            
            Npatch(ip) =  numel(unique(L2))-1;
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % expand habitat patch arrays, compare the old patch and new patch, update
            currentN = max(patchID(:));
            %  %  %  %  % - - - Case 1, examine born, loop through new patches
            L3 = L2;
            id2 = L2~=0;
            patchn = unique(L2(id2));
            patchn = patchn(patchn~=0);% remove patch id 0
            if numel(patchn)~=0
                for pati = 1:numel(patchn)
                    id3 = L2==patchn(pati);
                    if sum(L02(id3),'all') ==0&&sum(id3(:))>=crita
                        % update patch stats
                        incpatch = incpatch+1;
                        patchID(ip,originalN+incpatch) = originalN+incpatch;
                        patchAge(ip,originalN+incpatch) = patchAge(ip-1,originalN+incpatch)+1;
                        patchArea(ip,originalN+incpatch) = sum(id3(:));
                        patchCenterX(ip,originalN+incpatch) = mean(xq(id3));
                        patchCenterY(ip,originalN+incpatch) = mean(yq(id3));
                        patchhistory(ip,originalN+incpatch) = 1; % new born, no ancestors
                        % update patch matrix
                        L3(id3) = originalN+incpatch;
                        Nborn(ip) = Nborn(ip)+1;
                    end
                end
            end
            
            % Case others, loop through old patchs
            for kpi = 1:currentN
                
                id = L02==kpi;
                id2 = L2~=0;
                patchn = unique(L2(id&id2)); % the number of current patches overlap with the old patch
                patchn = patchn(patchn~=0); % remove patch id 0
                patcha = AREAobj(patchn); % current patches area
                
                %  %  %  %  % - - - Case 2, patch die
                if isempty(patchn)
                    % update patch stats
                    patchID(ip,kpi) = patchID(ip-1,kpi); % patch id froze
                    patchAge(ip,kpi) = patchAge(ip-1,kpi); % age froze
                    patchArea(ip,kpi) = 0; % patch area becomes zero
                    patchCenterX(ip,kpi) = nan; % make patch center nan
                    patchCenterY(ip,kpi) = nan; % make patch center nan
                    patchhistory(ip,kpi) = 2; % die
                    %  no need to update patch matrix L3
                    %  L3(id&id2) = 0;
                    Ndie(ip) = Ndie(ip)+1;
                end
                
                %  %  %  %  % - - - Case 3, patch split, newly splited patches are counted as
                % newborns but inherit the parent age
                if numel(patchn)>1
                    % calculate overlapping area of current patches
                    patchna = zeros(numel(patchn),1);
                    for pati = 1:numel(patchn)
                        id4 = L2==patchn(pati);
                        patchna(pati) = sum(id4&id,'all');
                    end
                    % sort the inheritors with ovelapping area
                    [~,I] = sort(patchna,'descend');
                    % the old ancestor patch evolves into the biggest overlapping inheritor
                    patchID(ip,kpi) = kpi; % patch id not change
                    patchAge(ip,kpi) = patchAge(ip-1,kpi)+1; % age counting continue
                    patchArea(ip,kpi) = patcha(I(1)); % patch area uses the biggest overlapping inheritor's area
                    id4 = L2==patchn(I(1));
                    patchCenterX(ip,kpi) = mean(xq(id4)); % update patch center x
                    patchCenterY(ip,kpi) = mean(yq(id4)); % update patch center y
                    patchhistory(ip,kpi) = 3; % inherit the old one from split
                    % update patch matrix L3
                    L3(id4) = patchID(ip,kpi);
                    
                    Nsplit(ip) = Nsplit(ip)+1;
                    splitEvent(Nsplit(ip),1) = kpi;
                    
                    % the other inheritors are treated as newborns with
                    % inherited age
                    for pati = 2:numel(patchn)
                        % update newborn patch stats
                        incpatch = incpatch+1;
                        id4 = L2==patchn(I(pati));
                        patchID(ip,originalN+incpatch) = originalN+incpatch;
                        patchAge(ip,originalN+incpatch) = patchAge(ip-1,kpi)+1;
                        patchArea(ip,originalN+incpatch) = sum(id4(:));
                        patchCenterX(ip,originalN+incpatch) = mean(xq(id4));
                        patchCenterY(ip,originalN+incpatch) = mean(yq(id4));
                        patchhistory(ip,originalN+incpatch) = 7;% small new borns from split
                        % update patch matrix L3
                        L3(id4) = originalN+incpatch;
                        %save children ids to array
                        splitEvent(Nsplit(ip), pati) = originalN+incpatch;
                    end
                    
                end
                
                
                if numel(patchn)==1
                    id5 = L2==patchn;
                    patchn2 = unique(L02(id5)); % the number of old patches overlap with the current patch
                    patchn2 = patchn2(patchn2~=0);  % remove patch id 0
                    
                    %  %  %  %  % - - - Case 4, patch migrate and deform (should distinguish with merge)
                    if numel(patchn2)==1 % the current patch overlap with only one old patch
                        % update patch stats
                        patchID(ip,kpi) = patchID(ip-1,kpi); % keep patch id
                        patchAge(ip,kpi) = patchAge(ip-1,kpi)+1;% update patch age
                        patchArea(ip,kpi) = sum(id&id2,'all');% update patch area
                        patchCenterX(ip,kpi) = mean(xq(id&id2));
                        patchCenterY(ip,kpi) = mean(yq(id&id2));
                        patchhistory(ip,kpi) = 4; % deform and migrate
                        % update patch matrix L3
                        L3(id5) = patchID(ip,kpi);
                        % update patch migration of current step
                        Nmigrate(ip) = Nmigrate(ip)+1;
                    end
                    
                    %  %  %  %  % - - - Case 5, patch merge, should distinguish with migrate
                    if numel(patchn2)>1 % the current patch overlap with more than one old patch
                        
                        % calculate overlapping area of ancestors
                        patchn2a = zeros(numel(patchn2),1);
                        for pati = 1:numel(patchn2)
                            id6 = L02==patchn2(pati);
                            patchn2a(pati) = sum(id6&id5,'all');
                        end
                        % sort ancestor patches with overlapping area
                        [~,I] = sort(patchn2a,'descend');
                        % the largest overlapping patch evolves into the newly
                        % merged patch
                        id4 = L2==patchn;
                        patchID(ip,patchn2(I(1))) = patchn2(I(1));
                        patchAge(ip,patchn2(I(1))) = patchAge(ip-1,patchn2(I(1)))+1; % inherit the largest ancestor' age
                        patchArea(ip,patchn2(I(1))) = sum(id4(:));
                        patchCenterX(ip,patchn2(I(1))) = mean(xq(id4));
                        patchCenterY(ip,patchn2(I(1))) = mean(yq(id4));
                        patchhistory(ip,patchn2(I(1))) = 5; % merge and dominate the merged patch
                        % update patch matrix L3
                        L3(id4) = patchn2(I(1));
                        
                        if ~ismember(patchn2(I(1)),mergeEvent(:,1))
                            % update merge events
                            Nmerge(ip) = Nmerge(ip)+1;
                            
                            % save merge list
                            mergeEvent(Nmerge(ip),1) = patchn2(I(1));
                            mergeEvent(Nmerge(ip),2:numel(patchn2)) = patchn2(I(2:numel(patchn2)));
                        end
                        % update stats of other old patches that merged
                        for pati = 2:numel(patchn2)
                            patchID(ip,patchn2(I(pati))) = patchn2(I(pati)); % patch id froze
                            patchAge(ip,patchn2(I(pati))) = patchAge(ip-1,patchn2(I(pati))); % age froze
                            patchCenterX(ip,patchn2(I(pati))) = nan; % old patch treated as die
                            patchCenterY(ip,patchn2(I(pati))) = nan; % old patch treated as die
                            
                            if patchhistory(ip,patchn2(I(pati)))==0
                                patchhistory(ip,patchn2(I(pati))) = 6; % history flag as merge with a bigger one,5+1
                                patchArea(ip,patchn2(I(pati))) = patchArea(ip-1,patchn2(I(pati))); % old patch area froze,
                            else
                                patchhistory(ip,patchn2(I(pati))) = 7; % history flag as merge and split
                                % area not change because it has been treated
                                % already
                            end
                            
                            
                        end
                        
                        
                    end
                end
                % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                
            end
            
            % update patch info arrays
            L02 = L3;
            
            save(mergename,'mergeEvent')
            save(splitname,'splitEvent')
            
        end
        
        
    end
    
    
    Ndie2 = [Ndie(1);diff(Ndie)];
    Ndie2(2) = Ndie2(2)- nobj0belowcrita; % those below crita patch of first step should not be counted as death
    
    filename = strcat(dir,'/patch_evolution_statistics_habitat_',num2str(targetH),'.mat');
    save(filename,'patchAge','patchArea','patchCenterX','patchCenterY','patchID','patchhistory','Nmerge','Nsplit','Nborn','Ndie2','Nmigrate','Npatch')
    
end






