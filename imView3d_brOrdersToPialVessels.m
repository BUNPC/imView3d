function inView3d_brOrdersToPialVessels(vesselType)

% vesselType: 1 = arteries - calc. branch order from pial arteries; 3 -
% veins - calc. br. order from pial veins

global im;
%[fname pname] = uigetfile('*.seed','Select seed file with Im structure');
%load([pname fname],'-mat'); % load 'im2' structure

if ~isfield(im,'grpStat'),
    disp('No info on group statistics!');
    disp('Run Calc. Segment Distances first!');
    return;
end;

grp = im.grpStat.selectedGroupNumber; 
segLst = im.grpStat.segLst;
numSegments = im.grpStat.numSegments;
lookuptable = im.grpStat.lookuptable;
numEndNodes = im.grpStat.numEndNodes;
dSP = im.grpStat.dSP;
nSP = im.grpStat.nSP;
dSPbr = im.grpStat.dSPbr;
nSPbr = im.grpStat.nSPbr;
segBranchOrder = im.grpStat.segBranchOrder;


% NOW, CREATE THE VESSEL MASK DEPENDING ON VESSELTYPE IN FUNCTION PARAMETER
% LIST

%%
% now calculate PO2 distribution by branching order from pial arteries
    

for ii = 1:length(im.PO2pts.PO2),
    po2Seg = im.PO2pts.seg(ii);
    po2_BRorder(ii) = segBranchOrder(find(segBranchOrder(:,1)==po2Seg),2);
    po2_PO2(ii) = im.PO2pts.PO2(ii);
end;
BRorderIdx = sort(unique(po2_BRorder));
for ii=1:length(BRorderIdx),
    foolst = find(po2_BRorder == BRorderIdx(ii));
    averArtPO2_BRorder(ii) = mean(po2_PO2(foolst));
    if length(foolst) > 1,
        averArtPO2err_BRorder(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
    else
        averArtPO2err_BRorder(ii)=0;
    end;
end;
figure; 
errorbar(BRorderIdx,averArtPO2_BRorder,averArtPO2err_BRorder,'.-');
xlabel('Branching Order from Pial Arteries');
ylabel('Average PO2 (mmHg)');

%%
% now calculate PO2 distribution by branching order from pial arteries
% BUT EXCLUDE SEGMENTS LABELED AS VEINS (EXCLUDE SEGMENTS LABELED WITH VESSEL TYPE = 3). Find segment diameters too
    
n=0;
clear po2_BRorder;
clear po2_PO2;
clear po2_SO2;
clear averArtPO2_BRorder;
clear averArtPO2err_BRorder;
clear averArtDiam_BRorder;
clear averArtDiamerr_BRorder;
clear averSO2;
clear averSO2err;

hillC = 2.7; p50 = 37; % Rat, Ellis C.G. et al., Am. J. Phys Heart Circ Phys 258, H1216-, 2002
hillC = 2.59; p50 = 40.2; % C57BL/6 mice, Uchida K. et al., Zoological Science 15, 703-706, 1998
for ii = 1:length(im.PO2pts.PO2),
    po2Seg = im.PO2pts.seg(ii);
    %if (segBranchOrder(find(segBranchOrder(:,1)==po2Seg),4) > 4), %try to exclude
    %venules, veins, everything labeled as venous side
    if( (im.segVesType(po2Seg) == 1) || (im.segVesType(po2Seg) ==2) ), % if arterial or cappilary
        n=n+1;
        po2_BRorder(n) = segBranchOrder(find(segBranchOrder(:,1)==po2Seg),2);
        po2_PO2(n) = im.PO2pts.PO2(ii);
        po2_SO2(n) = (po2_PO2(n))^hillC / ( (po2_PO2(n))^hillC + p50^hillC  ); 
        segDiam(n) = segBranchOrder(find(segBranchOrder(:,1)==po2Seg),6);
    end;
end;
BRorderIdx = sort(unique(po2_BRorder));
for ii=1:length(BRorderIdx),
    foolst = find(po2_BRorder == BRorderIdx(ii));
    averArtPO2_BRorder(ii) = mean(po2_PO2(foolst));
    averArtDiam_BRorder(ii) = mean(segDiam(foolst));
    averSO2(ii) = mean(po2_SO2(foolst));
    if length(foolst) > 1,
        averArtPO2err_BRorder(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
        averArtDiamerr_BRorder(ii) = std(segDiam(foolst))/sqrt(length(foolst)-1);
        averSO2err(ii) = std(po2_SO2(foolst))/sqrt(length(foolst)-1);
    else
        averArtPO2err_BRorder(ii)=0;
        averArtDiamerr_BRorder(ii) = 0;
        averSO2err(ii) = 0;
    end;
end;
fig = figure; 
[AX, H1, H2] = plotyy(BRorderIdx,averArtPO2_BRorder,BRorderIdx,averArtDiam_BRorder);

set(fig,'CurrentAxes',AX(1));
hold on;
errorbar(BRorderIdx,averArtPO2_BRorder,averArtPO2err_BRorder,'.-','Color','r');
axis([min(BRorderIdx) max(BRorderIdx) 0 120]);
foo = [BRorderIdx' averArtPO2_BRorder' averArtPO2err_BRorder'];
save('ArtPO2.dat','foo','-ascii');

set(fig,'CurrentAxes',AX(2));
hold on;
errorbar(BRorderIdx,averArtDiam_BRorder,averArtDiamerr_BRorder,'.-b');
axis([min(BRorderIdx) max(BRorderIdx) 5 35]);

set(AX(1),'YColor','r','YTick',0:20:120);
set(get(AX(1),'YLabel'),'String','pO2 (mmHg)','FontSize',15);
set(AX(2),'YColor','b');
set(get(AX(2),'YLabel'),'String','Vessel Diameter (um)','FontSize',15);
set(get(AX(2),'XLabel'),'String','Branch Order from Pial Arteries','FontSize',15);
grid on;
print(fig,'-dpsc2','PO2_vs_BranchOrder_ARTERIES.ps');
%print(fig,'-djpeg','PO2_vs_BranchOrder.jpg','-r300');


fig = figure; 
[AX, H1, H2] = plotyy(BRorderIdx,averSO2.*100,BRorderIdx,averArtDiam_BRorder);

set(fig,'CurrentAxes',AX(1));
hold on;
errorbar(BRorderIdx,averSO2.*100,averSO2err.*100,'.-','Color','r');
axis([min(BRorderIdx) max(BRorderIdx) 0 120]);
foo = [BRorderIdx' (averSO2.*100)' (averSO2err.*100)'];
save('ArtSO2.dat','foo','-ascii');

set(fig,'CurrentAxes',AX(2));
hold on;
errorbar(BRorderIdx,averArtDiam_BRorder,averArtDiamerr_BRorder,'.-b');
axis([min(BRorderIdx) max(BRorderIdx) 5 35]);
foo = [BRorderIdx' averArtDiam_BRorder' averArtDiamerr_BRorder'];
save('ArtDiam.dat','foo','-ascii');

set(AX(1),'YColor','r','YTick',0:20:120);
set(get(AX(1),'YLabel'),'String','SO2 (%)','FontSize',15);
set(AX(2),'YColor','b','YTick',5:5:35);
set(get(AX(2),'YLabel'),'String','Vessel Diameter (um)','FontSize',15);
set(get(AX(2),'XLabel'),'String','Branch Order from Pial Arteries','FontSize',15);
grid on;
print(fig,'-dpsc2','SO2_vs_BranchOrder_ARTERIES.ps');
%print(fig,'-djpeg','SO2_vs_BranchOrder.jpg','-r300');



%%
% now calculate PO2 distribution by branching order from pial veins
% BUT EXCLUDE SEGMENTS LABELED AS ARTERIES (EXCLUDE SEGMENTS LABELED WITH VESSEL TYPE = 1). Find segment diameters too
    
n=0;
clear po2_BRorder;
clear po2_PO2;
clear po2_SO2;
clear averVeinPO2_BRorder;
clear averVeinPO2err_BRorder;
clear averVeinDiam_BRorder;
clear averVeinDiamerr_BRorder;
clear averSO2;
clear averSO2err;

hillC = 2.7; p50 = 37; % Rat, Ellis C.G. et al., Am. J. Phys Heart Circ Phys 258, H1216-, 2002
hillC = 2.59; p50 = 40.2; % C57BL/6 mice, Uchida K. et al., Zoological Science 15, 703-706, 1998
for ii = 1:length(im.PO2pts.PO2),
    po2Seg = im.PO2pts.seg(ii);
    %if (segBranchOrder(find(segBranchOrder(:,1)==po2Seg),4) > 4), %try to exclude
    %venules, veins, everything labeled as venous side
    if( (im.segVesType(po2Seg) == 3) || (im.segVesType(po2Seg) ==2) ), % if vein or cappilary
        n=n+1;
        po2_BRorder(n) = segBranchOrder(find(segBranchOrder(:,1)==po2Seg),4);
        po2_PO2(n) = im.PO2pts.PO2(ii);
        po2_SO2(n) = (po2_PO2(n))^hillC / ( (po2_PO2(n))^hillC + p50^hillC  ); 
        segDiam(n) = segBranchOrder(find(segBranchOrder(:,1)==po2Seg),6);
    end;
end;
BRorderIdx = sort(unique(po2_BRorder));
for ii=1:length(BRorderIdx),
    foolst = find(po2_BRorder == BRorderIdx(ii));
    averVeinPO2_BRorder(ii) = mean(po2_PO2(foolst));
    averVeinDiam_BRorder(ii) = mean(segDiam(foolst));
    averSO2(ii) = mean(po2_SO2(foolst));
    if length(foolst) > 1,
        averVeinPO2err_BRorder(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
        averVeinDiamerr_BRorder(ii) = std(segDiam(foolst))/sqrt(length(foolst)-1);
        averSO2err(ii) = std(po2_SO2(foolst))/sqrt(length(foolst)-1);
    else
        averVeinPO2err_BRorder(ii)=0;
        averVeinDiamerr_BRorder(ii) = 0;
        averSO2err(ii) = 0;
    end;
end;
fig = figure; 
[AX, H1, H2] = plotyy(BRorderIdx,averVeinPO2_BRorder,BRorderIdx,averVeinDiam_BRorder);

set(fig,'CurrentAxes',AX(1));
hold on;
errorbar(BRorderIdx,averVeinPO2_BRorder,averVeinPO2err_BRorder,'.-','Color','r');
axis([min(BRorderIdx) max(BRorderIdx) 0 120]);
foo = [BRorderIdx' averVeinPO2_BRorder' averVeinPO2err_BRorder'];
save('VeinPO2.dat','foo','-ascii');

set(fig,'CurrentAxes',AX(2));
hold on;
errorbar(BRorderIdx,averVeinDiam_BRorder,averVeinDiamerr_BRorder,'.-b');
axis([min(BRorderIdx) max(BRorderIdx) 5 35]);

set(AX(1),'YColor','r','YTick',0:20:120);
set(get(AX(1),'YLabel'),'String','pO2 (mmHg)','FontSize',15);
set(AX(2),'YColor','b');
set(get(AX(2),'YLabel'),'String','Vessel Diameter (um)','FontSize',15);
set(get(AX(2),'XLabel'),'String','Branch Order from Pial Veins','FontSize',15);
grid on;
print(fig,'-dpsc2','PO2_vs_BranchOrder_VEINS.ps');
%print(fig,'-djpeg','PO2_vs_BranchOrder.jpg','-r300');


fig = figure; 
[AX, H1, H2] = plotyy(BRorderIdx,averSO2.*100,BRorderIdx,averVeinDiam_BRorder);

set(fig,'CurrentAxes',AX(1));
hold on;
errorbar(BRorderIdx,averSO2.*100,averSO2err.*100,'.-','Color','r');
axis([min(BRorderIdx) max(BRorderIdx) 0 120]);
foo = [BRorderIdx' (averSO2.*100)' (averSO2err.*100)'];
save('VeinSO2.dat','foo','-ascii');

set(fig,'CurrentAxes',AX(2));
hold on;
errorbar(BRorderIdx,averVeinDiam_BRorder,averVeinDiamerr_BRorder,'.-b');
axis([min(BRorderIdx) max(BRorderIdx) 5 35]);
foo = [BRorderIdx' averVeinDiam_BRorder' averVeinDiamerr_BRorder'];
save('VeinDiam.dat','foo','-ascii');

set(AX(1),'YColor','r','YTick',0:20:120);
set(get(AX(1),'YLabel'),'String','SO2 (%)','FontSize',15);
set(AX(2),'YColor','b','YTick',5:5:35);
set(get(AX(2),'YLabel'),'String','Vessel Diameter (um)','FontSize',15);
set(get(AX(2),'XLabel'),'String','Branch Order from Pial Veins','FontSize',15);
grid on;
print(fig,'-dpsc2','SO2_vs_BranchOrder_VEINS.ps');
%print(fig,'-djpeg','SO2_vs_BranchOrder.jpg','-r300');




%%
% now calculate PO2 distribution by distance from pial arteries
    

for ii = 1:length(im.PO2pts.PO2),
    po2Seg = im.PO2pts.seg(ii);
    po2_DIS(ii) = segBranchOrder(find(segBranchOrder(:,1)==po2Seg),3);
    po2_PO2(ii) = im.PO2pts.PO2(ii);
end;
% bin in 100 um bins...
dDIS = 100;
maxDIS = max(po2_DIS)+1;
minDIS = min(po2_DIS)-1;
nDIS = floor((maxDIS-minDIS)/dDIS);
for ii=1:nDIS,
    foolst = find( (po2_DIS<ii*dDIS) & (po2_DIS>=(ii-1)*dDIS) );
    DIST(ii) = (ii-0.5)*dDIS;
    if ~isempty(foolst),
        DIST_PO2(ii) = mean(po2_PO2(foolst));
        if length(foolst)>1,
            DIST_PO2err(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
        else
            DIST_PO2err(ii) = 0;
        end;
    else
        DIST_PO2(ii) = -1;
        DIST_PO2err(ii) = 0;
    end;
end;
figure; 
errorbar(DIST,DIST_PO2,DIST_PO2err,'.-');
xlabel('Distance from Pial Arteries (um)');
ylabel('Average PO2 (mmHg)');        

%%
% now calculate PO2 distribution by branching order from pial veins
    

for ii = 1:length(im.PO2pts.PO2),
    po2Seg = im.PO2pts.seg(ii);
    po2_BRorder(ii) = segBranchOrder(find(segBranchOrder(:,1)==po2Seg),4);
    po2_PO2(ii) = im.PO2pts.PO2(ii);
end;
BRorderIdx = sort(unique(po2_BRorder));
for ii=1:length(BRorderIdx),
    foolst = find(po2_BRorder == BRorderIdx(ii));
    averVeinPO2_BRorder(ii) = mean(po2_PO2(foolst));
    if length(foolst) > 1,
        averVeinPO2err_BRorder(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
    else
        averVeinPO2err_BRorder(ii)=0;
    end;
end;
figure; 
errorbar(BRorderIdx,averVeinPO2_BRorder,averVeinPO2err_BRorder,'.-');
xlabel('Branching Order from Pial Veins');
ylabel('Average PO2 (mmHg)');

%%
% now calculate PO2 distribution by distance from pial veins
    

for ii = 1:length(im.PO2pts.PO2),
    po2Seg = im.PO2pts.seg(ii);
    po2_DIS(ii) = segBranchOrder(find(segBranchOrder(:,1)==po2Seg),5);
    po2_PO2(ii) = im.PO2pts.PO2(ii);
end;
% bin in 100 um bins...
dDIS = 100;
maxDIS = max(po2_DIS)+1;
minDIS = min(po2_DIS)-1;
nDIS = floor((maxDIS-minDIS)/dDIS);
for ii=1:nDIS,
    foolst = find( (po2_DIS<ii*dDIS) & (po2_DIS>=(ii-1)*dDIS) );
    DIST(ii) = (ii-0.5)*dDIS;
    if ~isempty(foolst),
        DIST_PO2(ii) = mean(po2_PO2(foolst));
        if length(foolst)>1,
            DIST_PO2err(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
        else
            DIST_PO2err(ii) = 0;
        end;
    else
        DIST_PO2(ii) = -1;
        DIST_PO2err(ii) = 0;
    end;
end;
figure; 
errorbar(DIST,DIST_PO2,DIST_PO2err,'.-');
xlabel('Distance from Pial Veins (um)');
ylabel('Average PO2 (mmHg)');                



