% TO DO: 
% 1) Add to TOOLS menu the estimation of tortuosity and option to
%     straighten first 5, 10, 15, or 20 times for tortCurv. 
%     Don't calculate tortAng since it is not standard
% 2) Add pruning procedure to TOOLS menu with subsequent call to nodeGrps
% 3) Check out and verify behavior of removeNodes (see below)
%
% Strategy for pruning segments
% 1) Loop through end point nodes of segments with length < diameter
%    and query user if segment should be removed
% 2) Loop through segments with length < diameter of segments without
%    end point nodes and query user if segment should be removed
%
% NOT NEEDED BECAUSE ALL PARAMETERS UPDATED BY CALL TO nodeGrps
% removeNodes needs updating if it is to preserve segment information
% This includes: nodeSegN, edgeSegN, 
%      segLen, segDiam, segPos, tort, tortAng, tortCurv
% Also: segNedges, segVesType
% What should I do with: nodeGrp, nodeSegs
% I would still have to update nB
% Only nodeSegN, edgeSegN, and nodeGrp makes sense to update
%   the others are not clearly indicated by nodeFlag for removal in removeNodes
%   and therefore segment info needs to be updated with nodeGrps after
%   removing nodes! So, if I need to run nodeGrps then no need to update
%   nodeSegN, edgeSegN, and nodeGrp in removeNodes since they get reset in
%   nodeGrps
%
% What is: segNodeMap. Not needed so delete it. DONE
%           it was used in imView3d_flowCircuit but I think it is not
%           needed
%
% Not updated by removeNodes but should have been: nodeBC, nodeBCType, 
%           nodeType, nodeDiamThetaIdx
% Should these also be updated?    nodeVel, nodeDiamEst, nodeTypeUpdated
% I NEED TO LOOK AT THE WORK FLOW AND BETTER MANAGE THESE VARIABLES
% RENAME: segTort, segTortAng, segTortCurv
%         nodeSegs -> segEndNodes                  DONE

nNodes = size(im2.nodePos,1);
nB=zeros(1,nNodes);
for ii=1:nNodes
    nB(ii)=length(find(im2.nodeEdges(:,1)==ii | im2.nodeEdges(:,2)==ii));
end
im2.nB = nB;

im2 = nodeGrps(im2);

%%
% what segment info do we have
% segDiam
% segLen
% segPos
% tort, tortAng, tortCurv
% 
%
% I also want 
% nB
nSeg = length(im2.segDiam);
im2.segPos = squeeze(mean(reshape(im2.nodePos(im2.segEndNodes,:),[2 nSeg 3]),1));


%%
% find segments with length less than a fraction of the diameter
% which ones should be pruned? May have to manually view all of them
% 1) those that have an end point node. But what if it is a true end point
% 
lst2 = find((im2.segLen./max(im2.segDiam,1))<1 & im2.segLen>0);

% majority of segments with short lengths are at deep depths where we have
% lots of noise. Maybe we can view those that are more shallow and then manually
% prune.
% Also select segments without end point nodes
zThresh = 400;
lst4=im2.segEndNodes(lst2,:);
foo = find(im2.nodePos(lst4(:,1),3)<zThresh & im2.nodePos(lst4(:,2),3)<zThresh & im2.nB(lst4(:,1))'>1 & im2.nB(lst4(:,2))'>1);
lst2b = lst2(foo);
for ii = 1:length(lst2b)
    lst3 = find(im2.edgeSegN==lst2b(ii));
    eFlag = ones(size(im2.nodeEdges,1),1);
    eFlag(lst3) = 2;
    [pMin,pMax,Tree]=plotTreeSeg(im2.nodeEdges,im2.nodePos,im2.segEndNodes(lst2b(ii),1),15,[1 1],im2.I,[],[],[],eFlag);
    title( sprintf('Seg%d  diam=%.0f  len=%.0f',lst2b(ii),im2.segDiam(lst2b(ii)),im2.segLen(lst2b(ii))) )

    figure(1)
    z1 = floor(min(im2.nodePos(im2.segEndNodes(lst2b(ii),:),3)));
    z2 = ceil(max(im2.nodePos(im2.segEndNodes(lst2b(ii),:),3)));
    imagesc( max(im2.I(:,:,z1:z2),[],3) )
    hold on
    for jj=1:length(Tree)  
        n1 = im2.nodeEdges(Tree(jj),1);
        n2 = im2.nodeEdges(Tree(jj),2);
        hl=plot([im2.nodePos(n1,1) im2.nodePos(n2,1)],[im2.nodePos(n1,2) im2.nodePos(n2,2)],'y-');
        if eFlag(Tree(jj))==2
            set(hl,'color','r');
        end
%        plot([im2.nodePos(n1,1)],[im2.nodePos(n1,2)],'y.')
%        plot([im2.nodePos(n2,1)],[im2.nodePos(n2,2)],'y.')
    end
    hold off
%    xlim( [pMin(1) pMax(1)] )
%    ylim( [pMin(2) pMax(2)] )
    
    pause
end

%%
% view end point nodes
% with segLen < segDiam
% and query to manually remove
lst = find(im2.nB==1);
lstSeg = im2.nodeSegN(lst);
im2.segLen(1403)=0; im2.segDiam(1403)=0;  % check nodeGrps to see why nodeSegN has 
     %1403 segments when segDiam and segLen only have 1402
lst2 = find(im2.segLen(lstSeg)<im2.segDiam(lstSeg));
lst = lst(lst2);

nodeFlag = ones(size(im2.nodePos,1),1);
eFlag = ones(size(im2.nodeEdges,1),1);
for ii=1:length(lst)
    iSeg = im2.nodeSegN(lst(ii));
    lst3 = find(im2.edgeSegN==iSeg);

    ch = 3;
    nEdges = 15;
    while ch>2
        eFlag(lst3) = 2;
        [pMin,pMax,Tree]=plotTreeSeg(im2.nodeEdges,im2.nodePos,lst(ii),nEdges,[1 1],im2.I,[],[],[],eFlag);
        iSeg = im2.nodeSegN(lst(ii));
        title( sprintf('(%d of %d) Node%d  Seg%d  diam=%.0f  len=%.0f',ii,length(lst),lst(ii),iSeg,im2.segDiam(iSeg),im2.segLen(iSeg)) )

        figure(1)
        z1 = floor(min(im2.nodePos(im2.segEndNodes(iSeg,:),3)));
        z2 = ceil(max(im2.nodePos(im2.segEndNodes(iSeg,:),3)));
        imagesc( max(im2.I(:,:,z1:z2),[],3) )
        hold on
        for jj=1:length(Tree)
            n1 = im2.nodeEdges(Tree(jj),1);
            n2 = im2.nodeEdges(Tree(jj),2);
            hl=plot([im2.nodePos(n1,1) im2.nodePos(n2,1)],[im2.nodePos(n1,2) im2.nodePos(n2,2)],'r-');
            set(hl,'linewidth',2)
            if eFlag(Tree(jj))~=2
                set(hl,'color','y')
            end
        end
        plot([im2.nodePos(lst(ii),1)],[im2.nodePos(lst(ii),2)],'r*')
        hold off
        title( sprintf('Z = %d',im2.segPos(iSeg,3)) )
        colormap gray

        eFlag(lst3) = 1;
        
        ch = menu('What to do?','Keep','Delete','20 edges','25 edges','+5 edges','Quit');
        if ch==3
            nEdges = 20;
        elseif ch==4
            nEdges = 25;
        elseif ch==5
            nEdges = nEdges + 5;
        elseif ch==6
            break
        elseif ch==2
            eFlag(lst3) = 0;
            nLst = im2.nodeEdges(lst3,:);
            nLst = unique(nLst(:));
            foo = find(im2.nB(nLst)==1 | im2.nB(nLst)==2);
            nodeFlag(nLst(foo)) = 0;
        end
    end
    if ch==6
        break
    end
end

ch = menu('Remove deleted nodes and recalculate Segments','Yes','No');
drawnow
if ch==1
    [im2.nodePos,im2.nodeDiam,im2.nodeDiamThetaIdx,im2.nodeBC,im2.nodeBCType,im2.nodeType,...
        im2.nodeSegN,im2.nodeEdges,im2.edgeFlag] = removeNodes( nodeFlag, im2.nodePos, ...
        im2.nodeDiam, im2.nodeDiamThetaIdx, im2.nodeBC, im2.nodeBCType, im2.nodeType, im2.nodeSegN, im2.nodeEdges );
    im2 = nodeGrps(im2);
end





%%
% find segments with length less than a fraction of the diameter
% and less than a certain depth
zThresh = 400;
Len2DiamRatio = 2;
lst2 = find((im2.segLen./max(im2.segDiam,1))<Len2DiamRatio);
lst4=im2.segEndNodes(lst2,:);
foo = find(im2.nodePos(lst4(:,1),3)<zThresh & im2.nodePos(lst4(:,2),3)<zThresh);
lstBad = lst2(foo);
nSeg = length(im2.segLen);
nSegBad = length(lstBad);

lst2 = find((im2.segLen./max(im2.segDiam,1))>=Len2DiamRatio);
lst4=im2.segEndNodes(lst2,:);
foo = find(im2.nodePos(lst4(:,1),3)<zThresh & im2.nodePos(lst4(:,2),3)<zThresh);
lstGood = lst2(foo);

nSeg = length(im2.segLen);
nSegBad = length(lstBad);
nSegGood = length(lstGood);
[nSeg nSegBad+nSegGood nSegBad nSegGood]

% histogram depth of "bad" segments
nIdx = im2.segEndNodes(lstBad,:);
figure(1)                 
hist( im2.nodePos(nIdx(:),3), [25:50:425] )
title( 'Depth of "bad" segments' )

%%
% calculate tortuosity
nIdx = im2.segEndNodes(lstGood,:);
rsep = sqrt(sum( (im2.nodePos(nIdx(:,1),:) - im2.nodePos(nIdx(:,2),:)).^2, 2));
tortGood = im2.segLen(lstGood)'./rsep;

nIdx = im2.segEndNodes(lstBad,:);
rsep = sqrt(sum( (im2.nodePos(nIdx(:,1),:) - im2.nodePos(nIdx(:,2),:)).^2, 2));
tortBad = im2.segLen(lstBad)'./rsep;

tort = zeros(length(im2.segLen),1);
tort(lstGood) = tortGood;
tort(lstBad) = tortBad;

%%
% straighten the segments
im = im2;
for jj=1:15 % this is needed to minimize noise accumulation
    imView3d_CenterNodes( 2, 1, 0, 0, 2, [] ); % eventdata, centerstep, flag, flag, Ithresh         
end
im2.nodePosS = im.nodePos;

%%
% calculate tortuosity by summing angular change. NOT A RECOGNIZED METHOD
% this probably correlates very well with simpler metric
tortAng = zeros(length(im2.segLen),1);
for ii=1:length(lstGood)
    iSeg = lstGood(ii);
    Tree = find(im2.edgeSegN==iSeg);
    
    nStart = im2.segEndNodes(iSeg,1);
    nEnd = im2.segEndNodes(iSeg,2);
    lst = find(im2.nodeEdges(Tree,1)==nStart | im2.nodeEdges(Tree,2)==nStart);
    if length(lst)~=1
        error('This should be length 1')
    end
    n0 = setdiff(im2.nodeEdges(Tree(lst),:),nStart);
    nLst = nStart;
    aSum = 0;
    while n0~=nEnd
        lst = find(im2.nodeEdges(Tree,1)==n0 | im2.nodeEdges(Tree,2)==n0);
        if length(lst)~=2
            error('This should be length 2')
        end
        e1 = im.nodePos(im2.nodeEdges(Tree(lst(1)),1),:)-im.nodePos(im2.nodeEdges(Tree(lst(1)),2),:);
        e2 = im.nodePos(im2.nodeEdges(Tree(lst(2)),1),:)-im.nodePos(im2.nodeEdges(Tree(lst(2)),2),:);
        e1 = e1 / norm(e1);
        e2 = e2 / norm(e2);
        aSum = aSum + acos(abs(sum(e1.*e2)))*180/3.14159;
        nLst(end+1) = n0;
        foo = im2.nodeEdges(Tree(lst),:);
        n0 = setdiff(foo(:),nLst);
    end
    tortAng(iSeg) = aSum;
end    


%%
% calculate tortuosity by integral of square of derivative of curvature,
% divided by the length of a curve
tortCurv = zeros(length(im2.segLen),1);

FilterOrder = 3;
fs = 1;
lpf = 1/30;
[fb,fa]=butter(FilterOrder,lpf*2/fs,'low');

for ii=1:length(lstGood)
    iSeg = lstGood(ii);
    Tree = find(im2.edgeSegN==iSeg);
    
    nStart = im2.segEndNodes(iSeg,1);
    nEnd = im2.segEndNodes(iSeg,2);
    lst = find(im2.nodeEdges(Tree,1)==nStart | im2.nodeEdges(Tree,2)==nStart);
    if length(lst)~=1
        error('This should be length 1')
    end
    n0 = setdiff(im2.nodeEdges(Tree(lst),:),nStart);
    nLst = nStart;
    curv = [];
    curvLen = [];
    curvCumLen = [];
    while n0~=nEnd
        lst = find(im2.nodeEdges(Tree,1)==n0 | im2.nodeEdges(Tree,2)==n0);
        if length(lst)~=2
            error('This should be length 2')
        end

        foo = im2.nodeEdges(Tree(lst),:);
        n12 = setdiff(foo(:),n0);
        e1 = im.nodePos(n0,:) - im.nodePos(n12(1),:);
        e2 = im.nodePos(n12(2),:) - im.nodePos(n0,:);
        e3 = im.nodePos(n12(2),:) - im.nodePos(n12(1),:);
        b = norm(e1);
        c = norm(e2);
        e1 = e1 / norm(e1);
        e2 = e2 / norm(e2);
        e3 = e3 / norm(e3);
        epsilon = acos(sum(e1.*e3));
        curv(end+1) = 1 / (b * cos(epsilon) / sin(2*epsilon)); 
        curvLen(end+1) = 0.5*(b+c);
        curvCumLen(end+1) = sum(curvLen);
        
        nLst(end+1) = n0;
        foo = im2.nodeEdges(Tree(lst),:);
        n0 = setdiff(foo(:),nLst);
    end
    
    if length(curv)>=2
        curvTmp = interp1(curvCumLen,curv,curvLen(1):curvCumLen(end),'cubic');
        if length(curvTmp)>3*FilterOrder
            curvTmp2 = filtfilt(fb,fa,curvTmp);
        else
            curvTmp2 = curvTmp;
        end

        %    figure(20);
        %    plot(curvCumLen,curv,'.',curvCumLen(1):curvCumLen(end),curvTmp,'b-',curvCumLen(1):curvCumLen(end),curvTmp2,'r-')
        tortCurv(iSeg) = sum(diff(curvTmp2).^2) / (curvCumLen(end)-curvCumLen(1));
        %    title( sprintf('Tortuosity = %.2e',tortCurv(iSeg)) )
    end
    
%    pause
end    


%%
% histogram length and diameter of "good" segments
% plot diam vs length


figure(2)
subplot(1,3,1)
hist( im2.segLen(lstGood), [25:50:500] )
title( 'Segment length ("good" segments)' )
subplot(1,3,2)
edges = [0 10 20 30 40 50 100];
n=histc( im2.segDiam(lstGood), edges );
bar(edges,n,'histc')
title( 'Segment diameter ("good" segments)' )
subplot(1,3,3)
plot( im2.segLen(lstGood), im2.segDiam(lstGood), '.' )
xlabel( 'Length' )
ylabel( 'Diameter' )

figure(3)
clf
plot( im2.segLen(lstGood), im2.segDiam(lstGood), 'g.', im2.segLen(lstBad), im2.segDiam(lstBad), 'r.' )
xlabel( 'Length' )
ylabel( 'Diameter' )

% plot tortuosity
figure(4)
subplot(1,2,1)
plot( im2.segDiam(lstGood), tort, 'g.', im2.segDiam(lstBad), tortBad, 'r.')
xlabel( 'Diameter' )
ylabel( 'Tortuosity' )
subplot(1,2,2)
plot( im2.segLen(lstGood), tort, 'g.', im2.segLen(lstBad), tortBad, 'r.')
xlabel( 'Len' )
ylabel( 'Tortuosity' )


%%
% plot various parameters binned in X
lstX = [0:50:600];
xx = im2.segPos(lstGood,3);
%lstX = [0:10:60];
%xx = im2.segDiam(lstGood);

yy = im2.tortCurv(lstGood);
%yy = im2.segLen(lstGood);

x = [];
y = [];
y2 = [];
for ii=1:length(lstX)-1
    x(ii) = mean(lstX(ii+[0:1]));
    lst = find(xx>=lstX(ii) & xx<=lstX(ii+1));
    y(ii) = mean(yy(lst));
    y2(ii) = median(yy(lst));
end
figure(20)
[ax,h1,h2]=plotyy(x,y,x,y2);
set(get(ax(1),'ylabel'),'string','Mean')
set(get(ax(2),'ylabel'),'string','Median')

%%
% view segments with high tortuosity
%lst = find(tort<2 & tort>1 & tortAng>300);
%lst = find(im2.segLen>200);
%lst = find(tortCurv>1e-5);
lst = find(tortCurv<1e-7 & tortCurv>0);
lstS = lst;

for ii=1:length(lst)
    lst3 = find(im2.edgeSegN==lstS(ii));
    eFlag = ones(size(im2.nodeEdges,1),1);
    eFlag(lst3) = 2;

    iSeg = lstS(ii);
    
    if 0
        nIdx = im2.segEndNodes(lstS(ii),:);
        [pMin,pMax,Tree]=plotTreeSeg(im2.nodeEdges,im2.nodePos,nIdx(1),15,[1 1],im2.I,[],[],[],eFlag);
        title( sprintf('Node%d  Seg%d  diam=%.0f  len=%.0f  tort=%.1f',nIdx(1),iSeg,im2.segDiam(iSeg),im2.segLen(iSeg),tort(iSeg)) )
    end
    if 1
        figure(4)
        Tree = find(im2.edgeSegN==iSeg);
        lstN = find(im2.nodeSegN==iSeg);
        
        subplot(2,2,1)
        z1 = floor(max(min(im2.nodePos(lstN,3))-10,1));
        z2 = ceil(min(max(im2.nodePos(lstN,3))+10,size(im2.I,2)));
        imagesc( max(im2.I(:,:,z1:z2),[],3) )
        hold on
        for jj=1:length(Tree)
            n1 = im2.nodeEdges(Tree(jj),1);
            n2 = im2.nodeEdges(Tree(jj),2);
            hl=plot([im2.nodePos(n1,1) im2.nodePos(n2,1)],[im2.nodePos(n1,2) im2.nodePos(n2,2)],'r-');
            set(hl,'linewidth',2)
        end
        hold off
        title( sprintf('Node%d  Seg%d  diam=%.0f  len=%.0f  tort=%.1f  tortA=%.0f tortC=%.2f',nIdx(1),iSeg,im2.segDiam(iSeg),im2.segLen(iSeg),tort(iSeg),tortAng(iSeg),tortCurv(iSeg)*1e6 ) )
        
        subplot(2,2,2)
        z1 = floor(max(min(im2.nodePos(lstN,1))-10,1));
        z2 = ceil(min(max(im2.nodePos(lstN,1))+10,size(im2.I,2)));
        imagesc( squeeze(max(im2.I(:,z1:z2,:),[],2)) )
        hold on
        for jj=1:length(Tree)
            n1 = im2.nodeEdges(Tree(jj),1);
            n2 = im2.nodeEdges(Tree(jj),2);
            hl=plot([im2.nodePos(n1,3) im2.nodePos(n2,3)],[im2.nodePos(n1,2) im2.nodePos(n2,2)],'r-');
            set(hl,'linewidth',2)
        end
        hold off

        subplot(2,2,3)
        z1 = floor(max(min(im2.nodePos(lstN,2))-10,1));
        z2 = ceil(min(max(im2.nodePos(lstN,2))+10,size(im2.I,2)));
        imagesc( squeeze(max(im2.I(z1:z2,:,:),[],1))' )
        hold on
        for jj=1:length(Tree)
            n1 = im2.nodeEdges(Tree(jj),1);
            n2 = im2.nodeEdges(Tree(jj),2);
            hl=plot([im2.nodePos(n1,1) im2.nodePos(n2,1)],[im2.nodePos(n1,3) im2.nodePos(n2,3)],'r-');
            set(hl,'linewidth',2)
        end
        hold off
        
    end
    pause
end
