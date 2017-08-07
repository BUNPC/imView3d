function inView3d_AssignSegmentPositionInVascularThree

global im;
%[fname pname] = uigetfile('*.seed','Select seed file with Im structure');
%load([pname fname],'-mat'); % load 'im2' structure

%%

choice = questdlg('Do you want to recalculate segment distance matrices?','Long or quick calculation?','Yes, recalculate all','No, just load pial vessels','No, just load pial vessels');
if strcmp(choice, 'Yes, recalculate all'), % recalculate matrices...
    
    nGroups = length(unique(im.nodeGrp));
    disp(['Number of groups is ' num2str(nGroups) ]);
    disp('Grp. number      num. elements');
    for i=1:nGroups,
        grpElements(i) = length(find(im.nodeGrp==i));
        disp([num2str(i) '         ' num2str( grpElements(i) )  ]);
    end;
    [mg, mi] = max(grpElements);
    answer = inputdlg(['Select group 1 to ' num2str(nGroups)],'Select which group to process',1,{num2str(mi(1))});
    grp = str2num(answer{1});
    im.grpStat.selectedGroupNumber = grp;
    % find list of segments that belong to selected group
    totNumSegments = length(im.segDiam); % number of segments
    foo1 = im.nodeGrp(im.segEndNodes(:,1));
    foo2 = im.nodeGrp(im.segEndNodes(:,2));
    segLst = find(foo1==grp | foo2==grp); % indexes of segments which are in group grp
    numSegments = length(segLst);
    im.grpStat.segLst = segLst;
    im.grpStat.numSegments = numSegments;

    %%
    % start from largest artery (vessel types: A=1, C=2, V=3)
    %idxArt = find(im.segVesType==1);
    %[foo, maxArtIdx] = sort(im.segDiam(idxArt),'descend');
    %maxArtIdx = idxArt(maxArtIdx); % descending indexes of arterial segments
    %%
    segEndNodes = im.segEndNodes(segLst,:); % use only segments that belong to group grp
    endNodes = unique(segEndNodes(:));
    endNodes = sort(endNodes);
    numEndNodes = length(endNodes);
    lookuptable(1:numEndNodes,1)=endNodes;
    im.grpStat.lookuptable = lookuptable;
    im.grpStat.numEndNodes = numEndNodes;

    %%
    %         % create connectivity matrix
    %         segConnectivityMartix = zeros(numEndNodes);
    %         %segConnectivityMartix(numSegments,numSegments)=1;
    %         hp = waitbar(0,'Populating first order segment connectivity matrix'); 
    %         for ii=1:numEndNodes-1,
    %             waitbar(ii/(numEndNodes-1),hp);
    %             node1 = lookuptable(ii);
    %             %%segConnectivityMartix(ii,ii)=1;
    %             %%disp(num2str(ii));
    %             for jj=ii+1:numEndNodes,
    %                 node2 = lookuptable(jj);
    %                 lst = find(im.segEndNodes(:,1) == node1);
    %                 lst2 = find(im.segEndNodes(:,2) == node2);
    %                 c = intersect(lst, lst2);
    %                 if ~isempty(c),
    %                     segConnectivityMartix(ii,jj)=im.segLen_um(c(1));
    %                     segConnectivityMartix(jj,ii)=segConnectivityMartix(ii,jj);
    %                 else
    %                     lst = find(im.segEndNodes(:,1) == node2);
    %                     lst2 = find(im.segEndNodes(:,2) == node1);
    %                     c = intersect(lst, lst2);
    %                     if ~isempty(c),
    %                         segConnectivityMartix(ii,jj)=im.segLen_um(c(1));
    %                         segConnectivityMartix(jj,ii)=segConnectivityMartix(ii,jj);
    %                     end;
    %                 end;
    %             end;
    %         end;
    %         close(hp);
    %         segConnectivityMartix_1stOrder = segConnectivityMartix;

    %%
    % create matrix E(1:numSegments,1:3) with elements (node1, node2, segmentLEngth)
    segLen_um = im.segLen_um(segLst);
    E = zeros(2*numSegments,3); % 2 times for two directions of arrow
    for ii=1:numSegments,
        nodes = segEndNodes(ii,:);
        E(2*ii-1,1)=find(lookuptable == nodes(1));
        E(2*ii-1,2)=find(lookuptable == nodes(2));
        E(2*ii-1,3)=segLen_um(ii);

        E(2*ii,1)=E(2*ii-1,2);
        E(2*ii,2)=E(2*ii-1,1);
        E(2*ii,3)=E(2*ii-1,3);
    end;

    %%
    % calculate shortest distance between nodes
    m=2*numSegments; % number of oriented edges
    n = numEndNodes;
    %[m,n,E] = grValidation(E); % E data validation

    % ================ Initial values ===============
    dSP=ones(n)*inf; % initial distances
    dSP((E(:,2)-1)*n+E(:,1))=E(:,3);
    %dSP0=dSP;
    % ========= The main cycle of Floyd-Warshall algorithm =========
    hp = waitbar(0,'Calculating minimum distance for all graph nodes...'); 
    for j=1:n,
      waitbar(j/n,hp);
      i=setdiff((1:n),j);
      dSP(i,i)=min(dSP(i,i),repmat(dSP(i,j),1,n-1)+repmat(dSP(j,i),n-1,1));
    end
    close(hp);

    %%
    % find shortest path from node1 ('s') to node2 ('t')

    % populate matrix nSP (numEndNodes x numEndNodes), where nSP(i,j) is index of node
    % before 'j' in shortest path from i to j

    nSP = zeros(numEndNodes);
    dSP1=dSP;
    dSP1(1:n+1:n^2)=0; % modified dSP

    hp = waitbar(0,'Calculating shortest paths for all pairs of graph nodes...'); 
    for ii=1:numEndNodes-1,
        waitbar(ii/(numEndNodes-1),hp);
        s=ii;
        %s=1;
        for jj=ii+1:numEndNodes,
            t=jj;
            %t=915;

            sp=[];
            %s=s(1);
            %t=t(1);
            if isinf(dSP(s,t)), % t is not accessible from s
                disp('Grrrrr! This pair of nodes is not connected');
            else

                l=ones(m,1); % label for each arrow
                sp=t; % final vertex
                while ~(sp(1)==s),
                  nv=find((E(:,2)==sp(1))&l); % nv contains all labeled arrows (segments) pointing at node sp(1)
                  vnv=abs((dSP1(s,sp(1))-dSP1(s,E(nv,1)))'-E(nv,3))<eps*1e4; % valided arrows (I changed 1e3 to 1e4)
                  l(nv(~vnv))=0; % labels of not valided arrows
                  if all(~vnv), % invalided arrows
                    l(find((E(:,1)==sp(1))&(E(:,2)==sp(2))))=0; 
                    sp=sp(2:end); % one step back
                  else
                    nv=nv(vnv); % rested valided arrows
                    sp=[E(nv(1),1) sp]; % add one vertex to shortest path
                  end
                end;
            end;
            nSP(ii,jj) = sp(end-1);
            nSP(jj,ii) = sp(2);

        end;
    end;

    for ii=1:numEndNodes, % make sure that distances on diagonal are zero
        dSP(ii,ii)=0;
    end;
    im.grpStat.dSP = dSP;
    im.grpStat.nSP = nSP;

    close(hp);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Above was calculation based on distances. Now for branching order
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%
    % create matrix E(1:numSegments,1:3) with elements (node1, node2, segmentLEngth)

    %Ebr = zeros(2*numSegments,3); % 2 times for two directions of arrow
    Ebr = E;
    Ebr(:,3) = 1; % this will assure that we calculate shortest paths based on branching order

    %%
    % calculate shortest distance between nodes
    m=2*numSegments; % number of oriented edges
    n = numEndNodes;
    %[m,n,E] = grValidation(E); % E data validation

    % ================ Initial values ===============
    dSPbr=ones(n)*inf; % initial distances
    dSPbr((Ebr(:,2)-1)*n+Ebr(:,1))=Ebr(:,3);
    %dSP0=dSP;
    % ========= The main cycle of Floyd-Warshall algorithm =========
    hp = waitbar(0,'Calculating minimum distance for all graph nodes...'); 
    for j=1:n,
      waitbar(j/n,hp);
      i=setdiff((1:n),j);
      dSPbr(i,i)=min(dSPbr(i,i),repmat(dSPbr(i,j),1,n-1)+repmat(dSPbr(j,i),n-1,1));
    end
    close(hp);

    %%
    % find shortest path from node1 ('s') to node2 ('t')

    % populate matrix nSP (numEndNodes x numEndNodes), where nSP(i,j) is index of node
    % before 'j' in shortest path from i to j

    nSPbr = zeros(numEndNodes);
    dSP1br=dSPbr;
    dSP1br(1:n+1:n^2)=0; % modified dSP

    hp = waitbar(0,'Calculating shortest paths for all pairs of graph nodes...'); 
    for ii=1:numEndNodes-1,
        waitbar(ii/(numEndNodes-1),hp);
        s=ii;
        %s=1;
        for jj=ii+1:numEndNodes,
            t=jj;
            %t=915;

            sp=[];
            %s=s(1);
            %t=t(1);
            if isinf(dSPbr(s,t)), % t is not accessible from s
                disp('Grrrrr! This pair of nodes is not connected');
            else

                l=ones(m,1); % label for each arrow
                sp=t; % final vertex
                while ~(sp(1)==s),
                  nv=find((Ebr(:,2)==sp(1))&l); % nv contains all labeled arrows (segments) pointing at node sp(1)
                  vnv=abs((dSP1br(s,sp(1))-dSP1br(s,Ebr(nv,1)))'-Ebr(nv,3))<eps*1e4; % valided arrows (I changed 1e3 to 1e4)
                  l(nv(~vnv))=0; % labels of not valided arrows
                  if all(~vnv), % invalided arrows
                    l(find((Ebr(:,1)==sp(1))&(Ebr(:,2)==sp(2))))=0; 
                    sp=sp(2:end); % one step back
                  else
                    nv=nv(vnv); % rested valided arrows
                    sp=[Ebr(nv(1),1) sp]; % add one vertex to shortest path
                  end
                end;
            end;
            nSPbr(ii,jj) = sp(end-1);
            nSPbr(jj,ii) = sp(2);

        end;
    end;

    for ii=1:numEndNodes, % make sure that distances on diagonal are zero
        dSPbr(ii,ii)=0;
    end;
    im.grpStat.dSPbr = dSPbr;
    im.grpStat.nSPbr = nSPbr;

    close(hp);

end; % recalculate matrices...



% NOW, IMPORT NEW PIAL VESSEL BRANCH ORDER AND CALCULATE SEGMENT DISTANCES
% TO THESE PIAL SEGMENTS...

%[fname pname] = uigetfile('*.seed','Select seed file with Im structure');
%load([pname fname],'-mat'); % load 'im2' structure
%%
if ~isfield(im,'grpStat'),
    disp('No info on group statistics!');
    disp('Run Calc. Segment Distances first!');
    return;
end;
%%
grp = im.grpStat.selectedGroupNumber; 
segLst = im.grpStat.segLst;
numSegments = im.grpStat.numSegments;
lookuptable = im.grpStat.lookuptable;
numEndNodes = im.grpStat.numEndNodes;
dSP = im.grpStat.dSP;
nSP = im.grpStat.nSP;
dSPbr = im.grpStat.dSPbr;
nSPbr = im.grpStat.nSPbr;



%%

% First, supply segment numbers in pial arteries and veins together
% with branch order. Algorithm will calculate branching order for each 
% vessel segment in selected connected group of segments


%Bill Sep 10 2015. avoid making the text file with pial segment seeds by
%having these segments labeled in the GUI and saved in im.grpStat.pialseg
% only 3 col.

if isfield(im.grpStat, 'pialseg')
pialseg=im.grpStat.pialseg;
else
% load file with info about pial segments
[ff pp] = uigetfile('*.txt','Select File name with Pial Segments info');
pialseg = load([pp ff],'-ascii'); % 3 columns, segment number, branch order, vessel type (A=1, 2=C, 3=V)
end

for ii=1:size(pialseg,1), % find end nodes of pial segments and convert them based on lookuptable
    pialseg(ii,4) = find(lookuptable == im.segEndNodes(pialseg(ii,1),1));
    pialseg(ii,5) = find(lookuptable == im.segEndNodes(pialseg(ii,1),2));
end;
im.grpStat.pialseg = pialseg;



% first calculate branch order for each segment in selected graph group -
% this is not neccessary for PO2, but it is important for testing the
% algorithm
segBranchOrder = zeros(numSegments,6); 
% columns are: true segment number, branch order from A, length from A, BR
% order from V, length from V, diameter in 'um'

% BRANCHING ORDER FROM PIAL ARTERIES
maxBR = max(dSPbr(:));
for ii=1:numSegments,
    segBranchOrder(ii,1) = segLst(ii);
    segBranchOrder(ii,6) = im.segDiam(segLst(ii)); % segment diameters
    segnode1 = find(lookuptable == im.segEndNodes(segLst(ii),1));
    segnode2 = find(lookuptable == im.segEndNodes(segLst(ii),2));
    % find minimum branch distance to each node of pial arterial segments
    lstArt = find(pialseg(:,3)==1);
    minBR = maxBR+1;
    sameAsPial = 0;
    for jj = 1:length(lstArt),
        Anode1 = pialseg(lstArt(jj),4);
        Anode2 = pialseg(lstArt(jj),5);
        if ( ((Anode1 == segnode1) && (Anode2 == segnode2)) || ((Anode1 == segnode2) && (Anode2 == segnode1))  ),
            segBranchOrder(ii,2) = pialseg(lstArt(jj),2); % same branch order as pial artery
            %segBranchOrder(ii,3) = 0; % length from pial artery is 0 -  not good, it should look for '0' order arteries
            sameAsPial = 1;
            break;
        else            
            foo1 = dSPbr(segnode1,Anode1);
            foo2 = dSPbr(segnode1,Anode2);
            foo3 = dSPbr(segnode1,Anode1);
            foo4 = dSPbr(segnode1,Anode2);
            foo = min([foo1 foo2 foo3 foo4]);
            if foo<minBR,
                minBR = foo;
                segBRord = pialseg(lstArt(jj),2);
            end;
        end;
    end;
    if ~sameAsPial,
        segBranchOrder(ii,2) = minBR+segBRord+1; 
    end;
end;


% DISTANCE FROM SURFACE PIAL ARTERIES WITH BRANCHING ORDER ZERO
maxDIS = 1000000; % some large number of micrones
for ii=1:numSegments,
    segnode1 = find(lookuptable == im.segEndNodes(segLst(ii),1));
    segnode2 = find(lookuptable == im.segEndNodes(segLst(ii),2));
    % find minimum distance to each node of pial arterial segments
    lstArt = find( (pialseg(:,3)==1) & (pialseg(:,2) == 0) ); % find arterial pial vessels with branching order = 0
    minDIS = maxDIS+1;
    sameAsPial = 0;
    for jj = 1:length(lstArt),
        Anode1 = pialseg(lstArt(jj),4);
        Anode2 = pialseg(lstArt(jj),5);
        if ( ((Anode1 == segnode1) && (Anode2 == segnode2)) || ((Anode1 == segnode2) && (Anode2 == segnode1))  ),
            segBranchOrder(ii,3) = 0; % segment equal to pial atery with BR order zero, so length is 0
            sameAsPial = 1;
            break;
        else            
            foo1 = dSP(segnode1,Anode1);
            foo2 = dSP(segnode1,Anode2);
            foo3 = dSP(segnode1,Anode1);
            foo4 = dSP(segnode1,Anode2);
            foo = min([foo1 foo2 foo3 foo4]);
            if foo<minDIS,
                minDIS = foo;
                segDIS = dSP(Anode1, Anode2);
            end;
        end;
    end;
    if ~sameAsPial,
        segBranchOrder(ii,3) = minDIS+segDIS/2.0+dSP(segnode1,segnode2)/2.0; % add 1/2 of lengths of segments which are compared
    end;
end;


% BRANCHING ORDER FROM PIAL VEINS
maxBR = max(dSPbr(:));
for ii=1:numSegments,
    segnode1 = find(lookuptable == im.segEndNodes(segLst(ii),1));
    segnode2 = find(lookuptable == im.segEndNodes(segLst(ii),2));
    % find minimum branch distance to each node of pial arterial segments
    lstVein = find(pialseg(:,3)==3);
    minBR = maxBR+1;
    sameAsPial = 0;
    for jj = 1:length(lstVein),
        Vnode1 = pialseg(lstVein(jj),4);
        Vnode2 = pialseg(lstVein(jj),5);
        if ( ((Vnode1 == segnode1) && (Vnode2 == segnode2)) || ((Vnode1 == segnode2) && (Vnode2 == segnode1))  ),
            segBranchOrder(ii,4) = pialseg(lstVein(jj),2); % same branch order as pial vein
            sameAsPial = 1;
            break;
        else            
            foo1 = dSPbr(segnode1,Vnode1);
            foo2 = dSPbr(segnode1,Vnode2);
            foo3 = dSPbr(segnode1,Vnode1);
            foo4 = dSPbr(segnode1,Vnode2);
            foo = min([foo1 foo2 foo3 foo4]);
            if foo<minBR,
                minBR = foo;
                segBRord = pialseg(lstVein(jj),2);
            end;
        end;
    end;
    if ~sameAsPial,
        segBranchOrder(ii,4) = minBR+segBRord+1; 
    end;
end;


% DISTANCE FROM SURFACE PIAL VEINS WITH BRANCHING ORDER ZERO

% this needs improvement! For example, there may not be zero branches of
% veins in FOV. It should look to min distance to given vessels in pialseg,
% but then it should know distances of vessels inside pialseg from surface
% pial vessels!

maxDIS = 1000000; % some large number of micrones
for ii=1:numSegments,
    segnode1 = find(lookuptable == im.segEndNodes(segLst(ii),1));
    segnode2 = find(lookuptable == im.segEndNodes(segLst(ii),2));
    % find minimum distance to each node of pial venous segments
    lstVein = find( (pialseg(:,3)==3) & (pialseg(:,2) == 0) ); % find venous pial vessels with branching order = 0
    minDIS = maxDIS+1;
    sameAsPial = 0;
    for jj = 1:length(lstVein),
        Vnode1 = pialseg(lstVein(jj),4);
        Vnode2 = pialseg(lstVein(jj),5);
        if ( ((Vnode1 == segnode1) && (Vnode2 == segnode2)) || ((Vnode1 == segnode2) && (Vnode2 == segnode1))  ),
            segBranchOrder(ii,5) = 0; % segment equal to pial atery with BR order zero, so length is 0
            sameAsPial = 1;
            break;
        else            
            foo1 = dSP(segnode1,Vnode1);
            foo2 = dSP(segnode1,Vnode2);
            foo3 = dSP(segnode1,Vnode1);
            foo4 = dSP(segnode1,Vnode2);
            foo = min([foo1 foo2 foo3 foo4]);
            if foo<minDIS,
                minDIS = foo;
                segDIS = dSP(Vnode1, Vnode2);
            end;
        end;
    end;
    if ~sameAsPial,
        segBranchOrder(ii,5) = minDIS+segDIS/2.0+dSP(segnode1,segnode2)/2.0; % add 1/2 of lengths of segments which are compared
    end;
end;


im.grpStat.segBranchOrder = segBranchOrder;

% assign branching orders to all edges in the graph...
for ee=1:length(im.edgeSegN),
    [Lia, Locb] = ismember(im.edgeSegN(ee),im.grpStat.segLst);
    if Lia, % find if edge ee belongs to selected group
        im.edgeBRorderVeins(ee) = im.grpStat.segBranchOrder(Locb,4);
        im.edgeBRorderArt(ee) = im.grpStat.segBranchOrder(Locb,2);
    else % edge outside of desired group
        im.edgeBRorderVeins(ee) = 0;
        im.edgeBRorderArt(ee) = 0; % I may need better number for edges outside the group. This will make them branch order zero
    end; 
    %im.edgeBRorderVeins(ee) = segBranchOrder(im.edgeSegN(ee),4);
    %im.edgeBRorderArt(ee) = segBranchOrder(im.edgeSegN(ee),2);
end;


% %%
% 
% % calculate branch order from pial arteries for each vessel where we have
% % pO2 value measured. First, supply edge numbers in pial arteries together
% % with branch order of the segments with these edges. Algorithm will assign
% % segment numbers to these edges and calculate for each pO2 segment minimum
% % branching order.
% 
% % load file with info about pial segments
% [ff pp] = uigetfile('*.txt','Select File name with Pial Segments info');
% pialseg = load([pp ff],'-ascii'); % 3 columns, segment number, branch order, vessel type (A=1, 2=C, 3=V)
% 
% % make sure that dSPbr diagonal is zero
% for ii = 1:numEndNodes,
%     dSPbr(ii,ii)=0;
% end;
% 
% for ii=1:size(pialseg,1), % find end nodes of pial segments and convert them based on lookuptable
%     pialseg(ii,4) = find(lookuptable == im.segEndNodes(pialseg(ii,1),1));
%     pialseg(ii,5) = find(lookuptable == im.segEndNodes(pialseg(ii,1),2));
% end;
% %%
% % first calculate branch order for each segment in selected graph group -
% % this is not neccessary for PO2, but it is important for testing the
% % algorithm
% segBR_Art = zeros(numSegments,6); 
% % columns are: true segment number, branch order from A, length from A, BR
% % order from V, length from V, diameter
% 
% % BRANCHING ORDER FROM PIAL ARTERIES
% maxBR = max(dSPbr(:));
% for ii=1:numSegments,
%     segBR_Art(ii,1) = segLst(ii);
%     segBR_Art(ii,6) = im.segDiam(segLst(ii)); % segment diameters
%     segnode1 = find(lookuptable == im.segEndNodes(segLst(ii),1));
%     segnode2 = find(lookuptable == im.segEndNodes(segLst(ii),2));
%     % find minimum branch distance to each node of pial arterial segments
%     lstArt = find(pialseg(:,3)==1);
%     minBR = maxBR+1;
%     sameAsPial = 0;
%     for jj = 1:length(lstArt),
%         Anode1 = pialseg(lstArt(jj),4);
%         Anode2 = pialseg(lstArt(jj),5);
%         if ( ((Anode1 == segnode1) && (Anode2 == segnode2)) || ((Anode1 == segnode2) && (Anode2 == segnode1))  ),
%             segBR_Art(ii,2) = pialseg(lstArt(jj),2); % same branch order as pial artery
%             %segBR_Art(ii,3) = 0; % length from pial artery is 0 -  not good, it should look for '0' order arteries
%             sameAsPial = 1;
%             break;
%         else            
%             foo1 = dSPbr(segnode1,Anode1);
%             foo2 = dSPbr(segnode1,Anode2);
%             foo3 = dSPbr(segnode1,Anode1);
%             foo4 = dSPbr(segnode1,Anode2);
%             foo = min([foo1 foo2 foo3 foo4]);
%             if foo<minBR,
%                 minBR = foo;
%                 segBRord = pialseg(lstArt(jj),2);
%             end;
%         end;
%     end;
%     if ~sameAsPial,
%         segBR_Art(ii,2) = minBR+segBRord+1; 
%     end;
% end;
% 
% %%
% 
% % make sure that dSP diagonal is zero
% for ii = 1:numEndNodes,
%     dSP(ii,ii)=0;
% end;
% 
% % DISTANCE FROM SURFACE PIAL ARTERIES WITH BRANCHING ORDER ZERO
% maxDIS = 1000000; % some large number of micrones
% for ii=1:numSegments,
%     segnode1 = find(lookuptable == im.segEndNodes(segLst(ii),1));
%     segnode2 = find(lookuptable == im.segEndNodes(segLst(ii),2));
%     % find minimum distance to each node of pial arterial segments
%     lstArt = find( (pialseg(:,3)==1) & (pialseg(:,2) == 0) ); % find arterial pial vessels with branching order = 0
%     minDIS = maxDIS+1;
%     sameAsPial = 0;
%     for jj = 1:length(lstArt),
%         Anode1 = pialseg(lstArt(jj),4);
%         Anode2 = pialseg(lstArt(jj),5);
%         if ( ((Anode1 == segnode1) && (Anode2 == segnode2)) || ((Anode1 == segnode2) && (Anode2 == segnode1))  ),
%             segBR_Art(ii,3) = 0; % segment equal to pial atery with BR order zero, so length is 0
%             sameAsPial = 1;
%             break;
%         else            
%             foo1 = dSP(segnode1,Anode1);
%             foo2 = dSP(segnode1,Anode2);
%             foo3 = dSP(segnode1,Anode1);
%             foo4 = dSP(segnode1,Anode2);
%             foo = min([foo1 foo2 foo3 foo4]);
%             if foo<minDIS,
%                 minDIS = foo;
%                 segDIS = dSP(Anode1, Anode2);
%             end;
%         end;
%     end;
%     if ~sameAsPial,
%         segBR_Art(ii,3) = minDIS+segDIS/2.0+dSP(segnode1,segnode2)/2.0; % add 1/2 of lengths of segments which are compared
%     end;
% end;
% 
% %%
% % BRANCHING ORDER FROM PIAL VEINS
% maxBR = max(dSPbr(:));
% for ii=1:numSegments,
%     segnode1 = find(lookuptable == im.segEndNodes(segLst(ii),1));
%     segnode2 = find(lookuptable == im.segEndNodes(segLst(ii),2));
%     % find minimum branch distance to each node of pial arterial segments
%     lstVein = find(pialseg(:,3)==3);
%     minBR = maxBR+1;
%     sameAsPial = 0;
%     for jj = 1:length(lstVein),
%         Vnode1 = pialseg(lstVein(jj),4);
%         Vnode2 = pialseg(lstVein(jj),5);
%         if ( ((Vnode1 == segnode1) && (Vnode2 == segnode2)) || ((Vnode1 == segnode2) && (Vnode2 == segnode1))  ),
%             segBR_Art(ii,4) = pialseg(lstVein(jj),2); % same branch order as pial vein
%             sameAsPial = 1;
%             break;
%         else            
%             foo1 = dSPbr(segnode1,Vnode1);
%             foo2 = dSPbr(segnode1,Vnode2);
%             foo3 = dSPbr(segnode1,Vnode1);
%             foo4 = dSPbr(segnode1,Vnode2);
%             foo = min([foo1 foo2 foo3 foo4]);
%             if foo<minBR,
%                 minBR = foo;
%                 segBRord = pialseg(lstVein(jj),2);
%             end;
%         end;
%     end;
%     if ~sameAsPial,
%         segBR_Art(ii,4) = minBR+segBRord+1; 
%     end;
% end;
% 
% %%
% % DISTANCE FROM SURFACE PIAL VEINS WITH BRANCHING ORDER ZERO
% 
% % this needs improvement! For example, there may not be zero branches of
% % veins in FOV. It should look to min distance to given vessels in pialseg,
% % but then it should know distances of vessels inside pialseg from surface
% % pial vessels!
% 
% maxDIS = 1000000; % some large number of micrones
% for ii=1:numSegments,
%     segnode1 = find(lookuptable == im.segEndNodes(segLst(ii),1));
%     segnode2 = find(lookuptable == im.segEndNodes(segLst(ii),2));
%     % find minimum distance to each node of pial venous segments
%     lstVein = find( (pialseg(:,3)==3) & (pialseg(:,2) == 0) ); % find venous pial vessels with branching order = 0
%     minDIS = maxDIS+1;
%     sameAsPial = 0;
%     for jj = 1:length(lstVein),
%         Vnode1 = pialseg(lstVein(jj),4);
%         Vnode2 = pialseg(lstVein(jj),5);
%         if ( ((Vnode1 == segnode1) && (Vnode2 == segnode2)) || ((Vnode1 == segnode2) && (Vnode2 == segnode1))  ),
%             segBR_Art(ii,5) = 0; % segment equal to pial atery with BR order zero, so length is 0
%             sameAsPial = 1;
%             break;
%         else            
%             foo1 = dSP(segnode1,Vnode1);
%             foo2 = dSP(segnode1,Vnode2);
%             foo3 = dSP(segnode1,Vnode1);
%             foo4 = dSP(segnode1,Vnode2);
%             foo = min([foo1 foo2 foo3 foo4]);
%             if foo<minDIS,
%                 minDIS = foo;
%                 segDIS = dSP(Vnode1, Vnode2);
%             end;
%         end;
%     end;
%     if ~sameAsPial,
%         segBR_Art(ii,5) = minDIS+segDIS/2.0+dSP(segnode1,segnode2)/2.0; % add 1/2 of lengths of segments which are compared
%     end;
% end;
% 
% %%
% % now calculate PO2 distribution by branching order from pial arteries
%     
% 
% for ii = 1:length(im.PO2pts.PO2),
%     po2Seg = im.PO2pts.seg(ii);
%     po2_BRorder(ii) = segBR_Art(find(segBR_Art(:,1)==po2Seg),2);
%     po2_PO2(ii) = im.PO2pts.PO2(ii);
% end;
% BRorderIdx = sort(unique(po2_BRorder));
% for ii=1:length(BRorderIdx),
%     foolst = find(po2_BRorder == BRorderIdx(ii));
%     averArtPO2_BRorder(ii) = mean(po2_PO2(foolst));
%     if length(foolst) > 1,
%         averArtPO2err_BRorder(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
%     else
%         averArtPO2err_BRorder(ii)=0;
%     end;
% end;
% figure; 
% errorbar(BRorderIdx,averArtPO2_BRorder,averArtPO2err_BRorder,'.-');
% xlabel('Branching Order from Pial Arteries');
% ylabel('Average PO2 (mmHg)');
% 
% %%
% % now calculate PO2 distribution by branching order from pial arteries
% % BUT EXCLUDE SEGMENTS LABELED AS VEINS (EXCLUDE SEGMENTS LABELED WITH VESSEL TYPE = 3). Find segment diameters too
%     
% n=0;
% clear po2_BRorder;
% clear po2_PO2;
% clear po2_SO2;
% clear averArtPO2_BRorder;
% clear averArtPO2err_BRorder;
% clear averArtDiam_BRorder;
% clear averArtDiamerr_BRorder;
% clear averSO2;
% clear averSO2err;
% 
% hillC = 2.7; p50 = 37; % Rat, Ellis C.G. et al., Am. J. Phys Heart Circ Phys 258, H1216-, 2002
% hillC = 2.59; p50 = 40.2; % C57BL/6 mice, Uchida K. et al., Zoological Science 15, 703-706, 1998
% for ii = 1:length(im.PO2pts.PO2),
%     po2Seg = im.PO2pts.seg(ii);
%     %if (segBR_Art(find(segBR_Art(:,1)==po2Seg),4) > 4), %try to exclude
%     %venules, veins, everything labeled as venous side
%     if( (im.segVesType(po2Seg) == 1) || (im.segVesType(po2Seg) ==2) ), % if arterial or cappilary
%         n=n+1;
%         po2_BRorder(n) = segBR_Art(find(segBR_Art(:,1)==po2Seg),2);
%         po2_PO2(n) = im.PO2pts.PO2(ii);
%         po2_SO2(n) = (po2_PO2(n))^hillC / ( (po2_PO2(n))^hillC + p50^hillC  ); 
%         segDiam(n) = segBR_Art(find(segBR_Art(:,1)==po2Seg),6);
%     end;
% end;
% BRorderIdx = sort(unique(po2_BRorder));
% for ii=1:length(BRorderIdx),
%     foolst = find(po2_BRorder == BRorderIdx(ii));
%     averArtPO2_BRorder(ii) = mean(po2_PO2(foolst));
%     averArtDiam_BRorder(ii) = mean(segDiam(foolst));
%     averSO2(ii) = mean(po2_SO2(foolst));
%     if length(foolst) > 1,
%         averArtPO2err_BRorder(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
%         averArtDiamerr_BRorder(ii) = std(segDiam(foolst))/sqrt(length(foolst)-1);
%         averSO2err(ii) = std(po2_SO2(foolst))/sqrt(length(foolst)-1);
%     else
%         averArtPO2err_BRorder(ii)=0;
%         averArtDiamerr_BRorder(ii) = 0;
%         averSO2err(ii) = 0;
%     end;
% end;
% fig = figure; 
% [AX, H1, H2] = plotyy(BRorderIdx,averArtPO2_BRorder,BRorderIdx,averArtDiam_BRorder);
% 
% set(fig,'CurrentAxes',AX(1));
% hold on;
% errorbar(BRorderIdx,averArtPO2_BRorder,averArtPO2err_BRorder,'.-','Color','r');
% axis([min(BRorderIdx) max(BRorderIdx) 0 120]);
% 
% set(fig,'CurrentAxes',AX(2));
% hold on;
% errorbar(BRorderIdx,averArtDiam_BRorder,averArtDiamerr_BRorder,'.-b');
% axis([min(BRorderIdx) max(BRorderIdx) 5 35]);
% 
% set(AX(1),'YColor','r','YTick',0:20:120);
% set(get(AX(1),'YLabel'),'String','pO2 (mmHg)','FontSize',15);
% set(AX(2),'YColor','b');
% set(get(AX(2),'YLabel'),'String','Vessel Diameter (um)','FontSize',15);
% set(get(AX(2),'XLabel'),'String','Branch Order from Pial Arteries','FontSize',15);
% grid on;
% print(fig,'-dpsc2','PO2_vs_BranchOrder_ARTERIES.ps');
% %print(fig,'-djpeg','PO2_vs_BranchOrder.jpg','-r300');
% 
% 
% fig = figure; 
% [AX, H1, H2] = plotyy(BRorderIdx,averSO2.*100,BRorderIdx,averArtDiam_BRorder);
% 
% set(fig,'CurrentAxes',AX(1));
% hold on;
% errorbar(BRorderIdx,averSO2.*100,averSO2err.*100,'.-','Color','r');
% axis([min(BRorderIdx) max(BRorderIdx) 0 120]);
% 
% set(fig,'CurrentAxes',AX(2));
% hold on;
% errorbar(BRorderIdx,averArtDiam_BRorder,averArtDiamerr_BRorder,'.-b');
% axis([min(BRorderIdx) max(BRorderIdx) 5 35]);
% 
% set(AX(1),'YColor','r','YTick',0:20:120);
% set(get(AX(1),'YLabel'),'String','SO2 (%)','FontSize',15);
% set(AX(2),'YColor','b','YTick',5:5:35);
% set(get(AX(2),'YLabel'),'String','Vessel Diameter (um)','FontSize',15);
% set(get(AX(2),'XLabel'),'String','Branch Order from Pial Arteries','FontSize',15);
% grid on;
% print(fig,'-dpsc2','SO2_vs_BranchOrder_ARTERIES.ps');
% %print(fig,'-djpeg','SO2_vs_BranchOrder.jpg','-r300');
% 
% 
% 
% %%
% % now calculate PO2 distribution by branching order from pial veins
% % BUT EXCLUDE SEGMENTS LABELED AS ARTERIES (EXCLUDE SEGMENTS LABELED WITH VESSEL TYPE = 1). Find segment diameters too
%     
% n=0;
% clear po2_BRorder;
% clear po2_PO2;
% clear po2_SO2;
% clear averVeinPO2_BRorder;
% clear averVeinPO2err_BRorder;
% clear averVeinDiam_BRorder;
% clear averVeinDiamerr_BRorder;
% clear averSO2;
% clear averSO2err;
% 
% hillC = 2.7; p50 = 37; % Rat, Ellis C.G. et al., Am. J. Phys Heart Circ Phys 258, H1216-, 2002
% hillC = 2.59; p50 = 40.2; % C57BL/6 mice, Uchida K. et al., Zoological Science 15, 703-706, 1998
% for ii = 1:length(im.PO2pts.PO2),
%     po2Seg = im.PO2pts.seg(ii);
%     %if (segBR_Art(find(segBR_Art(:,1)==po2Seg),4) > 4), %try to exclude
%     %venules, veins, everything labeled as venous side
%     if( (im.segVesType(po2Seg) == 3) || (im.segVesType(po2Seg) ==2) ), % if vein or cappilary
%         n=n+1;
%         po2_BRorder(n) = segBR_Art(find(segBR_Art(:,1)==po2Seg),2);
%         po2_PO2(n) = im.PO2pts.PO2(ii);
%         po2_SO2(n) = (po2_PO2(n))^hillC / ( (po2_PO2(n))^hillC + p50^hillC  ); 
%         segDiam(n) = segBR_Art(find(segBR_Art(:,1)==po2Seg),6);
%     end;
% end;
% BRorderIdx = sort(unique(po2_BRorder));
% for ii=1:length(BRorderIdx),
%     foolst = find(po2_BRorder == BRorderIdx(ii));
%     averVeinPO2_BRorder(ii) = mean(po2_PO2(foolst));
%     averVeinDiam_BRorder(ii) = mean(segDiam(foolst));
%     averSO2(ii) = mean(po2_SO2(foolst));
%     if length(foolst) > 1,
%         averVeinPO2err_BRorder(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
%         averVeinDiamerr_BRorder(ii) = std(segDiam(foolst))/sqrt(length(foolst)-1);
%         averSO2err(ii) = std(po2_SO2(foolst))/sqrt(length(foolst)-1);
%     else
%         averVeinPO2err_BRorder(ii)=0;
%         averVeinDiamerr_BRorder(ii) = 0;
%         averSO2err(ii) = 0;
%     end;
% end;
% fig = figure; 
% [AX, H1, H2] = plotyy(BRorderIdx,averVeinPO2_BRorder,BRorderIdx,averVeinDiam_BRorder);
% 
% set(fig,'CurrentAxes',AX(1));
% hold on;
% errorbar(BRorderIdx,averVeinPO2_BRorder,averVeinPO2err_BRorder,'.-','Color','r');
% axis([min(BRorderIdx) max(BRorderIdx) 0 120]);
% 
% set(fig,'CurrentAxes',AX(2));
% hold on;
% errorbar(BRorderIdx,averVeinDiam_BRorder,averVeinDiamerr_BRorder,'.-b');
% axis([min(BRorderIdx) max(BRorderIdx) 5 35]);
% 
% set(AX(1),'YColor','r','YTick',0:20:120);
% set(get(AX(1),'YLabel'),'String','pO2 (mmHg)','FontSize',15);
% set(AX(2),'YColor','b');
% set(get(AX(2),'YLabel'),'String','Vessel Diameter (um)','FontSize',15);
% set(get(AX(2),'XLabel'),'String','Branch Order from Pial Veins','FontSize',15);
% grid on;
% print(fig,'-dpsc2','PO2_vs_BranchOrder_VEINS.ps');
% %print(fig,'-djpeg','PO2_vs_BranchOrder.jpg','-r300');
% 
% 
% fig = figure; 
% [AX, H1, H2] = plotyy(BRorderIdx,averSO2.*100,BRorderIdx,averVeinDiam_BRorder);
% 
% set(fig,'CurrentAxes',AX(1));
% hold on;
% errorbar(BRorderIdx,averSO2.*100,averSO2err.*100,'.-','Color','r');
% axis([min(BRorderIdx) max(BRorderIdx) 0 120]);
% 
% set(fig,'CurrentAxes',AX(2));
% hold on;
% errorbar(BRorderIdx,averVeinDiam_BRorder,averVeinDiamerr_BRorder,'.-b');
% axis([min(BRorderIdx) max(BRorderIdx) 5 35]);
% 
% set(AX(1),'YColor','r','YTick',0:20:120);
% set(get(AX(1),'YLabel'),'String','SO2 (%)','FontSize',15);
% set(AX(2),'YColor','b','YTick',5:5:35);
% set(get(AX(2),'YLabel'),'String','Vessel Diameter (um)','FontSize',15);
% set(get(AX(2),'XLabel'),'String','Branch Order from Pial Veins','FontSize',15);
% grid on;
% print(fig,'-dpsc2','SO2_vs_BranchOrder_VEINS.ps');
% %print(fig,'-djpeg','SO2_vs_BranchOrder.jpg','-r300');
% 
% 
% 
% 
% %%
% % now calculate PO2 distribution by distance from pial arteries
%     
% 
% for ii = 1:length(im.PO2pts.PO2),
%     po2Seg = im.PO2pts.seg(ii);
%     po2_DIS(ii) = segBR_Art(find(segBR_Art(:,1)==po2Seg),3);
%     po2_PO2(ii) = im.PO2pts.PO2(ii);
% end;
% % bin in 100 um bins...
% dDIS = 100;
% maxDIS = max(po2_DIS)+1;
% minDIS = min(po2_DIS)-1;
% nDIS = floor((maxDIS-minDIS)/dDIS);
% for ii=1:nDIS,
%     foolst = find( (po2_DIS<ii*dDIS) & (po2_DIS>=(ii-1)*dDIS) );
%     DIST(ii) = (ii-0.5)*dDIS;
%     if ~isempty(foolst),
%         DIST_PO2(ii) = mean(po2_PO2(foolst));
%         if length(foolst)>1,
%             DIST_PO2err(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
%         else
%             DIST_PO2err(ii) = 0;
%         end;
%     else
%         DIST_PO2(ii) = -1;
%         DIST_PO2err(ii) = 0;
%     end;
% end;
% figure; 
% errorbar(DIST,DIST_PO2,DIST_PO2err,'.-');
% xlabel('Distance from Pial Arteries (um)');
% ylabel('Average PO2 (mmHg)');        
% 
% %%
% % now calculate PO2 distribution by branching order from pial veins
%     
% 
% for ii = 1:length(im.PO2pts.PO2),
%     po2Seg = im.PO2pts.seg(ii);
%     po2_BRorder(ii) = segBR_Art(find(segBR_Art(:,1)==po2Seg),4);
%     po2_PO2(ii) = im.PO2pts.PO2(ii);
% end;
% BRorderIdx = sort(unique(po2_BRorder));
% for ii=1:length(BRorderIdx),
%     foolst = find(po2_BRorder == BRorderIdx(ii));
%     averVeinPO2_BRorder(ii) = mean(po2_PO2(foolst));
%     if length(foolst) > 1,
%         averVeinPO2err_BRorder(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
%     else
%         averVeinPO2err_BRorder(ii)=0;
%     end;
% end;
% figure; 
% errorbar(BRorderIdx,averVeinPO2_BRorder,averVeinPO2err_BRorder,'.-');
% xlabel('Branching Order from Pial Veins');
% ylabel('Average PO2 (mmHg)');
% 
% %%
% % now calculate PO2 distribution by distance from pial veins
%     
% 
% for ii = 1:length(im.PO2pts.PO2),
%     po2Seg = im.PO2pts.seg(ii);
%     po2_DIS(ii) = segBR_Art(find(segBR_Art(:,1)==po2Seg),5);
%     po2_PO2(ii) = im.PO2pts.PO2(ii);
% end;
% % bin in 100 um bins...
% dDIS = 100;
% maxDIS = max(po2_DIS)+1;
% minDIS = min(po2_DIS)-1;
% nDIS = floor((maxDIS-minDIS)/dDIS);
% for ii=1:nDIS,
%     foolst = find( (po2_DIS<ii*dDIS) & (po2_DIS>=(ii-1)*dDIS) );
%     DIST(ii) = (ii-0.5)*dDIS;
%     if ~isempty(foolst),
%         DIST_PO2(ii) = mean(po2_PO2(foolst));
%         if length(foolst)>1,
%             DIST_PO2err(ii) = std(po2_PO2(foolst))/sqrt(length(foolst)-1);
%         else
%             DIST_PO2err(ii) = 0;
%         end;
%     else
%         DIST_PO2(ii) = -1;
%         DIST_PO2err(ii) = 0;
%     end;
% end;
% figure; 
% errorbar(DIST,DIST_PO2,DIST_PO2err,'.-');
% xlabel('Distance from Pial Veins (um)');
% ylabel('Average PO2 (mmHg)');                
% 
% 
% 
