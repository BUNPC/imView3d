% TO DO
% Implement velocity BC. If velocity not specified then use Lipowsky's data
% on velocity versus diameter for Arterioles and venules. My polynomial fit
% gives 
%    vel (mm/s) = -5.28e-6 d^4 + 2.11e-4 d^3 + 0.0113 d^2 - 0.405 c + 4.70
% where  d=8 um for capillaries, d>8 is venules
% d<8 is arterioles but the true diameter is d_art = (8-d) + 8
%
% Mr and M are sparse. COuld I define them as sparse initially rather than 
% converting them from huge matrices to sparse?


function imView3d_flowCircuitEq( )
global im

nodeEdges = im.nodeEdges;
nodePos = im.nodePos;
nodeDiam = im.nodeDiam;
nodeBC = im.nodeBC;
nodeBCType = im.nodeBCType;
nodeType = im.nodeType;

nNodes = size(nodePos,1);
nEdges = size(nodeEdges,1);

nSegs = size(im.segEndNodes,1);
nSegNodes = max(im.segEndNodes(:));
segEndNodes = im.segEndNodes;
segLen = im.segLen;
segDiam = im.segDiam;
%segNodeMap = im.segNodeMap; 
segNodeMap = 1:nNodes;
segVesType = im.segVesType;


[nB,im] = nBupdate( im );


% Viscosity is ~ 2cP
% convert from cP to mmHg s
% 1 mmHg = 133.3 Pa
% 1 Pa s = 10 Poise = 1000 cP
% 1 Poise = 100 cP
% THEREFORE 1 cP = 0.01 Poise = 1e-3 Pa s = 7.5e-6 mmHg s



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATING RESISTANCE MATRIX (Mr matrix)
%

if ~im.flagUseSegments
    
    % USE EDGES
    
    Mr = zeros(nEdges,nNodes);
    for iE = 1:nEdges
        p1 = nodeEdges(iE,1);
        p2 = nodeEdges(iE,2);
        d1 = nodeDiam(p1);
        d2 = nodeDiam(p2);
        diam = (d1+d2)/2;
        len = norm(nodePos(p1,:)-nodePos(p2,:));
        R = 128 * 2*7.5e-6 * len / (3.14159 * diam.^4); % the 2 in the numerator is estimated viscosity
        eR(iE) = R;
        Mr(iE,p1) = 1/R;
        Mr(iE,p2) = -1/R;
    end

else
    
    % USE SEGMENTS
    
    Mr = zeros(nSegs,nSegNodes);
    for iS = 1:nSegs
        p1 = segEndNodes(iS,1);
        p2 = segEndNodes(iS,2);
        R = 128 * 2 * 7.5e-6 * segLen(iS) / (3.14159 * segDiam(iS)^4);
        eR(iS) = R;
        Mr(iS,p1) = 1/R;
        Mr(iS,p2) = -1/R;
    end
    
    nSegB = zeros(nSegNodes,1);
    for ii=1:nSegNodes
        nSegB(ii) = length(find(segEndNodes==ii));
    end

end
Mr = sparse(Mr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUNDARY CONDITION (M matrix, conservation of flow, etc)
% implement boundary conditions
% to construct M p = y
% where p is the pressure at each node point
% we implement flow conservation at all node points
% except at end points where we can specify pressure
% or specify velocity




hwait = waitbar( 0, 'Calculating flow : creating matrix' );
if ~im.flagUseSegments
    
    % USE EDGES
    
    M = zeros(nNodes,nNodes);
    y = zeros(nNodes,1);
    for ii=1:nNodes
        waitbar(ii/nNodes,hwait);
        [lstR, lstC] = find(nodeEdges==ii);

        if nB(ii)>1   % flow conservation
            for jj=1:length(lstR)
                p1 = nodeEdges(lstR(jj),lstC(jj));
                p2 = nodeEdges(lstR(jj),mod(lstC(jj),2)+1);
                M(ii,p1) = M(ii,p1) + 1/eR(lstR(jj));
                M(ii,p2) = M(ii,p2) - 1/eR(lstR(jj));
            end

        elseif nB(ii)==1  % use a BC
            p1 = nodeEdges(lstR,lstC);
            p2 = nodeEdges(lstR,mod(lstC,2)+1);

            if nodeBCType(ii)==1 || nodeBCType(ii)==3  % pressure BC

                if nodeBCType(ii)==1
                    M(ii,ii) = 1;
                    y(ii) = nodeBC(ii);

                elseif nodeBCType(ii)==3
                    M(ii,ii) = 1;
                    y(ii) = getPressure( nodeDiam(ii), nodeType(ii) );
                end

            elseif nodeBCType(ii)==2 || nodeBCType(ii)==4   % velocity BC

                if nodeBCType(ii)==2

                    M(ii,p1) = M(ii,p1) + 1/eR(lstR);
                    M(ii,p2) = M(ii,p2) - 1/eR(lstR);

                    vel = nodeBC(ii);
                    flow = vel * 3.14159 * (nodeDiam(ii)/2)^2;
                    if lstC==1
                        y(ii) = -flow;
                    else
                        y(ii) = flow;
                    end

                elseif nodeBCType(ii)==4

                    M(ii,p1) = M(ii,p1) + 1/eR(lstR);
                    M(ii,p2) = M(ii,p2) - 1/eR(lstR);

                    d1 = nodeDiam(p1);
                    d2 = nodeDiam(p2);
                    diam = (d1+d2)/2;
                    vel = -5.28e-6*diam^4 + 2.11e-4*diam^3 + 0.0113*diam^2 ...
                        - 0.405*diam + 4.70;
                    flow = vel * 3.14159 * (nodeDiam(ii)/2)^2;
                    if lstC==1  % This needs to consider if it is a vein or artery
                        y(ii) = -flow;
                    else
                        y(ii) = flow;
                    end
                    
                end

            else % no BC specified so assume vel = 0
                M(ii,p1) = M(ii,p1) + 1/eR(lstR);
                M(ii,p2) = M(ii,p2) - 1/eR(lstR);

                vel=0;

                if lstC==1
                    y(ii) = -vel;
                else
                    y(ii) = vel;
                end
            end

        end  % End of use a BC
    end % End of loop over nodes
    
else
    
    % USE SEGMENTS
    
    M = zeros(nSegNodes,nSegNodes);
    y = zeros(nSegNodes,1);
    for ii=1:nSegNodes
        waitbar(ii/nSegNodes,hwait);
        [lstR, lstC] = find(segEndNodes==ii);

        if nSegB(ii)>1   % flow conservation
            for jj=1:length(lstR)
                p1 = segEndNodes(lstR(jj),lstC(jj));
                p2 = segEndNodes(lstR(jj),mod(lstC(jj),2)+1);
                M(ii,p1) = M(ii,p1) + 1/eR(lstR(jj));
                M(ii,p2) = M(ii,p2) - 1/eR(lstR(jj));
            end

        elseif nSegB(ii)==1  % use a BC
            p1 = segEndNodes(lstR,lstC);
            p2 = segEndNodes(lstR,mod(lstC,2)+1);

            if nodeBCType(segNodeMap(ii))==1 || nodeBCType(segNodeMap(ii))==3  % pressure BC

                if nodeBCType(segNodeMap(ii))==1 %user defined
                    M(ii,ii) = 1;
                    y(ii) = nodeBC(segNodeMap(ii));

                elseif nodeBCType(segNodeMap(ii))==3 %from literature
                    M(ii,ii) = 1;
                    y(ii) = getPressure( segDiam(lstR), segVesType(lstR) ); %added 1/10 need to remove in future
                end

            elseif nodeBCType(segNodeMap(ii))==2 || nodeBCType(segNodeMap(ii))==4   % velocity BC

                if nodeBCType(segNodeMap(ii))==2 %user defined

                    M(ii,p1) = M(ii,p1) + 1/eR(lstR);
                    M(ii,p2) = M(ii,p2) - 1/eR(lstR);

                    vel = nodeBC(segNodeMap(ii));
                    flow = vel * 3.14159 * (segDiam(lstR)/2)^2;
                    
                   
                    if lstC==1
                       y(ii) = -1000*flow;
                    else
                       y(ii) = 1000*flow;
                    end

                    
                elseif nodeBCType(segNodeMap(ii))==4 %from litterature. How do we set the sign ?????

                    M(ii,p1) = M(ii,p1) + 1/eR(lstR);
                    M(ii,p2) = M(ii,p2) - 1/eR(lstR);

%                    d1 = segDiam(p1);
%                    d2 = segDiam(p2);
                    diam = segDiam(lstR);%(d1+d2)/2;
                    vel = -5.28e-6*diam^4 + 2.11e-4*diam^3 + 0.0113*diam^2 ...
                        - 0.405*diam + 4.70;
                    flow = vel * 3.14159 * (segDiam(lstR)/2)^2;

                end

            else % no BC specified so assume vel = 0
                M(ii,p1) = M(ii,p1) + 1/eR(lstR);
                M(ii,p2) = M(ii,p2) - 1/eR(lstR);

                vel=0;

                if lstC==1
                    y(ii) = -vel;
                else
                    y(ii) = vel;
                end
            end

        end  % End of use a BC
    end % End of loop over nodes
        
end
M = sparse(M);
close(hwait)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE PRESSURE from BC (M*p = BC)
%
% row normalize before inversion
hwait = waitbar(0,'Calculating flow : matrix inversion');
scl = max(M,[],2);
lstSegMap = find(scl~=0);
M = M(lstSegMap,lstSegMap);
scl = scl(lstSegMap);
y = y(lstSegMap);
M = spdiags(1./scl,0,length(scl),length(scl)) * M;
%M ./ (scl*ones(1,size(M,2)));
y = y ./ scl;
P(lstSegMap) = (M'*M)\M'*y;

close(hwait)

hwait = waitbar(0,'Calculating flow : map pressure');
if ~im.flagUseSegments
    im.nodePressure = P;
else
    im.nodePressure = zeros(nNodes,1);
    im.segPressure = P(lstSegMap);
    for iSeg = 1:nSegs
        waitbar(iSeg/nSegs,hwait);
        lstE = find(im.edgeSegN==iSeg);
        i1 = im.segEndNodes(iSeg,1);
        i2 = im.segEndNodes(iSeg,2);
        dP = P(i2) - P(i1);
        nSteps = length(lstE);
        dPstep = dP/nSteps;
        n1 = segNodeMap(i1);
        n2 = segNodeMap(i2);
        im.nodePressure(n1) = P(i1);
        nn = n1;
        Po = P(i1);
        nold = [];
        while nn~=n2
            lstE1 = find(im.nodeEdges(lstE,1)==nn | im.nodeEdges(lstE,2)==nn);
            lstN = im.nodeEdges(lstE(lstE1),:);
            nold(end+1) = nn;
            nn = setdiff( lstN(:), nold );
            Po = Po + dPstep;
            im.nodePressure(nn) = Po;
        end
    end
end
close(hwait)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE FLOW (F = M*p)
%
hwait = waitbar(0,'Calculating flow : calculate velocity step 1 of 3');
F = Mr(:,lstSegMap) * P(lstSegMap)';

%map segment flows to edge flows  (problem here, most edge Vel's are zero LG 1/25/2012)
if im.flagUseSegments
    Fseg = F;
    F = zeros(size(im.nodeEdges,1),1);
    for iSeg = 1:nSegs
        lstE = find(im.edgeSegN==iSeg);
        n1 = segNodeMap( im.segEndNodes(iSeg,1) );
        n2 = segNodeMap( im.segEndNodes(iSeg,2) );
        nn = n1;
        nold = [];
        eold = [];
        while nn~=n2     
            lstE1 = find(im.nodeEdges(lstE,1)==nn | im.nodeEdges(lstE,2)==nn);
            eIdx = setdiff( lstE(lstE1), eold );
            eold(end+1) = eIdx;
            if im.nodeEdges(eIdx,1)==nn
                F(eIdx) = Fseg(iSeg);
            else
                F(eIdx) = -Fseg(iSeg);
            end
            
            lstN = im.nodeEdges(lstE(lstE1),:);
            nold(end+1) = nn;
            nn = setdiff( lstN(:), nold );
        end
    end
%    F = Fseg(im.edgeSegN);   % Need to correct sign based on column order of edges
    im.segFlow = Fseg;
    im.segVel = Fseg ./ (3.14159*(im.segDiam(:)/2).^2);
end
Fedges = F;

%map edge flow to node flow (node flow = averaged abs(flow) of adjacent edges)
Fnode = zeros(nNodes,1);
for ii=1:nNodes
    waitbar(ii/nNodes,hwait);
    [lstR,lstC] = find(nodeEdges==ii);
    %    Fnode(ii) = mean(F(lstC).*((-1).^(lstC+1)));
    Fnode(ii) = mean(abs(F(lstR)));
    %    F(lstR)'
    %    nodeEdges(lstR,:)
    %    pause
end
close(hwait)

if ~im.flagUseSegments

    % USE EDGES

    im.nodeVel = zeros(nNodes,1);
    hwait = waitbar(0,'Calculating flow : calculate velocity step 2 of 3');
    for iN=1:nNodes
        waitbar(iN/nNodes,hwait);
        im.nodeVel(iN) = Fnode(iN) / (3.14159*nodeDiam(iN)^2/4);
    end
    close(hwait)
    im.edgeFlow = Fedges;
    im.edgeVel = zeros(nEdges,1);
    hwait = waitbar(0,'Calculating flow : calculate velocity step 3 of 3');
    for iE=1:size(nodeEdges,1)
        waitbar(iE/nEdges,hwait);
        rad = mean(nodeDiam(nodeEdges(iE,:)))/2;
        im.edgeVel(iE) = Fedges(iE) / (3.14159*rad^2);
    end
    close(hwait)

else

    % USE SEGMENTS
    %compute node velocities
    im.nodeVel = zeros(nNodes,1);
    hwait = waitbar(0,'Calculating flow : calculate velocity step 2 of 3');
    for iN=1:nNodes
        waitbar(iN/nNodes,hwait);
        if nB(iN)==2
            im.nodeVel(iN) = Fnode(iN) / (3.14159*segDiam(im.nodeSegN(iN))^2/4);
        else
            im.nodeVel(iN) = 0;
        end
    end
    close(hwait)
    
    %compute edge velocities
    im.edgeFlow = Fedges; 
    im.edgeVel = zeros(nEdges,1);
    hwait = waitbar(0,'Calculating flow : calculate velocity step 3 of 3');
    for iE=1:size(nodeEdges,1)
        waitbar(iE/nEdges,hwait);
%        rad = mean(segDiam(im.nodeSegN(nodeEdges(iE,:))))/2;
        rad = max(segDiam(im.nodeSegN(nodeEdges(iE,:))))/2;
        im.edgeVel(iE) = Fedges(iE) / (3.14159*rad^2);
    end
    close(hwait)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quiver plot of flow
%%
% xq = [];
% yq = [];
% zq = [];
% uq = [];
% vq = [];
% wq = [];
% for ii=1:size(nodeEdges,1)
%     xq(end+1) = mean(nodePos(nodeEdges(ii,:),1));
%     yq(end+1) = mean(nodePos(nodeEdges(ii,:),2));
%     zq(end+1) = mean(nodePos(nodeEdges(ii,:),3));
%     rq = nodePos(nodeEdges(ii,2),:) - nodePos(nodeEdges(ii,1),:);
%     rq = rq / norm(rq);
%     uq(end+1) = F(ii)*rq(1);
%     vq(end+1) = F(ii)*rq(2);
%     wq(end+1) = F(ii)*rq(3);
% end
% figure(10);
% quiver3(xq,yq,zq,uq,vq,wq,1,'linewidth',1.5,'maxheadsize',1)
% 


