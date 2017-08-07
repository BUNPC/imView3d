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
%
% 08/17/09
% RETHINKING FLOW OF UPDATE EQUATIONS TO HANDLE DYNAMICS
% 1) Given resistance everywhere calculate Flow using pressure BC and flow
% concervation (utilize dV/dt from previous time step). For baseline, the
% pressure BC is only at the end points. Afterwards, pressure is defined
% everywhere using the Windkessel Compliance relating pressure with volume.
% The capillaries and veins have a fixed compliance term and pressure and volume passively
% respond to dynamics induced by active arterioles. Arterioles actively
% modulated the compliance term. Specifically, dilation follows from: an increase
% in compliance drops pressure for a fixed volume and thus in-flow
% increases and out-flow decreases resulting in a positive dV/dt.
% 2) Given flow everywhere, calculate dV/dt for each segment
% 3) update volume for each segment, and pressure for each segment based on
% Windkessel Compliance
% BUT, AFTER TALKING WITH BUXTON I REALIZE IT IS BETTER TO 
% 0) Calculate Ao, Ro, Vo for each segment
% 1) Calculate Resistance for all segments based on segment diameter/volume
%    for use in (2) and (4). Is this most efficiently done as Ro (Vo/V)^2?
% 2) Calculate Pressure: using pressure BC at end points, flow conservation
%    at bifurcations, and flow conservation within segments considering
%    dV/dt. dV/dt in arterioles given by active dilation and constriction.
%    dV/dt in capillaries and venules given by first order temporal
%    response to pressure difference  dV/dt=(Ao P(t)^(1/beta) - V(t))/tau_comp
% 3) Update segment volumes based on dV/dt from (2) 
% 4) Calculate Flow: given pressure and resistance
% 5) go back to (1) for next time step
%
% Allow viscosity to vary in segments by varying Hct in segments
%

function im=flowCircuitEq_dyn( im, options )

if ~exist('options')
    options = [0 0 0];  % options(1) - map pressure from seg node to graph node
                        % options(2) - map edge flow to node flow
end

nodeEdges = im.nodeEdges;
nodePos = im.nodePos;
nodeDiam = im.nodeDiam;
nodeBC = im.nodeBC;
nodeBCType = im.nodeBCType;
nodeType = im.nodeType;

nNodes = size(nodePos,1);
nEdges = size(nodeEdges,1);

nSegs = size(im.segEndNodes,1);
segNodes = setdiff(unique(im.segEndNodes),0);
nSegNodes = length(segNodes);
%nSegNodes = max(im.segEndNodes(:));
[foo,segEndNodes] = ismember(im.segEndNodes,segNodes);
segLen = im.segLen';
segDiam = max(im.segDiam',5);
segVolo = 3.14159*segLen.*(segDiam/2).^2;
segVol = segVolo;
segBeta = 1*ones(nSegs,1);
segTau_v = .1 * ones(nSegs,1);
%segNodeMap = im.segNodeMap; 
%segNodeMap = 1:nNodes;
segNodeMap = segNodes;
segVesType = im.segVesType;


[nB,im] = nBupdate( im );


% Viscosity is ~ 2cP
% convert from cP to mmHg s
% 1 mmHg = 133.3 Pa
% 1 Pa s = 10 Poise = 1000 cP
% 1 Poise = 100 cP
% THEREFORE 1 cP = 0.01 Poise = 1e-3 Pa s = 7.5e-6 mmHg s


% Flow = Mr * Pressure
% Mr is presently defined as nSegs x nSegNodes (specifically seg end nodes)
% This only gives us flow in a segment. It won't give us in-flow and
% out-flow. For this I need pressure defined at segment half nodes.
% My Mr will be constructed as 2*nSegs x (nSegNodes+nSegs)
%    where we have F = [Fin; Fout] and P = [SegEndNodes; SegMidNodes]

%hwait = waitbar(0,'Calculating flow : creating resistance matrix');    
nSSN = nSegs+nSegNodes;
nSegs2 = 2*nSegs;
for iS = 1:nSegs
    p1 = segEndNodes(iS,1);
    p2 = segEndNodes(iS,2);
    lstMr(iS,1) = iS + (p1-1)*nSegs2;
    lstMr(iS,2) = iS + (nSegNodes+iS-1)*nSegs2;
    lstMr(iS,3) = nSegs+iS + (nSegNodes+iS-1)*nSegs2;
    lstMr(iS,4) = nSegs+iS + (p2-1)*nSegs2;
end
%close(hwait)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE PRESSURE
% implement boundary conditions
% to construct M p = y
% where p is the pressure at each node point
% we implement flow conservation at all node points
% except at end points where we can specify pressure
% or specify velocity
% RECALL
%    F = Mr P
%    Mr will be constructed as 2*nSegs x (nSegNodes+nSegs)
%    where we have F = [Fin; Fout] and P = [SegEndNodes; SegMidNodes]
% SO
%    M P = y
%    M is (nSegNodes+nSegs)x(nSegNodes+nSegs) where the bottom nSegs
%      simply comes from Mr where y is dV/dt
%    y is dV/dt for segments mid nodes when nSegB==2
%         0 for segment end nodes when nSegB>2
%         and BC when nSegB==1

%hwait = waitbar( 0, 'Calculating flow : creating pressure matrix' );
nSegB = zeros(nSegNodes,1);
for ii=1:nSegNodes
    nSegB(ii) = length(find(segEndNodes==ii));
end
maxNsegB = max(nSegB);

% USE SEGMENTS
lstM1a = [];
%lstM1aa = [];
lstM1b = [];
lstM1seg = [];
lstM2 = [];
lstM3 = [];
lsty2 = [];
lsty3 = [];
lstyRow2 = [];
lstyRow3 = [];

% LOOP OVER SEGMENT END NODES
for ii=1:nSegNodes
    [lstR, lstC] = find(segEndNodes==ii);
    
    if nSegB(ii)>1   % flow conservation
                     % I leave this as 1 rather than 2 in case some segment
                     % end nodes exist with no bifurcation. dV/dt is
                     % handled in the next loop below over nSegs
        for jj=1:length(lstR)
            p1 = segEndNodes(lstR(jj),lstC(jj)); % this is just ii
%            p2 = segEndNodes(lstR(jj),mod(lstC(jj),2)+1);
            lstM1a(end+1) = ii + (p1-1)*nSSN;
%            lstM1aa(end+1) = ii + (jj-1)*nSSN + (p1-1)*nSSN*maxNsegB;
            lstM1b(end+1) = ii + (nSegNodes+lstR(jj)-1)*nSSN;
            lstM1seg(end+1) = lstR(jj);
        end
        
    elseif nSegB(ii)==1  % use a BC
        p1 = segEndNodes(lstR,lstC);
%        p2 = segEndNodes(lstR,mod(lstC,2)+1);
        
        if nodeBCType(segNodeMap(ii))==1 || nodeBCType(segNodeMap(ii))==3  % pressure BC            
            if nodeBCType(segNodeMap(ii))==1
                lstM2(end+1) = ii + (ii-1)*nSSN;
                lstyRow2(end+1) = ii;
                lsty2(end+1) = segNodeMap(ii);
%                M(ii,ii) = 1;
%                y(ii) = nodeBC(segNodeMap(ii));                
            elseif nodeBCType(segNodeMap(ii))==3
                lstM3(end+1) = ii + (ii-1)*nSSN;
                lstyRow3(end+1) = ii;
                lsty3(end+1) = lstR;
%                M(ii,ii) = 1;
%                y(ii) = getPressure( segDiam(lstR), segVesType(lstR) );
            end

% HANDLE LATER            
%         elseif nodeBCType(segNodeMap(ii))==2 || nodeBCType(segNodeMap(ii))==4   % velocity BC            
%             if nodeBCType(segNodeMap(ii))==2                
%                 M(ii,p1) = M(ii,p1) + 1/eR(lstR);
%                 M(ii,p2) = M(ii,p2) - 1/eR(lstR);                
%                 vel = nodeBC(segNodeMap(ii));
%                 flow = vel * 3.14159 * (segDiam(lstR)/2)^2;
%                 if lstC==1
%                     y(ii) = -flow;
%                 else
%                     y(ii) = flow;
%                 end                
%             elseif nodeBCType(segNodeMap(ii))==4                
%                 M(ii,p1) = M(ii,p1) + 1/eR(lstR);
%                 M(ii,p2) = M(ii,p2) - 1/eR(lstR);
%                 diam = segDiam(lstR);
%                 vel = -5.28e-6*diam^4 + 2.11e-4*diam^3 + 0.0113*diam^2 ...
%                     - 0.405*diam + 4.70;
%                 flow = vel * 3.14159 * (segDiam(lstR)/2)^2;
%                 
%             end
            
        else % no BC specified so assume vel = 0
            
            error( 'we should not get here if only pressure BC employed' )
%             M(ii,p1) = M(ii,p1) + 1/eR(lstR);
%             M(ii,p2) = M(ii,p2) - 1/eR(lstR);            
%             vel=0;            
%             if lstC==1
%                 y(ii) = -vel;
%             else
%                 y(ii) = vel;
%             end
        end
        
    end  % End of use a BC
end % End of loop over nodes
%close(hwait)


for ii=1:length(lstyRow3)   % Maybe I want this fixed during iterations
    yRow3(ii) = getPressure(segDiam(lsty3(ii)),floor(segVesType(lsty3(ii))));
end


% get list of end segs
lst = find(nSegB==1);
lstEndSegA = [];
lstEndSegC = [];
lstEndSegV = [];
for ii=1:length(lst)
    [lstR,lstC] = find(segEndNodes==lst(ii));
    if floor(segVesType(lstR))==1
        lstEndSegA(end+1) = lstR;
    elseif segVesType(lstR)==2
        lstEndSegC(end+1) = lstR;
    elseif segVesType(lst)==3
        lstEndSegV(end+1) = lstR;
    end
end

% find end seg Art with largest diam
%[foo,iAmax]=max(segDiam(lstEndSegA))
%iAmax = lstEndSegA(iAmax);
lstDilateArt = find(segVesType==1.1);
%lstDilateArt = lstDilateArt(3); %3,4

lstSegCV = find(segVesType==2 | segVesType==3);
lstSegA = find(floor(segVesType)==1);

% SOLVE

% initialize for first step
% R = 8 viscosity * len / (pi radius^4)
% divide segLen by 2 to split resistance between in-flow and out-flow
Rsego = (128 * 2*7.5e-6 * (segLen/2) ./ (3.14159 * segDiam.^4));
Rseg = Rsego;
segDiamo = segDiam;

Mr = sparse(nSegs*2, nSegNodes+nSegs);
M = sparse(nSSN,nSSN);
y = zeros(nSSN,1);
Po = zeros(nSSN,1);

seg_dVdt = zeros(nSegs,1);

% iterate
dt = 10e-3;
nT = 100;
Ptime = zeros(nSSN,nT);
VelTime = zeros(nSegs,nT);
FTime = zeros(nSegs*2,nT);
segVolTime = zeros(nSegs,nT);

for iT = 1:nT
    disp( sprintf('Iter %d of %d  (%.3f sec)',iT,nT,iT*dt) )
    drawnow
    Mr(lstMr(:)) = [1./Rseg; -1./Rseg; 1./Rseg; -1./Rseg];
    M = sparse(nSSN,nSSN);
    for ii=1:length(lstM1a) % This can be vectorized
        M(lstM1a(ii)) = M(lstM1a(ii)) + 1/Rseg(lstM1seg(ii));
        M(lstM1b(ii)) = M(lstM1b(ii)) - 1/Rseg(lstM1seg(ii));
    end
    M(lstM2) = 1;
    M(lstM3) = 1;
    M(nSegNodes+[1:nSegs],:) = Mr(1:nSegs,:)-Mr(nSegs+[1:nSegs],:);
    
    y(lstyRow2) = nodeBC(lsty2);
    y(lstyRow3) = yRow3;
    y(nSegNodes+[1:nSegs]) = seg_dVdt;
    
    scl=sum(M.^2,2).^0.5;
    M = spdiags(1./scl,0,length(scl),length(scl)) * M;
    y = y ./ scl;
%    tic 
%    if iT<2
        P = (M'*M)\M'*y;
%    else
%        [P,flag,relres,iter,resvec] = qmr(M'*M,M'*y,1e-12,1000,[],[],Po);
%        [P,flag,relres,iter,resvec] = bicg(M'*M,M'*y,1e-12,1000,[],[],Po);
%        [P,flag,relres,iter,resvec] = bicgstab(M'*M,M'*y,1e-12,1000,[],[],Po);
%        [P,flag,relres,iter,resvec] = gmres(M'*M,M'*y,[],1e-12,1000,[],[],Po);
%    end
%    toc

    if iT==1
        segAo = (segVolo.^segBeta) ./ P(nSegNodes+[1:nSegs]);
    end
    
    % update arteriole dVdt separately from capillary and venule
    % capillaries and venules
    foo = ( (segAo.*P(nSegNodes+[1:nSegs])).^(1./segBeta) - ...
                 segVol )./segTau_v;
    seg_dVdt(lstSegCV) = foo(lstSegCV);
    % arterioles   ...  It seems that this allows some arterioles to be
    % compliantly defined by the above equation
    seg_dVdt(lstSegA) = 0;
    if iT<20
        seg_dVdt(lstDilateArt) = 1e-1 * segVol(lstDilateArt);
    else
        seg_dVdt(lstDilateArt) = 0;
    end
    
    segVol = segVol + seg_dVdt*dt;
    Rseg = Rsego .* (segVolo./segVol).^2;
    segDiam = segDiamo .* (segVol./segVolo).^(0.5);
    Po = P;
    Ptime(:,iT) = P;
    segVolTime(:,iT) = segVol';
    FTime(:,iT) = Mr * P;
    
    if mod(iT,1)==0
        if 1
            % Pressure
            % need to use code from imView3d_flowCircuitEq if I want smooth map
            % of P to nodes
%            Pnode(find(im.nodeSegN~=0))=P(nSegNodes+im.nodeSegN(find(im.nodeSegN~=0)));
            Pnode(find(im.nodeSegN~=0))=Ptime(nSegNodes+im.nodeSegN(find(im.nodeSegN~=0)),iT) - ...
                Ptime(nSegNodes+im.nodeSegN(find(im.nodeSegN~=0)),1);
            Pnode(find(im.nodeSegN~=0))=Ptime(nSegNodes+im.nodeSegN(find(im.nodeSegN~=0)),iT);
        elseif 0
            % Velocity or Flow
            Fseg2 = FTime(:,iT);
            VelSeg = Fseg2([1:nSegs]) ./ (3.14159*(segDiam/2).^2); % Need to use updated Diameter
%            Fnode(find(im.nodeSegN~=0)) = Fseg(im.nodeSegN(find(im.nodeSegN~=0)))
            VelNode2(find(im.nodeSegN~=0)) = VelSeg(im.nodeSegN(find(im.nodeSegN~=0)));
            Fnode2(find(im.nodeSegN~=0)) = abs(Fseg2(im.nodeSegN(find(im.nodeSegN~=0))));
%            Pnode = VelNode;
%            max(VelSeg)
            VelTime(:,iT) = VelSeg;

            Fseg1 = FTime(:,1);
            VelSeg = Fseg1([1:nSegs]) ./ (3.14159*(segDiam/2).^2);
            VelNode1(find(im.nodeSegN~=0)) = VelSeg(im.nodeSegN(find(im.nodeSegN~=0)));
            Fnode1(find(im.nodeSegN~=0)) = abs(Fseg1(im.nodeSegN(find(im.nodeSegN~=0))));
            
%            Pnode = min(abs(VelNode1)/1e3,10);
            Pnode = min((abs(VelNode2) - abs(VelNode1))/1e3,5);
%            Pnode = min((Fnode2 - Fnode1)/1e6,1);
%            Pnode = min((abs(Fnode2)/1e6),15);
        else
            %Seg Vol
            Pnode(find(im.nodeSegN~=0))=(segVolTime(im.nodeSegN(find(im.nodeSegN~=0)),iT) - ...
                segVolTime(im.nodeSegN(find(im.nodeSegN~=0)),1)) ./ ...
                segVolTime(im.nodeSegN(find(im.nodeSegN~=0)),1);
            Pnode = min(Pnode,1e-3);
        end
        
        foo = zeros(size(im.Mesh.node,1),1);
        lst = find(im.VesFlux.gfMap>0);
        boo = Pnode;
        foo(im.Mesh.vesWallnode) = boo(im.VesFlux.gfMap(lst));
        figure(12)
        trisurf( im.Mesh.boundary(find(im.Mesh.boundary(:,end)==0),1:3), im.Mesh.node(:,1), im.Mesh.node(:,2), im.Mesh.node(:,3), foo, 'linestyle','none','facecolor','flat');
        colorbar
        axis image
        colormap(jet(64))
        
        theta = 120;
        campos([4000*cos(theta*pi/180) 400*sin(theta*pi/180) -800]);
        camup([0 0 -1]);
        camva(12);
        if iT>1
            pause(0.1)
        end
    end

end



return















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Velocity
Fseg = Mr * P;

% map flow in segments to flow in edges
% make sure we get the sign right
Fsign = sign(im.nodePressure(im.nodeEdges(:,1)) - im.nodePressure(im.nodeEdges(:,2)))';
F = Fsign .* abs(Fseg(im.edgeSegN));
Fnode(find(im.nodeSegN~=0)) = Fseg(im.nodeSegN(find(im.nodeSegN~=0)))
nodeVel = Fnode ./ (3.14159*(im.segDiam/2).^2);



im.segFlow = Fseg;
im.segVel = Fseg ./ (3.14159*(im.segDiam(:)/2).^2);
im.edgeFlow = F;
im.edgeVel = im.segVel(im.edgeSegN);



% This is wrong, check imView3d_flowCircuitEq
if im.flagUseSegments
    Fseg = F;
    F = Fseg(im.edgeSegN);
    
    im.segFlow = Fseg;
    im.segVel = Fseg ./ (3.14159*(im.segDiam(:)/2).^2);
    im.edgeFlow = F;
    im.edgeVel = im.segVel(im.edgeSegN);
end
Fedges = F;


return




Fnode = zeros(nNodes,1);
if options(2)
    hwait = waitbar(0,'Calculating flow : mapping edge flow to node flow');
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
end


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


