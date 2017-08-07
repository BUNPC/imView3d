function imView3d_vesselMask( handles )
global im

if ~isfield(im,'maskMaxRad')
    im.maskMaxRad = 3;
end
if ~isfield(im,'maskCircleFlag')
    im.maskCircleFlag = 1;
end
if ~isfield(im,'maskUseSegDiam')
    im.maskUseSegDiam = 0;
end
if ~isfield(im,'maskSphereSpacing')
    im.maskSphereSpacing = 0.7;
end

maxR = im.maskMaxRad;  % max radius for mask. This speeds up the generation of the mask
% useful if you want to just view the vessel type label, etc

if im.maskUseSegDiam
    im = nodeGrps( im );
end



nx = im.nX;
ny = im.nY;
nz = im.nZ;

if ~isfield(im,'Hvox')
    im.Hvox = 1;
    hwait = waitbar(0,'Please enter the voxel size!');
    while length(im.Hvox~=3)
        im.Hvox = str2num(input('Enter voxel size [x y z]?','s'));
        if length(im.Hvox~=3)
            disp('Please enter 3 lengths!');
        end
    end
    close(hwait);
end

if isfield(im,'edgeVel')
    edgeVel = abs(im.edgeVel);
    maxEdgeVel = max(edgeVel);
    edgeFlow = abs(im.edgeFlow);
    maxEdgeFlow = max(edgeFlow);
    edgeFlowLog = log10(abs(im.edgeFlow));
    maxEdgeFlowLog = max(edgeFlowLog);
    nodePressure = im.nodePressure;
    maxPressure = max(im.nodePressure(:));
    
else
    edgeVel = zeros(size(im.nodeEdges,1),1);
    maxEdgeVel = 1;
    edgeFlow = zeros(size(im.nodeEdges,1),1);
    maxEdgeFlow = 1;
    edgeFlowLog = zeros(size(im.nodeEdges,1),1);
    maxEdgeFlowLog = 1;
    nodePressure = zeros(size(im.nodePos,1),1);
    maxPressure = 1;
end

flagPO2 = 0;
if isfield(im,'PO2pts')
    if isfield(im.PO2pts,'nodePO2')
        flagPO2 = 1;
        nodePO2 = im.PO2pts.nodePO2;
        maxnodePO2 = max(nodePO2);
        minnodePO2 = min(nodePO2(find(nodePO2>0)));
        minnodePO2 = minnodePO2 - (maxnodePO2 - minnodePO2)/32;
    end
end
if ~flagPO2
    nodePO2 = zeros(size(im.nodePos,1),1);
    maxnodePO2 = 1;
    minnodePO2 = 0;
end

if isfield(im,'edgeBRorderArt')
    edgeBRorderArt = im.edgeBRorderArt;
    maxBRorderArt = max(edgeBRorderArt);
    edgeBRorderVeins = im.edgeBRorderVeins;
    maxBRorderVeins = max(edgeBRorderVeins);
else
    edgeBRorderArt = zeros(size(im.nodeEdges,1),1);
    maxBRorderArt = 1;
    edgeBRorderVeins = zeros(size(im.nodeEdges,1),1);
    maxBRorderVeins = 1;
end

updateFlag = 0;
if ~isfield(im,'nodeTypeUpdated')
    im.nodeTypeUpdated = ones(length(im.nodeType),1);
elseif ~isempty(find(im.nodeTypeUpdated==1))
    ch = menu('Updated mask for all nodes or just modified node?','All','Modified');
    if ch==1
        im.nodeTypeUpdated = ones(length(im.nodeType),1);
    else
        updateFlag = 1;
    end
else
    im.nodeTypeUpdated = ones(length(im.nodeType),1);
end


nEdges = size(im.nodeEdges,1);
nNodes = size(im.nodePos,1);

hwait = waitbar(0,'Allocating memory for mask');
Vm = uint8(zeros(ny,nx,nz));
VmDiam = single(zeros(nx,ny,nz));
waitbar(.25,hwait);
if isfield(im,'edgeVel')
    Vvel = uint8(zeros(ny,nx,nz));
    waitbar(.5,hwait);
    Vflow = uint8(zeros(ny,nx,nz));
    VflowLog = uint8(zeros(ny,nx,nz));
    waitbar(.75,hwait);
    Vpres = uint8(zeros(ny,nx,nz));
end
if flagPO2
    Vpo2 = uint8(zeros(ny,nx,nz));
end
if isfield(im,'edgeBRorderArt')
    VbrA = uint8(zeros(ny,nx,nz));
    waitbar(.5,hwait);
    VbrV = uint8(zeros(ny,nx,nz));
    waitbar(.75,hwait);
end

close(hwait)

if updateFlag==1
    Vm = im.Vm;
    Vupdate = logical(zeros(ny,nx,nz));
    if isfield(im,'edgeVel')
        Vvel = im.Vvel;
        Vflow = im.Vflow;
        VflowLog = im.VflowLog;
        Vpres = im.Vpres;
    end
    if flagPO2 & isfield(im,'Vpo2')
        Vpo2 = im.Vpo2;
    end
    if isfield(im,'edgeBRorderArt')
        VbrA = im.VbrA;
        VbrV = im.VbrV;
    end;
end


% Raj 2015   make this poition about 6x faster. #########################
[r,w]=unix('hostname');
empt = isempty(gcp('nocreate')); %% edit Raj 2016 matlabpool not a function anymore
if strcmp (w(1:5), 'ralph') && ~isequal(empt,1)
    
    hWait = waitbar( 0, sprintf('Masking vessels in paralell...'));
    tic
    nodeMasked = zeros(nNodes,1);
    if updateFlag
        lstUpdate = find(im.nodeTypeUpdated==1);
        foo1 = ismember(im.nodeEdges(:,1),lstUpdate);
        foo2 = ismember(im.nodeEdges(:,2),lstUpdate);
        lstEdges = find(foo1==1 | foo2==1);
    else
        lstEdges = 1:nEdges;
    end
    %break up the work into 12 equal sized chunks.
    %this should only be done on a high mem system
    
    
    
    % [Mm Nn Oo]=size(Vm);
    p_Vm=cell(10,1);
    p_VmDiam=cell(10,1);
    p_Vupdate=cell(10,1);
    p_Vvel=cell(10,1);
    p_Vflow=cell(10,1);
    p_VflowL=cell(10,1);
    p_Vpres=cell(10,1);
    p_Vpo2=cell(10,1);
    p_VbrA=cell(10,1);
    p_VbrV=cell(10,1);
    
    for ii=1:10
        
        p_Vm{ii}=Vm;
        p_VmDiam{ii}=VmDiam;
        if  exist('Vupdate', 'var')
            p_Vupdate{ii}=Vupdate;
        end
        if  exist('Vvel', 'var')
            p_Vvel{ii}=Vvel;
        end
        if exist('Vflow', 'var')
            p_Vflow{ii}=Vflow;
        end
        if exist('VflowL', 'var')
            p_VflowL{ii}=VflowL;
        end
        if exist('Vpres', 'var')
            p_Vpres{ii}=Vpres;
        end
        if exist('Vpo2', 'var')
            p_Vpo2{ii}=Vpo2;
        end
        if exist('VbrA', 'var')
            p_VbrA{ii}= VbrA;
        end
        if exist('VbrV', 'var')
            p_VbrV{ii}= VbrV;
        end
        
    end %for
    
    
    p_lstEdges=cell(10,1);
    p_nodeMasked=cell(10,1);
    Pp= round (length(lstEdges)/10);
    for ii=1:9
        p_lstEdges{ii}=lstEdges((ii*Pp-Pp+1): ii*Pp);
        p_nodeMasked{ii}=zeros(Pp);
    end
    ii=10;
    p_lstEdges{ii}=lstEdges((ii*Pp-Pp+1): end);
    p_nodeMasked{ii}=zeros(Pp);
    
    parfor p_index=1:10
        
        for iiE=1:length(p_lstEdges{p_index})
            iE = p_lstEdges{p_index}(iiE);
            %     if isequal(rem(iiE,1000),0)
            %     waitbar(iiE/length(lstEdges),hWait); %waitbars take a long time.
            %     end
            % for debugging:
            %     global i1   %making things global may take a long time.
            %     global i2
            
            i1 = im.nodeEdges(iE,1);
            i2 = im.nodeEdges(iE,2);
            %nodeMasked becomes p_nodeMasked{p_index}
            if (im.nodeTypeUpdated(i1) | im.nodeTypeUpdated(i2)) & (p_nodeMasked{p_index}(i1)==0 |  p_nodeMasked{p_index}(i2)==0)
                if ~im.maskUseSegDiam
                    r1 = im.nodeDiam(i1)/2;
                    r1 = min(max(r1,1),maxR);
                    r2 = im.nodeDiam(i2)/2;
                    r2 = min(max(r2,1),maxR);
                else
                    if ~im.nodeSegN(i1)==0 || ~im.nodeSegN(i2)==0
                        foo = max(im.segDiam(im.nodeSegN(i1)),im.segDiam(im.nodeSegN(i2)) )/2;
                        foo = min(max(foo,1),maxR);
                        r1 = foo;
                        r2 = foo;
                    end
                end
                p1 = round(im.nodePos(i1,:));
                p2 = round(im.nodePos(i2,:));
                d12 = norm(p2-p1);
                dxyz = (p2-p1);%/d12;
                rd = (r2-r1);%/d12;
                nSteps = 1;
                stepLen = max(r1*im.maskSphereSpacing,1);
                if stepLen<d12
                    dxyz = (p2-p1)*stepLen/d12;
                    rd = (r2-r1)*stepLen/d12;
                    nSteps = floor(d12/stepLen)+1;
                end
                
                if im.nodeType(i1)>0
                    nType = im.nodeType(i1)+1;
                else
                    nType = 1;
                end
                
                % this should consider im.Hvox
                %        warning( 'this should consider im.Hvox' )
                lst = find(sum((im.nodePos-ones(nNodes,1)*p1).^2,2).^0.5<(r1*im.maskSphereSpacing) );
                p_nodeMasked{p_index}(lst) = 1;
                lst = find(sum((im.nodePos-ones(nNodes,1)*p2).^2,2).^0.5<(r2*im.maskSphereSpacing) );
                p_nodeMasked{p_index}(lst) = 1;
                
                p = p1;
                r = r1;
                flag = 1;
                while flag
                    if norm(round(p)-p2)==0
                        flag = 0;
                    end
                    pr = round(p);
                    rTmp = min(r,maxR);
                    rx = ceil(rTmp/im.Hvox(1));
                    ry = ceil(rTmp/im.Hvox(2));
                    rz = ceil(rTmp/im.Hvox(3));
                    if im.maskCircleFlag
                        rz = 0;
                    end
                    for iX = -rx:+rx
                        for iY = -ry:+ry
                            for iZ = -rz:+rz
                                if norm([iX*im.Hvox(1) iY*im.Hvox(2) iZ*im.Hvox(3)])<=r
                                    iix = min(max(pr(1)+iX,1),nx);
                                    iiy = min(max(pr(2)+iY,1),ny);
                                    iiz = min(max(pr(3)+iZ,1),nz);
                                    p_Vm{p_index}(iiy,iix,iiz) = 1 + nType;
                                    if nType == 4, % if vein
                                        p_VmDiam{p_index}(iiy,iix,iiz) = -single(2*r);
                                    else
                                        p_VmDiam{p_index}(iiy,iix,iiz) = single(2*r);
                                    end;
                                    if updateFlag
                                        p_Vupdate{p_index}(iiy,iix,iiz) = 1;
                                    end
                                    if isfield(im,'edgeVel')
                                        p_Vvel{p_index}(iiy,iix,iiz) = max(round(32*min(edgeVel(iE)/maxEdgeVel,1)),1);%add that to avoid black spot
                                        p_Vflow{p_index}(iiy,iix,iiz) = max(round(32*min(edgeFlow(iE)/maxEdgeFlow,1)),1);%add that to avoid black spot
                                        p_VflowLog{p_index}(iiy,iix,iiz) = max(round(32*min(edgeFlowLog(iE)/maxEdgeFlowLog,1)),1);%add that to avoid black spot
                                        p_Vpres{p_index}(iiy,iix,iiz) = round(32*min(nodePressure(i1)/maxPressure,1));
                                    end
                                    if flagPO2
                                        p_Vpo2{p_index}(iiy,iix,iiz) = round(32*min((nodePO2(i1)-minnodePO2)/(maxnodePO2-minnodePO2),1));
                                    end
                                    if isfield(im,'edgeBRorderArt')
                                        p_VbrA{p_index}(iiy,iix,iiz) = round(32*min(edgeBRorderArt(iE)/maxBRorderArt,1));
                                        p_VbrV{p_index}(iiy,iix,iiz) = round(32*min(edgeBRorderVeins(iE)/maxBRorderVeins,1));
                                    end
                                end
                            end
                        end
                    end
                    
                    %         for iX = max(pr(1)-rx,1):min(pr(1)+rx,nx)
                    %             for iY = max(pr(2)-ry,1):min(pr(2)+ry,ny)
                    %                 for iZ = max(pr(3)-rz,1):min(pr(3)+rz,nz)
                    %                     if norm([iX-pr(1) iY-pr(2) iZ-pr(3)])<=r
                    %                         Vm(iY,iX,iZ) = 1 + nType;
                    %                     end
                    %                 end
                    %             end
                    %         end
                    if flag
                        p = p + dxyz;
                        r = r + rd;
                    end
                    nSteps = nSteps - 1;
                    if nSteps==0
                        flag = 0;
                    end
                    if norm(round(p)-p2)==0
                        flag = 0;
                    end
                    
                end
            end % end of check if node type updated
        end
    end
    [Mm, Nn, Oo]=size(Vm);
    
    
    
    %bring everything back togetehr after paralell
    for ii=1:10
        if ii==1;
            p_Vm{1}=reshape(p_Vm{1},1,Mm*Nn*Oo);
        end
        p_Vm{ii}=reshape(p_Vm{ii},1,Mm*Nn*Oo);
        [xx, yy, vv]=find(p_Vm{ii}); %get all the non zero points in p_Vm
        p_Vm{1}(yy)=p_Vm{ii}(yy);
        
        
        if ii==1;
            p_VmDiam{1}=reshape(p_VmDiam{1},1,Mm*Nn*Oo);
        end
        p_VmDiam{ii}=reshape(p_VmDiam{ii},1,Mm*Nn*Oo);
        [xx, yy, vv]=find(p_VmDiam{ii}); %get all the non zero points in p_Vm
        p_VmDiam{1}(yy)=p_VmDiam{ii}(yy);
        
        
        if  exist('Vupdate', 'var')
            if ii==1;
                p_Vupdate{1}=reshape(p_Vupdate{1},1,Mm*Nn*Oo);
            end
            p_Vupdate{ii}=reshape(p_Vupdate{ii},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_Vupdate{ii}); %get all the non zero points in p_Vm
            p_Vupdate{1}(yy)=p_Vupdate{ii}(yy);
        end
        
        if  exist('Vvel', 'var')
            if ii==1;
                p_Vvel{1}=reshape(p_Vvel{1},1,Mm*Nn*Oo);
            end
            p_Vvel{ii}=reshape(p_Vvel{ii},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_Vvel{ii}); %get all the non zero points in p_Vm
            p_Vvel{1}(yy)=p_Vvel{ii}(yy);
        end
        
        if exist('Vflow', 'var')
            if ii==1;
                p_Vflow{1}=reshape(p_Vflow{1},1,Mm*Nn*Oo);
            end
            p_Vflow{ii}=reshape(p_Vflow{ii},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_Vflow{ii}); %get all the non zero points in p_Vm
            p_Vflow{1}(yy)=p_Vflow{ii}(yy);
            
        end
        if exist('VflowL', 'var')
            if ii==1;
                p_VflowL{1}=reshape(p_VflowL{1},1,Mm*Nn*Oo);
            end
            p_VflowL{ii}=reshape(p_VflowL{ii},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_VflowL{ii}); %get all the non zero points in p_Vm
            p_VflowL{1}(yy)=p_VflowL{ii}(yy);
        end
        
        if exist('Vpres', 'var')
            if ii==1;
                p_Vpres{1}=reshape(p_Vpres{1},1,Mm*Nn*Oo);
            end
            p_Vpres{ii}=reshape(p_Vpres{ii},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_Vpres{ii}); %get all the non zero points in p_Vm
            p_Vpres{1}(yy)=p_Vpres{ii}(yy);
        end
        
        if exist('Vpo2', 'var')
            if ii==1;
                p_Vpo2{1}=reshape(p_Vpo2{1},1,Mm*Nn*Oo);
            end
            p_Vpo2{ii}=reshape(p_Vpo2{ii},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_Vpo2{ii}); %get all the non zero points in p_Vm
            p_Vpo2{1}(yy)=p_Vpo2{ii}(yy);
        end
        
        if exist('VbrA', 'var')
            if ii==1;
                p_VbrA{1}=reshape(p_VbrA{1},1,Mm*Nn*Oo);
            end
            p_VbrA{ii}=reshape(p_VbrA{ii},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_VbrA{ii}); %get all the non zero points in p_Vm
            p_VbrA{1}(yy)=p_VbrA{ii}(yy);
        end
        
        if exist('VbrV', 'var')
            if ii==1;
                p_VbrV{1}=reshape(p_VbrV{1},1,Mm*Nn*Oo);
            end
            p_VbrV{ii}=reshape(p_VbrV{ii},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_VbrV{ii}); %get all the non zero points in p_Vm
            p_VbrV{1}(yy)=p_VbrV{ii}(yy);
        end
        
        
    end
    if ~updateFlag %then Vm started as zeros
        %     Vm=p_Vm{1};
        Vm=reshape(p_Vm{1},Mm,Nn,Oo);
        
        VmDiam=p_VmDiam{1};
        
        if  exist('Vupdate', 'var')
            Vupdate=p_Vupdate{1};
        end
        if  exist('Vvel', 'var')
            Vvel=p_Vvel{1};
        end
        if exist('Vflow', 'var')
            Vflow=p_Vflow{1};
        end
        if exist('VflowL', 'var')
            VflowL=p_VflowL{1};
        end
        if exist('Vpres', 'var')
            Vpres=p_Vpres{1};
        end
        if exist('Vpo2', 'var')
            Vpo2=p_Vpo2{1};
        end
        if exist('VbrA', 'var')
            VbrA=p_VbrA{1};
        end
        if exist('VbrV', 'var')
            VbrV=p_VbrV{1};
        end
        
    else
        %only update what has changed
        [Mm, Nn ,Oo]=size(Vm);
        
        Vm=reshape(Vm,1,Mm*Nn*Oo);
        p_Vm{1}=reshape(p_Vm{1},1,Mm*Nn*Oo);
        [xx, yy, vv]=find(p_Vm{1}); %get all the non zero points in p_Vm
        Vm(yy)=p_Vm{1}(yy);
        Vm=reshape(Vm,Mm,Nn,Oo);
        
        VmDiam=reshape(VmDiam,1,Mm*Nn*Oo);
        p_VmDiam{1}=reshape(p_VmDiam{1},1,Mm*Nn*Oo);
        [xx ,yy ,vv]=find(p_VmDiam{1}); %get all the non zero points in p_VmDiam
        VmDiam(yy)=p_VmDiam{1}(yy);
        VmDiam=reshape(VmDiam,Mm,Nn,Oo);
        
        
        if  exist('Vupdate', 'var')
            Vupdate=reshape(Vupdate,1,Mm*Nn*Oo);
            p_Vupdate{1}=reshape(p_Vupdate{1},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_Vupdate{1}); %get all the non zero points in p_Vupdate
            Vupdate(yy)=p_Vupdate{1}(yy);
            Vupdate=reshape(Vupdate,Mm,Nn,Oo);
        end
        if  exist('Vvel', 'var')
            Vvel=reshape(Vvel,1,Mm*Nn*Oo);
            p_Vvel{1}=reshape(p_Vvel{1},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_Vvel{1}); %get all the non zero points in p_Vvel
            Vvel(yy)=p_Vvel{1}(yy);
            Vvel=reshape(Vvel,Mm,Nn,Oo);
        end
        if exist('Vflow', 'var')
            Vflow=reshape(Vflow,1,Mm*Nn*Oo);
            p_Vflow{1}=reshape(p_Vflow{1},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_Vflow{1}); %get all the non zero points in p_Vflow
            Vflow(yy)=p_Vflow{1}(yy);
            Vflow=reshape(Vflow,Mm,Nn,Oo);
        end
        if exist('VflowL', 'var')
            VflowL=reshape(VflowL,1,Mm*Nn*Oo);
            p_VflowL{1}=reshape(p_VflowL{1},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_VflowL{1}); %get all the non zero points in p_VflowL
            VflowL(yy)=p_VflowL{1}(yy);
            VflowL=reshape(VflowL,Mm,Nn,Oo);
        end
        if exist('Vpres', 'var')
            Vpres=reshape(Vpres,1,Mm*Nn*Oo);
            p_Vpres{1}=reshape(p_Vpres{1},1,Mm*Nn*Oo);
            [xx, yy ,vv]=find(p_Vpres{1}); %get all the non zero points in p_Vpres
            Vpres(yy)=p_Vpres{1}(yy);
            Vpres=reshape(Vpres,Mm,Nn,Oo);
        end
        if exist('Vpo2', 'var')
            Vpo2=reshape(Vpo2,1,Mm*Nn*Oo);
            p_Vpo2{1}=reshape(p_Vpo2{1},1,Mm*Nn*Oo);
            [xx ,yy, vv]=find(p_Vpo2{1}); %get all the non zero points in p_Vpo2
            Vpo2(yy)=p_Vpo2{1}(yy);
            Vpo2=reshape(Vpo2,Mm,Nn,Oo);
        end
        if exist('VbrA', 'var')
            VbrA=reshape(VbrA,1,Mm*Nn*Oo);
            p_VbrA{1}=reshape(p_VbrA{1},1,Mm*Nn*Oo);
            [xx, yy, vv]=find(p_VbrA{1}); %get all the non zero points in p_VbrA
            VbrA(yy)=p_VbrA{1}(yy);
            VbrA=reshape(VbrA,Mm,Nn,Oo);
        end
        if exist('VbrV', 'var')
            VbrV=reshape(VbrV,1,Mm*Nn*Oo);
            p_VbrV{1}=reshape(p_VbrV{1},1,Mm*Nn*Oo);
            [xx, yy ,vv]=find(p_VbrV{1}); %get all the non zero points in p_VbrV
            VbrV(yy)=p_VbrV{1}(yy);
            VbrV=reshape(VbrV,Mm,Nn,Oo);
        end
        
    end
    
    
    toc
else %if statement about paralell
    % ####### END Raj 2015 make this poition about 6x faster. #############
    
    
    
    
    
    
    %hWait = waitbar( 0, sprintf('Masking vessels...\nMake sure diameters\nand velocities\nhave been estimated.'));
    % ######################################################################
    %
    %             ORIGINAL code commented out Raj 2015
    %
    % ######################################################################
    
    hWait = waitbar( 0, sprintf('Masking vessels...'));
    tic
    nodeMasked = zeros(nNodes,1);
    if updateFlag
        lstUpdate = find(im.nodeTypeUpdated==1);
        foo1 = ismember(im.nodeEdges(:,1),lstUpdate);
        foo2 = ismember(im.nodeEdges(:,2),lstUpdate);
        lstEdges = find(foo1==1 | foo2==1);
    else
        lstEdges = 1:nEdges;
    end
    
    for iiE=1:length(lstEdges)
        iE = lstEdges(iiE);
        %     if isequal(rem(iiE,1000),0)
        %     waitbar(iiE/length(lstEdges),hWait); %waitbars take a long time.
        %     end
        % for debugging:
        %     global i1   %making things global may take a long time.
        %     global i2
        
        i1 = im.nodeEdges(iE,1);
        i2 = im.nodeEdges(iE,2);
        
        if (im.nodeTypeUpdated(i1) | im.nodeTypeUpdated(i2)) & (nodeMasked(i1)==0 | nodeMasked(i2)==0)
            if ~im.maskUseSegDiam
                r1 = im.nodeDiam(i1)/2;
                r1 = min(max(r1,1),maxR);
                r2 = im.nodeDiam(i2)/2;
                r2 = min(max(r2,1),maxR);
            else
                if ~im.nodeSegN(i1)==0 || ~im.nodeSegN(i2)==0
                    foo = max(im.segDiam(im.nodeSegN(i1)),im.segDiam(im.nodeSegN(i2)) )/2;
                    foo = min(max(foo,1),maxR);
                    r1 = foo;
                    r2 = foo;
                end
            end
            p1 = round(im.nodePos(i1,:));
            p2 = round(im.nodePos(i2,:));
            d12 = norm(p2-p1);
            dxyz = (p2-p1);%/d12;
            rd = (r2-r1);%/d12;
            nSteps = 1;
            stepLen = max(r1*im.maskSphereSpacing,1);
            if stepLen<d12
                dxyz = (p2-p1)*stepLen/d12;
                rd = (r2-r1)*stepLen/d12;
                nSteps = floor(d12/stepLen)+1;
            end
            
            if im.nodeType(i1)>0
                nType = im.nodeType(i1)+1;
            else
                nType = 1;
            end
            
            % this should consider im.Hvox
            %        warning( 'this should consider im.Hvox' )
            lst = find(sum((im.nodePos-ones(nNodes,1)*p1).^2,2).^0.5<(r1*im.maskSphereSpacing) );
            nodeMasked(lst) = 1;
            lst = find(sum((im.nodePos-ones(nNodes,1)*p2).^2,2).^0.5<(r2*im.maskSphereSpacing) );
            nodeMasked(lst) = 1;
            
            p = p1;
            r = r1;
            flag = 1;
            while flag
                if norm(round(p)-p2)==0
                    flag = 0;
                end
                pr = round(p);
                rTmp = min(r,maxR);
                rx = ceil(rTmp/im.Hvox(1));
                ry = ceil(rTmp/im.Hvox(2));
                rz = ceil(rTmp/im.Hvox(3));
                if im.maskCircleFlag
                    rz = 0;
                end
                for iX = -rx:+rx
                    for iY = -ry:+ry
                        for iZ = -rz:+rz
                            if norm([iX*im.Hvox(1) iY*im.Hvox(2) iZ*im.Hvox(3)])<=r
                                iix = min(max(pr(1)+iX,1),nx);
                                iiy = min(max(pr(2)+iY,1),ny);
                                iiz = min(max(pr(3)+iZ,1),nz);
                                Vm(iiy,iix,iiz) = 1 + nType;
                                if nType == 4, % if vein
                                    VmDiam(iiy,iix,iiz) = -single(2*r);
                                else
                                    VmDiam(iiy,iix,iiz) = single(2*r);
                                end;
                                if updateFlag
                                    Vupdate(iiy,iix,iiz) = 1;
                                end
                                if isfield(im,'edgeVel')
                                    Vvel(iiy,iix,iiz) = max(round(32*min(edgeVel(iE)/maxEdgeVel,1)),1);%add that to avoid black spot
                                    Vflow(iiy,iix,iiz) = max(round(32*min(edgeFlow(iE)/maxEdgeFlow,1)),1);%add that to avoid black spot
                                    VflowLog(iiy,iix,iiz) = max(round(32*min(edgeFlowLog(iE)/maxEdgeFlowLog,1)),1);%add that to avoid black spot
                                    Vpres(iiy,iix,iiz) = round(32*min(nodePressure(i1)/maxPressure,1));
                                end
                                if flagPO2
                                    Vpo2(iiy,iix,iiz) = round(32*min((nodePO2(i1)-minnodePO2)/(maxnodePO2-minnodePO2),1));
                                end
                                if isfield(im,'edgeBRorderArt')
                                    VbrA(iiy,iix,iiz) = round(32*min(edgeBRorderArt(iE)/maxBRorderArt,1));
                                    VbrV(iiy,iix,iiz) = round(32*min(edgeBRorderVeins(iE)/maxBRorderVeins,1));
                                end
                            end
                        end
                    end
                end
                
                %         for iX = max(pr(1)-rx,1):min(pr(1)+rx,nx)
                %             for iY = max(pr(2)-ry,1):min(pr(2)+ry,ny)
                %                 for iZ = max(pr(3)-rz,1):min(pr(3)+rz,nz)
                %                     if norm([iX-pr(1) iY-pr(2) iZ-pr(3)])<=r
                %                         Vm(iY,iX,iZ) = 1 + nType;
                %                     end
                %                 end
                %             end
                %         end
                if flag
                    p = p + dxyz;
                    r = r + rd;
                end
                nSteps = nSteps - 1;
                if nSteps==0
                    flag = 0;
                end
                if norm(round(p)-p2)==0
                    flag = 0;
                end
                
            end
        end % end of check if node type updated
    end
    toc
    
end %if statement about parallell
% ######################################################################
%
%             end  ORIGINAL code commented out Raj 2015
%
% ######################################################################







if updateFlag % update the image volume
    lst = find(Vupdate==1);
    % overlay Vessel Mask or Velocity
    %    if get(handles.radiobuttonOverlayVesselType,'value');
    im.III(lst) = mod(im.III(lst),32)+(Vm(lst)-1)*32;
    %    elseif get(handles.radiobuttonOverlayVel,'value');
    %        im.III(lst) = 160 + Vvel(lst);
    %    elseif get(handles.radiobuttonOverlayFlow,'value');
    %        im.III(lst) = 160 + Vflow(lst);
    %    elseif get(handles.radiobuttonOverlayPressure,'value');
    %        im.III(lst) = 160 + Vpres(lst);
    %    end
    
    % redraw the portion of the graph updated
    if isfield(im,'edgeFlag')
        edgeFlag = im.edgeFlag;
    else
        edgeFlag = zeros(size(nodeEdges,1),1);
    end
    
    if get(handles.checkboxDisplayGrps,'value') & isfield(im,'nodeGrp')
        if ~get(handles.checkboxHighlightGrp,'value')
            grp = im.nodeGrp;
        else
            grp = 2*ones(nNodes,1);
            grp(find(im.nodeGrp==str2num(get(handles.editHighlightGrp,'string')))) = 1;
        end
    else
        grp = ones(nNodes,1);
    end
    for iiE=1:length(lstEdges)
        ii = lstEdges(iiE);
        pos0 = max(im.nodePos(im.nodeEdges(ii,1),:),1);
        pos1 = max(im.nodePos(im.nodeEdges(ii,2),:),1);
        rsep = norm(pos1-pos0);
        if rsep>0
            cxyz = (pos1-pos0) / rsep;
            rstep = 0;
            pos = pos0;
            while rstep<rsep
                im.III(round(pos(2)),round(pos(1)),max(round(pos(3)),1)) = min(250 - edgeFlag(ii) + grp(im.nodeEdges(ii,1)),254);
                pos = pos + cxyz*0.5;
                if pos(1)<2 & pos(2)<2
                    keyboard
                end
                rstep = rstep + 0.5;
            end
        end
        im.III(round(pos0(2)),round(pos0(1)),round(pos0(3))) = 255;
        im.III(round(pos1(2)),round(pos1(1)),round(pos1(3))) = 255;
    end
end

close(hWait);

im.nodeTypeUpdated = zeros(length(im.nodeType),1);

im.Vm = uint8(Vm);
im.VmDiam = VmDiam;
if isfield(im,'edgeVel')
    im.Vvel = uint8(Vvel);
    im.Vflow = uint8(Vflow);
    im.VflowLog = uint8(VflowLog);
    im.Vpres = uint8(Vpres);
end
if flagPO2
    im.Vpo2 = uint8(Vpo2);
end
if isfield(im,'edgeBRorderArt')
    im.VbrA = uint8(VbrA);
    im.VbrV = uint8(VbrV);
end;
