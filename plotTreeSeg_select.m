function plotTreeSeg_select( edgeSel )
global im

E = im.nodeEdges;
V = im.nodePos;

% hlSel = [];
% for ii=1:length(im.Tree) 
%     jj = im.Tree(ii);
%     hl=plot3( V(E(jj,:),2), V(E(jj,:),1), V(E(jj,:),3), '.-' );
%     set(hl,'ButtonDownFcn',sprintf('plotTreeSeg_select(%d)',jj));
%     set(hl,'linewidth',3);
%     if jj == edgeSel
%         set(hl,'color','r');
%         hlSel = hl;
%     end
% end


%check if there is already a segment selected
if isfield(im,'multiDeleteFlag') 
    edgeSel=[im.deletion_list; edgeSel];
end

foo = get(gca,'children');
hlSel = [];
for ii=1:length(foo);
    for jj=1:length(edgeSel)
        if get(foo(ii),'UserData')==edgeSel(jj)
            hlSel(end+1) = foo(ii);
        end
    end
end
    
set(hlSel,'color','r')

%check for join flag and if flag "ON" then join segments
if isfield(im,'joinflag')
   if im.joinflag==1;
       
       %identify previous segment handle
       previous_hlSel=[];
        for ii=1:length(foo);
            if get(foo(ii),'UserData')==im.previous_edgeSel
                previous_hlSel(end+1) = foo(ii);
            end
        end
       
       woo=menu('Join these two segments?','Yes','No');
       if woo==1;
           
           %%%%%%% PUT CODE FOR JOIN SEGMENTS HERE %%%%%%%%%%%% 
            
                
            %identify the 4 edges of the two segments
            nodeIdx1=E(im.previous_edgeSel,:);
            nodeIdx2=E(edgeSel,:);
              
%             %debug
%             disp('previous selection')
%             nodeIdx1
%             V(nodeIdx1(1),:)
%             V(nodeIdx1(2),:)
%             disp('new selection')
%             nodeIdx2
%             V(nodeIdx2(1),:)
%             V(nodeIdx2(2),:)
            
            node1=V(nodeIdx1(1),:);
            node2=V(nodeIdx1(2),:);
            node3=V(nodeIdx2(1),:);
            node4=V(nodeIdx2(2),:);
            
            %test which combination of edge gives shortest distance (1-3,
            %1-4, 2-3 or 2-4) (node 1 & 2 belongs to the same segment as well as 3 & 4)
            dist13=( (node1(1)-node3(1)).^2+(node1(2)-node3(2)).^2+(node1(3)-node3(3)).^2 )^(1/2);
            dist14=( (node1(1)-node4(1)).^2+(node1(2)-node4(2)).^2+(node1(3)-node4(3)).^2 )^(1/2);
            dist23=( (node2(1)-node3(1)).^2+(node2(2)-node3(2)).^2+(node2(3)-node3(3)).^2 )^(1/2);
            dist24=( (node2(1)-node4(1)).^2+(node2(2)-node4(2)).^2+(node2(3)-node4(3)).^2 )^(1/2);
     
            lowDist=sort([dist13 dist14 dist23 dist24]);
            lowDist=lowDist(1);
            
            %identify the two relevant edges
            if lowDist==dist13
                E1=nodeIdx1(1);E2=nodeIdx2(1);
            elseif lowDist==dist14;
                E1=nodeIdx1(1);E2=nodeIdx2(2);
            elseif lowDist==dist23;
                E1=nodeIdx1(2);E2=nodeIdx2(1);
            elseif lowDist==dist24;
                E1=nodeIdx1(2);E2=nodeIdx2(2);
            end
            
            %update the list
            nodeEdges=im.nodeEdges;
            nodeEdges=[nodeEdges; E1 E2];
            im.nodeEdges=nodeEdges;
            
%             %debug
%             disp('we added')
%             [E1 E2]
%             V(E1,:)
%             V(E2,:)
            
            %finally, remove redundant edges (edges can be duplicated because of the plotting routine ...)
            % point edges
            nodeEdges = im.nodeEdges;
            nodeEdges = nodeEdges(find(nodeEdges(:,1)~=nodeEdges(:,2)),:);
            % redundant edges
            sE = cell(size(nodeEdges,1),1);
            for ii=1:length(nodeEdges)
                if nodeEdges(ii,1)<nodeEdges(ii,2)
                    sE{ii} = sprintf('%05d%05d',nodeEdges(ii,1),nodeEdges(ii,2));
                else
                    sE{ii} = sprintf('%05d%05d',nodeEdges(ii,2),nodeEdges(ii,1));
                end
            end
            [b,i,j]=unique(sE);
            im.nodeEdges = nodeEdges(sort(i),:);        

            im.edgeFlag = zeros(size(im.nodeEdges,1),1);

            %then re-plot
            I = im.III;
            lst = find(I>=65);
            I(lst) = 0;
            [pMin,pMax,im.Tree] = plotTreeSeg(im.nodeEdges,im.nodePos,im.nodeSelected,im.TreeNedges,[1 im.TreeThresh],I);            
       
       else % if woo==1
       
       
           %put back segment in blue
           if ~isempty(hlSel)
               set(hlSel,'color','b');
           end
           
           %put back previous segment in blue
           set(previous_hlSel,'color','b');

       end %if woo==1
       
       %remove the flag
       im=rmfield(im,'joinflag');
       
       %refresh display
       hVesGraph = findobj('tag','vesselGraph');
       if ~isempty(hVesGraph)
           h = guidata(hVesGraph);
           set(h.pushbuttonImageGraph,'enable','on');
       end
       
   end
else
    
if isfield(im,'multiDeleteFlag')
    chm = menu('Action?','Delete all selected','Select for future deletion','Cancel');
    if chm==1
        ch=1;
    elseif chm==2
        ch=2;
    else
        ch=6;
    end
else    
    ch = menu('Action?','Delete all selected','Select for future deletion','Recenter','Collapse','Join (select a 2nd segment)','Cancel');
end

if ch==1
    % DELETE EDGE
    edgeFlag = ones(size(E,1),1);
    edgeFlag(edgeSel) = 0;%simple deletion
    lstEdges = find(edgeFlag==1);
    im.nodeEdges = E(lstEdges,:);
    im.edgeFlag = im.edgeFlag(lstEdges);
    
%     %to remove later LG
%     disp('we removed')
%     foo=E(edgeSel,:)
%     V(foo(1),:)
%     V(foo(2),:)
    
    % check for and remove abandoned nodes
    nodeDel = [];
    for kk=1:length(edgeSel)
        if length(find(E==E(edgeSel(kk),1)))==1
            nodeDel = [nodeDel; E(edgeSel(kk),1)];
        elseif length(find(E==E(edgeSel(kk),2)))==1
            nodeDel = [nodeDel; E(edgeSel(kk),2)];
        end
    end
    
    if ~isempty(nodeDel)
        edgeFlag = ones(size(V,1),1);
        edgeFlag(nodeDel) = 0;
        
        nNodes = 0;
        nodePosTmp = [];
        nodeMap = zeros(size(im.nodePos,1),1);
        for iN = 1:size(im.nodePos,1)
            if edgeFlag(iN)==1
                nNodes = nNodes + 1;
                nodePosTmp(nNodes,:) = im.nodePos(iN,:);
                nodeDiamTmp(nNodes) = im.nodeDiam(iN);
                nodeMap(iN) = nNodes;
            end
        end
        im.nodePos = nodePosTmp;
        im.nodeDiam = nodeDiamTmp;
        im.nodeEdges = nodeMap(im.nodeEdges);
        im.nBflag = 1;
    end
    
    % remove redundant edges
    % point edges
    nodeEdges = im.nodeEdges;
    nodeEdges = nodeEdges(find(nodeEdges(:,1)~=nodeEdges(:,2)),:);
    % redundant edges
    sE = cell(size(nodeEdges,1),1);
    for ii=1:length(nodeEdges)
        if nodeEdges(ii,1)<nodeEdges(ii,2)
            sE{ii} = sprintf('%05d%05d',nodeEdges(ii,1),nodeEdges(ii,2));
        else
            sE{ii} = sprintf('%05d%05d',nodeEdges(ii,2),nodeEdges(ii,1));
        end
    end
    [b,i,j]=unique(sE);
    im.nodeEdges = nodeEdges(sort(i),:);
    
    im.edgeFlag = zeros(size(im.nodeEdges,1),1);

    I = im.III;
    lst = find(I>=65);
    I(lst) = 0;

    %remove the the flag for multiple deletion if there is one
    if isfield(im,'multiDeleteFlag')
        im=rmfield(im,'multiDeleteFlag');
        im=rmfield(im,'deletion_list');
        
    end
    
    %re-plot
    [pMin,pMax,im.Tree] = plotTreeSeg(im.nodeEdges,im.nodePos,im.nodeSelected,im.TreeNedges,[1 im.TreeThresh],I);    
    hVesGraph = findobj('tag','vesselGraph');
    if ~isempty(hVesGraph)
        h = guidata(hVesGraph);
        set(h.pushbuttonImageGraph,'enable','on');
    end

elseif ch==2
    %Select for future deletion (just put a flag)
    
    if isfield(im,'multiDeleteFlag')   
        deletion_list=im.deletion_list;
        deletion_list=[deletion_list; edgeSel];
        im.deletion_list=deletion_list;
    else
        im.multiDeleteFlag=1;
        im.deletion_list=edgeSel;
    end
    
elseif ch==3
    % RECENTER TREE    
    I = im.III;
    lst = find(I>=65);
    I(lst) = 0;

    im.nodeSelected = E(edgeSel,1);
    [pMin,pMax,im.Tree] = plotTreeSeg(im.nodeEdges,im.nodePos,im.nodeSelected,im.TreeNedges,[1 im.TreeThresh],I);

elseif ch==4
    % COLLAPSE NODES
    % collapse around node with most branches, if equal then choose 1
    pts = E(edgeSel,:);
    if length(find(E==pts(2))) > length(find(E==pts(1)))
        pt = pts(2);
    else
        pt = pts(1);
    end
    
    nNodes = size(V,1);
    rho = sum( (V - ones(nNodes,1)*V(pt,:) ).^2,2 ).^0.5; 
    lst = find(rho<im.nodeDiam(pt));
    
    for ii=1:length(lst)
        figure(10)
        plot3(V(lst(ii),2),V(lst(ii),1),V(lst(ii),3),'r*')
    end
    plot3(V(pt,2),V(pt,1),V(pt,3),'g*');
    
    ch = menu('Collapse these points?','Yes','No');
    if ch==1
        % collapse edges to selected node
        for ii=1:length(lst)
            if lst(ii)~=pt
                lst2 = find(E==lst(ii));
                E(lst2) = pt;
            end
        end
        % remove abandoned nodes
        edgeFlag = ones(size(V,1),1);
        edgeFlag(lst) = 0;
        edgeFlag(pt) = 1;
        nNodes = 0;
        nodePosTmp = [];
        nodeMap = zeros(size(V,1),1);
        for iN = 1:size(V,1)
            if edgeFlag(iN)==1
                nNodes = nNodes + 1;
                nodePosTmp(nNodes,:) = V(iN,:);
                nodeDiamTmp(nNodes) = im.nodeDiam(iN);
                nodeMap(iN) = nNodes;
                if iN==pt
                    im.nodeSelected = nNodes;
                end
            end
        end
        nodeEdges = nodeMap(E);
        im.nodePos = nodePosTmp;
        im.nodeDiam = nodeDiamTmp;

        % remove redundant edges
        % point edges
        nodeEdges = nodeEdges(find(nodeEdges(:,1)~=nodeEdges(:,2)),:);
        % redundant edges
        sE = cell(size(nodeEdges,1),1);
        for ii=1:length(nodeEdges)
            if nodeEdges(ii,1)<nodeEdges(ii,2)
                sE{ii} = sprintf('%05d%05d',nodeEdges(ii,1),nodeEdges(ii,2));
            else
                sE{ii} = sprintf('%05d%05d',nodeEdges(ii,2),nodeEdges(ii,1));
            end
        end
        [b,i,j]=unique(sE);
        im.nodeEdges = nodeEdges(sort(i),:);        

        im.edgeFlag = zeros(size(im.nodeEdges,1),1);

    end
    
    % RECENTER TREE    
    I = im.III;
    lst = find(I>=65);
    I(lst) = 0;

    [pMin,pMax,im.Tree] = plotTreeSeg(im.nodeEdges,im.nodePos,im.nodeSelected,im.TreeNedges,[1 im.TreeThresh],I);
    
    hVesGraph = findobj('tag','vesselGraph');
    if ~isempty(hVesGraph)
        h = guidata(hVesGraph);
        set(h.pushbuttonImageGraph,'enable','on');
    end

elseif ch==5
    %PUT THE FLAG ON FOR JOINING (user need to select another segment)
    im.joinflag=1;
    im.previous_edgeSel=edgeSel;
    
%     %debug
%     disp('first selection')
%     E(edgeSel,:)
%     V(E(edgeSel,1),:)
%     V(E(edgeSel,2),:)
else
    if ~isempty(hlSel)
        set(hlSel,'color','b');
    end
    
    %remove any flag we might have put
    if isfield(im,'multiDeleteFlag')
        im=rmfield(im,'multiDeleteFlag');
    end
    if isfield(im,'deletion_list')
        im=rmfield(im,'deletion_list');
    end
    
end
end%if no join to do

