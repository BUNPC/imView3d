function AddXYZcoordinatesToBRorderTable
% load .seed file and .txt file with the pial vessel Branching Orders
% adds XYZ coordinate of mid-point node in each pial segment to .txt table,
% so that procedures that calculate BR order do not need unreliable
% segment numbers of pial vessels as a starting point...

[BRfname, BRpname] = uigetfile('*.txt','Select TEXT file with pial vessel Brancing Orders');
oldBRordTable = load([BRpname BRfname],'-ascii'); % nx3 matrix: [segment #, BR order, Vessel Type]
%%
[seedFname, seedPname] = uigetfile('*.seed','Select .seed file with graph structure in sync with BR order table');
load([seedPname seedFname],'-mat'); % loads im2 structure...
%%
nSeg = size(oldBRordTable,1);
newBRordTable = zeros(nSeg,6);
newBRordTable(1:nSeg,1:3) = oldBRordTable(1:nSeg,1:3);
for ii=1:nSeg,
    seg = oldBRordTable(ii,1);
    endNodes = im2.segEndNodes(seg,:);
    nodes = find(im2.nodeSegN == seg);
    lst = find(nodes==endNodes(1) | nodes==endNodes(2));
    nodes(lst)=[];
    if ~isempty(nodes),
        nodes = sort(nodes);
        pos = im2.nodePos(nodes(round((length(nodes)+1)/2)),:);
    else % segment has only two nodes!!! which one to use? This doesn't solve the problem, just pick one...
        pos = im2.nodePos(endNodes(1),:);
    end;
    newBRordTable(ii,4:6) = pos(1:3);
end
%%
[fname, pname] = uiputfile('*.txt','Type the name for a new Pial BR order table',[BRpname BRfname(1:end-4) '_withSegXYZ.txt']);
save([pname fname],'newBRordTable','-ascii');