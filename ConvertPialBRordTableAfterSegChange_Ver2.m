function ConvertPialBRordTableAfterSegChange
% after segment numbers are changed in .seed file, we need to convert
% segment numbers in Pial Vessel Branching Order Table, based on the
% segment mid-point XYZ coordinates

[BRfname, BRpname] = uigetfile('*.txt','Select TEXT file with pial vessel Brancing Orders and Segment XYZ ');
oldBRordTable = load([BRpname BRfname],'-ascii'); % nx6 matrix: [segment #, BR order, Vessel Type, node X, Y, Z]
%%
[seedFname, seedPname] = uigetfile('*.seed','Select .seed file with new graph structure and changed seg. numbers');
load([seedPname seedFname],'-mat'); % loads im2 structure...
%%
nSeg = size(oldBRordTable,1);
newBRordTable = oldBRordTable;
h=50;
delLst = [];
delLstCount = 0;
for ii=1:nSeg,
    nodePos = oldBRordTable(ii,4:6);
    lst = find( abs(nodePos(1)-im2.nodePos(:,1))<h & abs(nodePos(2)-im2.nodePos(:,2))<h & ...
        abs(nodePos(3)-im2.nodePos(:,3))<h );
    if ~isempty(lst)
        rsep = sum( (ones(length(lst),1)*nodePos - im2.nodePos(lst,:)).^2, 2).^0.5;
        [foo,idx] = min(rsep);
        node = lst(idx(1));
        seg = im2.nodeSegN(node);
        pos = im2.nodePos(node,:);
        newBRordTable(ii,1) = seg; % new segment number
        newBRordTable(ii,4:6) = pos;
    else % if close node was not found in new seed file
        delLstCount = delLstCount+1; % remember to delete 
        delLst(delLstCount) = ii;
    end;
end;

if delLstCount,
    newBRordTable(delLst,:)=[];
end;

%%
[fname, pname] = uiputfile('*.txt','Type the name for a new Pial BR order table',[BRpname BRfname(1:end-4) '_new.txt']);
save([pname fname],'newBRordTable','-ascii');