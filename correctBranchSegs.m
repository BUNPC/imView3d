function correctBranchSegs

[seedname seedpath] = uigetfile( '*.seed','Load seed file');
load([seedpath seedname],'-mat','im2');
im = im2;
clear im2;

[ff pp] = uigetfile('*.txt','Select File name with Pial Segments info');
pialseg = load([pp ff],'-ascii'); % 3 columns, segment number, branch order, vessel type (A=1, 2=C, 3=V)

numberOfSegs=size(pialseg,1);

for ii=1:numberOfSegs
    nodeIdx = find(pialseg(ii,1)==im.nodeSegN);
    foo = im.nB(nodeIdx);
    fooIdx = find(foo==2);
    segNodeIdx = nodeIdx(fooIdx(1));
    newTable(ii,:) = [pialseg(ii,1) pialseg(ii,2) pialseg(ii,3) im.nodePos(segNodeIdx,1) im.nodePos(segNodeIdx,2) im.nodePos(segNodeIdx,3)];
end

[newFile newPath] = uiputfile('*.txt','Save new text file',sprintf(ff(1:(end-4))));

dlmwrite([newPath newFile],newTable,'delimiter','\t');
