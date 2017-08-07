function ExtractAVCmaskForMCandEdwardCode
%%
% load seed file with processed AVC mask and segment and node diameters...
[fname pname] = uigetfile('*.seed','Select .seed file with AVC mask');
load([pname fname],'-mat');
%%
% select segment diameters or node diameters
diamFlag = menu('Assign diameter based on','whole segment','individual node');

foo = zeros(size(im2.nodePos,1),1);
VmNew = 0.*im2.Vm;
[Nx, Ny, Nz] = size(im2.Vm);
for iX = 1:Nx,
    iX
    for iY = 1:Ny,
        for iZ = 1:Nz,
            if im2.Vm(iX,iY,iZ)>0,
                for ii = 1:size(im2.nodePos,1),
                    foo(ii) = norm(im2.nodePos(ii,:)-[iX iY iZ]);
                end;
                [mval, idx] = min(foo);
                if diamFlag == 2, % diameter of individual node
                    VmNew(iX,iY,iZ)=im2.nodeDiam(idx);
                    if im2.Vm(iX,iY,iZ) == 5,
                        VmNew(iX,iY,iZ) = -VmNew(iX,iY,iZ); %veins get negative values of diameters
                    end;
                else % dimeter based on segment diameter
                    % not trivial for now...
                end;
                
            end;
        end;
    end;
end;

                
