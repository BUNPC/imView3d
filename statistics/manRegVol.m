function varargout = manRegVol(varargin)
% MANREGVOL M-file for manRegVol.fig
%      MANREGVOL, by itself, creates a new MANREGVOL or raises the existing
%      singleton*.
%
%      H = MANREGVOL returns the handle to a new MANREGVOL or the handle to
%      the existing singleton*.
%
%      MANREGVOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANREGVOL.M with the given input arguments.
%
%      MANREGVOL('Property','Value',...) creates a new MANREGVOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before manRegVol_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to manRegVol_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help manRegVol

% Last Modified by GUIDE v2.5 25-Oct-2010 15:49:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @manRegVol_OpeningFcn, ...
                   'gui_OutputFcn',  @manRegVol_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before manRegVol is made visible.
function manRegVol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manRegVol (see VARARGIN)
global im

% Choose default command line output for manRegVol
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if length(varargin)~=2 & length(varargin)~=4 & length(varargin)~=1
    disp( 'USAGE: manRegVol target_image moveable_image [target_voxel_size] [moveable_voxel_size]' );
    disp( '   OR  manRegVol file' );
    close;
end

if exist('im')
    clear im
    global im
end

if length(varargin)==1
    
    if ~exist( varargin{1}, 'file' )
        disp( sprintf( 'File %s does not exist.',varargin{1}) );
        close;
    end
    hwait = waitbar(0,'Loading file');
    load( varargin{1} );
    close(hwait)
    if ~exist('im')
        disp( sprintf(' File %s does not contain structure im',varargi{1}));
        close;
    end
    set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) )
    set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) )
    set(handles.editVoxelSizeTarget,'string',sprintf('%.1f ',im.ht));
    set(handles.editVoxelSizeMoveable,'string',sprintf('%.1f ',im.hm));
    set(handles.pushbuttonSave,'visible','off');
    set(handles.editWid,'string',sprintf('%d',im.h))
    
    im.flagView = 3;
    im.flagButton=1;
    im.MarkerPoint = 0;
    
else
    
    im.It = varargin{1};
    im.Im = varargin{2};
    im.Imo = im.Im;
    
    im.CrangeM = [min(im.Im(:)) max(im.Im(:))];
    im.CrangeT = [min(im.It(:)) max(im.It(:))];
    
    set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) )
    set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) )
    
    if length(varargin)==4
        ht = varargin{3};
        hm = varargin{4};
        if length(ht)==3
            set(handles.editVoxelSizeTarget,'string',sprintf('%.1f ',ht));
        end
        if length(hm)==3
            set(handles.editVoxelSizeMoveable,'string',sprintf('%.1f ',hm));
        end
    end
    
    im.pos = [1 1 1];
    im.h = str2num( get(handles.editWid,'string') );
    
    im.ht = str2num( get(handles.editVoxelSizeTarget,'string') );
    im.hm = str2num( get(handles.editVoxelSizeMoveable,'string') );
    
    [nyt nxt nzt] = size(im.It);
    [nym nxm nzm] = size(im.Im);
    im.posMax = [max(nxt*im.ht(1),nxm*im.hm(1)) max(nyt*im.ht(2),nym*im.hm(2)) max(nzt*im.ht(3),nzm*im.hm(3))];
    
    im.MarkerNum = 0;
    im.MarkerState = 0;  % 0 - no active marker
    % 1 - marker placed in axes 1, waiting for
    %     corresponding marker in axes 2
    % 2 - marker placed in axes 2, waiting for
    %     corresponding marker in axes 1
    im.MarkerPos1 = [];
    im.MarkerPos2 = [];
    im.T = eye(4,4);
    
    im.flagButton = 1;
    im.flagView = 3;
    im.MarkerPoint = 0;
    
end


updateImages( handles );



% --- Outputs from this function are returned to the command line.
function varargout = manRegVol_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonUp.
function pushbuttonUp_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; if fv==1, fv=2;elseif fv==2,fv=1; end
im.pos(fv)=max(im.pos(fv)-1,1);
updateImages( handles );


% --- Executes on button press in pushbuttonDown10.
function pushbuttonUp10_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; if fv==1, fv=2;elseif fv==2,fv=1; end
im.pos(fv)=max(im.pos(fv)-10,1);
updateImages( handles );


% --- Executes on button press in pushbuttonDown.
function pushbuttonDown_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; if fv==1, fv=2;elseif fv==2,fv=1; end
im.pos(fv)=min(im.pos(fv)+1,im.posMax(fv));
updateImages( handles );


% --- Executes on button press in pushbuttonUp10.
function pushbuttonDown10_Callback(hObject, eventdata, handles)
global im
fv = im.flagView; if fv==1, fv=2;elseif fv==2,fv=1; end
im.pos(fv)=min(im.pos(fv)+10,im.posMax(fv));
updateImages( handles );



function editWid_Callback(hObject, eventdata, handles)
global im
foo = str2num( get(handles.editWid,'string') );
if isempty(foo)
    set(handles.editWid,'string',num2str(im.h));
elseif foo==0
    set(handles.editWid,'string',num2str(im.h));
else
    im.h = foo;
    updateImages( handles );
end




function editVoxelSizeTarget_Callback(hObject, eventdata, handles)
global im
foo = str2num( get(handles.editVoxelSizeTarget,'string') );
if isempty(foo)
    set(handles.editWid,'string',num2str(im.h));
elseif length(foo)~=3 || min(foo)<=0
    set(handles.editWid,'string',num2str(im.h));
else
    im.ht = foo;
    [nyt nxt nzt] = size(im.It);
    [nym nxm nzm] = size(im.Im);
    im.posMax = [max(nxt*im.ht(1),nxm*im.hm(1)) max(nyt*im.ht(2),nym*im.hm(2)) max(nzt*im.ht(3),nzm*im.hm(3))];
    updateImages( handles );
end


function editVoxelSizeMoveable_Callback(hObject, eventdata, handles)
global im
foo = str2num( get(handles.editVoxelSizeMoveable,'string') );
if isempty(foo)
    set(handles.editWid,'string',num2str(im.h));
elseif length(foo)~=3 || min(foo)<=0
    set(handles.editWid,'string',num2str(im.h));
else
    im.hm = foo;
    [nyt nxt nzt] = size(im.It);
    [nym nxm nzm] = size(im.Im);
    im.posMax = [max(nxt*im.ht(1),nxm*im.hm(1)) max(nyt*im.ht(2),nym*im.hm(2)) max(nzt*im.ht(3),nzm*im.hm(3))];
    updateImages( handles );
end




% --- Executes when selected object is changed in uipanelProjection1.
function uipanelProjection1_SelectionChangeFcn(hObject, eventdata, handles)
updateImages( handles )


% --- Executes when selected object is changed in uipanelProjection2.
function uipanelProjection2_SelectionChangeFcn(hObject, eventdata, handles)
updateImages( handles )





% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
global im

pos = get(gca,'CurrentPoint');
pos = pos(1,:);

fv = im.flagView;
if fv==2
    pos = [im.pos(1)  pos(1)  pos(2)];
elseif fv==1
    pos = [pos(1)  im.pos(2)  pos(2)];
elseif fv==3
    pos = [pos(1)  pos(2)  im.pos(3)];
end

if im.flagButton==2
    if im.MarkerState==0 | im.MarkerState==2
        pos1 = im.MarkerPos1;
        r = sum( (ones(size(pos1,1),1)*pos-pos1).^2,2).^0.5;
        [foo,idx] = min(r);
        pos = pos1(idx,:);
        im.MarkerPoint = idx;
        im.MarkerState = 1;
        set(handles.textMarkerStatus,'string', sprintf('Move Marker %d in Target Axes',im.MarkerPoint))
        set(handles.pushbuttonMarkerCancel,'visible','on')
        set(handles.pushbuttonMarkerDelete,'visible','on')
        
    elseif im.MarkerState==1
        im.MarkerPos1(im.MarkerPoint,:) = pos;
        im.MarkerState = 0;
        im.MarkerPoint = 0;
        set(handles.textMarkerStatus,'string','');
        set(handles.pushbuttonSave,'visible','on')
        set(handles.pushbuttonMarkerCancel,'visible','off')
        set(handles.pushbuttonMarkerDelete,'visible','off')
    end

elseif im.MarkerState==0 & im.flagButton~=4
    % Add marker and wait for corresponding marker in axes 2
    im.MarkerState = 1;
%    im.MarkerPos1(end+1,1:3) = [pos(1)*im.ht(1) pos(2)*im.ht(2) im.pos(3)+im.h/2  ];
    im.MarkerPos1(end+1,1:3) = pos;
    set(handles.textMarkerStatus,'string', sprintf('Place Marker %d in Moveable Axes',im.MarkerNum+1))
    set(handles.pushbuttonMarkerDelete,'visible','on')
    
elseif im.MarkerState==2 & im.flagButton~=4
    % Marker pair completed
    im.MarkerState = 0;
    im.MarkerNum = im.MarkerNum + 1;
%    im.MarkerPos1(end+1,1:3) = [pos(1)*im.ht(1) pos(2)*im.ht(2) im.pos(3)+im.h/2  ];
    im.MarkerPos1(end+1,1:3) = pos;
    set(handles.textMarkerStatus,'string','');
    set(handles.pushbuttonMarkerDelete,'visible','off')
    set(handles.pushbuttonSave,'visible','on')
end

im.pos = pos;
updateImages( handles )



% --- Executes on mouse press over axes background.
 function axes2_ButtonDownFcn(hObject, eventdata, handles)
global im

pos = get(gca,'CurrentPoint');
%pos = [im.T*[pos(1,:) 1]']';
pos = pos(1,:);
posC = im.pos;

fv = im.flagView;
if fv==2
    pos = [im.T*[posC(1)  pos(1)  pos(2) 1]']';
elseif fv==1
    pos = [im.T*[pos(1)  posC(2)  pos(2) 1]']';
elseif fv==3
    pos = [im.T*[pos(1)  pos(2)  posC(3) 1]']';
end
pos = pos(1:3);

if im.flagButton==2
    if im.MarkerState==0 | im.MarkerState==1
        pos2 = im.MarkerPos2;
        r = sum( (ones(size(pos2,1),1)*pos-pos2).^2,2).^0.5;
        [foo,idx] = min(r);
        pos = pos2(idx,:);
        im.MarkerPoint = -idx;
        im.MarkerState = 2;
        set(handles.textMarkerStatus,'string', sprintf('Move Marker %d in Target Axes',im.MarkerPoint))
        set(handles.pushbuttonMarkerCancel,'visible','on')
        set(handles.pushbuttonMarkerDelete,'visible','on')
        pos = [im.T\[pos 1]']';
        pos = pos(1:3);
        
    elseif im.MarkerState==2
        im.MarkerPos2(-im.MarkerPoint,:) = pos;
        im.MarkerState = 0;
        im.MarkerPoint = 0;
        set(handles.textMarkerStatus,'string','');
        set(handles.pushbuttonSave,'visible','on')
        set(handles.pushbuttonMarkerCancel,'visible','off')
        set(handles.pushbuttonMarkerDelete,'visible','off')
        pos = [im.T\[pos 1]']';
        pos = pos(1:3);
    end
    
elseif im.MarkerState==0 & im.flagButton~=4
    % Add marker and wait for corresponding marker in axes 2
    im.MarkerState = 2;
%    im.MarkerPos2(end+1,1:3) = [pos(1)*im.hm(1) pos(2)*im.hm(2) im.pos(3)+im.h/2  ];
    im.MarkerPos2(end+1,1:3) = pos;
    set(handles.textMarkerStatus,'string', sprintf('Place Marker %d in Target Axes',im.MarkerNum+1))
    set(handles.pushbuttonMarkerDelete,'visible','on')
    pos = [im.T\[pos 1]']';
    pos = pos(1:3);

elseif im.MarkerState==1 & im.flagButton~=4
    % Marker pair completed
    im.MarkerState = 0;
    im.MarkerNum = im.MarkerNum + 1;
%    im.MarkerPos2(end+1,1:3) = [pos(1)*im.hm(1) pos(2)*im.hm(2) im.pos(3)+im.h/2  ];
    im.MarkerPos2(end+1,1:3) = pos;
    set(handles.textMarkerStatus,'string','');
    set(handles.pushbuttonMarkerDelete,'visible','off')
    set(handles.pushbuttonSave,'visible','on')
    pos = [im.T\[pos 1]']';
    pos = pos(1:3);
    
else
    pos = [im.T\[pos 1]']';
    pos = pos(1:3);
    
end

%pos = [im.T\[pos 1]']';
im.pos = pos(1:3);
updateImages( handles )


% --- Executes on button press in pushbuttonMarkerDelete.
function pushbuttonMarkerDelete_Callback(hObject, eventdata, handles)
global im

if im.MarkerPoint == 0
    if im.MarkerState==1
        im.MarkerPos1(end,:) = [];
    elseif im.MarkerState==2
        im.MarkerPos2(end,:) = [];
    end
else
    im.MarkerPos1(abs(im.MarkerPoint),:) = [];
    im.MarkerPos2(abs(im.MarkerPoint),:) = [];
    im.MarkerPoint = 0;
    set(handles.pushbuttonSave,'visible','on')
end


im.MarkerState=0;
set(handles.textMarkerStatus,'string','');
set(handles.pushbuttonMarkerDelete,'visible','off')
set(handles.pushbuttonMarkerCancel,'visible','off')
updateImages( handles )


% --- Executes on button press in pushbuttonRegisterXY.
function pushbuttonRegisterXY_Callback(hObject, eventdata, handles)
global im

pT = im.MarkerPos1;
pM = im.MarkerPos2;


% Get XY transformation
A = [pT(:,1:2) ones(size(pT,1),1)];
Ttmp = [A\pM(:,1) A\pM(:,2) [0 0 1]']';
T = eye(4,4);
T(1:2,1:2) = Ttmp(1:2,1:2);
T(1:2,4) = Ttmp(1:2,3);

% transform the moveable image
im.Im = vol2vol(im.Imo,im.Imo,T);
im.T = T;

updateImages( handles )
set(handles.pushbuttonSave,'visible','on');


% --- Executes on button press in pushbuttonRegisterXYZ.
function pushbuttonRegisterXYZ_Callback(hObject, eventdata, handles)
global im

pT = im.MarkerPos1;
pM = im.MarkerPos2;


% Get XYZ transformation
A = [pT(:,1:3) ones(size(pT,1),1)];
Ttmp = [A\pM(:,1) A\pM(:,2) A\pM(:,3)]';
T = eye(4,4);
T(1:3,1:4) = Ttmp(1:3,1:4);

% transform the moveable image
im.Im = vol2vol(im.Imo,im.Imo,T);
im.T = T;

updateImages( handles )
set(handles.pushbuttonSave,'visible','on');

    

% --- Executes on button press in pushbuttonAxes2Linc.
function pushbuttonAxes2Linc_Callback(hObject, eventdata, handles)
global im
im.CrangeM(1) = min(im.CrangeM(1)*2,im.CrangeM(2));
if im.CrangeM(1)==0
    im.CrangeM(1) = 1;
end
set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
updateImages( handles )


% --- Executes on button press in pushbuttonAxes2Udec.
function pushbuttonAxes2Udec_Callback(hObject, eventdata, handles)
global im
im.CrangeM(2) = max(im.CrangeM(2)/2,im.CrangeM(1));
set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
updateImages( handles )


% --- Executes on button press in pushbuttonAxes2Uinc.
function pushbuttonAxes2Uinc_Callback(hObject, eventdata, handles)
global im
im.CrangeM(2) = im.CrangeM(2)*2;
set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
updateImages( handles )


% --- Executes on button press in pushbuttonAxes2Ldec.
function pushbuttonAxes2Ldec_Callback(hObject, eventdata, handles)
global im
im.CrangeM(1) = max(im.CrangeM(1)/2,0);
if im.CrangeM(1)<=1
    im.CrangeM(1) = 0;
end
set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
updateImages( handles )



function editAxes2Crange_Callback(hObject, eventdata, handles)
global im

crange = str2num(get(handles.editAxes2Crange,'string'));
if length(crange)~=2
    set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
elseif crange(1)<0 | crange(1)>crange(2)
    set(handles.editAxes2Crange,'string',sprintf('%d ',im.CrangeM) );
else
    im.CrangeM = crange;
    updateImage( handles );
end



function editAxes1Crange_Callback(hObject, eventdata, handles)
global im

crange = str2num(get(handles.editAxes1Crange,'string'));
if length(crange)~=2
    set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
elseif crange(1)<0 | crange(1)>crange(2)
    set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
else
    im.CrangeT = crange;
    updateImage( handles );
end


% --- Executes on button press in pushbuttonAxes1Ldec.
function pushbuttonAxes1Ldec_Callback(hObject, eventdata, handles)
global im
if im.CrangeT(1)<=1
    im.CrangeT(1) = 0;
end
im.CrangeT(1) = max(im.CrangeT(1)/2,0);
set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
updateImages( handles )


% --- Executes on button press in pushbuttonAxes1Uinc.
function pushbuttonAxes1Uinc_Callback(hObject, eventdata, handles)
global im
im.CrangeT(2) = im.CrangeT(2)*2;
set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
updateImages( handles )



% --- Executes on button press in pushbuttonAxes1Udec.
function pushbuttonAxes1Udec_Callback(hObject, eventdata, handles)
global im
im.CrangeT(2) = max(im.CrangeT(2)/2,im.CrangeT(1));
set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
updateImages( handles )



% --- Executes on button press in pushbuttonAxes1Linc.
function pushbuttonAxes1Linc_Callback(hObject, eventdata, handles)
global im
im.CrangeT(1) = min(im.CrangeT(1)*2,im.CrangeT(2));
if im.CrangeT(1)==0
    im.CrangeT(1) = 1;
end
set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
updateImages( handles )





%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE IMAGES
%%%%%%%%%%%%%%%%%%%%%%%
function updateImages( handles )
global im

fv = im.flagView;
if fv==1
    fv0 = 2;
    fv1 = 1;
    fv2 = 3;
    set(handles.textRange,'string',sprintf('Y [%.1f %.1f]',im.pos(fv0)-im.h/2,im.pos(fv0)+im.h/2))
elseif fv==2
    fv0 = 1;
    fv1 = 2;
    fv2 = 3;
    set(handles.textRange,'string',sprintf('X [%.1f %.1f]',im.pos(fv0)-im.h/2,im.pos(fv0)+im.h/2))
elseif fv==3
    fv0 = 3;
    fv1 = 1;
    fv2 = 2;
    set(handles.textRange,'string',sprintf('Z [%.1f %.1f]',im.pos(fv0)-im.h/2,im.pos(fv0)+im.h/2))
end

% get overlay options
flagOverlay = get(handles.checkboxOverlay,'value');
overlayRight = 0;
if strcmpi(get(handles.togglebuttonLeftRight,'string'),'Right')
    overlayRight = 1;
end
overlayInvertRight = get(handles.checkboxInvertRight,'value');
overlayInvertLeft = get(handles.checkboxInvertLeft,'value');
overlayChColor = get(handles.popupmenuChColor,'value');

% Axes 1
axes(handles.axes1)
pos = [1:size(im.It,fv)]*im.ht(fv0);
lst = find(pos>=(im.pos(fv0)-im.h/2) & pos<=(im.pos(fv0)+im.h/2));
if get(handles.radiobuttonMin1,'value')
    if fv==1
        foo = squeeze(min(im.It(lst,:,:),[],1))';
    elseif fv==2
        foo = squeeze(min(im.It(:,lst,:),[],2))';
    elseif fv==3
        foo = min(im.It(:,:,lst),[],3);
    end
else
    if fv==1
        foo = squeeze(max(im.It(lst,:,:),[],1))';
    elseif fv==2
        foo = squeeze(max(im.It(:,lst,:,:),[],2))';
    elseif fv==3
        foo = max(im.It(:,:,lst),[],3);
    end
end
imagesc([1:size(foo,2)]*im.ht(fv1),[1:size(foo,1)]*im.ht(fv2),foo, im.CrangeT)
axis off
axis image
colormap gray

hold on
for ii=1:size(im.MarkerPos1,1)
    if im.MarkerPos1(ii,fv0)>=(im.pos(fv0)-im.h/2) & im.MarkerPos1(ii,fv0)<=(im.pos(fv0)+im.h/2)
%        ht = text(im.MarkerPos1(ii,1)/im.ht(1),im.MarkerPos1(ii,2)/im.ht(2),sprintf('%d',ii));
        ht = text(im.MarkerPos1(ii,fv1),im.MarkerPos1(ii,fv2),sprintf('%d',ii));
        set(ht,'color','r');
        if im.MarkerPoint==ii
            set(ht,'color','g');
        end
    end
end
if size(im.MarkerPos1,1)<size(im.MarkerPos2,1)
    pM = [im.T\[im.MarkerPos2(end,:) 1]']';
    ht = text(pM(fv1),pM(fv2),sprintf('%d',size(im.MarkerPos2,1)));
    set(ht,'color',[1 1 0]);
end
if im.flagButton==4 | im.flagButton==3
    pM = im.pos;
    ht = text(pM(fv1),pM(fv2),'X');
    set(ht,'color',[1 0 1]);
end
hold off

% plot spots 
if isfield(im,'Spot')
    if isfield(im.Spot,'flagDisplayCentroids')
        if im.Spot.flagDisplayCentroids==1
            posC = im.Spot.stats.posC;
            rad = im.Spot.stats.Rad;
            lst = find(posC(:,fv0)>=(im.pos(fv0)-rad) & posC(:,fv0)<=(im.pos(fv0)+rad));
            for ii=1:length(lst)
                ht = text(posC(lst(ii),fv1),posC(lst(ii),fv2),'x');
                set(ht,'color',[1 0 0]);
                set(ht,'HorizontalAlignment','center')
                set(ht,'clipping','on')
            end
        end
    end
end





% Axes 2
axes(handles.axes2)
pos = [1:size(im.Im,fv)]*im.hm(fv0);
lst = find(pos>=(im.pos(fv0)-im.h/2) & pos<=(im.pos(fv0)+im.h/2));
if get(handles.radiobuttonMin2,'value')
    if fv==1
        foo2 = squeeze(min(im.Im(lst,:,:),[],1))';
    elseif fv==2
        foo2 = squeeze(min(im.Im(:,lst,:,:),[],2))';
    elseif fv==3
        foo2 = min(im.Im(:,:,lst),[],3);
    end
else
    if fv==1
        foo2 = squeeze(max(im.Im(lst,:,:),[],1))';
    elseif fv==2
        foo2 = squeeze(max(im.Im(:,lst,:,:),[],2))';
    elseif fv==3
        foo2 = max(im.Im(:,:,lst),[],3);
    end
end
if flagOverlay==0 | overlayRight==0
    imagesc([1:size(foo2,fv2)]*im.hm(fv1),[1:size(foo2,fv1)]*im.hm(fv2),foo2, im.CrangeM)
    axis off
    axis image
    colormap gray
    
    if ~isempty(im.MarkerPos2)
        hold on
        pM = [im.T\[im.MarkerPos2 ones(size(im.MarkerPos2,1),1)]']';
        for ii=1:size(im.MarkerPos2,1)
            if pM(ii,fv0)>=(im.pos(fv0)-im.h/2) & pM(ii,fv0)<=(im.pos(fv0)+im.h/2)
                %            ht = text(pM(ii,1)/im.hm(1),pM(ii,2)/im.hm(2),sprintf('%d',ii));
                ht = text(pM(ii,fv1),pM(ii,fv2),sprintf('%d',ii));
                set(ht,'color','r');
                if im.MarkerPoint==-ii
                    set(ht,'color','g');
                end
            end
        end
        if size(im.MarkerPos2,1)<size(im.MarkerPos1,1)
            pT = im.MarkerPos1(end,:);
            ht = text(pT(fv1),pT(fv2),sprintf('%d',size(im.MarkerPos1,1)));
            set(ht,'color',[0 1 1]);
        end
        if im.flagButton==4 | im.flagButton==3
            %        pM = im.T\[im.pos 1]';
            pM = im.pos;
            ht = text(pM(fv1),pM(fv2),'X');
            set(ht,'color',[1 0 1]);
        end
        hold off
    end
    
    % plot spots
    if isfield(im,'Spot')
        if isfield(im.Spot,'flagDisplayCentroids')
            if im.Spot.flagDisplayCentroids==1
                posC = im.Spot.stats.posC;
                rad = im.Spot.stats.Rad;
                lst = find(posC(:,fv0)>=(im.pos(fv0)-rad) & posC(:,fv0)<=(im.pos(fv0)+rad));
                for ii=1:length(lst)
                    ht = text(posC(lst(ii),fv1),posC(lst(ii),fv2),'x');
                    set(ht,'color',[1 0 0]);
                    set(ht,'HorizontalAlignment','center')
                    set(ht,'clipping','on')
                end
            end
        end
    end
end


% check for overlay
if flagOverlay==1
    foo = max(min((foo-im.CrangeT(1)) / (im.CrangeT(2) - im.CrangeT(1)),1),0);
    boo = im.CrangeM(1);
    foo2 = max(min((foo2-double(boo)) / double(im.CrangeM(2) - im.CrangeM(1)),1),0);
    
    if overlayInvertLeft==1
        foo = 1-foo;
    end
    if overlayInvertRight==1
        foo2 = 1-foo2;
    end
    
    [ny,nx] = size(foo);
    [ny2,nx2] = size(foo2);
    
    fooO = zeros(max(ny,ny2),max(nx,nx2),3);
    fooO(1:ny,1:nx,2) = foo;
    fooO(1:ny2,1:nx2,1) = foo2;
    fooO(:,:,3) = 0;
    
    if overlayRight==1
        axes(handles.axes2)
    else
        axes(handles.axes1)
    end
    imagesc( fooO );
end


% set up button down callbacks
axes(handles.axes1)
if im.flagButton==1 | im.flagButton==2 | im.flagButton==4
    fhandle1='manRegVol(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
    set(handles.axes1,'ButtonDownFcn',fhandle1);
    set(get(handles.axes1,'children'), 'ButtonDownFcn',fhandle1);
end

axes(handles.axes2)
if im.flagButton==1 | im.flagButton==2 | im.flagButton==4
    zoom off
    fhandle2='manRegVol(''axes2_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
    set(handles.axes2,'ButtonDownFcn',fhandle2);
    set(get(handles.axes2,'children'), 'ButtonDownFcn',fhandle2);
elseif im.flagButton == 3
    h = zoom;
    set(h,'enable','on');
    linkaxes([handles.axes1 handles.axes2],'xy')
end



% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
global im

[fn,pn] = uiputfile( '*.mat', 'Save file');
if fn~=0
    hwait=waitbar(0,'Saving');
    save([pn fn],'im');
    close(hwait)
    set(handles.pushbuttonSave,'visible','off')
end


% --- Executes when selected object is changed in uipanelAxesButton.
function uipanelAxesButton_SelectionChangeFcn(hObject, eventdata, handles)
global im

if get(handles.radiobuttonMarker,'value')
    im.flagButton = 1; % Marker
elseif get(handles.radiobuttonMarkerMove,'value')
    im.flagButton = 2; % move
elseif get(handles.radiobuttonZoom,'value')
    im.flagButton = 3; % Zoom
elseif get(handles.radiobuttonPoint,'value')
    im.flagButton = 4; % point
end

im.MarkerPoint = 0;
    
if im.MarkerState==1
    im.MarkerPos1(end,:) = [];
elseif im.MarkerState==2
    im.MarkerPos2(end,:) = [];
end

im.MarkerState=0;
set(handles.textMarkerStatus,'string','');
set(handles.pushbuttonMarkerDelete,'visible','off')
set(handles.pushbuttonMarkerCancel,'visible','off')
updateImages( handles )


% --- Executes when selected object is changed in uipanelView.
function uipanelView_SelectionChangeFcn(hObject, eventdata, handles)
global im

if get(handles.radiobuttonXY,'value')
    im.flagView=3;
elseif get(handles.radiobuttonXZ,'value')
    im.flagView=1;
elseif get(handles.radiobuttonYZ,'value')
    im.flagView=2;
end    
updateImages( handles )


% --- Executes on button press in pushbuttonMarkerCancel.
function pushbuttonMarkerCancel_Callback(hObject, eventdata, handles)
global im

im.MarkerPoint = 0;
im.MarkerState=0;
set(handles.textMarkerStatus,'string','');
set(handles.pushbuttonMarkerCancel,'visible','off')
updateImages( handles )


% --- Executes on button press in checkboxOverlay.
function checkboxOverlay_Callback(hObject, eventdata, handles)
updateImages( handles )


% --- Executes on button press in togglebuttonLeftRight.
function togglebuttonLeftRight_Callback(hObject, eventdata, handles)
if strcmpi(get(hObject,'string'),'right')
    set(hObject,'string','Left');
else
    set(hObject,'string','Right');
end
    

% --- Executes on button press in checkboxInvertLeft.
function checkboxInvertLeft_Callback(hObject, eventdata, handles)
updateImages( handles )


% --- Executes on button press in checkboxInvertRight.
function checkboxInvertRight_Callback(hObject, eventdata, handles)
updateImages( handles )


% --- Executes on selection change in popupmenuChColor.
function popupmenuChColor_Callback(hObject, eventdata, handles)
updateImages( handles )


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_Spots_Callback(hObject, eventdata, handles)

spotGUI();

    
    
