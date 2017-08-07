function varargout = doStatistics(varargin)
% DOSTATISTICS M-file for doStatistics.fig
%      DOSTATISTICS, by itself, creates a new DOSTATISTICS or raises the existing
%      singleton*.
%
%      H = DOSTATISTICS returns the handle to a new DOSTATISTICS or the handle to
%      the existing singleton*.
%
%      DOSTATISTICS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DOSTATISTICS.M with the given input arguments.
%
%      DOSTATISTICS('Property','Value',...) creates a new DOSTATISTICS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before doStatistics_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to doStatistics_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help doStatistics

% Last Modified by GUIDE v2.5 12-May-2011 13:46:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @doStatistics_OpeningFcn, ...
                   'gui_OutputFcn',  @doStatistics_OutputFcn, ...
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


% --- Executes just before doStatistics is made visible.
function doStatistics_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to doStatistics (see VARARGIN)
global im

% Choose default command line output for doStatistics
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

[fname pname] = uigetfile('*.seed','Select Vascular Graph .seed File with AVC mask');

hwait = waitbar(0,'Loading file');
load([pname fname],'-mat' );
close(hwait)
if (~exist('im2') && isempty(im)),
    msgbox('No im2 structure in file','','error');
    close;
end

if ~isempty(im), % data was perviously already processed by this GUI
    set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
    set(handles.editVoxelSizeTarget,'string',sprintf('%.2f ',im.ht));
    im.flagButton = 1;
    im.flagView = 3;
    im.h = str2num( get(handles.editWid,'string') );
    [nyt nxt nzt] = size(im.It);
    im.posMax = [nxt*im.ht(1)  nyt*im.ht(2) nzt*im.ht(3)];
    im.MarkerNum = 0;
    im.MarkerState = 0;  % 0 - no active marker
        % 1 - marker placed in axes 1, waiting for
        %     corresponding marker in axes 2
        % 2 - marker placed in axes 2, waiting for
        %     corresponding marker in axes 1
    im.MarkerPos1 = [];
    im.MarkerPos2 = [];
    im.T = eye(4,4);
    im.MarkerPoint = 0;
    %set(handles.pushbuttonSave,'visible','off');
else
    im.It = im2.I; % im is new structure with parameters needed to select brain surface and inspect statistical results
    im.CrangeT = [min(im.It(:)) max(im.It(:))];
    set(handles.editAxes1Crange,'string',sprintf('%d ',im.CrangeT) );
    im.ht = im2.Hvox;
    set(handles.editVoxelSizeTarget,'string',sprintf('%.2f ',im.ht));
    %set(handles.pushbuttonSave,'visible','off');
    im.pos = [1 1 1];
    im.h = str2num( get(handles.editWid,'string') );
    [nyt nxt nzt] = size(im.It);
    im.posMax = [nxt*im.ht(1)  nyt*im.ht(2) nzt*im.ht(3)];
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

    im.im2 = im2; % here are all data about graph from imView...
    clear im2;
end;

updateImages( handles );



% --- Outputs from this function are returned to the command line.
function varargout = doStatistics_OutputFcn(hObject, eventdata, handles) 
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
    im.posMax = [nxt*im.ht(1)  nyt*im.ht(2)   nzt*im.ht(3)];
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
    im.posMax = [nxt*im.ht(1)   nyt*im.ht(2)    nzt*im.ht(3)];
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
        %set(handles.pushbuttonSave,'visible','on')
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
    %set(handles.pushbuttonSave,'visible','on')
end

im.pos = pos;
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
    %set(handles.pushbuttonSave,'visible','on')
end


im.MarkerState=0;
set(handles.textMarkerStatus,'string','');
set(handles.pushbuttonMarkerDelete,'visible','off')
set(handles.pushbuttonMarkerCancel,'visible','off')
updateImages( handles )




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

% flagOverlay = get(handles.checkboxOverlay,'value');
% overlayRight = 0;
% if strcmpi(get(handles.togglebuttonLeftRight,'string'),'Right')
%     overlayRight = 1;
% end
% overlayInvertRight = get(handles.checkboxInvertRight,'value');
% overlayInvertLeft = get(handles.checkboxInvertLeft,'value');
% overlayChColor = get(handles.popupmenuChColor,'value');

% Axes 1
axes(handles.axes1);

I = im.It;
if get(handles.checkbox_viewSurface,'value') == 1,
    if isfield(im,'Isurf'),
        I = im.Isurf;
    else
        msgbox('No surface fitted yet!','','error');
    end;
end;

pos = [1:size(im.It,fv)]*im.ht(fv0);
lst = find(pos>=(im.pos(fv0)-im.h/2) & pos<=(im.pos(fv0)+im.h/2));
if get(handles.radiobuttonMin1,'value')
    if fv==1
        foo = squeeze(min(I(lst,:,:),[],1))';
    elseif fv==2
        foo = squeeze(min(I(:,lst,:),[],2))';
    elseif fv==3
        foo = min(I(:,:,lst),[],3);
    end
else
    if fv==1
        foo = squeeze(max(I(lst,:,:),[],1))';
    elseif fv==2
        foo = squeeze(max(I(:,lst,:,:),[],2))';
    elseif fv==3
        foo = max(I(:,:,lst),[],3);
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

% % plot spots 
% if isfield(im,'Spot')
%     if isfield(im.Spot,'flagDisplayCentroids')
%         if im.Spot.flagDisplayCentroids==1
%             posC = im.Spot.stats.posC;
%             rad = im.Spot.stats.Rad;
%             lst = find(posC(:,fv0)>=(im.pos(fv0)-rad) & posC(:,fv0)<=(im.pos(fv0)+rad));
%             for ii=1:length(lst)
%                 ht = text(posC(lst(ii),fv1),posC(lst(ii),fv2),'x');
%                 set(ht,'color',[1 0 0]);
%                 set(ht,'HorizontalAlignment','center')
%                 set(ht,'clipping','on')
%             end
%         end
%     end
% end








% check for overlay
% if flagOverlay==1
%     foo = max(min((foo-im.CrangeT(1)) / (im.CrangeT(2) - im.CrangeT(1)),1),0);
%     boo = im.CrangeM(1);
%     foo2 = max(min((foo2-double(boo)) / double(im.CrangeM(2) - im.CrangeM(1)),1),0);
%     
%     if overlayInvertLeft==1
%         foo = 1-foo;
%     end
%     if overlayInvertRight==1
%         foo2 = 1-foo2;
%     end
%     
%     [ny,nx] = size(foo);
%     [ny2,nx2] = size(foo2);
%     
%     fooO = zeros(max(ny,ny2),max(nx,nx2),3);
%     fooO(1:ny,1:nx,2) = foo;
%     fooO(1:ny2,1:nx2,1) = foo2;
%     fooO(:,:,3) = 0;
%     
%     if overlayRight==1
%         axes(handles.axes2)
%     else
%         axes(handles.axes1)
%     end
%     imagesc( fooO );
% end


% set up button down callbacks
axes(handles.axes1)
if im.flagButton==1 | im.flagButton==2 | im.flagButton==4
    fhandle1='manRegVol(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))';
    set(handles.axes1,'ButtonDownFcn',fhandle1);
    set(get(handles.axes1,'children'), 'ButtonDownFcn',fhandle1);
end





% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
global im

[fn,pn] = uiputfile( '*.seed', 'Save file');
if fn~=0
    hwait=waitbar(0,'Saving');
    save([pn fn],'im');
    close(hwait)
    %set(handles.pushbuttonSave,'visible','off')
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


% --- Executes on button press in pushbutton_SelSurfacePts.
function pushbutton_SelSurfacePts_Callback(hObject, eventdata, handles)
global im;

if (im.flagView > 2),
    msgbox('View must be XZ or YZ','','error');
else
    [xfoo zfoo] = ginput; % if flagView = 2, xfoo is Y coordinate
    nPts = length(xfoo);
    if nPts > 0,
        if ~isfield(im,'surface'),
            im.surface.nPts = nPts;
        elseif ~isfield(im.surface,'nPts'),
            im.surface.nPts = nPts;
        else
            im.surface.nPts = im.surface.nPts + nPts;
        end;
        if im.flagView == 1, % XZ projection
            ypos = round(im.pos(2)); % it looks that im.pos has wrong order as xyz instead yxz
            oldNpts = im.surface.nPts - nPts;
            im.surface.YXZ(oldNpts+1:oldNpts+nPts,1) = ypos;
            im.surface.YXZ(oldNpts+1:oldNpts+nPts,2) = round(xfoo);
            im.surface.YXZ(oldNpts+1:oldNpts+nPts,3) = round(zfoo);
        else  % projection is YZ
            xpos = round(im.pos(1)); % it looks that im.pos has wrong order as xyz instead yxz
            oldNpts = im.surface.nPts - nPts;
            im.surface.YXZ(oldNpts+1:oldNpts+nPts,2) = xpos;
            im.surface.YXZ(oldNpts+1:oldNpts+nPts,1) = round(xfoo);
            im.surface.YXZ(oldNpts+1:oldNpts+nPts,3) = round(zfoo);
        end;
    updateImages( handles );
    end;
end;


% --- Executes on button press in pushbutton_Fit2Dsurface.
function pushbutton_Fit2Dsurface_Callback(hObject, eventdata, handles)
global im;

if ~isfield(im,'surface'),
    msgbox('Points must be selected first','','error');
elseif ~isfield(im.surface,'YXZ'),
    msgbox('More than 6 points must be selected first','','error');
elseif im.surface.nPts < 6,
    msgbox('More than 6 points must be selected first','','error');
else
    y = im.surface.YXZ(:,1);
    x = im.surface.YXZ(:,2);
    z = im.surface.YXZ(:,3);
    Zprofile = @(c,y,x) (c(1).*y.*y + c(2).*x.*x + c(3).*y.*x + c(4).*y + c(5).*x +c(6) );
    im.Zprofile = Zprofile;
    Zprofile_error = @(c,y,x,z) (z - (c(1).*y.*y + c(2).*x.*x + c(3).*y.*x + c(4).*y + c(5).*x +c(6) ) ) ;
    MyOptions = optimset('Display','off','LevenbergMarquardt','on');
    c0 = [1 1 1 1 1 0];
    [cAll,resnorm,resid,exitflag,output,lambda,jacobian_c] = lsqnonlin(Zprofile_error, c0, [],[],MyOptions,y,x,z);
    im.cAll = cAll;
    im.Isurf = im.It;
    maxC = max(im.It(:));
    
    [ny, nx, nz] = size(im.It);
    Zsurf = zeros(ny,nx);
    hw = waitbar(0,'Calc. Z surface...');
   
    for y = 1:ny,
        waitbar(y/ny);
        for x = 1:nx;
            z = cAll(1)*y*y + cAll(2)*x*x + cAll(3)*y*x + cAll(4)*y + cAll(5)*x +cAll(6) ;
            z = abs(round(z));
            Zsurf(y,x) = z;
            if z == 0, z = 1; end;
            im.Isurf(y,x,z) = maxC;
        end;
    end;
    im.Zsurf = Zsurf;
    close(hw);
end;
    
    


% --- Executes on button press in checkbox_viewSurface.
function checkbox_viewSurface_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox_viewSurface




% --- Executes on button press in radiobutton_SelectArterials.
function radiobutton_SelectArterials_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton_SelectArterials
set(handles.radiobutton_SelectArterials,'Value',1);
set(handles.radiobutton_SelectVeins,'value',0);
set(handles.radiobutton_SelectCapillaries,'value',0);

% --- Executes on button press in radiobutton_SelectVeins.
function radiobutton_SelectVeins_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton_SelectVeins
set(handles.radiobutton_SelectArterials,'value',0);
set(handles.radiobutton_SelectVeins,'Value',1);
set(handles.radiobutton_SelectCapillaries,'value',0);

% --- Executes on button press in radiobutton_SelectCapillaries.
function radiobutton_SelectCapillaries_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton_SelectCapillaries
set(handles.radiobutton_SelectArterials,'value',0);
set(handles.radiobutton_SelectVeins,'value',0);
set(handles.radiobutton_SelectCapillaries,'Value',1);



function edit_TopOfRegion_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TopOfRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TopOfRegion as text
%        str2double(get(hObject,'String')) returns contents of edit_TopOfRegion as a double
global im;

foo = get(handles.edit_TopOfRegion,'String');
num = str2num(foo);
if num < 0, num = 0; end; 
if num > str2num(get(handles.edit_SelectBottom,'String')), num = 0; end; % if num > bottom border
if num > size(im.It,3)*im.ht(3), % if num > max depth in stack
    num = size(im.It,3)*im.ht(3);
end;
set(handles.edit_TopOfRegion,'String',sprintf('%.1f',num));


% --- Executes during object creation, after setting all properties.
function edit_TopOfRegion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TopOfRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    set(hObject,'String','0.0');
end


% --- Executes on button press in pushbutton_ForctTopSurface.
function pushbutton_ForctTopSurface_Callback(hObject, eventdata, handles)
set(handles.edit_TopOfRegion,'String','0.0');



function edit_SelectBottom_Callback(hObject, eventdata, handles)
global im;

foo = get(handles.edit_SelectBottom,'String');
num = str2num(foo);
if num < 0, num = 0; end;
if num < str2num(get(handles.edit_TopOfRegion,'String')), num = size(im.It,3)*im.ht(3); end; % if num > bottom border
if num > size(im.It,3)*im.ht(3), % if num > max depth in stack
    num = size(im.It,3)*im.ht(3);
end;
set(handles.edit_SelectBottom,'String',sprintf('%.1f',num));

% --- Executes during object creation, after setting all properties.
function edit_SelectBottom_CreateFcn(hObject, eventdata, handles)
global im;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    num = size(im.It,3)*im.ht(3);
    set(hObject,'String',sprintf('%.1f',num));
end


% --- Executes on button press in pushbutton_ForceBottomSelection.
function pushbutton_ForceBottomSelection_Callback(hObject, eventdata, handles)
global im;
num = size(im.It,3)*im.ht(3);
set(handles.edit_SelectBottom,'String',sprintf('%.1f',num));




% --- Executes on button press in pushbutton_ClcNodeGrps.
function pushbutton_ClcNodeGrps_Callback(hObject, eventdata, handles)
global im;

% nB is 1D array, #elements is nNodes, values are number of edges from each
% node
nNodes = size(im.im2.nodePos,1);
nB=zeros(1,nNodes);
for ii=1:nNodes
    nB(ii)=length(find(im.im2.nodeEdges(:,1)==ii | im.im2.nodeEdges(:,2)==ii));
end
im.im2.nB = nB;
im.im2 = nodeGrps(im.im2);

% ?
nSeg = length(im.im2.segDiam);
im.im2.segPos = squeeze(mean(reshape(im.im2.nodePos(im.im2.segEndNodes,:),[2 nSeg 3]),1));

% find segments with length less than a fraction of the diameter
% and less than a certain depth
zThresh = size(im.It,3)*im.ht(3); %400;
Len2DiamRatio = 2;
lst2 = find((im.im2.segLen./max(im.im2.segDiam,1))<Len2DiamRatio);
lst4=im.im2.segEndNodes(lst2,:);
foo = find(im.im2.nodePos(lst4(:,1),3)<zThresh & im.im2.nodePos(lst4(:,2),3)<zThresh);
lstBad = lst2(foo);
nSeg = length(im.im2.segLen);
nSegBad = length(lstBad);

lst2 = find((im.im2.segLen./max(im.im2.segDiam,1))>=Len2DiamRatio);
lst4=im.im2.segEndNodes(lst2,:);
foo = find(im.im2.nodePos(lst4(:,1),3)<zThresh & im.im2.nodePos(lst4(:,2),3)<zThresh);
lstGood = lst2(foo);

% remove loops from good list
nIdx = im.im2.segEndNodes(lstGood,:);
rsep = sqrt(sum( (im.im2.nodePos(nIdx(:,1),:) - im.im2.nodePos(nIdx(:,2),:)).^2, 2));
% search for instances of loops (rsep = 0)
foo = find(rsep == 0);  % bad ones are lstGood(foo)
if ~isempty(foo),
    lstBad = [lstBad lstGood(foo)];
    lstGood(foo) = [];
end;

nSeg = length(im.im2.segLen);
nSegBad = length(lstBad);
nSegGood = length(lstGood);
[nSeg nSegBad+nSegGood nSegBad nSegGood]

im.im2.lstBad = lstBad;
im.im2.lstGood = lstGood;

% histogram depth of "bad" segments
nIdx = im.im2.segEndNodes(lstBad,:);
figure(1)                 
hist( im.im2.nodePos(nIdx(:),3), [zThresh/20:zThresh/10:zThresh] )
title([ 'Depth of "bad" segments (Lenght / Diam < ' num2str(Len2DiamRatio) ')'] );

% Assign each segment AVC
% already done. im.im2.segVesType 1D array
% 2 - capillary
% 3 - vein
% 1 - artery
% 0 - not assigned

% Assign each segment a depth
im.im2.segDepth_um = zeros(1,nSeg);
for i = 1:nSeg,
    pos = im.im2.segPos(i,:); % pos is (x,y,z) oriented. Zprofile needs (y,x) input.
    %im.Zprofile = @(c,y,x) (c(1).*y.*y + c(2).*x.*x + c(3).*y.*x + c(4).*y + c(5).*x +c(6) );
    im.im2.segDepth_um(i) = im.ht(3) * abs ( pos(3) - im.Zprofile(im.cAll, pos(2), pos(1)) ); % in um
    % depth is in um!
end;

% Assign each segment a volume
im.im2.segVol_um3 = zeros(1,nSeg);
for i=1:nSeg,
    im.im2.segVol_um3(i) = (im.im2.segDiam(i))^2 * (pi/4) * im.im2.segLen_um(i);
end;






% --- Executes on button press in pushbutton_CalcTortuosity.
function pushbutton_CalcTortuosity_Callback(hObject, eventdata, handles)
% use segLen_um and rsep must also be calculated in um. voxels are not
% cubes!
global im;

hxy = im.ht(1);
hz = im.ht(3);

nIdx = im.im2.segEndNodes(im.im2.lstGood,:);
rsep = sqrt(sum( (im.im2.nodePos_um(nIdx(:,1),:) - im.im2.nodePos_um(nIdx(:,2),:)).^2, 2));
tortGood = im.im2.segLen_um(im.im2.lstGood)'./rsep;

nIdx = im.im2.segEndNodes(im.im2.lstBad,:);
rsep = sqrt(sum( (im.im2.nodePos_um(nIdx(:,1),:) - im.im2.nodePos_um(nIdx(:,2),:)).^2, 2));
tortBad = im.im2.segLen_um(im.im2.lstBad)'./rsep;

tort = zeros(length(im.im2.segLen_um),1);
tort(im.im2.lstGood) = tortGood;
tort(im.im2.lstBad) = tortBad;
im.im2.tort = tort;

% Find Tortuosity based on AVC and depth region

% extract Z-region boundaries and AVC flag
topZ_um = str2num(get(handles.edit_TopOfRegion,'string'));
bottomZ_um = str2num(get(handles.edit_SelectBottom,'string'));
if get(handles.radiobutton_SelectArterials,'Value')==1,
    AVCflag = 1; % 1 = Artery
    AVCflagStr = 'Arterials';
elseif get(handles.radiobutton_SelectVeins,'Value')==1,
    AVCflag = 3; % 3 = Veins
    AVCflagStr = 'Venules';
else
    AVCflag = 2; % 2 = capillary
    AVCflagStr = 'Capillaries';
end;

%im.im2.segDepth_um - list of segment depths
%im.im2.segVesType - list of segment-vessel types
%im.im2.lstGood - list of good segments

lstDepth = find(im.im2.segDepth_um >= topZ_um  &  im.im2.segDepth_um <= bottomZ_um);
lstVesType = find(im.im2.segVesType == AVCflag);
% find segments inside XY_ROI
if isfield(im,'XY_ROIs'),
    if im.XY_ROIs.flagUseROIs, % use XY ROI...
        bw = im.XY_ROIs.ROI; % mask where 1 is at good pixels and 0 at rejected pixels ariund large vessels
        foo = round(squeeze(im.im2.segPos(:,1:2)));
        cntr = 0;
        lstXYroi = [];
        for iii = 1:size(foo,1),
            if bw(foo(iii,1),foo(iii,2)), % CHECK IF X,Y ARE IN CORRECT ORDER HERE!!!
                cntr=cntr+1;
                lstXYroi(cntr)=iii;
            end;
        end;
        lstResult = intersect(intersect(lstXYroi,intersect(lstDepth,lstVesType)),im.im2.lstGood);
    else
        lstResult = intersect(intersect(lstDepth,lstVesType),im.im2.lstGood); % list of segments with selected properties
    end;
else
    lstResult = intersect(intersect(lstDepth,lstVesType),im.im2.lstGood); % list of segments with selected properties
end;





% get number of result segments , mean and error of tortuosity, total
% volume between selected topZ and bottomZ

% calc volume between Ztop and Zbottom
[ny, nx, nz] = size(im.It); % # voxels

Zsurf = im.Zsurf;
Zsurf = Zsurf(:); % Zsurf is in voxels
if isfield(im,'XY_ROIs'),
    if im.XY_ROIs.flagUseROIs, % use XY ROI...
        Zsurf = Zsurf(find(bw(:)>0.5)); %CHECK IF ZSURF(X,Y) or (Y,X) AND WHAT I USE HERE IS CORRECT!!!
    end;
end;
Z2 = Zsurf + bottomZ_um/hz;
Z2(find(Z2>nz))=nz;
Z1 = Zsurf + topZ_um/hz;
Z1(find(Z1>nz))=nz;
vol_um3 = sum(Z2-Z1)*hz*hxy^2;
vol_mm3 = vol_um3/1e9;

% get remaining statistics
prntStr(1) = {['Segment type     : ' AVCflagStr]};
prntStr(2) = {['Volume top       : ' num2str(topZ_um) ' um']};
prntStr(3) = {['Volume bottom    : ' num2str(bottomZ_um) ' um']};
prntStr(4) = {['Volume size      : ' num2str(vol_um3) ' um^3']};

numResSegments = length(lstResult);
prntStr(5) = {['Number of segments: ' num2str(numResSegments)]};
prntStr(6) = {'  '};
if numResSegments > 0,
 
    meanTort = mean(im.im2.tort(lstResult));
    stdErrTort = std(im.im2.tort(lstResult))/sqrt(numResSegments);
    prntStr(7) = {['Mean tortuosity   : ' num2str(meanTort)]};
    prntStr(8) = {['StdErr tortuosity : ' num2str(stdErrTort)]};
    prntStr(9) = {'  '};
    
    meanSegLength_um = mean(im.im2.segLen_um(lstResult));
    stdErrSegLength_um = std(im.im2.segLen_um(lstResult))/sqrt(numResSegments);
    totSegLength_um = sum(im.im2.segLen_um(lstResult));
    fracSegLength_um_mm3 = totSegLength_um / vol_mm3;
    prntStr(10)  = {['Mean Segment Length     : ' num2str(meanSegLength_um) ' um']};
    prntStr(11)  = {['StdErr Segment Length   : ' num2str(stdErrSegLength_um) ' um']};
    prntStr(12)  = {['Total Segment Length    : ' num2str(totSegLength_um) ' um']};
    prntStr(13)  = {['Fract Segment Length    : ' num2str(fracSegLength_um_mm3) ' um/mm^3']};
    prntStr(14)  = {'  '};
    
    meanSegVol_um3 = mean(im.im2.segVol_um3(lstResult));
    stdErrSegVol_um3 = std(im.im2.segVol_um3(lstResult))/sqrt(numResSegments);
    totSegVol_um3 = sum(im.im2.segVol_um3(lstResult));
    fracSegVol = totSegVol_um3 / vol_um3;
    fracSegVolErr = stdErrSegVol_um3 / vol_um3;
    
    prntStr(15)  = {['Mean Segment Volume     : ' num2str(meanSegVol_um3) ' um^3']};
    prntStr(16)  = {['stdErr Segment Volume   : ' num2str(stdErrSegVol_um3) ' um^3']};
    prntStr(17)  = {['Sum of Segments Volume  : ' num2str(totSegVol_um3) ' um^3']};
    prntStr(18)  = {['Fract. Segments Volume  : ' num2str(fracSegVol) ]};
    prntStr(19)  = {['Fract. Segments Volume Err. : ' num2str(fracSegVolErr) ]};
    
    resFigH = figure; 
    set(gca,'Visible','off');
    text(0.1,0.6,prntStr,'FontSize',12);
    
end;


% --- Executes on button press in radiobutton_UseROI.
function radiobutton_UseROI_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_UseROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_UseROI
global im;
if im.flagView ~= 3,
    msgbox('View must be XY','','error');
    foo = get(hObject,'Value');
    set(hObject,'Value',1-foo);
else
    if get(hObject,'Value'), % radiobutton set to 1 => use ROIs
        
        if isfield(im,'XY_ROIs'), % ask to select ROIs or to use existing ones
            choice = questdlg('Select new rejection ROI(s) or use existing one(s)?','New ROI(s)?','Yes, select new','No, use existing','No, use existing');
            if strcmp(choice,'Yes, select new'),
                flagResetROIs = 1;
            else
                flagResetROIs = 0;
                im.XY_ROIs.flagUseROIs = 1;
            end;
        else
            flagResetROIs = 1;
        end;
        
        if flagResetROIs, % select ROIs and set flagUseROIs to 1
            numROIs = 0;
            im.XY_ROIs.ROI = [];
            answer = inputdlg({'How many ROIs do you want?'},'Number ROIs',1,{'1'} );
            numROIs = round(str2num(answer{1}));
            if numROIs < 1,
                numROIs = 1; % force at leat one ROI ifsomeone inputs 0 or negative number
            end;
            for ii = 1:numROIs,
                [bw, xCoor yCoor] = roipoly;
                hold on;
                line(round(xCoor),round(yCoor));
                hold off;
                im.XY_ROIs.ROI{ii} = bw;
            end;
            bw = im.XY_ROIs.ROI{1};
            for ii = 1:numROIs,
                bw = bw | im.XY_ROIs.ROI{ii};
            end;
            im.XY_ROIs.ROI = [];
            im.XY_ROIs.ROI = 1-double(bw); % this mask will be 0 when pixel is in rejected ROI and 1 if we need this pixel
            im.XY_ROIs.flagUseROIs = 1;
        end;
        
    else % radiobutton set to 0
        im.XY_ROIs.flagUseROIs = 0;
    end;
    
end;
            
  


% --- Executes on button press in pushbutton_PlotROI.
function pushbutton_PlotROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_PlotROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global im;
if im.flagView ~= 3,
    msgbox('View must be XY','','error');
else
    if isfield(im,'XY_ROIs'), % if ROI was selected before...
        if isfield(im.XY_ROIs,'ROI'), % if ROI was selected before
            foo=xor(im.XY_ROIs.ROI, imerode(im.XY_ROIs.ROI,strel('disk',1)));
            [row, col]=find(foo>0.5);
            hold on; 
            plot(col.*im.ht(1),row.*im.ht(1),'.');
            hold off;
        else
            msgbox('ROI was never selected','','error');
        end;
    else
        msgbox('ROI was never selected','','error');
    end;
end;
        
    
