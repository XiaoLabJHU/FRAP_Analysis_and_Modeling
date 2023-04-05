function varargout = FRAP_ring_circleFH(varargin)
% FRAP_RING_CIRCLEFH MATLAB code for FRAP_ring_circleFH.fig
%      FRAP_RING_CIRCLEFH, by itself, creates a new FRAP_RING_CIRCLEFH or raises the existing
%      singleton*.
%
%      H = FRAP_RING_CIRCLEFH returns the handle to a new FRAP_RING_CIRCLEFH or the handle to
%      the existing singleton*.
%
%      FRAP_RING_CIRCLEFH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FRAP_RING_CIRCLEFH.M with the given input arguments.
%
%      FRAP_RING_CIRCLEFH('Property','Value',...) creates a new FRAP_RING_CIRCLEFH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FRAP_ring_circleFH_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FRAP_ring_circleFH_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FRAP_ring_circleFH

% Last Modified by GUIDE v2.5 27-Feb-2023 13:04:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FRAP_ring_circleFH_OpeningFcn, ...
                   'gui_OutputFcn',  @FRAP_ring_circleFH_OutputFcn, ...
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


% --- Executes just before FRAP_ring_circleFH is made visible.
function FRAP_ring_circleFH_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FRAP_ring_circleFH (see VARARGIN)

% Choose default command line output for FRAP_ring_circleFH
handles.output = hObject;

%---set up your parameters for ROI cropping----%
% the radius of circle in pixel

%Old radius
%handles.R = 5;

%New radius
handles.R = 3;

% the scaling factor to expand the image, has to be a integer bigger than 1
% the bigger the more accurate, but slower
% The algorithm of expanding is 'bicubic' now
handles.Scale = 12;
handles.X = 49;
handles.Y = 50;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FRAP_ring_circleFH wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FRAP_ring_circleFH_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in BFinput.
function BFinput_Callback(hObject, eventdata, handles)
% hObject    handle to BFinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename pathname] = uigetfile('.tif','Input Bright field image.');
handles.pathname = pathname;
handles.BF_filename = filename;
handles.BFimage = imread([pathname filename],'tif');

handles.LowI=min(min(handles.BFimage));
handles.HighI=max(max(handles.BFimage));

axes(handles.BF_image);
hold off
imshow(handles.BFimage,[]);

set(handles.filename,'String',[pathname filename]);
guidata(hObject, handles);

% --- Executes on button press in FL_input.
function FL_input_Callback(hObject, eventdata, handles)
% hObject    handle to FL_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname]=uigetfile('*.tif','select the fluorescence image');
InfoFL=imfinfo([pathname filename],'tif');
FrameNum=length(InfoFL); % the number of frames the image has
handles.FrameFL=FrameNum;
for idxF=1:FrameNum
    ImFL(:,:,idxF)=imread([pathname filename],'tif',idxF);
end
% show this image in the figure
axes(handles.FL_image);
hold off
imshow(ImFL(:,:,1),[]);
handles.curr_image=ImFL(:,:,1);
handles.LowF=min(min(ImFL(:,:,1)));
handles.HighF=max(max(ImFL(:,:,1)));
% axes(handles.hist_image);
% hold off
% plothist(handles.curr_image);
handles.curr_frame=1;
% save variable
handles.FLimage=ImFL;
handles.FLpath=pathname;
handles.FLfile=filename;
set(handles.FL_slider,'SliderStep',[1/(FrameNum-1) 0]);
%set(handles.FL_slider,'Min',1/(FrameNum-1));
set(handles.FL_slider, 'Min', 0);
set(handles.FL_slider,'Value',1);
% handles.imagetype='radioFL';
% handles.output = hObject;
% update the data
guidata(hObject, handles);

% --- Executes on slider movement.
function FL_slider_Callback(hObject, eventdata, handles)
% hObject    handle to FL_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% get the current value of the max and min
    frame_p = 1-get(handles.FL_slider,'Value');
    Frame = floor(frame_p*handles.FrameFL)+1;
    axes(handles.FL_image);
    hold off
    imshow(handles.FLimage(:,:,Frame),[handles.LowF handles.HighF]);
    handles.curr_image = handles.FLimage(:,:,Frame);
    % axes(handles.hist_image);
    % hold off
    % plothist(handles.curr_image);
    handles.curr_frame = Frame;
    guidata(hObject, handles);
    if isfield(handles,'ROIring')
        ROI = handles.ROIring;
        Angles = [0:0.1:2*pi];
        Vx = ROI(3) * cos(Angles) + ROI(1);
        Vy = ROI(3) * sin(Angles) + ROI(2);
        axes(handles.FL_image);
        hold on
        plot(Vx,Vy,'-y');
    end

% --- Executes during object creation, after setting all properties.
function FL_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FL_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in ROI_ring.
function ROI_ring_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_ring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load the parameters
R = handles.R;
Scale = handles.Scale;

% calculate and give the meshgrid of the image
[SY,SX] = size(handles.FLimage(:,:,1));
[MeshX MeshY] = meshgrid([1:SX*Scale],[1:SY*Scale]);

% select center of the circle ROI
% X = handles.X;
% Y = handles.Y;
[X,Y] = ginput(1);

% plot the region on the image
Angles = [0:0.1:2*pi];
Vx = R * cos(Angles) + X;
Vy = R * sin(Angles) + Y;
ROI=[X,Y,R];

% ROI_FL=handles.FLimage(ROI(2,2):ROI(3,2),ROI(1,1):ROI(2,1),:);
axes(handles.FL_image);
hold on
plot(Vx,Vy,'-y');
axes(handles.BF_image);
hold on
plot(Vx,Vy,'-y');

% calculate the meanintensity in the ROI
XCenter = round(X * Scale - Scale/2);
YCenter = round(Y * Scale - Scale/2);
DistanceM = sqrt((MeshX - XCenter).^2 + (MeshY - YCenter).^2);
IndexM = find(DistanceM <= R * Scale);
for idx = 1 : handles.FrameFL
    FLimC = handles.FLimage(:,:,idx); % iterate the image
    FLimC_rescale = imresize(FLimC,Scale); % expand the image
    RingI(:,idx)=sum(sum(FLimC_rescale(IndexM)))/Scale^2;
end
handles.RingI = RingI;
handles.AreaRing = pi*R^2;
handles.ROIRing = ROI;

% update data
guidata(hObject, handles);

% --- Executes on button press in ROI_cell.
function ROI_cell_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_cell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% select rectangle ROI
[X,Y] = ginput(2);
% get the vertex of the rectangle
Vx=[X(1);X(2);X(2);X(1);X(1)];
Vy=[Y(1);Y(1);Y(2);Y(2);Y(1)];
ROI=round(cat(2,Vx,Vy));
ROI_FL=handles.FLimage(ROI(2,2):ROI(3,2),ROI(1,1):ROI(2,1),:);
axes(handles.FL_image);
hold on
plot(Vx,Vy,'-w');
axes(handles.BF_image);
hold on
plot(Vx,Vy,'-w');
response=questdlg('Use the current ROI?', ...
    'ROI selection', 'Yes' , 'No','Yes');
if  strcmp(response,'Yes')
    CellI(:,1)=sum(sum(ROI_FL,1),2);
    handles.CellI = CellI;
    handles.AreaCell = (abs(ROI(2,2)-ROI(3,2))+1).*(abs(ROI(2,1)-ROI(1,1))+1);
    handles.ROICell = ROI;
    handles.ROI_FL = ROI_FL;
end

% update data
guidata(hObject, handles);

% --- Executes on button press in ROI_background.
function ROI_background_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[X,Y] = ginput(2);
% get the vertex of the rectangle
Vx=[X(1);X(2);X(2);X(1);X(1)];
Vy=[Y(1);Y(1);Y(2);Y(2);Y(1)];
ROI=round(cat(2,Vx,Vy));
ROI_FL=handles.FLimage(ROI(2,2):ROI(3,2),ROI(1,1):ROI(2,1),:);
axes(handles.BF_image);
hold on
plot(Vx,Vy,'-r');
axes(handles.FL_image);
hold on
plot(Vx,Vy,'-r');
response=questdlg('Use the current ROI?', ...
    'ROI selection', 'Yes' , 'No','Yes');
if  strcmp(response,'Yes')
    BackI(:,1)=sum(sum(ROI_FL,1),2);
    handles.BackI = BackI;
    handles.AreaBack = (abs(ROI(2,2)-ROI(3,2))+1).*(abs(ROI(2,1)-ROI(1,1))+1);
    handles.ROIBack = ROI;
end

% update data
guidata(hObject, handles);


% --- Executes on button press in BF_rotate.
function BF_rotate_Callback(hObject, eventdata, handles)
% hObject    handle to BF_rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[X,Y]=ginput(2);
% get the angle
Angle=180*(atan((Y(2)-Y(1))/(X(2)-X(1)))+pi/2)/pi;

%This line makes the cell horizontal instead of vertical
Angle = Angle + 90;

BFimage=imrotate(handles.BFimage,Angle,'bicubic','crop');

axes(handles.BF_image);
hold off
imshow(BFimage,[handles.LowI handles.HighI]);
response=questdlg('The rotation Okay?', ...
    'Rotate the Cell', 'Yes' , 'No','Yes');
if strcmp(response,'Yes');
    newXY = rotate_point([handles.X handles.Y],Angle,[50.5,50.5]);
    handles.X = newXY(1);
    handles.Y = newXY(2);
    handles.BFimage=BFimage;
    for idxf=1:size(handles.FLimage,3)
        FLimage(:,:,idxf)=imrotate(handles.FLimage(:,:,idxf),Angle,'bicubic','crop');
    end
    handles.FLimage=FLimage;
    axes(handles.FL_image);
    hold off
    imshow(FLimage(:,:,1),[handles.LowF handles.HighF]);
    set(handles.FL_slider,'Value',1);
else
    axes(handles.BF_image);
    hold off
    imshow(handles.BFimage,[handles.LowI handles.HighI]);
end
guidata(hObject,handles);

% --- Executes on button press in line_long.
function line_long_Callback(hObject, eventdata, handles)
% hObject    handle to line_long (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[X_long Y_long] = ginput(2);
dist_long = sum(sqrt((X_long-Y_long).^2));
n = round(dist_long*5);
handles.FLimageav = mean(handles.FLimage,3);
[cx,cy,c] = improfile(handles.BFimage,X_long,Y_long,n,'bilinear');
handles.BFprof_long = [cx,cy,c];

[cx1,cy1,c1] = improfile(handles.FLimageav,X_long,Y_long,n,'bilinear');
handles.FLprof_long = [cx1,cy1,c1];
handles.index_long = [X_long, Y_long];

axes(handles.BF_image);
hold on
plot(X_long, Y_long,'y');

axes(handles.FL_image);
hold on
plot(X_long, Y_long,'y');

% figure
% plot(cx,c)
% figure
% plot(cx1,c1)
guidata(hObject,handles);


% --- Executes on button press in line_short.
function line_short_Callback(hObject, eventdata, handles)
% hObject    handle to line_short (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[X_short Y_short] = ginput(2);
dist_short = sum(sqrt((X_short-Y_short).^2));
n = round(dist_short*5);
handles.FLimageav = mean(handles.FLimage,3);
[cx,cy,c] = improfile(handles.BFimage,X_short,Y_short,n,'bilinear');
handles.BFprof_short = [cx,cy,c];

[cx1,cy1,c1] = improfile(handles.FLimageav,X_short,Y_short,n,'bilinear');
handles.FLprof_short = [cx1,cy1,c1];
handles.index_short = [X_short, Y_short];

axes(handles.BF_image);
hold on
plot(X_short, Y_short,'y');

axes(handles.FL_image);
hold on
plot(X_short, Y_short,'y');

guidata(hObject,handles);


% --- Executes on button press in Data_save.
function Data_save_Callback(hObject, eventdata, handles)
% hObject    handle to Data_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'pathname')
    errordlg('You need to input a BF image!','Input Error');
    return;
end

if ~isfield(handles,'FrameFL')
    errordlg('You need to input a fluorescence time lapse!','Input Error');
    return;
end

if isfield(handles,'FLprof_long')
    FRAP_data.FLprof_long = handles.FLprof_long;
    FRAP_data.BFprof_long = handles.BFprof_long;
    FRAP_data.index_long = handles.index_long;
end

if isfield(handles,'FLprof_short')
    FRAP_data.FLprof_short = handles.FLprof_short;
    FRAP_data.BFprof_short = handles.BFprof_short;
    FRAP_data.index_short = handles.index_short;
end

if ~isfield(handles,'RingI')
    errordlg('You need to choose the ROI of the ring!','Input Error');
    return;
else
FRAP_data.RingI = handles.RingI;
FRAP_data.AreaRing = handles.AreaRing;
FRAP_data.ROIRing = handles.ROIRing;
end

if ~isfield(handles,'CellI')
    errordlg('You need to choose the ROI of the cell!','Input Error');
    return;
else
FRAP_data.CellI = handles.CellI;
FRAP_data.AreaCell = handles.AreaCell;
FRAP_data.ROICell = handles.ROICell;
FRAP_data.ROI_FL = handles.ROI_FL;
FRAP_data.CellFL = handles.ROI_FL;

end

if ~isfield(handles,'BackI')
    errordlg('You need to choose the ROI of the background!','Input Error');
    return;
else
FRAP_data.BackI = handles.BackI;
FRAP_data.AreaBack = handles.AreaBack;
FRAP_data.ROIBack = handles.ROIBack;
end

%filename = [handles.BF_filename(1:find(handles.BF_filename=='_')-1) '_FRAP_data.mat'];
%save([handles.pathname 'radius5_FRAP_data'],'FRAP_data');

%file name is built from BF filename and radius
filename = [handles.BF_filename(1:find(handles.BF_filename=='B')-1) 'r' num2str(handles.R) '_data.mat'];
save([handles.pathname filename],'FRAP_data');



% --- Executes on button press in FL_delete.
function FL_delete_Callback(hObject, eventdata, handles)
% hObject    handle to FL_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
response=questdlg('Delete the current fluorescence image?', ...
    'Deletion', 'Yes' , 'No','Yes');
 FrameNum=handles.curr_frame;
 if strcmp(response,'Yes')
     handles.FLimage(:,:,FrameNum)=[];
     handles.curr_frame = 1;
 end
axes(handles.FL_image);
hold off
imshow(handles.FLimage(:,:,1),[]);
set(handles.FL_slider,'Value',1);
handles.FrameFL = handles.FrameFL - 1;
% update the data
guidata(hObject, handles);
