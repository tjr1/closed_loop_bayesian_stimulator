function varargout = on_demand_BPO_v1(varargin)
% ON_DEMAND_BPO_V1 MATLAB code for on_demand_BPO_v1.fig
%      ON_DEMAND_BPO_V1, by itself, creates a new ON_DEMAND_BPO_V1 or raises the existing
%      singleton*.
%
%      H = ON_DEMAND_BPO_V1 returns the handle to a new ON_DEMAND_BPO_V1 or the handle to
%      the existing singleton*.
%
%      ON_DEMAND_BPO_V1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ON_DEMAND_BPO_V1.M with the given input arguments.
%
%      ON_DEMAND_BPO_V1('Property','Value',...) creates a new ON_DEMAND_BPO_V1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before on_demand_BPO_v1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to on_demand_BPO_v1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help on_demand_BPO_v1

% Last Modified by GUIDE v2.5 19-Apr-2018 10:38:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @on_demand_BPO_v1_OpeningFcn, ...
                   'gui_OutputFcn',  @on_demand_BPO_v1_OutputFcn, ...
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


% --- Executes just before on_demand_BPO_v1 is made visible.
function on_demand_BPO_v1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to on_demand_BPO_v1 (see VARARGIN)

% Choose default command line output for on_demand_BPO_v1

clc

daqreset
imaqreset

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes on_demand_BPO_v1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = on_demand_BPO_v1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ET_load_path_Callback(hObject, eventdata, handles)
% hObject    handle to ET_load_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_load_path as text
%        str2double(get(hObject,'String')) returns contents of ET_load_path as a double


% --- Executes during object creation, after setting all properties.
function ET_load_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_load_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PB_LoadSettings.
function PB_LoadSettings_Callback(hObject, eventdata, handles)
% hObject    handle to PB_LoadSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load(get(handles.ET_load_path,'String'));

try
    set(handles.ET_n_ch_in,'String',struc.ET_n_ch_in);
end
try
    set(handles.ET_n_ch_out,'String',struc.ET_n_ch_out);
end
try
    set(handles.ET_seizure_detection_channels,'String',struc.ET_seizure_detection_channels);
end
try
    set(handles.ET_n_cams,'String',struc.ET_n_cams);
end
try
    set(handles.ET_fs,'String',struc.ET_fs);
end
try
    set(handles.ET_MonoOrBiPhasic,'String',struc.ET_MonoOrBiPhasic);
end
try
    set(handles.ET_AmplitudeRange,'String',struc.ET_AmplitudeRange);
end
try
    set(handles.ET_PulseWidthRange,'String',struc.ET_PulseWidthRange);
end
try
    set(handles.ET_FrequencyRange,'String',struc.ET_FrequencyRange);
end
try
    set(handles.ET_TrainDuration,'String',struc.ET_TrainDuration);
end
try
    set(handles.ET_SaveFolder,'String',struc.ET_SaveFolder);
end
try
    set(handles.ET_SaveName,'String',struc.ET_SaveName);
end
try
    set(handles.ET_save_path,'String',struc.ET_save_path);
end

guidata(hObject,handles);



function ET_save_path_Callback(hObject, eventdata, handles)
% hObject    handle to ET_save_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_save_path as text
%        str2double(get(hObject,'String')) returns contents of ET_save_path as a double


% --- Executes during object creation, after setting all properties.
function ET_save_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_save_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PB_SaveSettings.
function PB_SaveSettings_Callback(hObject, eventdata, handles)
% hObject    handle to PB_SaveSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.ET_load_path,'String',get(handles.ET_save_path,'String'));

struc.ET_n_ch_in=get(handles.ET_n_ch_in,'String');
struc.ET_n_ch_out=get(handles.ET_n_ch_out,'String');
struc.ET_seizure_detection_channels=get(handles.ET_seizure_detection_channels,'String');
struc.ET_n_cams=get(handles.ET_n_cams,'String');
struc.ET_fs=get(handles.ET_fs,'String');
struc.ET_MonoOrBiPhasic=get(handles.ET_MonoOrBiPhasic,'String');
struc.ET_AmplitudeRange=get(handles.ET_AmplitudeRange,'String');
struc.ET_PulseWidthRange=get(handles.ET_PulseWidthRange,'String');
struc.ET_FrequencyRange=get(handles.ET_FrequencyRange,'String');
struc.ET_TrainDuration=get(handles.ET_TrainDuration,'String');
struc.ET_SaveFolder=get(handles.ET_SaveFolder,'String');
struc.ET_SaveName=get(handles.ET_SaveName,'String');
struc.ET_save_path=get(handles.ET_save_path,'String');
struc.ET_load_path=get(handles.ET_load_path,'String');

save(struc.ET_save_path, 'struc');

guidata(hObject,handles);


% --- Executes on button press in PB_Go.
function PB_Go_Callback(hObject, eventdata, handles)
% hObject    handle to PB_Go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

daqreset
imaqreset

%% create save folder

save_folder_base = get(handles.ET_SaveFolder,'String');
save_name = get(handles.ET_SaveName,'String');
date_stamp = datestr(datetime,'mm-dd-yyyy_HH-MM-SS');
save_folder = [save_folder_base '\' date_stamp];

if not(exist(save_folder,'dir'))
    mkdir(save_folder)
else
    warning('folder already exists')
end

%% setup video cameras
n_cams = str2num(get(handles.ET_n_cams,'String'));
[vid, src] = setup_video(n_cams);

%% setup daq
devices = daq.getDevices;

for i = 1:length(devices)
    if strcmp(devices(i).Model,'PCI-6251')
        PCI_6251_dev = devices(i);
    end
    if strcmp(devices(i).Model,'USB-6229 (BNC)')
        USB_6229_dev = devices(i);
    end
end

s = daq.createSession('ni');

global fs
fs = str2num(get(handles.ET_fs,'String'));
dt = 1/fs;

s.Rate = fs;
s.IsContinuous = true;

n_ch_out = str2num(get(handles.ET_n_ch_out,'String')); % hard coded number of outputs

seizure_detection_ch_vec = str2num(get(handles.ET_seizure_detection_channels,'String'));% detect seizures on ACH0 and ACH2, these are BNC-2090 numbers

addAnalogOutputChannel(s, PCI_6251_dev.ID, (1:n_ch_out)-1, 'Voltage');

global stim_flag
stim_flag = zeros(1, n_ch_out,'logical');

n_ch_in = str2num(get(handles.ET_n_ch_in,'String')); % hard coded 1 input channel
for i_ch = 1:n_ch_in % add ephys recording channels
    ch(i_ch) = addAnalogInputChannel(s,PCI_6251_dev.ID, i_ch-1, 'Voltage');
end

%% initialize matfile
data=zeros(1,1+n_ch_in); % column 1 is time stamps, next n_ch_in columns are the input channels, and the last n_ch_out columns are the output channels

save_mat_path = [save_folder '\' save_name '.mat']
save(save_mat_path,'data','fs','-v7.3','-nocompression')

for i_cam = 1:n_cams
    eval(['time_cam_' num2str(i_cam) ' = 0;'])
    eval(['meta_time_cam_' num2str(i_cam) ' = zeros(1,9);'])
end

for i_cam = 1:n_cams
    save(save_mat_path,['time_cam_' num2str(i_cam)],'-nocompression','-append')
    save(save_mat_path,['meta_time_cam_' num2str(i_cam)],'-nocompression','-append')
end

clear data
mf = matfile(save_mat_path,'Writable',true);

%% add listeners
lh1 = addlistener(s,'DataAvailable', @(src, event) Process_Plot_Save(src,event,mf, handles.A_MainPlot,n_cams,vid,seizure_detection_ch_vec) );
lh2 = addlistener(s,'DataRequired', @Generate_Stim_Vec);

s.NotifyWhenDataAvailableExceeds = fix(fs/2); % input buffer treshold, hard coded
s.NotifyWhenScansQueuedBelow = fix(fs/3); % output buffer threshold, hard coded
% s.NotifyWhenDataAvailableExceeds = fix(fs/5); % input buffer treshold, hard coded
% s.NotifyWhenScansQueuedBelow = fix(fs/6); % output buffer threshold, hard coded


handles.s = s;

%% buffer 2 seconds of zeros to the daq output to start with
stim_vec_init = zeros(fix(2*fs),n_ch_out);
queueOutputData(s, [stim_vec_init]);

%% setup video save_file and start the video
for i_cam = 1:n_cams
    logfile(i_cam) = VideoWriter([save_folder '\' save_name '_cam_' num2str(i_cam) '.mp4'], 'MPEG-4');
    logfile(i_cam).FrameRate = 30;
    vid(i_cam).DiskLogger = logfile(i_cam);
    start(vid(i_cam))
end

%% trigger cameras to start
pause(.1)
disp('starting cam')
for i_cam = 1:n_cams
    trigger(vid(i_cam))
end
pause(0.6)

%% start daq
disp('starting daq')
daq_start_now = now;
daq_start_datetime = datetime;

startBackground(s); % start the daq

%% record daq start time info
[Y,M,D,H,MN,S] = datevec(daq_start_datetime);
daq_date_vector = [Y,M,D,H,MN,S];
save(save_mat_path,'daq_date_vector','daq_start_now','daq_start_datetime','-nocompression','-append')

handles.vid = vid;
handles.n_cams = n_cams;

guidata(hObject,handles);


% --- Executes on button press in PB_Stop.
function PB_Stop_Callback(hObject, eventdata, handles)
% hObject    handle to PB_Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s = handles.s;
s.stop()

vid = handles.vid;
n_cams = handles.n_cams;

%% stop video
for i_cam = 1:n_cams
    stop(vid(i_cam));
end

for i_cam = 1:n_cams
    numAvail(i_cam) = vid(i_cam).FramesAvailable;
    [~, time{i_cam,1}, metadata] = getdata(vid(i_cam),numAvail(i_cam));
end

% save([save_folder '\' 'video_time_train_' num2str(i_train) '_' date_stamp '.mat' '.mat'],'time','rep_start_time','rep_start_now', 'metadata','-v7.3')

clear time numAvail



function ET_fs_Callback(hObject, eventdata, handles)
% hObject    handle to ET_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_fs as text
%        str2double(get(hObject,'String')) returns contents of ET_fs as a double


% --- Executes during object creation, after setting all properties.
function ET_fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_AmplitudeRange_Callback(hObject, eventdata, handles)
% hObject    handle to ET_AmplitudeRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_AmplitudeRange as text
%        str2double(get(hObject,'String')) returns contents of ET_AmplitudeRange as a double


% --- Executes during object creation, after setting all properties.
function ET_AmplitudeRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_AmplitudeRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_PulseWidthRange_Callback(hObject, eventdata, handles)
% hObject    handle to ET_PulseWidthRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_PulseWidthRange as text
%        str2double(get(hObject,'String')) returns contents of ET_PulseWidthRange as a double


% --- Executes during object creation, after setting all properties.
function ET_PulseWidthRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_PulseWidthRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_MonoOrBiPhasic_Callback(hObject, eventdata, handles)
% hObject    handle to ET_MonoOrBiPhasic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_MonoOrBiPhasic as text
%        str2double(get(hObject,'String')) returns contents of ET_MonoOrBiPhasic as a double


% --- Executes during object creation, after setting all properties.
function ET_MonoOrBiPhasic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_MonoOrBiPhasic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_FrequencyRange_Callback(hObject, eventdata, handles)
% hObject    handle to ET_FrequencyRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_FrequencyRange as text
%        str2double(get(hObject,'String')) returns contents of ET_FrequencyRange as a double


% --- Executes during object creation, after setting all properties.
function ET_FrequencyRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_FrequencyRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_TrainDuration_Callback(hObject, eventdata, handles)
% hObject    handle to ET_TrainDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_TrainDuration as text
%        str2double(get(hObject,'String')) returns contents of ET_TrainDuration as a double


% --- Executes during object creation, after setting all properties.
function ET_TrainDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_TrainDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_SaveFolder_Callback(hObject, eventdata, handles)
% hObject    handle to ET_SaveFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_SaveFolder as text
%        str2double(get(hObject,'String')) returns contents of ET_SaveFolder as a double


% --- Executes during object creation, after setting all properties.
function ET_SaveFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_SaveFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_SaveName_Callback(hObject, eventdata, handles)
% hObject    handle to ET_SaveName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_SaveName as text
%        str2double(get(hObject,'String')) returns contents of ET_SaveName as a double


% --- Executes during object creation, after setting all properties.
function ET_SaveName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_SaveName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [vid, src] = setup_video(n_cams)
    %% setup video

low_res_flag = 1;
    
for i_cam = 1:n_cams

    if low_res_flag == 1
        vid(i_cam) = videoinput('winvideo', i_cam, 'MJPG_320x240'); % could use 640x480
    else
        vid(i_cam) = videoinput('winvideo', i_cam, 'MJPG_640x480'); % could use 640x480
    end

    src(i_cam) = getselectedsource(vid(i_cam));
    src(i_cam).BacklightCompensation =  'on';
    src(i_cam).ExposureMode='manual';
    src(i_cam).Exposure = -6;

    preview(vid(i_cam))
end

choice = menu('Preview Control','End preview and continue');

for i_cam = 1:n_cams
    closepreview(vid(i_cam))
end

for i_cam = 1:n_cams
    vid(i_cam).LoggingMode = 'disk&memory';
    vid(i_cam).FramesPerTrigger = Inf;
    triggerconfig(vid(i_cam), 'manual')
end

function Process_Plot_Save(src,event,mf,plot_handle,n_cams,vid, seizure_detection_ch_vec) 

% disp('Process_Plot_Save')

global stim_flag


%% log the frame timestamps
% save metadata too!
for i_cam = 1:n_cams
    numAvail(i_cam) = vid(i_cam).FramesAvailable;
    [~, time{i_cam,1}, metadata] = getdata(vid(i_cam),numAvail(i_cam));

    [ro_frame, ~] = size(mf,['time_cam_' num2str(i_cam)]);
    [ro_frame_meta, ~] = size(mf,['meta_time_cam_' num2str(i_cam)]);
    
    n_new_frames = size(time{i_cam,1},1);
    
    meta_mat = table2array(struct2table(metadata));
    n_new_frames_meta = size(meta_mat,1);
    
    if ro_frame ~= ro_frame_meta
        warning('ro_frame ~= ro_frame_meta')
    end
    if n_new_frames ~= n_new_frames_meta
        warning('n_new_frames ~= n_new_frames_meta')
    end
    
    eval(['mf.time_cam_' num2str(i_cam) '(ro_frame+(1:n_new_frames),1) = time{i_cam,1};'])
    
    eval(['mf.meta_time_cam_' num2str(i_cam) '(ro_frame_meta+(1:n_new_frames),:) = meta_mat;'])
end
%%

% also log the data
deci = 5;

% MainPlot_handle = findobj('Tag', 'A_MainPlot');

data_deci = event.Data(1:deci:end,:);
time_deci = event.TimeStamps(1:deci:end);
n_t_deci = length(time_deci);
n_ch = size(data_deci,2);

channel_spacing = 5;

plot(plot_handle, time_deci, data_deci+repmat(channel_spacing*(1:n_ch), n_t_deci,1))
ylim(plot_handle,[0 channel_spacing*(n_ch+1)])

data = event.Data;

[ro, co] = size(mf,'data');

[r, c] = size(event.Data);

if co~=c+1
    error('data dims not equal in columns')
end

mf.data(ro+(1:r), :) = [event.TimeStamps event.Data];

%% seizure detection code


for i_ch = 1:length(seizure_detection_ch_vec)
    diff_data = diff(event.Data(:,seizure_detection_ch_vec(i_ch)+1));
    if any(diff_data(:,1)>.1)
        if stim_flag(i_ch) == 1
            stim_flag(i_ch) = 0;
        elseif stim_flag(i_ch) == 0
            stim_flag(i_ch) = 1;
        else
            error('stim_flag must be 1 or 0')
        end
    end
end

function Generate_Stim_Vec(src, event)

% disp('Generate_Stim_Vec')

global stim_flag
global fs

n_ch_out = length(stim_flag);

out_chunk = 1/2; % hard coded half a second
% out_chunk = 1/6; % hard coded half a second

dt = 1/fs;

n_t = fix(fs*out_chunk);

data=zeros(n_t, n_ch_out);

disp(num2str(stim_flag))

ind_stim = find(stim_flag);

data(:,ind_stim) = 3;

queueOutputData(src,data);



function ET_n_cams_Callback(hObject, eventdata, handles)
% hObject    handle to ET_n_cams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_n_cams as text
%        str2double(get(hObject,'String')) returns contents of ET_n_cams as a double


% --- Executes during object creation, after setting all properties.
function ET_n_cams_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_n_cams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_n_ch_out_Callback(hObject, eventdata, handles)
% hObject    handle to ET_n_ch_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_n_ch_out as text
%        str2double(get(hObject,'String')) returns contents of ET_n_ch_out as a double


% --- Executes during object creation, after setting all properties.
function ET_n_ch_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_n_ch_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_n_ch_in_Callback(hObject, eventdata, handles)
% hObject    handle to ET_n_ch_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_n_ch_in as text
%        str2double(get(hObject,'String')) returns contents of ET_n_ch_in as a double


% --- Executes during object creation, after setting all properties.
function ET_n_ch_in_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_n_ch_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_seizure_detection_channels_Callback(hObject, eventdata, handles)
% hObject    handle to ET_seizure_detection_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_seizure_detection_channels as text
%        str2double(get(hObject,'String')) returns contents of ET_seizure_detection_channels as a double


% --- Executes during object creation, after setting all properties.
function ET_seizure_detection_channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_seizure_detection_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
