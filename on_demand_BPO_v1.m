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

% Last Modified by GUIDE v2.5 23-Apr-2018 16:37:00

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
    set(handles.ET_high_pass_filter,'String',struc.ET_high_pass_filter);
end
try
    set(handles.ET_low_pass_filter,'String',struc.ET_low_pass_filter);
end
try
    set(handles.ET_fast_slow_ratio_thresh,'String',struc.ET_fast_slow_ratio_thresh);
end
try
    set(handles.ET_G_fast,'String',struc.ET_G_fast);
end
try
    set(handles.ET_G_slow,'String',struc.ET_G_slow);
end
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


struc.ET_high_pass_filter=get(handles.ET_high_pass_filter,'String');
struc.ET_low_pass_filter=get(handles.ET_low_pass_filter,'String');
struc.ET_fast_slow_ratio_thresh=get(handles.ET_fast_slow_ratio_thresh,'String');
struc.ET_g_fast=get(handles.ET_G_fast,'String');
struc.ET_g_slow=get(handles.ET_G_slow,'String');
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

%% 

global q

p = gcp();
q{1,1} = parallel.pool.PollableDataQueue;
q{1,2} = parallel.pool.PollableDataQueue;

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

global out_chunk in_chunk
out_chunk = 1; % these are very long, 1/2 second would be better, but computation time of bayes opt is too slow to fit inside the loop
in_chunk = 1;

n_ch_out = str2num(get(handles.ET_n_ch_out,'String'));

seizure_detection_ch_vec = str2num(get(handles.ET_seizure_detection_channels,'String'));% detect seizures on ACH0 and ACH2, these are BNC-2090 numbers

if length(seizure_detection_ch_vec) ~= n_ch_out
    error('length(seizure_detection_ch_vec) ~= n_ch_out')
end

addAnalogOutputChannel(s, PCI_6251_dev.ID, (1:n_ch_out)-1, 'Voltage');

global stim_flag
stim_flag = zeros(1, n_ch_out,'logical');

n_ch_in = str2num(get(handles.ET_n_ch_in,'String')); % hard coded 1 input channel
for i_ch = 1:n_ch_in % add ephys recording channels
    ch(i_ch) = addAnalogInputChannel(s,PCI_6251_dev.ID, i_ch-1, 'Voltage');
end

%%
%%
global n_read fast_int slow_int Seizure_On Seizure_Off Seizure_Count Seizure_Duration stim_amp stim_freq Seizure_Start_Ind next_freq next_amp
n_read = 1;
fast_int = zeros(1,n_ch_out);
slow_int = zeros(1,n_ch_out);
Seizure_On = zeros(1,n_ch_out);
Seizure_Off = zeros(1,n_ch_out);
Seizure_Count = zeros(1,n_ch_out);
Seizure_Duration = zeros(1,n_ch_out);
Seizure_Start_Ind = zeros(1,n_ch_out);

next_freq = 10*ones(1,n_ch_out); % initialize arbitrarily to 10 Hz
next_amp = 1*ones(1,n_ch_out); % initialize arbitrarily to 1 unit

set(handles.T_NextFreq,'String',num2str(next_freq));
set(handles.T_NextAmp,'String',num2str(next_amp));

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

save(save_mat_path,'fast_int','slow_int','Seizure_On','Seizure_Off','stim_flag','Seizure_Duration','-nocompression','-append')

clear data
mf = matfile(save_mat_path,'Writable',true);

%% specify optimized variables
opt_freq = optimizableVariable('frequency',[1 300]); % Hz
opt_amp = optimizableVariable('amplitude',[-4 4]); % Volts?

%% add listeners
lh1 = addlistener(s,'DataAvailable', @(src, event) Process_Plot_Save(src,event,mf, handles.A_MainPlot,n_cams,vid,seizure_detection_ch_vec,n_ch_out,opt_freq,opt_amp, p, handles) );
lh2 = addlistener(s,'DataRequired', @Generate_Stim_Vec);

s.NotifyWhenDataAvailableExceeds = fix(fs*in_chunk); % input buffer treshold, hard coded
s.NotifyWhenScansQueuedBelow = fix(fs*out_chunk); % output buffer threshold, hard coded
% s.NotifyWhenDataAvailableExceeds = fix(fs/5); % input buffer treshold, hard coded
% s.NotifyWhenScansQueuedBelow = fix(fs/6); % output buffer threshold, hard coded


handles.s = s;

%% buffer 2 seconds of zeros to the daq output to start with
stim_vec_init = zeros(fix(5*fs),n_ch_out);
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

function Process_Plot_Save(src,event,mf,plot_handle,n_cams,vid, seizure_detection_ch_vec, n_ch_out,opt_freq,opt_amp, p, handles) 

% disp('Process_Plot_Save')

global stim_flag

global q

global n_read fast_int slow_int Seizure_On Seizure_Off Seizure_Count Seizure_Duration stim_amp stim_freq Seizure_Start_Ind next_freq next_amp

n_read = n_read+1;


%% log the frame timestamps
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
%% plot decimated data
deci = 5;

% MainPlot_handle = findobj('Tag', 'A_MainPlot');

data_deci = event.Data(1:deci:end,:);
time_deci = event.TimeStamps(1:deci:end);
n_t_deci = length(time_deci);
n_ch = size(data_deci,2);

channel_spacing = 5;

plot(plot_handle, time_deci, data_deci+repmat(channel_spacing*(1:n_ch), n_t_deci,1))
hold(plot_handle, 'on')
ylim(plot_handle,[0 channel_spacing*(n_ch+1)])

%% save raw data
data = event.Data;

[ro, co] = size(mf,'data');

[r, c] = size(event.Data);

if co~=c+1
    error('data dims not equal in columns')
end

mf.data(ro+(1:r), :) = [event.TimeStamps event.Data];

%% seizure detection code
% decimate and remove mean

detect_deci = 10;
fs = str2double(get(handles.ET_fs,'String'));
fs_deci = fs/detect_deci;
dt_deci = 1/fs_deci;

detect_data = event.Data(1:detect_deci:end,seizure_detection_ch_vec+1); % select channels and decimate
detect_data = detect_data - repmat(mean(detect_data,1), size(detect_data,1),1); % remove mean

n_t_detect_deci = size(detect_data,1);
detect_time_deci = event.TimeStamps(1:detect_deci:end);

%% filter data
hp_f = str2double(get(handles.ET_high_pass_filter,'String'));
lp_f = str2double(get(handles.ET_low_pass_filter,'String'));

if hp_f ~= 0
    [bH,aH] = butter(3,hp_f/(fs_deci/2),'high');
    detect_data = filtfilt(bH,aH,detect_data);
end
if lp_f ~= 0
    [bL,aL] = butter(3,lp_f/(fs_deci/2),'low');
    detect_data = filtfilt(bL,aL,detect_data);
end

plot(plot_handle, detect_time_deci, detect_data+repmat(channel_spacing*(seizure_detection_ch_vec+1), n_t_detect_deci,1))
hold(plot_handle, 'off')

%% fast/slow stdev detector
detect_data_std = std(detect_data,0,1);
G_fast = str2double(get(handles.ET_G_fast,'String'));
G_slow = str2double(get(handles.ET_G_slow,'String'));
if n_read == 2
    fast_int(1,:) = detect_data_std;
    slow_int(1,:) = detect_data_std;
    mf.fast_int(1,:) = fast_int(1,:);
    mf.slow_int(1,:) = slow_int(1,:);
end
fast_int(n_read,:) = (1-G_fast)*detect_data_std + G_fast*fast_int(n_read-1,:);
slow_int(n_read,:) = (1-G_slow)*detect_data_std + G_slow*slow_int(n_read-1,:);

mf.fast_int(n_read,:) = fast_int(n_read,:);
mf.slow_int(n_read,:) = slow_int(n_read,:);

disp('.')
disp(['fast std = ' num2str(fast_int(n_read,:))])
disp(['slow std = ' num2str(slow_int(n_read,:))])

fast_slow_ratio_thresh = str2double(get(handles.ET_fast_slow_ratio_thresh,'String'));

fast_slow_ratio = fast_int(n_read,:)./slow_int(n_read,:);

disp(['fast/slow ratio = ' num2str(fast_slow_ratio)])

fast_slow_ratio_trigger = fast_slow_ratio > fast_slow_ratio_thresh;

%% seizue starts
Seizure_On(n_read,:) = fast_slow_ratio_trigger; % add additional logic of triggers here 
mf.Seizure_On(n_read,:) = Seizure_On(n_read,:);

disp(['Seizure_On = ' num2str(Seizure_On(n_read,:))])

%% stim when seizure detected
stim_flag = and(Seizure_On(n_read,:)==1, Seizure_On(n_read-1,:)==0); % only the up thresholds
mf.stim_flag(n_read,:) = stim_flag;
Seizure_Start_Ind(stim_flag) = n_read;

%% seizure ends
Seizure_Off(n_read,:) = and(Seizure_On(n_read,:)==0, Seizure_On(n_read-1,:)==1);
mf.Seizure_Off(n_read,:) = Seizure_Off(n_read,:);

disp(['Seizure_Off = ' num2str(Seizure_Off(n_read,:))])

%% count seizures
Seizure_Count = Seizure_Count + Seizure_Off(n_read,:);

disp(['Seizure_Count = ' num2str(Seizure_Count)])

%% determine seizure duration

time_per_read = 0.5; % hard coded half second, should be connected to listenerj

global f

% Seizure_Duration stim_amp stim_freq
for i_ch_out = 1:n_ch_out
    if Seizure_Off(n_read,i_ch_out)==1
        % determine seizure length
        Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out) = time_per_read * (n_read- Seizure_Start_Ind(1,i_ch_out));
        mf.Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out) = Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out);
        
        stim_freq(Seizure_Count(1,i_ch_out),i_ch_out) = next_freq(1,i_ch_out);
        stim_amp(Seizure_Count(1,i_ch_out),i_ch_out) = next_amp(1,i_ch_out);
        
        %% Bayes Opt for next stimulation parameters
        InitialObjective = Seizure_Duration(1:Seizure_Count(1,i_ch_out),i_ch_out);
        freq = stim_freq(1:Seizure_Count(1,i_ch_out),i_ch_out);
        amp = stim_amp(1:Seizure_Count(1,i_ch_out),i_ch_out);
        InitialX = table(freq, amp);
        
%         res = bayesopt(@place_holder_fcn,[opt_freq, opt_amp],'InitialX',InitialX,'InitialObjective',InitialObjective, 'MaxObjectiveEvaluations', 1, 'PlotFcn', [], 'AcquisitionFunctionName', 'expected-improvement-plus', 'ExplorationRatio', .4);
        % might need to do this in a seperate instance of matlab so that it
        % can be asynchronous, https://www.mathworks.com/matlabcentral/answers/83051-how-to-make-two-instances-of-matlab-communicate
        
%         f{1,i_ch_out} = parfeval(@BO_wrapper,0,opt_amp, InitialX, InitialObjective, q);
        f(1,i_ch_out) = parfeval(@BO_wrapper,0,opt_freq, opt_amp, InitialX, InitialObjective, q{1,i_ch_out});
        
    end
end

% for i_ch_out = 1:n_ch_out
for i_ch_out = 1:n_ch_out
    [res, gotMsg] = poll(q{1,i_ch_out}, .05);
    gotMsg=gotMsg

    if gotMsg
%         close all
        tic
        figure(i_ch_out)
        plot(res,@plotObjectiveModel) %  A_BPO_vis, would be good to make these into a video
        toc

        next_freq(1,i_ch_out) = res.NextPoint{1,1};
        next_amp(1,i_ch_out) = res.NextPoint{1,2};

        set(handles.T_NextFreq,'String',num2str(next_freq));
        set(handles.T_NextAmp,'String',num2str(next_amp));
    end
end

disp('Seizure_Duration = ')
disp(num2str(Seizure_Duration))

function BO_wrapper(opt_freq, opt_amp, InitialX, InitialObjective, que)

res = bayesopt(@place_holder_fcn,[opt_freq, opt_amp],'InitialX',InitialX,'InitialObjective',InitialObjective, 'MaxObjectiveEvaluations', 1, 'PlotFcn', [], 'AcquisitionFunctionName', 'expected-improvement-plus', 'ExplorationRatio', .4);
send(que,res)

function Generate_Stim_Vec(src, event)

% disp('Generate_Stim_Vec')

global stim_flag
global fs
global out_chunk

n_ch_out = length(stim_flag);

dt = 1/fs;

n_t = fix(fs*out_chunk);

data=zeros(n_t, n_ch_out);

disp(['stim flag = ' num2str(stim_flag)])

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



function ET_G_fast_Callback(hObject, eventdata, handles)
% hObject    handle to ET_G_fast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_G_fast as text
%        str2double(get(hObject,'String')) returns contents of ET_G_fast as a double


% --- Executes during object creation, after setting all properties.
function ET_G_fast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_G_fast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_G_slow_Callback(hObject, eventdata, handles)
% hObject    handle to ET_G_slow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_G_slow as text
%        str2double(get(hObject,'String')) returns contents of ET_G_slow as a double


% --- Executes during object creation, after setting all properties.
function ET_G_slow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_G_slow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_fast_slow_ratio_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to ET_fast_slow_ratio_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_fast_slow_ratio_thresh as text
%        str2double(get(hObject,'String')) returns contents of ET_fast_slow_ratio_thresh as a double


% --- Executes during object creation, after setting all properties.
function ET_fast_slow_ratio_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_fast_slow_ratio_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_high_pass_filter_Callback(hObject, eventdata, handles)
% hObject    handle to ET_high_pass_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_high_pass_filter as text
%        str2double(get(hObject,'String')) returns contents of ET_high_pass_filter as a double


% --- Executes during object creation, after setting all properties.
function ET_high_pass_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_high_pass_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ET_low_pass_filter_Callback(hObject, eventdata, handles)
% hObject    handle to ET_low_pass_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ET_low_pass_filter as text
%        str2double(get(hObject,'String')) returns contents of ET_low_pass_filter as a double


% --- Executes during object creation, after setting all properties.
function ET_low_pass_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_low_pass_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [z] = place_holder_fcn(x, y)
z = x+y; % not actually used.
