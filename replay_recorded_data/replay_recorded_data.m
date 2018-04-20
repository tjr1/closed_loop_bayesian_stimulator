close all
clear
clc

% daqreset

path_to_data = 'C:\Users\Scientist\Desktop\seizure_data\good\20151119_DV001_DV002_DV003_DV004__2015_11_19___12_53_51.mat';

load(path_to_data,'fs','sbuf')

channels = [1, 3];

data = sbuf(:, channels);

clear sbuf

up_sample_factor = 1;

if up_sample_factor ~= 1
    data = resample(data,up_sample_factor,1);
    fs = up_sample_factor*fs;
end

data(data>10) = 10;
data(data<-10) = -10;

dt = 1/fs;
time = (0:size(data,1)-1)*dt;

% figure(1)
% plot(time,[data(:,1) data(:,2)+20])

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

s.Rate = fs;
s.IsContinuous = true;

addAnalogOutputChannel(s, USB_6229_dev.ID, (1:length(channels))-1, 'Voltage');

queueOutputData(s, data);

startBackground(s)

choice = menu('stop?','Stop');
if choice == 1
    s.stop()
end