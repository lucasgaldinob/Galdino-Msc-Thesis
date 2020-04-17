%% Data parameters

close all
clear
clc

Fs = 1000;
Experiment_Time = 600; %600 seconds = 10 minutes
dt = 1/Fs;
t = dt:dt:600; %time_vector

salPrl = 10;
salNac = 20;
salVta = 30;
salBla = 40;
salVh = 50;

gbrPrl = 60;
gbrNac = 70;
gbrVta = 80;
gbrBla = 90;
gbrVh = 100; 

nanimais_SAL = 5;
nanimais_GBR = 7;

Baseline_Sal = [];
for i=1:nanimais_SAL
    Baseline_Sal(i).AD01 =  (sin(2*pi*salPrl*t)+0.5*randn(size(t)))';
    Baseline_Sal(i).AD02 =  (sin(2*pi*salNac*t)+0.5*randn(size(t)))';
    Baseline_Sal(i).AD03 =  (sin(2*pi*salVta*t)+0.5*randn(size(t)))';
    Baseline_Sal(i).AD04 =  (sin(2*pi*salBla*t)+0.5*randn(size(t)))';
    Baseline_Sal(i).AD05 =  (sin(2*pi*salVh*t)+0.5*randn(size(t)))';
end

% Filtering signal
for i=1:nanimais_SAL
    Baseline_Sal(i).AD01=eegfilt(Baseline_Sal(i).AD01',Fs,0,100);
    Baseline_Sal(i).AD01=Baseline_Sal(i).AD01';
    Baseline_Sal(i).AD02=eegfilt(Baseline_Sal(i).AD02',Fs,0,100);
    Baseline_Sal(i).AD02=Baseline_Sal(i).AD02';
    Baseline_Sal(i).AD03=eegfilt(Baseline_Sal(i).AD03',Fs,0,100);
    Baseline_Sal(i).AD03=Baseline_Sal(i).AD03';
    Baseline_Sal(i).AD04=eegfilt(Baseline_Sal(i).AD04',Fs,0,100);
    Baseline_Sal(i).AD04=Baseline_Sal(i).AD04';
    Baseline_Sal(i).AD05=eegfilt(Baseline_Sal(i).AD05',Fs,0,100);
    Baseline_Sal(i).AD05=Baseline_Sal(i).AD05';
end


%Loop to create 300 windows containing 2s data
for i=1:nanimais_SAL
    Start=1;
    End=2000;
    j=2000;
    for k=1:(Experiment_Time/2) %Where Experiment_Time is the duration of 
        %the electrophysiological data and 2 is the window determined to windowing
        %Each line is a 2s window
        SAL_Baseline_Segmented(i).PRL(k,:)=Baseline_Sal(i).AD01(Start:End,:);
        SAL_Baseline_Segmented(i).NAC(k,:)=Baseline_Sal(i).AD02(Start:End,:);
        SAL_Baseline_Segmented(i).VTA(k,:)=Baseline_Sal(i).AD03(Start:End,:);
        SAL_Baseline_Segmented(i).BLA(k,:)=Baseline_Sal(i).AD04(Start:End,:);
        SAL_Baseline_Segmented(i).HV(k,:)=Baseline_Sal(i).AD05(Start:End,:);
        Start=Start+j;
        End=End+j;
    end
end
clear Start End i j k 

for i=1:nanimais_SAL
    for j=1:(Experiment_Time/2)
        if max((SAL_Baseline_Segmented(i).PRL(j,:)) <= (1.4*std(Baseline_Sal(i).AD01)))  
            SAL_Baseline_Seg{i,1}(j,:)=SAL_Baseline_Segmented(i).PRL(j,:);
        end
        if max((SAL_Baseline_Segmented(i).NAC(j,:)) <= (1.4*std(Baseline_Sal(i).AD02)))  
            SAL_Baseline_Seg{i,2}(j,:)=SAL_Baseline_Segmented(i).NAC(j,:);
        end
        if max((SAL_Baseline_Segmented(i).VTA(j,:)) <= (1.4*std(Baseline_Sal(i).AD03)))  
            SAL_Baseline_Seg{i,3}(j,:)=SAL_Baseline_Segmented(i).VTA(j,:);
        end  
        if max((SAL_Baseline_Segmented(i).BLA(j,:)) <= (1.4*std(Baseline_Sal(i).AD04)))  
            SAL_Baseline_Seg{i,4}(j,:)=SAL_Baseline_Segmented(i).BLA(j,:);
        end
        
        if max((SAL_Baseline_Segmented(i).HV(j,:)) <= (1.4*std(Baseline_Sal(i).AD05)))
            SAL_Baseline_Seg{i,5}(j,:)=SAL_Baseline_Segmented(i).HV(j,:); 
        else
            j=j+1;
        end
            
    end
end

clear Baseline_Data i j SAL_Baseline_Segmented
clc


WINDOW=2000; %Window length
Fs=1000; %Sampling rate
NFFT=4096; %NFFT is the number of points in a window (2s = 1000 NFFT) // no zero padding


for i=1:nanimais_SAL
    Pxx=[];
    for j=1:size(SAL_Baseline_Seg,2)
        Pxx=[];
        for z=1:size(SAL_Baseline_Seg{i, j},1)
            [Pxx F] = pwelch(SAL_Baseline_Seg{i, j}(z,:),WINDOW,[],NFFT,Fs);
            Baseline_SAL_All{i, j}(z,:)=Pxx;
        end
    end
end


CatSAL = {};
for i=1:5
    A = [];
    for j=1:nanimais_SAL
        A = [A; Baseline_SAL_All{j, i}];
    end
    CatSAL {1,i} = mean(A,1);
    clear A
end

%% GBR
Baseline_Gbr = [];
for i=1:nanimais_GBR
    Baseline_Gbr(i).AD01 =  (sin(2*pi*gbrPrl*t)+0.5*randn(size(t)))';
    Baseline_Gbr(i).AD02 =  (sin(2*pi*gbrNac*t)+0.5*randn(size(t)))';
    Baseline_Gbr(i).AD03 =  (sin(2*pi*gbrVta*t)+0.5*randn(size(t)))';
    Baseline_Gbr(i).AD04 =  (sin(2*pi*gbrBla*t)+0.5*randn(size(t)))';
    Baseline_Gbr(i).AD05 =  (sin(2*pi*gbrVh*t)+0.5*randn(size(t)))';
end

% Filtering signal
for i=1:nanimais_GBR
    Baseline_Gbr(i).AD01=eegfilt(Baseline_Gbr(i).AD01',Fs,0,100);
    Baseline_Gbr(i).AD01=Baseline_Gbr(i).AD01';
    Baseline_Gbr(i).AD02=eegfilt(Baseline_Gbr(i).AD02',Fs,0,100);
    Baseline_Gbr(i).AD02=Baseline_Gbr(i).AD02';
    Baseline_Gbr(i).AD03=eegfilt(Baseline_Gbr(i).AD03',Fs,0,100);
    Baseline_Gbr(i).AD03=Baseline_Gbr(i).AD03';
    Baseline_Gbr(i).AD04=eegfilt(Baseline_Gbr(i).AD04',Fs,0,100);
    Baseline_Gbr(i).AD04=Baseline_Gbr(i).AD04';
    Baseline_Gbr(i).AD05=eegfilt(Baseline_Gbr(i).AD05',Fs,0,100);
    Baseline_Gbr(i).AD05=Baseline_Gbr(i).AD05';
end


%Loop to create 300 windows containing 2s data
for i=1:nanimais_GBR
    Start=1;
    End=2000;
    j=2000;
    for k=1:(Experiment_Time/2) %Where Experiment_Time is the duration of 
        %the electrophysiological data and 2 is the window determined to windowing
        %Each line is a 2s window
        GBR_Baseline_Segmented(i).PRL(k,:)=Baseline_Gbr(i).AD01(Start:End,:);
        GBR_Baseline_Segmented(i).NAC(k,:)=Baseline_Gbr(i).AD02(Start:End,:);
        GBR_Baseline_Segmented(i).VTA(k,:)=Baseline_Gbr(i).AD03(Start:End,:);
        GBR_Baseline_Segmented(i).BLA(k,:)=Baseline_Gbr(i).AD04(Start:End,:);
        GBR_Baseline_Segmented(i).HV(k,:)=Baseline_Gbr(i).AD05(Start:End,:);
        Start=Start+j;
        End=End+j;
    end
end
clear Start End i j k 

for i=1:nanimais_GBR
    for j=1:(Experiment_Time/2)
        if max((GBR_Baseline_Segmented(i).PRL(j,:)) <= (1.4*std(Baseline_Gbr(i).AD01)))  
            GBR_Baseline_Seg{i,1}(j,:)=GBR_Baseline_Segmented(i).PRL(j,:);
        end
        if max((GBR_Baseline_Segmented(i).NAC(j,:)) <= (1.4*std(Baseline_Gbr(i).AD02)))  
            GBR_Baseline_Seg{i,2}(j,:)=GBR_Baseline_Segmented(i).NAC(j,:);
        end
        if max((GBR_Baseline_Segmented(i).VTA(j,:)) <= (1.4*std(Baseline_Gbr(i).AD03)))  
            GBR_Baseline_Seg{i,3}(j,:)=GBR_Baseline_Segmented(i).VTA(j,:);
        end  
        if max((GBR_Baseline_Segmented(i).BLA(j,:)) <= (1.4*std(Baseline_Gbr(i).AD04)))  
            GBR_Baseline_Seg{i,4}(j,:)=GBR_Baseline_Segmented(i).BLA(j,:);
        end
        
        if max((GBR_Baseline_Segmented(i).HV(j,:)) <= (1.4*std(Baseline_Gbr(i).AD05)))
            GBR_Baseline_Seg{i,5}(j,:)=GBR_Baseline_Segmented(i).HV(j,:); 
        else
            j=j+1;
        end
            
    end
end


for i=1:nanimais_GBR
    Pxx=[];
    for j=1:size(GBR_Baseline_Seg,2)
        Pxx=[];
        for z=1:size(GBR_Baseline_Seg{i, j},1)
            [Pxx F] = pwelch(GBR_Baseline_Seg{i, j}(z,:),WINDOW,[],NFFT,Fs);
            Baseline_GBR_All{i, j}(z,:)=Pxx;
        end
    end
end


CatGBR = {};
for i=1:5
    A = [];
    for j=1:nanimais_GBR
        A = [A; Baseline_GBR_All{j, i}];
    end
    CatGBR {1,i} = mean(A,1);
    clear A
end


%% Plot 

plot(F, CatGBR{1, 1}, 'linewidth', 2)
hold on
plot(F, CatGBR{1, 2}, 'linewidth', 2)
hold on
plot(F, CatGBR{1, 3}, 'linewidth', 2)
hold on
plot(F, CatGBR{1, 4}, 'linewidth', 2)
hold on
plot(F, CatGBR{1, 5}, 'linewidth', 2)
hold on
xlim ([1 120])
legend('PRL','NAc','VTA', 'BLA', 'VH')
xlabel('Frequency (Hz)','fontsize',16);
ylabel('Power Spectrum Density (W/Hz)','fontsize',16);
title('GBR group','fontsize',20);

plot(F, CatSAL{1, 1}, 'linewidth', 2)
hold on
plot(F, CatSAL{1, 2}, 'linewidth', 2)
hold on
plot(F, CatSAL{1, 3}, 'linewidth', 2)
hold on
plot(F, CatSAL{1, 4}, 'linewidth', 2)
hold on
plot(F, CatSAL{1, 5}, 'linewidth', 2)
hold on
xlim ([1 120])
legend('PRL','NAc','VTA', 'BLA', 'VH')
xlabel('Frequency (Hz)','fontsize',16);
ylabel('Power Spectrum Density (W/Hz)','fontsize',16);
title('SAL group','fontsize',20);

