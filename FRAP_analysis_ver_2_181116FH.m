% This script uses Jie's analysis code, which has the proper normalization
% calculation, to analyze the data. This script requires you to first have
% each dataset analyzed with Xinxing's code, either "FRAP_ring" or
% "FRAP_ring_circle".
%pathname = '/Users/bettridgeke/Dropbox/Data/RNAP_GFP_M9';
pathname = 'D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging';
cd(pathname)
%cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_M9')
%cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\23Feb23\cell_15')
%folders = dir('*Chlor*');
folders = dir('*20Mar23-FRAP-Chlor*');
%folders = dir('*15Mar23-FRAP-WT*');
addpath('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\MATLAB')

for ff = 1:length(folders)
    cd(folders(ff).name)
    vars = dir('cell*');
    
    for ii=1:length(vars)
        string = strsplit(vars(ii).name,'_');
        cd(vars(ii).name);
        num = string{2};
        %filename = sprintf('rpoC-GFPuv_I1309A_FRAP_cell_%s_30mW_r3_data.mat',num);
        %filename = sprintf('rpoC-GFPuv_WT_FRAP_cell_%s_30mW_r3_data.mat',num);
        filename = sprintf('rpoC-GFPuv_Chlor_FRAP_cell_%s_30mW_r3_data.mat',num);
        load(filename)
        FRAP_analysis_jie_KBmod_181116;
        save('FRAP_data_wCellROIFH.mat','CellI_aB','CellI_bB','RingI_bB',...
            'RingI_maxR','x','y_cell','y_ring','y_ringExp','y_ringNor',...
            'Elowitz_D','Elowitz_ampl','xx','yy','Brutton_D','ft','ft2',...
            'Brutton_mobile_fraction','Brutton_immobile_fraction',...
            'xxx','yyy','Brutton_gof','Fourier_gof','Brutton_output',...
            'Fourier_output','Fourier_fit','Brutton_fit');
       
            ii
           
        cd('..')

    end
    
    total_FRAP = [];
    k=1;
    for k=1:length(vars)
        cd(vars(k).name);
        load('FRAP_data_wCellROIFH.mat');
        total_FRAP(k).time = x;
        total_FRAP(k).y_ringNor = y_ringNor;
        total_FRAP(k).Fourier_fit = Fourier_fit;
        total_FRAP(k).Fourier_gof = Fourier_gof;
        total_FRAP(k).Fourier_output = Fourier_output;
        total_FRAP(k).Elowitz_D = Elowitz_D;
        total_FRAP(k).Elowitz_ampl = Elowitz_ampl;
        total_FRAP(k).xx = xx;
        total_FRAP(k).yy = yy;
        total_FRAP(k).Brutton_fit = Brutton_fit;
        total_FRAP(k).Brutton_gof = Brutton_gof;
        total_FRAP(k).Brutton_output = Brutton_output;
        total_FRAP(k).Brutton_D = Brutton_D;
        total_FRAP(k).ft = ft;
        total_FRAP(k).time = x';
        total_FRAP(k).ft2 = ft2;
        total_FRAP(k).Brutton_mobile_fraction = Brutton_mobile_fraction;
        total_FRAP(k).xxx = xxx;
        total_FRAP(k).yyy = yyy;
        cd('..')
    end
    average_FRAP.time = mean(vertcat(total_FRAP.time));
    average_FRAP.y_ringNor = mean(vertcat(total_FRAP.y_ringNor));
    average_FRAP.Elowitz_D = mean(vertcat(total_FRAP.Elowitz_D));
    average_FRAP.Brutton_D = mean(vertcat(total_FRAP.Brutton_D));
    average_FRAP.ft = mean(vertcat(total_FRAP.ft));
    average_FRAP.ft2 = mean(vertcat(total_FRAP.ft2));
    average_FRAP.Brutton_mobile_fraction = mean(vertcat(total_FRAP.Brutton_mobile_fraction));
    
    save('total_FRAPFH.mat','total_FRAP');
    save('average_FRAPFH.mat', 'average_FRAP');
    cd('..')
end
%% Combine Normalized Data for all Days
pathname = 'D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging';
cd(pathname)
%cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_M9')
%cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\23Feb23\cell_15')
folders = dir('*FRAP-Chlor');
%folders = dir('*Feb23-FRAP-Chlor');
addpath('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\MATLAB')
cd(folders(1).name)
load('total_FRAPFH.mat','total_FRAP');
load('average_FRAPFH.mat', 'average_FRAP');
total_FRAP_all_cells = total_FRAP;
avg_FRAP_all_cells = average_FRAP;
cd('..')
for ff = 2:length(folders)
    cd(folders(ff).name)
    load('total_FRAPFH.mat','total_FRAP');
    load('average_FRAPFH.mat', 'average_FRAP');
    total_FRAP_all_cells = cat(1,[total_FRAP_all_cells,total_FRAP]);
    avg_FRAP_all_cells = cat(1,[avg_FRAP_all_cells,average_FRAP]);
    cd('..')
end
save('total_FRAP_all_cells_FH.mat','total_FRAP_all_cells');
save('average_FRAP_all_cells_FH.mat', 'avg_FRAP_all_cells');

time = 1:601;
%avg_FRAP30 = mean(vertcat(total_FRAP_all_cells.ft2));
%std_FRAP = std(vertcat(total_FRAP_all_cells.ft2));
%avg_FRAP40 = mean(vertcat(total_FRAP_all_cells.ft2));
%std_FRAP = std(vertcat(total_FRAP_all_cells.ft2));
avg_FRAP = mean(vertcat(total_FRAP_all_cells.ft2));
std_FRAP = std(vertcat(total_FRAP_all_cells.ft2));
%total_FRAP_all_cells_30mw = total_FRAP_all_cells(1:5)
save('average_FRAP_FH.mat','time','avg_FRAP','std_FRAP')
FRAP_fit = fit(time(:), avg_FRAP(:), ...
   '1-a*exp(-koff1*x)-b*exp(-koff2*x)-C','StartPoint',...
   [0.45 0.2 0.012 0.5 0.15],'Lower',[0 0 0 0 0])

save('FRAP_fit.mat','FRAP_fit')

%%

%load('total_FRAP_all_cells_FH.mat')
figure
hold on
for ii=1:length(total_FRAP_all_cells)
    plot(total_FRAP_all_cells(ii).time(1:end), total_FRAP_all_cells(ii).ft2(1:end), 'Color', [rand rand rand])
    ii
    total_FRAP_all_cells(ii).ft2(326)
end
%axis([10 112 0.35 1])
axis([0 600 -2 5])
title('Example FRAP, cell 15 23Feb23','FontSize',14)
xlabel('time (s)', 'FontSize', 14)
ylabel('normalized FRAP signal', 'FontSize', 14)
set(gca,'FontSize',14)
hold off

%% combine different days
%mutants
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\16Mar23-FRAP-I1309A\total_FRAPFH.mat')
%total_FRAP1 = total_FRAP;
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\15Mar23-FRAP-I1309A\total_FRAPFH.mat')
%total_FRAP2= total_FRAP;
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\19Mar23-FRAP-I1309A\total_FRAPFH.mat')
%total_FRAP3= total_FRAP;
%WT
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\14Mar23-FRAP-WT\total_FRAPFH.mat')
%total_FRAP1 = total_FRAP;
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\15Mar23-FRAP-WT\total_FRAPFH.mat')
%total_FRAP2= total_FRAP;
%60mW Chlor
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\23Feb23-FRAP-Chlor\total_FRAPFH.mat')
%total_FRAP1 = total_FRAP;
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\25Feb23-FRAP-Chlor\total_FRAPFH.mat')
%total_FRAP2= total_FRAP;
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\26Feb23-FRAP-Chlor\total_FRAPFH.mat')
%total_FRAP3 = total_FRAP;
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\28Feb23-FRAP-Chlor\total_FRAPFH.mat')
%total_FRAP4= total_FRAP;
%30mW Chlor
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\18Mar23-FRAP-Chlor\total_FRAPFH.mat')
total_FRAP1 = total_FRAP;
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\19Mar23-FRAP-Chlor\total_FRAPFH.mat')
total_FRAP2= total_FRAP;
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\20Mar23-FRAP-Chlor\total_FRAPFH.mat')
total_FRAP3 = total_FRAP;



total_FRAP_all_cells = [total_FRAP1, total_FRAP2,total_FRAP3];
time = 1:601;
avg_FRAP = mean(vertcat(total_FRAP_all_cells.ft2));
std_FRAP = std(vertcat(total_FRAP_all_cells.ft2));

%%

figure
hold on
for ii=1:length(total_FRAP_all_cells)
    plot(total_FRAP_all_cells(ii).time(1:end), total_FRAP_all_cells(ii).ft2(1:end), 'Color', [rand rand rand])
    ii
  
end
figure
plot(total_FRAP_all_cells(1).time(1:end), avg_FRAP)

%%
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM-Rif\Analysis\FRAP_analysis_181119')
load('average_FRAP_all_cells.mat')
Rifavg = avg_FRAP;
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM\new_data_collection_method\Analysis\updated_analysis')
load('best_ezrdm_model.mat')
load("average_FRAP_all_cells.mat")
load("total_FRAP_all_cells.mat")
WTavg = avg_FRAP;WTmodel = ezrdm_model.FRAP;
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\avgFRAP30mW.mat')
avg30mW = avg_FRAP30mW;
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\avgFRAPmutant.mat')
avgmutant = avg_FRAPmutant;
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\18Mar23-FRAP-Chlor\total_FRAPFH.mat')
%FRAP_fit = fit(time(:), avg_FRAP(:),'1-a*exp(-koff1*x)-b*exp(-koff2*x)','StartPoint',[0.3 0.1 0.001 0.1])
%%
time = 1:601;
figure
hold on
plot(time, avg_FRAP)
plot(time, avg_FRAP30mW)
plot(time, avg_FRAPmutant)
plot(time, avg_FRAPchlor)
%plot(time, mean(vertcat(total_FRAP.ft2)))
%plot(time,avg_FRAPchlor60)
%plot(time, Rifavg)
%plot(total_FRAP(1).time(1:end), mean(vertcat(total_FRAP([3,5,7:9,11:14]).ft2)))
%legend('60mW',"50mW","30mW", "Kelsey", 'Location', 'Southeast')
legend("Kelsey", "30mW","mutant","chlor", 'Location', 'Southeast')
axis([0 600 0 1.5]) 
xlabel('time (s)','FontSize',16)
ylabel('average FRAP recovery','FontSize',16)
title('rpoC-GFP, EZRDM, N=24','FontSize',16)
set(gca,'FontSize',16)
%%
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\60mWFRAPavg.mat')
avg60mW = avg_FRAP2;
figure
hold on
plot(time(1:end), avg60mW(1:end))
plot(time, avg_FRAP)
legend('Chlor',"WT")
%legend('60mW',"50mW","40mW","30mW", "Kelsey 5", 'Location', 'Southeast')
axis([0 600 0 1.5])
xlabel('time (s)','FontSize',16)
ylabel('average FRAP recovery','FontSize',16)
title('rpoC-GFP, EZRDM, N=24','FontSize',16)
set(gca,'FontSize',16)

%save('average_FRAP.mat','time', 'avg_FRAP', 'std_FRAP')
%%
FRAP_fit = fit(time(:), avg_FRAP(:),'1-a*exp(-koff1*x)-b*exp(-koff2*x)','StartPoint',[0.3 0.1 0.001 0.1])
save('FRAP_fit.mat','FRAP_fit')
%%
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\20Mar23-FRAP-Chlor\total_FRAPFH.mat')
figure
hold on
for ii=1:length(total_FRAP)
    plot(total_FRAP(ii).time(1:end), total_FRAP(ii).ft2(1:end), 'Color', [rand rand rand])
    ii
    total_FRAP(ii).ft2(203)
end

%axis([10 112 0.35 1])
axis([0 600 -2 5])
title('Example FRAP, cell 15 23Feb23','FontSize',14)
xlabel('time (s)', 'FontSize', 14)
ylabel('normalized FRAP signal', 'FontSize', 14)
set(gca,'FontSize',14)
hold off
figure
plot(total_FRAP(1).time(1:end), mean(vertcat(total_FRAP.ft2)))