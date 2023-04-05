% figure script

%% Figure 1A -- example BF image

% load the relevant data
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_PAmCherry_EZRDM\5ms_continuous_data\140625_SMT_rpoC-PAmCherry_EZRDM_5ms\cell_1')
load('movie1traces.mat')
load('movie2traces.mat')
state1_xCoords = movie2traces.TracksROI(31).Coordinates(:,2);
state1_yCoords = movie2traces.TracksROI(31).Coordinates(:,3);
state3_xCoords = movie1traces.TracksROI(11).Coordinates(:,2);
state3_yCoords = movie1traces.TracksROI(11).Coordinates(:,3);
switch_xCoords = movie1traces.TracksROI(24).Coordinates(:,2);
switch_yCoords = movie1traces.TracksROI(24).Coordinates(:,3);

% display the BF image
t = Tiff('rpoC_PAmCherry_SMT_5ms_streaming_cell_1bf-01_FH_fromjpg.tif');
imgData = read(t);
imshow(imgData);
adjusted = imadjust(imgData);
imshow(adjusted);
I = imtophat(adjusted,ones(12, 12));
imshow(I);
% add traces
hold on
plot(state1_xCoords, state1_yCoords, '-', 'Color', [ 0.4941    0.7569    0.8706], 'LineWidth',2.0)
plot(state3_xCoords, state3_yCoords, '-','Color',[ 0.1137    0.2941    0.3686], 'LineWidth',2.0)
plot(switch_xCoords, switch_yCoords, '-', 'Color', [1 1 0.2] ,'LineWidth',2.0)
hold off;
%Scale bar added in FIJI, for some reason I had to add it to the tif, save
%it as a jpg and then save again as a tif. I think Tiff() didn't know how to
%interpret the scale bar the way it was saved by FIJI

%% Figure 1B -- single step displacement w/ Gaussian fit

% script to get the SFD with single and three Gaussian fits
% structure remade by D:\Xiao Lab Dropbox\Lab
% Members\Alumni\Bettridge_Kelsey\RNAP_PAmCherry_EZRDM\5ms_continuous_data\current_HMM_analysis\redone_HMM_160517\RNAP_Dest_distributions_FH.m
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\remadeStructures')
%cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_PAmCherry_EZRDM\5ms_continuous_data\current_HMM_analysis');
load('Dest_separated_trajs.mat')
figure
hold on
bar(xx_RNAP_state3, yy_RNAP_state3, 'FaceColor', [ 0.1137    0.2941    0.3686], 'FaceAlpha', 0.5, 'EdgeColor', [0.25 0.25 0.25])
bar(xx_RNAP_state2, yy_RNAP_state2, 'FaceColor', [ 0.2275    0.5647    0.7098], 'FaceAlpha', 0.5, 'EdgeColor', [0.25 0.25 0.25])
bar(xx_RNAP_state1, yy_RNAP_state1, 'FaceColor', [ 0.4941    0.7569    0.8706], 'FaceAlpha', 0.5, 'EdgeColor', [0.25 0.25 0.25])
plot(0:0.02:2, RNAP_state3_gaussfit(0:0.02:2), '--', 'Color', [ 0.1137    0.2941    0.3686], 'LineWidth', 2)
plot(0:0.02:2, RNAP_state2_gaussfit(0:0.02:2), '--', 'Color', [ 0.2275    0.5647    0.7098], 'LineWidth', 2)
plot(0:0.02:2, RNAP_state1_gaussfit(0:0.02:2), '--', 'Color', [ 0.4941    0.7569    0.8706], 'LineWidth', 2)
axis([-0.05 2 0 0.6])
legend('State III', "State II", "State I", 'EdgeColor', [1 1 1],'FontSize',24)
ylabel('Fraction of Molecules','FontSize',24)
xlabel('Individual Trajectory D_e_s_t (\mum^2/s)','FontSize',24)
set(gca, 'FontSize', 24)


%% Figure 1C -- HMM Pie Chart

%% FIGURE 1D -- squared displacements with fits separated by HMM States
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_kusumi_results.mat', 'rpoC_ezrdm')

figure
hold on
plot(rpoC_ezrdm(1).kusumi_25frames.xx, rpoC_ezrdm(1).kusumi_25frames.F, '-', 'Color', [0.4941 0.7569 0.8706], 'LineWidth', 1)
plot(rpoC_ezrdm(2).kusumi_25frames.xx, rpoC_ezrdm(2).kusumi_25frames.F, '-', 'Color', [0.2275 0.5647 0.7098], 'LineWidth', 1)
plot(rpoC_ezrdm(3).kusumi_15frames.xx, rpoC_ezrdm(3).kusumi_15frames.F, '-', 'Color', [0.1137 0.2941 0.3686], 'LineWidth', 1)
errorbar(rpoC_ezrdm(1).kusumi_25frames.D_all(1:10, 1), rpoC_ezrdm(1).kusumi_25frames.D_all(1:10,2), rpoC_ezrdm(1).kusumi_25frames.D_all(1:10,3), '.', 'MarkerSize', 15,'Color', [0.4941 0.7569 0.8706])
errorbar(rpoC_ezrdm(2).kusumi_25frames.D_all(1:10, 1), rpoC_ezrdm(2).kusumi_25frames.D_all(1:10, 2), rpoC_ezrdm(2).kusumi_25frames.D_all(1:10, 3), '.', 'MarkerSize', 15,'Color', [0.2275 0.5647 0.7098])
errorbar(rpoC_ezrdm(3).kusumi_15frames.D_all(1:10, 1), rpoC_ezrdm(3).kusumi_15frames.D_all(1:10, 2), rpoC_ezrdm(3).kusumi_15frames.D_all(1:10, 3), '.','MarkerSize', 15, 'Color', [0.1137 0.2941 0.3686])
axis([0 10.2*0.00674 0 0.06])
%legend('State I', "State II", "State III", 'EdgeColor', [1 1 1])
xlabel('time (s)', 'FontSize', 24, 'FontName', 'Helvetica')
ylabel('MSD (\mum^2)', 'FontSize', 24, 'FontName', 'Helvetica')
set(gca, 'FontSize', 24, 'FontName', 'Helvetica')

%% Figure 1E - FRAP Data
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM\new_data_collection_method\Analysis\updated_analysis')
load('best_ezrdm_model.mat')
load("average_FRAP_all_cells.mat")
%Params [Kon koff 0 kterm kinit 0]
plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,lower(end:-1:1)],color),'EdgeColor',color);
figure
hold on
plot_variance(time, avg_FRAP-std_FRAP, avg_FRAP+std_FRAP, [0.7529    0.8706    0.9216]);
plot(time,avg_FRAP, 'Color',[0.1843    0.4745    0.6000],'LineWidth',1.5)
plot(time, ezrdm_model.FRAP, '-','Color', 'k','LineWidth',1.5)
alpha(.25)
set(gca, 'FontSize', 24)
xlabel('time (s)', 'FontSize', 24)
ylabel('fraction RNAP recovered', 'FontSize', 24)
axis([0 600 0 1.5])

%% Find example images

%% Figure 1F Kinetic Model

%% Figure 2a Mutant Pie Chart

%% Figure 2b Dwell Time State 1 WT vs Mutant vs Chlor
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_vbSPT_values.mat')
%sample sizes and number of cells from powerpoint summary in google drive,
%confirmed they matched hmm pngs listed in the comments. Chlor was the only
%number of cells not explicitly listed, so i counted cells in data folder
%that were not "moved"
%mutant HMM from D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_I1309A_PAmCherry_EZRDM\HMM\170131
%chlor HMM from D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_PAmCherry_EZRDM-chlor\HMM\171130_test
figure
hold on
bar(1, rpoC_ezrdm(1).dwelltime*1000, 'FaceColor',[ 0.4941    0.7569    0.8706])
bar(2, rpoC_I1309A(1).dwelltime*1000, 'FaceColor',[ 0.9020    0.9020    0.9020])
bar(3, rpoC_chlor(1).dwelltime*1000, 'FaceColor',[0.8549    0.9804    0.6941])
errorbar(1, rpoC_ezrdm(1).dwelltime*1000, rpoC_ezrdm(1).std_dwelltime*1000,'.','MarkerSize', 15, 'Color','k')
errorbar(2, rpoC_I1309A(1).dwelltime*1000, rpoC_I1309A(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
errorbar(3, rpoC_chlor(1).dwelltime*1000, rpoC_chlor(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
ylabel('State I Dwell Time (ms)', 'FontSize', 24)
set(gca, 'FontSize', 24, 'XTick', [])
axis([0.5 3.5 0 350])
xticks([1 2 3])
xticklabels({'WT', 'RNAP I1309A', 'Chloramphenicol'})
figure 
hold on
bar(1, rpoC_ezrdm(2).dwelltime*1000, 'FaceColor',[0.2275 0.5647 0.7098])
errorbar(1, rpoC_ezrdm(2).dwelltime*1000, rpoC_ezrdm(2).std_dwelltime*1000,'.','MarkerSize', 15, 'Color','k')
bar(2, rpoC_I1309A(2).dwelltime*1000, 'FaceColor',[ 0.5020    0.5020    0.5020])
errorbar(2, rpoC_I1309A(2).dwelltime*1000, rpoC_I1309A(2).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
bar(3, rpoC_chlor(2).dwelltime*1000, 'FaceColor',[0.4667    0.6745    0.1882])
errorbar(3, rpoC_chlor(2).dwelltime*1000, rpoC_chlor(2).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
ylabel('State II Dwell Time (ms)', 'FontSize', 24)
set(gca, 'FontSize', 24, 'XTick', [])
axis([0.5 3.5 0 120])
xticks([1 2 3])
xticklabels({'WT', 'RNAP I1309A', 'Chloramphenicol'})

%% Figure 2D MSD plots for Chlor and mutant State 2
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_kusumi_results.mat')
figure
hold on
plot(rpoC_ezrdm(1).kusumi_25frames.xx, rpoC_ezrdm(1).kusumi_25frames.F, '-', 'Color', [0.2275 0.5647 0.7098], 'LineWidth', 1)
plot(rpoC_I1309A(1).kusumi_25frames.xx, rpoC_I1309A(1).kusumi_25frames.F, '-', 'Color',[ 0.5020    0.5020    0.5020], 'LineWidth', 1)
plot(rpoC_chlor(1).kusumi_25frames.xx, rpoC_chlor(1).kusumi_25frames.F, '-', 'Color', [0.4667    0.6745    0.1882], 'LineWidth', 1)
errorbar(rpoC_ezrdm(1).kusumi_25frames.D_all(1:12, 1), rpoC_ezrdm(1).kusumi_25frames.D_all(1:12,2), rpoC_ezrdm(1).kusumi_25frames.D_all(1:12,3), '.','MarkerSize', 15, 'Color', [0.2275 0.5647 0.7098])
errorbar(rpoC_I1309A(1).kusumi_25frames.D_all(1:12, 1), rpoC_I1309A(1).kusumi_25frames.D_all(1:12, 2), rpoC_I1309A(1).kusumi_25frames.D_all(1:12, 3), '.', 'MarkerSize', 15, 'Color', [ 0.5020    0.5020    0.5020])
errorbar(rpoC_chlor(1).kusumi_25frames.D_all(1:12, 1), rpoC_chlor(1).kusumi_25frames.D_all(1:12, 2), rpoC_chlor(1).kusumi_25frames.D_all(1:12, 3),'.', 'MarkerSize', 15, 'Color', [0.4667    0.6745    0.1882])
axis([0 12.2*0.00674 0 0.01])
legend("WT", 'Mutant','Chloramphenicol', 'EdgeColor', [1 1 1], 'FontSize', 24,'Location','Northwest')
xlabel('time (s)', 'FontSize', 24)
ylabel('State I MSD (\mum^2)', 'FontSize', 24)
set(gca, 'FontSize', 24)

figure
hold on
plot(rpoC_ezrdm(2).kusumi_25frames.xx, rpoC_ezrdm(2).kusumi_25frames.F, '-', 'Color', [0.2275 0.5647 0.7098], 'LineWidth', 1)
plot(rpoC_I1309A(2).kusumi_15frames.xx, rpoC_I1309A(2).kusumi_15frames.F, '-', 'Color',[ 0.5020    0.5020    0.5020] , 'LineWidth', 1)
plot(rpoC_chlor(2).kusumi_25frames.xx, rpoC_chlor(2).kusumi_25frames.F, '-', 'Color', [0.4667    0.6745    0.1882], 'LineWidth', 1)
errorbar(rpoC_ezrdm(2).kusumi_25frames.D_all(1:12, 1), rpoC_ezrdm(2).kusumi_25frames.D_all(1:12,2), rpoC_ezrdm(2).kusumi_25frames.D_all(1:12,3), '.','MarkerSize', 15, 'Color', [0.2275 0.5647 0.7098])
errorbar(rpoC_I1309A(2).kusumi_15frames.D_all(1:12, 1), rpoC_I1309A(2).kusumi_15frames.D_all(1:12, 2), rpoC_I1309A(2).kusumi_15frames.D_all(1:12, 3), '.', 'MarkerSize', 15, 'Color', [0.502 0.502 0.502])
errorbar(rpoC_chlor(2).kusumi_25frames.D_all(1:12, 1), rpoC_chlor(2).kusumi_25frames.D_all(1:12, 2), rpoC_chlor(2).kusumi_25frames.D_all(1:12, 3),'.', 'MarkerSize', 15, 'Color',  [0.4667 0.6745 0.1882])
axis([0 12.2*0.00674 0 0.05])
%legend("WT", 'Mutant','Chloramphenicol', 'EdgeColor', [1 1 1], 'FontSize', 24,'Location','Northwest')
xlabel('time (s)', 'FontSize', 24)
ylabel('State II MSD (\mum^2)', 'FontSize', 24)
set(gca, 'FontSize', 24)
%% Figure 2E Chlor and mutant FRAP curves
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM\new_data_collection_method\Analysis\updated_analysis')
load('best_ezrdm_model.mat')
load("average_FRAP_all_cells.mat")
WTavg = avg_FRAP; WTmodel = ezrdm_model.FRAP;
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\average_FRAP_Chlor.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\bestModel6chlor.mat')
chloravg = avg_FRAP; chlormodel = best.FRAP;
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\average_FRAP_mutant.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\model2mutant.mat')
mutavg = avg_FRAP; mutmodel = model2.FRAP.FRAP;
figure
hold on
%plot(time,WTavg, 'Color',[0.1843    0.4745    0.6000],'LineWidth',1.5)
%plot(time,WTmodel, '-','Color', [0.25 0.25 0.25],'LineWidth',1.5)
plot(time,WTavg, 'Color',[0.5176    0.6510    0.7098],'LineWidth',1.5)
plot(time,WTmodel, '-','Color', [0.45 0.45 0.45],'LineWidth',1.5)
plot(time,chloravg, 'Color',[0.3137    0.4784    0.0941],'LineWidth',1.5)
plot(time, chlormodel, '-','Color', [0.25 0.25 0.25],'LineWidth',1.5)
alpha(.25)
set(gca, 'FontSize', 24)
xlabel('time (s)', 'FontSize', 24)
ylabel('fraction RNAP recovered', 'FontSize', 24)
axis([0 600 0 1.2])
figure
hold on
%plot(time,WTavg, 'Color',[0.1843    0.4745    0.6000],'LineWidth',1.5)
%plot(time,WTmodel, '-','Color', [0.25 0.25 0.25],'LineWidth',1.5)
plot(time,WTavg, 'Color',[0.5176    0.6510    0.7098],'LineWidth',1.5)
plot(time,WTmodel, '-','Color', [0.45 0.45 0.45],'LineWidth',1.5)
plot(time,mutavg, 'Color',[ 0.5020    0.5020    0.5020] ,'LineWidth',1.5)
plot(time, mutmodel, '-','Color', [0.15 0.15 0.15],'LineWidth',1.5)
alpha(.25)
set(gca, 'FontSize', 24)
xlabel('time (s)', 'FontSize', 24)
ylabel('fraction RNAP recovered', 'FontSize', 24)
axis([0 600 0 1.2])

%% Figure 3a M9 Pie Chart
%% Figure 3b M9 FRAP Curve
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM\new_data_collection_method\Analysis\updated_analysis')
load('best_ezrdm_model.mat')
load("average_FRAP_all_cells.mat")
WTavg = avg_FRAP; WTmodel = ezrdm_model.FRAP;
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_M9\Analysis\160229_analysis')
load('average_FRAP_all_cells.mat')
load('model6.mat')
M9avg = average_FRAP; M9model = model6(2).FRAP;
figure
hold on
%plot(time,WTavg, 'Color',[0.1843    0.4745    0.6000],'LineWidth',1.5)
%plot(time,WTmodel, '-','Color', [0.25 0.25 0.25],'LineWidth',1.5)
plot(time,WTavg, 'Color',[0.5176    0.6510    0.7098],'LineWidth',1.5)
plot(time,WTmodel, '-','Color', [0.45 0.45 0.45],'LineWidth',1.5)
plot(time,M9avg, 'Color',[0.9882    0.3922    0.1373],'LineWidth',1.5)
plot(time, M9model, '-','Color', [0.25 0.25 0.25],'LineWidth',1.5)
alpha(.25)
set(gca, 'FontSize', 24)
xlabel('time (s)', 'FontSize', 24)
ylabel('fraction RNAP recovered', 'FontSize', 24)
axis([0 600 0 1.2])
%% Figure 3c M9 Kinetic Scheme
%% Figure 3d AsiA Pie Chart
%% Figure 3e AsiA FRAP Curve
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM\new_data_collection_method\Analysis\updated_analysis')
load('best_ezrdm_model.mat')
load("average_FRAP_all_cells.mat")
WTavg = avg_FRAP;
WTmodel = ezrdm_model.FRAP;
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\AsiA\Analysis\FRAP_analysis_181121')
load('best_asia_model.mat')
load('average_FRAP_all_cells.mat')
Asiaavg = avg_FRAP;
Asiamodel = asia_model.FRAP;
figure
hold on
plot(time,WTavg, 'Color',[0.5176    0.6510    0.7098],'LineWidth',1.5)
plot(time,WTmodel, '-','Color', [0.45 0.45 0.45],'LineWidth',1.5)
plot(time,Asiaavg, 'Color',[0.4353    0.0706    0.5020],'LineWidth',1.5)
plot(time, Asiamodel, '-','Color', [0.25 0.25 0.25],'LineWidth',1.5)
alpha(.25)
set(gca, 'FontSize', 24)
xlabel('time (s)', 'FontSize', 24)
ylabel('fraction RNAP recovered', 'FontSize', 24)
axis([0 600 0 1.2])

%% Figure 3f AsiA Kinetic Scheme

%% Figure 4a Rifampicin Pie Chart

%% Figure 4b Rifampicin FRAP Curve
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM\new_data_collection_method\Analysis\updated_analysis')
load('best_ezrdm_model.mat')
load("average_FRAP_all_cells.mat")
WTavg = avg_FRAP;WTmodel = ezrdm_model.FRAP;
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM-Rif\Analysis\FRAP_analysis_181119')
load('average_FRAP_all_cells.mat')
load('best_rif_model.mat')
Rifavg = avg_FRAP;Rifmodel = rif_model.FRAP.FRAP;
figure
hold on
plot(time,WTavg, 'Color',[0.5176    0.6510    0.7098],'LineWidth',1.5)
plot(time,WTmodel, '-','Color', [0.45 0.45 0.45],'LineWidth',1.5)
plot(time,Rifavg, 'Color',[0.8392    0.1765    0.5843],'LineWidth',1.5)
plot(time, Rifmodel, '-','Color', [0.25 0.25 0.25],'LineWidth',1.5)
alpha(.25)
set(gca, 'FontSize', 24)
xlabel('time (s)', 'FontSize', 24)
ylabel('fraction RNAP recovered', 'FontSize', 24)

axis([0 600 0 1.2])

%% Figure 4c Rifampicin Kinetic Scheme

%% Figure 3C Rif and Chlor Dwell Times for State I
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_vbSPT_values.mat')
 
figure
hold on
bar(1, rpoC_rif(1).dwelltime*1000, 'FaceColor',[0.9608    0.5843    0.8157])
errorbar(1, rpoC_rif(1).dwelltime*1000, rpoC_rif(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
bar(2, rpoC_ezrdm(1).dwelltime*1000, 'FaceColor',[ 0.4941    0.7569    0.8706])
errorbar(2, rpoC_ezrdm(1).dwelltime*1000, rpoC_ezrdm(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
bar(3, rpoC_chlor(1).dwelltime*1000, 'FaceColor', [0.8549    0.9804    0.6941])
errorbar(3, rpoC_chlor(1).dwelltime*1000, rpoC_chlor(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
ylabel('State I Dwell Time (ms)', 'FontSize', 15, 'FontName', 'Helvetica')
set(gca, 'FontSize', 15, 'FontName', 'Helvetica', 'XTick', [])
axis([0.5 3.5 0 350])
xticks([1 2 3])
xticklabels({'Rifampicin','WT', 'Chloramphenicol'})

%% Figure 3D MSD plots for Chlor and Rif State 2
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_kusumi_results.mat')

figure
hold on
plot(rpoC_chlor(2).kusumi_25frames.xx, rpoC_chlor(2).kusumi_25frames.F, '-', 'Color', [0.4667    0.6745    0.1882], 'LineWidth', 1)
plot(rpoC_ezrdm(2).kusumi_25frames.xx, rpoC_ezrdm(2).kusumi_25frames.F, '-', 'Color', [0.2275 0.5647 0.7098], 'LineWidth', 1)
plot(rpoC_rif(2).kusumi_12frames.xx, rpoC_rif(2).kusumi_12frames.F, '-', 'Color', [0.7882    0.2510    0.5843], 'LineWidth', 1)
errorbar(rpoC_chlor(2).kusumi_25frames.D_all(1:12, 1), rpoC_chlor(2).kusumi_25frames.D_all(1:12,2), rpoC_chlor(2).kusumi_25frames.D_all(1:12,3), '.','MarkerSize', 15, 'Color', [0.4667    0.6745    0.1882])
errorbar(rpoC_ezrdm(2).kusumi_25frames.D_all(1:12, 1), rpoC_ezrdm(2).kusumi_25frames.D_all(1:12, 2), rpoC_ezrdm(2).kusumi_25frames.D_all(1:12, 3), '.', 'MarkerSize', 15, 'Color', [0.2275 0.5647 0.7098])
errorbar(rpoC_rif(2).kusumi_12frames.D_all(1:12, 1), rpoC_rif(2).kusumi_12frames.D_all(1:12, 2), rpoC_rif(2).kusumi_12frames.D_all(1:12, 3),'.', 'MarkerSize', 15, 'Color',  [0.7882    0.2510    0.5843])
axis([0 12.2*0.00674 0 0.05])
legend('Chloramphenicol', "WT", "Rifampicin", 'EdgeColor', [1 1 1], 'FontSize', 16,'Location','Northwest')
xlabel('time (s)', 'FontSize', 24)
ylabel('State II MSD (\mum^2)', 'FontSize', 24)
set(gca, 'FontSize', 24)

figure
hold on
plot(rpoC_chlor(2).kusumi_25frames.xx, rpoC_chlor(2).kusumi_25frames.F, '-', 'Color', [0.4667    0.6745    0.1882], 'LineWidth', 1)
plot(rpoC_ezrdm(2).kusumi_25frames.xx, rpoC_ezrdm(2).kusumi_25frames.F, '-', 'Color', [0.2275 0.5647 0.7098], 'LineWidth', 1)
plot(rpoC_rif(2).kusumi_12frames.xx, rpoC_rif(2).kusumi_12frames.F, '-', 'Color', [0.7882    0.2510    0.5843], 'LineWidth', 1)
errorbar(rpoC_chlor(2).kusumi_25frames.D_all(1:12, 1), rpoC_chlor(2).kusumi_25frames.D_all(1:12,2), rpoC_chlor(2).kusumi_25frames.D_all(1:12,3), '.','MarkerSize', 15, 'Color', [0.4667    0.6745    0.1882])
errorbar(rpoC_ezrdm(2).kusumi_25frames.D_all(1:12, 1), rpoC_ezrdm(2).kusumi_25frames.D_all(1:12, 2), rpoC_ezrdm(2).kusumi_25frames.D_all(1:12, 3), '.', 'MarkerSize', 15, 'Color', [0.2275 0.5647 0.7098])
errorbar(rpoC_rif(2).kusumi_12frames.D_all(1:12, 1), rpoC_rif(2).kusumi_12frames.D_all(1:12, 2), rpoC_rif(2).kusumi_12frames.D_all(1:12, 3),'.', 'MarkerSize', 15, 'Color',  [0.7882    0.2510    0.5843])
axis([0 12.2*0.00674 0 0.05])
legend('Chloramphenicol', "WT", "Rifampicin", 'EdgeColor', [1 1 1], 'FontSize', 16,'Location','Northwest')
xlabel('time (s)', 'FontSize', 24)
ylabel('State II MSD (\mum^2)', 'FontSize', 24)
set(gca, 'FontSize', 24)


%% Figure 3a Asia and M9 Pie Chart


%% Figure 3d Asia and M9 FRAP data
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM\new_data_collection_method\Analysis\updated_analysis')
load('best_ezrdm_model.mat')
load("average_FRAP_all_cells.mat")
WTavg = avg_FRAP;WTmodel = ezrdm_model.FRAP;
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\AsiA\Analysis\FRAP_analysis_181121')
load('best_asia_model.mat')
load('average_FRAP_all_cells.mat')
Asiaavg = avg_FRAP;Asiamodel = asia_model.FRAP;
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_M9\Analysis\160229_analysis')
load('average_FRAP_all_cells.mat')
load('model6.mat')
M9avg = average_FRAP;M9model = model6(2).FRAP;
figure
hold on
plot(time,M9avg, 'Color',[0.9882    0.3922    0.1373],'LineWidth',1.5)
plot(time,WTavg, 'Color',[0.1843    0.4745    0.6000],'LineWidth',1.5)
plot(time,Asiaavg, 'Color',[0.4353    0.0706    0.5020],'LineWidth',1.5)
plot(time,WTmodel, '-','Color', [0.25 0.25 0.25],'LineWidth',1.5)
plot(time, M9model, '-','Color', [0.25 0.25 0.25],'LineWidth',1.5)
plot(time, Asiamodel, '-','Color', [0.25 0.25 0.25],'LineWidth',1.5)
alpha(.25)
set(gca, 'FontSize', 24)
xlabel('time (s)', 'FontSize', 24)
ylabel('fraction RNAP recovered', 'FontSize', 24)
legend('M9', 'WT-EZRDM', ['AsiA' char(10) 'Overexpression'], 'EdgeColor', [1 1 1], 'FontSize', 14,'Location','Southeast')

axis([0 600 0 1.25])

%% Supplemental Figures

%% SF1a population fitsclear
%cd('/Users/kbettridge/Dropbox/Data/RNAP_PAmCherry_EZRDM/5ms_continuous_data/current_HMM_analysis/redone_HMM_160517')

cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_PAmCherry_EZRDM\5ms_continuous_data\current_HMM_analysis\redone_HMM_160517')
load('rpoCp_tot_rotated_hmm.mat')

cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\MATLAB\')
traj_results1 = SingleMolD(trajectories, 1, 0.001, 0.00674);
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_PAmCherry_EZRDM\5ms_continuous_data\current_HMM_analysis\redone_HMM_160517')

n = 1; % number of steps
Displ_n = [];
Dx = [];
Dy = [];
D2 = [];

traj_results = traj_results1;

for ii = 1: length(traj_results)
    
    Displ_n = cat(1, Displ_n, traj_results{ii}.Displacements{n});
    Dx = cat (1, Dx, traj_results{ii}.Dx);
    Dy = cat (1, Dy, traj_results{ii}.Dy);
    D2 = cat (1, D2, traj_results{ii}.D2);
end

D(n).Displ_n = Displ_n;
D(n).Dx = Dx;
D(n).Dy = Dy;
D(n).D2 = D2;

n = 1;
Displ = D(n).Displ_n;
[y x] = hist(Displ, 100);
y = y/trapz(x, y);

% single Gaussian fit
mu = mean(Displ(:, 1));
sigma = std(Displ(:, 1));
y_fit = normpdf(x, mu, sigma);
sing_resid = y(:) - y_fit(:);

% three Gaussian fit
Displ_1_3G_fit = fit(x, y, 'gauss3', 'StartPoint', [5 0 0.03 1 0 0.03 1 0 0.03]);
f = Displ_1_3G_fit;
para = coeffvalues(f);

sigma1 = para(3)/sqrt(2);
sigma2 = para(6)/sqrt(2);
sigma3 = para(9)/sqrt(2);

t = 0.00674;

p1 = para(1)*para(3)*sqrt(pi);
D1 = para(3)^2/(4*n*t);

p2 = para(4)*para(6)*sqrt(pi);
D2 = para(6)^2/(4*n*t);

p3 = para(7)*para(9)*sqrt(pi);
D3 = para(9)^2/(4*n*t);

Displ_1_3G_para = [D1 p1 D2 p2 D3 p3];

g1 = para(1)*exp(-((x-para(2))/para(3)).^2);
g2 = para(4)*exp(-((x-para(5))/para(6)).^2);
g3 = para(7)*exp(-((x-para(8))/para(9)).^2);
g = g1 + g2 + g3;

tri_resid = y(:) - g(:);

figure
subplot(4,1,1:3)
hold on
bar(x, y, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.5 0.5 0.5])
%plot(x, y, '.k', 'MarkerSize', 24)
plot(x, y_fit, '-r', 'LineWidth', 1.5)
plot(x, g, '-k', 'LineWidth', 1.5)
set(gca, 'FontSize', 18, 'FontName', 'Helvetica', 'XTick', [-0.2 -0.1 0 0.1 0.2], 'XTickLabel', {})
axis([-0.2 0.2 0 8])

subplot(4, 1, 4)
hold on
plot(x, sing_resid, '--r')
plot(x, tri_resid, '--k')
set(gca, 'FontSize', 18, 'FontName', 'Helvetica')
maxy = max(abs(sing_resid))+0.5;
axis([-0.2 0.2 -maxy maxy])
box on
%legend('displacements', 'single Gaussian fit', 'triple Gaussian fit')
hold off


%% Fixed Cell vs WT Survival Times
% clear
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_PAmCherry_fixed\150127_SMT_rpoC-PAmCherry_fixed_cells_5ms')
load('rpoCp_fixed_e5_150127.mat')
cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_PAmCherry_fixed\150128_SMT_rpoC-PAmCherry_fixed_cells_5ms')
load('rpoCp_fixed_150128.mat')
total_rpoCp_fixed = [rpoCp_fixed_150128,rpoCp_fixed_e5_150127];
trajlength = zeros(length(total_rpoCp_fixed),1);
traj2 = [];
j = 1;
for i = 1:length(total_rpoCp_fixed)
   trajlength(i) = length(total_rpoCp_fixed(i).frames);
   if length(total_rpoCp_fixed(i).frames) > 1
       traj2(j) = length(total_rpoCp_fixed(i).frames);
       j = j+1;
   end
end

%cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\remadeStructures\')
%save('total_rpoCp_fixed_traj.mat','total_rpoCp_fixed')

cd('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_PAmCherry_EZRDM\5ms_continuous_data\')
load('all_combined_consc_traj_wcelldata.mat')

trajlengthEzrdm = zeros(length(s_consc),1);
traj2Ezrdm = [];
j = 1;
for i = 1:length(s_consc)
   trajlengthEzrdm(i) = length(s_consc(i).frames);
    if length(s_consc(i).frames) > 1
       traj2Ezrdm(j) = length(s_consc(i).frames);
       j = j+1;
   end
end

figure
hold on
h1= histogram(trajlength,BinWidth=1, FaceColor="red", FaceAlpha=0.75)
h2 = histogram(trajlengthEzrdm, BinWidth=1, FaceColor="cyan", FaceAlpha=0.5)
set(gca, "YScale", "log")
legend("Fixed Cell Trajectories", "WT-EZRDM Trajectories",'EdgeColor', [1 1 1],"FontSize",24)
title("Trajectory Length Distribution","FontSize",24)
ylabel("Number of Trajectories","FontSize",24)
xlabel("Trajectory Length (Frames)","FontSize",24)
%%
figure
hold on
h1values= h1.Values/sum(h1.Values);
h2values= h2.Values/sum(h2.Values);
bar(h1values, FaceColor="red", FaceAlpha=1)
bar(h2values,  FaceColor="cyan", FaceAlpha=1)
set(gca, "YScale", "log")
legend("Fixed Cell Trajectories", "WT-EZRDM Trajectories",'EdgeColor', [1 1 1],"FontSize",26)
title("Trajectory Length Distribution", "FontSize",30)
ylabel("Fraction of Trajectories","FontSize",30)
xlabel("Trajectory Length (Frames)","FontSize",30)
%%
figure 
hold on
bar(h1values, FaceColor="red", FaceAlpha=0.75)
bar(h2values,  FaceColor="cyan", FaceAlpha=0.5)
set(gca, "YScale", "linear")
legend("Fixed Cell Trajectories", "WT-EZRDM Trajectories",'EdgeColor', [1 1 1],"FontSize",24)
title("Trajectory Length Distribution", "FontSize",24)
ylabel("Fraction of Trajectories","FontSize",24)
xlabel("Trajectory Length (Frames)","FontSize",24)
%% Figure S4 Asia and M9 Dwell Time Comparison
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_vbSPT_values.mat')
%M9 HMM from D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_PAmCherry_M9\5ms_data\HMM_160603
%AsiA HMM from
figure
hold on
bar(1, rpoC_m9(1).dwelltime*1000, 'FaceColor', [0.9686    0.6392    0.4941])
bar(2, rpoC_ezrdm(1).dwelltime*1000, 'FaceColor',[ 0.4941    0.7569    0.8706])
bar(3, rpoC_asia(1).dwelltime*1000, 'FaceColor', [0.8784    0.7098    0.9098])
bar(5, rpoC_m9(2).dwelltime*1000, 'FaceColor', [1.0000    0.3961    0.1412])
bar(6, rpoC_ezrdm(2).dwelltime*1000, 'FaceColor',[ 0.2275    0.5647    0.7098])
bar(7, rpoC_asia(2).dwelltime*1000, 'FaceColor',[0.4941    0.1843    0.5569])
errorbar(1, rpoC_m9(1).dwelltime*1000, rpoC_m9(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
errorbar(2, rpoC_ezrdm(1).dwelltime*1000, rpoC_ezrdm(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
errorbar(3, rpoC_asia(1).dwelltime*1000, rpoC_asia(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
errorbar(5, rpoC_m9(2).dwelltime*1000, rpoC_m9(2).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
errorbar(6, rpoC_ezrdm(2).dwelltime*1000, rpoC_ezrdm(2).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
errorbar(7, rpoC_asia(2).dwelltime*1000, rpoC_asia(2).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
ylabel('Dwell Time (ms)', 'FontSize', 24)
set(gca, 'FontSize', 24, 'XTick', [])
legend('M9','WT-EZRDM', ['AsiA' char(10) 'Overexpression'], 'EdgeColor', [1 1 1], 'FontSize', 24,'Location','Northeast')
axis([0.5 7.5 0 350])
xticks([2 6])
xticklabels({'State I','State II'})

%% Figure S4 diffusive domain comparison
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_kusumi_results.mat')
figure
hold on
bar(1, rpoC_m9(1).kusumi_25frames.L*1000, 'FaceColor', [0.9686    0.6392    0.4941])
errorbar(1, rpoC_m9(1).kusumi_25frames.L*1000, rpoC_m9(1).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
bar(2, rpoC_ezrdm(1).kusumi_25frames.L*1000, 'FaceColor',[ 0.4941    0.7569    0.8706])
errorbar(2, rpoC_ezrdm(1).kusumi_25frames.L*1000, rpoC_ezrdm(1).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
bar(3, rpoC_asia(1).kusumi_25frames.L*1000, 'FaceColor', [0.8784    0.7098    0.9098])
errorbar(3, rpoC_asia(1).kusumi_25frames.L*1000, rpoC_asia(1).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
bar(5, rpoC_m9(2).kusumi_25frames.L*1000, 'FaceColor', [1.0000    0.3961    0.1412])
errorbar(5, rpoC_m9(2).kusumi_25frames.L*1000, rpoC_m9(2).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
bar(6, rpoC_ezrdm(2).kusumi_25frames.L*1000, 'FaceColor',[ 0.2275    0.5647    0.7098])
errorbar(6, rpoC_ezrdm(2).kusumi_25frames.L*1000, rpoC_ezrdm(2).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
bar(7, rpoC_asia(2).kusumi_20frames.L*1000, 'FaceColor',[0.4941    0.1843    0.5569])
errorbar(7, rpoC_asia(2).kusumi_20frames.L*1000, rpoC_asia(2).kusumi_20frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
ylabel('Diffusive Domain Size (nm)', 'FontSize', 24)
set(gca, 'FontSize', 24,'XTick', [])
axis([0.5 7.5 0 350])
xticks([2 6])
xticklabels({'State I','State II'})

%% Rif Supplemental Figures
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_vbSPT_values.mat')
figure
hold on
bar(1, rpoC_ezrdm(1).dwelltime*1000, 'FaceColor',[ 0.4941    0.7569    0.8706])
bar(2, rpoC_rif(1).dwelltime*1000, 'FaceColor',[0.9608    0.5843    0.8157])
errorbar(1, rpoC_ezrdm(1).dwelltime*1000, rpoC_ezrdm(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
errorbar(2, rpoC_rif(1).dwelltime*1000, rpoC_rif(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
bar(4, rpoC_ezrdm(2).dwelltime*1000, 'FaceColor',[0.2275 0.5647 0.7098])
errorbar(4, rpoC_ezrdm(2).dwelltime*1000, rpoC_ezrdm(2).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
bar(5, rpoC_rif(2).dwelltime*1000, 'FaceColor',[0.7882    0.2510    0.5843])
errorbar(5, rpoC_rif(2).dwelltime*1000, rpoC_rif(2).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
ylabel('Dwell Time (ms)', 'FontSize', 24)
set(gca, 'FontSize', 24, 'XTick', [])
axis([0.5 5.5 0 350])
xticks([1.5 4.5])
xticklabels({'State I', 'State II'})
legend("WT", "Rifampicin",'EdgeColor', [1 1 1],'FontSize', 24)
clear

%%
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_kusumi_results.mat')
%%
figure
hold on
plot(rpoC_ezrdm(1).kusumi_25frames.xx, rpoC_ezrdm(1).kusumi_25frames.F, '-', 'Color', [ 0.4941    0.7569    0.8706], 'LineWidth', 1)
plot(rpoC_rif(1).kusumi_15frames.xx, rpoC_rif(1).kusumi_15frames.F, '-', 'Color', [0.9608    0.5843    0.8157], 'LineWidth', 1)
errorbar(rpoC_ezrdm(1).kusumi_25frames.D_all(1:12, 1), rpoC_ezrdm(1).kusumi_25frames.D_all(1:12, 2), rpoC_ezrdm(1).kusumi_25frames.D_all(1:12, 3), '.', 'MarkerSize', 15, 'Color', [ 0.4941    0.7569    0.8706])
errorbar(rpoC_rif(1).kusumi_15frames.D_all(1:12, 1), rpoC_rif(1).kusumi_15frames.D_all(1:12, 2), rpoC_rif(1).kusumi_15frames.D_all(1:12, 3),'.', 'MarkerSize', 15, 'Color', [0.9608    0.5843    0.8157])
axis([0 12.2*0.00674 0 0.01])
legend( "WT", "Rifampicin", 'EdgeColor', [1 1 1], 'FontSize', 24,'Location','Northwest')
xlabel('time (s)', 'FontSize', 24)
ylabel('State I MSD (\mum^2)', 'FontSize', 24)
set(gca, 'FontSize', 24)

figure
hold on
plot(rpoC_ezrdm(2).kusumi_25frames.xx, rpoC_ezrdm(2).kusumi_25frames.F, '-', 'Color', [0.2275 0.5647 0.7098], 'LineWidth', 1)
plot(rpoC_rif(2).kusumi_12frames.xx, rpoC_rif(2).kusumi_12frames.F, '-', 'Color', [0.7882    0.2510    0.5843], 'LineWidth', 1)
errorbar(rpoC_ezrdm(2).kusumi_25frames.D_all(1:12, 1), rpoC_ezrdm(2).kusumi_25frames.D_all(1:12, 2), rpoC_ezrdm(2).kusumi_25frames.D_all(1:12, 3), '.', 'MarkerSize', 15, 'Color', [0.2275 0.5647 0.7098])
errorbar(rpoC_rif(2).kusumi_12frames.D_all(1:12, 1), rpoC_rif(2).kusumi_12frames.D_all(1:12, 2), rpoC_rif(2).kusumi_12frames.D_all(1:12, 3),'.', 'MarkerSize', 15, 'Color',  [0.7882    0.2510    0.5843])
axis([0 12.2*0.00674 0 0.06])
legend( "WT", "Rifampicin", 'EdgeColor', [1 1 1], 'FontSize', 24,'Location','Northwest')
xlabel('time (s)', 'FontSize', 24)
ylabel('State II MSD (\mum^2)', 'FontSize', 24)
set(gca, 'FontSize', 24)


%% Figure 4 C diffusive domain comparison
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_kusumi_results.mat')
figure
hold on
bar(1, rpoC_m9(1).kusumi_25frames.L*1000, 'FaceColor', [0.9686    0.6392    0.4941])
errorbar(1, rpoC_m9(1).kusumi_25frames.L*1000, rpoC_m9(1).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
bar(2, rpoC_ezrdm(1).kusumi_25frames.L*1000, 'FaceColor',[ 0.4941    0.7569    0.8706])
errorbar(2, rpoC_ezrdm(1).kusumi_25frames.L*1000, rpoC_ezrdm(1).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
bar(4, rpoC_m9(2).kusumi_25frames.L*1000, 'FaceColor', [1.0000    0.3961    0.1412])
errorbar(4, rpoC_m9(2).kusumi_25frames.L*1000, rpoC_m9(2).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
bar(5, rpoC_ezrdm(2).kusumi_25frames.L*1000, 'FaceColor',[ 0.2275    0.5647    0.7098])
errorbar(5, rpoC_ezrdm(2).kusumi_25frames.L*1000, rpoC_ezrdm(2).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
ylabel('Diffusive Domain Size (nm)', 'FontSize', 15, 'FontName', 'Helvetica')
set(gca, 'FontSize', 15, 'FontName', 'Helvetica', 'XTick', [])
axis([0.5 5.5 0 350])
xticks([1.5 4.5])
xticklabels({'State I','State II'})
%% Asia only
clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_vbSPT_values.mat')
figure
hold on
bar(1, rpoC_ezrdm(1).dwelltime*1000, 'FaceColor',[ 0.4941    0.7569    0.8706])
bar(2, rpoC_asia(1).dwelltime*1000, 'FaceColor', [0.8784    0.7098    0.9098])
bar(4, rpoC_ezrdm(2).dwelltime*1000, 'FaceColor',[ 0.2275    0.5647    0.7098])
bar(5, rpoC_asia(2).dwelltime*1000, 'FaceColor',[0.4941    0.1843    0.5569])
errorbar(1, rpoC_ezrdm(1).dwelltime*1000, rpoC_ezrdm(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
errorbar(2, rpoC_asia(1).dwelltime*1000, rpoC_asia(1).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
errorbar(4, rpoC_ezrdm(2).dwelltime*1000, rpoC_ezrdm(2).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
errorbar(5, rpoC_asia(2).dwelltime*1000, rpoC_asia(2).std_dwelltime*1000, '.','MarkerSize', 15, 'Color','k')
ylabel('Dwell Time (ms)', 'FontSize', 15, 'FontName', 'Helvetica')
set(gca, 'FontSize', 15, 'FontName', 'Helvetica', 'XTick', [])
legend('WT-EZRDM', ['AsiA' char(10) 'Overexpression'], 'EdgeColor', [1 1 1], 'FontSize',24,'Location','Northeast')
axis([0.5 5.5 0 350])
xticks([1.5 4.5])
xticklabels({'State I','State II'})

clear
cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\Scripts\StructuresFromGoogleDrive')
load('all_kusumi_results.mat')
figure
hold on
bar(1, rpoC_ezrdm(1).kusumi_25frames.L*1000, 'FaceColor',[ 0.4941    0.7569    0.8706])
errorbar(1, rpoC_ezrdm(1).kusumi_25frames.L*1000, rpoC_ezrdm(1).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
bar(2, rpoC_asia(1).kusumi_25frames.L*1000, 'FaceColor', [0.8784    0.7098    0.9098])
errorbar(2, rpoC_asia(1).kusumi_25frames.L*1000, rpoC_asia(1).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
bar(4, rpoC_ezrdm(2).kusumi_25frames.L*1000, 'FaceColor',[ 0.2275    0.5647    0.7098])
errorbar(4, rpoC_ezrdm(2).kusumi_25frames.L*1000, rpoC_ezrdm(2).kusumi_25frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
bar(5, rpoC_asia(2).kusumi_20frames.L*1000, 'FaceColor',[0.4941    0.1843    0.5569])
errorbar(5, rpoC_asia(2).kusumi_20frames.L*1000, rpoC_asia(2).kusumi_20frames.std_L*1000, '.','MarkerSize', 15, 'Color','k')
ylabel('Diffusive Domain Size (nm)', 'FontSize', 15, 'FontName', 'Helvetica')
set(gca, 'FontSize', 15, 'FontName', 'Helvetica', 'XTick', [])
axis([0.5 5.5 0 350])
xticks([1.5 4.5])
xticklabels({'State I','State II'})
