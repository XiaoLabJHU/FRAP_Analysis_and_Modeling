%% Create Bootstrapped Data
%load ('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\total_FRAP_all_cells_Chlor.mat')
%load('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM\new_data_collection_method\Analysis\updated_analysis\total_FRAP_all_cells.mat')
%load('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_M9\total_FRAP_all_cells_FH.mat')
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\total_FRAP_all_cells_Chlor.mat');
%load('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM-Rif\Analysis\FRAP_analysis_181119\total_FRAP_all_cells.mat')
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\totalFRAP_all_cells_mutant.mat');
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\AsiA\Analysis\FRAP_analysis_181121\total_FRAP_all_cells.mat')

avg_FRAP = zeros(100,601);
std_FRAP = zeros(100,601);
bootPclust = zeros(100,1);
time = 1:601;
rng(1);
for i= 1:100
    bootFRAP = randsample(total_FRAP_all_cells, 41, false);
    bootFRAP(42:46) = randsample(bootFRAP,5,true);
    avg_FRAP(i,:) = mean(vertcat(bootFRAP.ft2));
    std_FRAP(i,:) = std(vertcat(bootFRAP.ft2));
    bootFRAP_fit = fit(time(:), avg_FRAP(i,:)','1-a*exp(-koff1*x)-b*exp(-koff2*x)-C','StartPoint',[0.45 0.2 0.012 0.5 0.15],'Lower',[0 0 0 0 0]);
    values = coeffvalues(bootFRAP_fit);
    bootPclust(i) = values(1);
end

cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\')
save('average_FRAP_bootM92.mat','time','avg_FRAP','std_FRAP',"bootPclust")