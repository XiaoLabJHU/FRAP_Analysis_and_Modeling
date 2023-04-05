clear
%load ('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\total_FRAP_all_cells_Chlor.mat')
%load('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM\new_data_collection_method\Analysis\updated_analysis\total_FRAP_all_cells.mat')
%load('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_M9\total_FRAP_all_cells_FH.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\total_FRAP_all_cells_Chlor.mat');
%load('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM-Rif\Analysis\FRAP_analysis_181119\total_FRAP_all_cells.mat')
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\totalFRAP_all_cells_mutant.mat');
%load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\AsiA\Analysis\FRAP_analysis_181121\total_FRAP_all_cells.mat')

avg_FRAP = zeros(100,601);
std_FRAP = zeros(100,601);
bootPclust = zeros(100,1);
time = 1:601;
rng(1);
parfor i= 1:100
    bootFRAP = randsample(total_FRAP_all_cells, 40, false);
    bootFRAP(41:45) = randsample(bootFRAP,5,true);
    avg_FRAP(i,:) = mean(vertcat(bootFRAP.ft2));
    std_FRAP(i,:) = std(vertcat(bootFRAP.ft2));
    bootFRAP_fit = fit(time(:), avg_FRAP(i,:)','1-a*exp(-koff1*x)-b*exp(-koff2*x)-C','StartPoint',[0.45 0.2 0.012 0.5 0.15],'Lower',[0 0 0 0 0]);
    values = coeffvalues(bootFRAP_fit);
    bootPclust(i) = values(1);
end

cd('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\')
save('average_FRAP_bootChlor.mat','time','avg_FRAP','std_FRAP',"bootPclust")

%%
figure
hold on
for i = 1:100
    plot(time, avg_FRAP(i,:))
end
%%
clear
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\average_FRAP_bootChlor.mat');

addpath('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\MATLAB\')
addpath('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\MATLAB\FRAP_simulation')
tic
%f = waitbar(0,'Initiliazing...','Name','LSQNONLIN takes forever',...
%    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
nrounds=0;
known = [NaN NaN 0 0 0 0];
N = 1;

for ii = 1:100
    flag = 0;
     while flag == 0
        nrounds = nrounds + 1;
        try
            test = FRAP_lsm(time, avg_FRAP(ii,:), std_FRAP(ii,:), 0.25, known);
            populationTest = test.FRAP.Pelong >0 && test.FRAP.Pfree >0;
            paramTest = test.params(1) > 0 && test.params(2) >0;
            endFRAP = mean(test.FRAP.FRAP(end-40:end));
            if populationTest == 1 && paramTest == 1 && endFRAP <= 1
                flag = 1;
                model2 = test;
            end
        catch
        end
     end
    bootmodel2(ii) = test;
end
%delete(f);
toc
save("bootModel2Chlor.mat","bootmodel2")
%%
allModel2Params = horzcat(bootmodel2.params);

allModel2Pfree = zeros(100,1);
allModel2Pelong = zeros(100,1);
for i = 1:100
    allModel2Pfree(i) = bootmodel2(i).FRAP.Pfree;
    allModel2Pelong(i) = bootmodel2(i).FRAP.Pelong;
end
figure
hold on
subplot(2,2,1)
histogram(allModel2Params(1,:))
title('k_o_n')
subplot(2,2,2)
histogram(allModel2Params(2,:))
title('k_o_f_f')
subplot(2,2,3)
histogram(allModel2Pfree)
title('P_f_r_e_e')
subplot(2,2,4)
histogram(allModel2Pelong)
title('P_b_o_u_n_d')
fprintf('k_on mean: %1.4f, std: %1.6f\n', mean(allModel2Params(1,:)),std(allModel2Params(1,:)))
fprintf("k_off mean: %1.4f, std: %1.6f\n", mean(allModel2Params(2,:)),std(allModel2Params(2,:)))
fprintf("P_free mean: %1.4f, std: %1.6f\n", mean(allModel2Pfree),std(allModel2Pfree))
fprintf("P_bound mean: %1.4f, std: %1.6f\n", mean(allModel2Pelong),std(allModel2Pelong))
