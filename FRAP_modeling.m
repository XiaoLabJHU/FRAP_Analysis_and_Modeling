% FRAP modeling 
%% One-Step Modeling
clear

load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\average_FRAP_mutant.mat');
% file includes 1x601 double avg_FRAP (normalized FRAP signal)
% 1x601 double time (1:601 here)
% 1x601 std_FRAP (used to weight fitting)
addpath('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\MATLAB\')
addpath('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\MATLAB\FRAP_simulation')

model1 = FRAP_lsm(time, avg_FRAP, std_FRAP, 0.25, [NaN NaN 0 0 0 0]);

%% Two-Step Modeling -- assume k3 = 0, k6 = 0, float all other parameters, input Ptranscribing
% essentially modified Singer model
% Rfree <---Rpromoter (k2)
% Rfree --->Rpromoter (k1)
% Rpromoter ---> Relongation (k5)
% Relongation --> Rfree (k4)
clear
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\average_FRAP_mutant.mat');
% file includes 1x601 double avg_FRAP (normalized FRAP signal)
% 1x601 double time (1:601 here)
% 1x601 std_FRAP (used to weight fitting)
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\FRAP_fit_mutant.mat');
%cfit object from code below
%FRAP_fit = fit(time(:), avg_FRAP(:), ...
   %'1-a*exp(-koff1*x)-b*exp(-koff2*x)-C','StartPoint',...
   %[0.45 0.2 0.012 0.5 0.15],'Lower',[0 0 0 0 0])
addpath('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\MATLAB\')
addpath('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\MATLAB\FRAP_simulation')


f = waitbar(0,'Initiliazing...','Name','LSQNONLIN takes forever',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
known = [NaN NaN 0 NaN NaN 0];
values = coeffvalues(FRAP_fit);
Ptranscribing = values(1);

N = 100;
idx = 0;
nrounds = 0;
while idx < N
    nrounds = nrounds + 1;
    try
        test = frap_lsm2(time, avg_FRAP(1,:), std_FRAP(1,:), 0.25, known, Ptranscribing);
        residual_threshold = 2;
        sqsum_threshold = 50;
        g = test.resnorm;
        passes_resids = sum(find(test.residual>2));
        endFRAP = mean(test.FRAP(end-40:end));
        if passes_resids == 0 && g <= sqsum_threshold && endFRAP <= 1
            idx = idx + 1;
            model6(idx) = test;
            waitbar(idx/N,f,sprintf('On round %i of %i', idx, N))
        end
    catch
    end
end

delete(f);


%% calculate AIC
%AIC = 2k+nln(sigma^2)
%k = m = parameters, n = observations, sigma^2 = RSS/(n-m)
%2 step model
% n = number of cells, m = 6
%3 step model
% n = number of cells, m = 12
clear
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\boot6WTresults.mat')
load("D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM\new_data_collection_method\Analysis\updated_analysis\average_FRAP_all_cells.mat")
bootavg_FRAP = mean(boot6FRAP,1);
resnormWT6 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICWT6 = 2*12 + 45*log((resnormWT6/(45-12))^2);
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\boot6rifresults.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM-Rif\Analysis\FRAP_analysis_181119\average_FRAP_all_cells.mat')
bootavg_FRAP = mean(boot6FRAP,1);
resnormRif6 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICRif6 = 2*12 + 45*log((resnormRif6/(45-12))^2);
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\average_FRAP_Chlor.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\boot6chlorresults.mat')
bootavg_FRAP = mean(boot6FRAP,1);
resnormChlor6 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICChlor6 = 2*12 + 45*log((resnormChlor6/(45-12))^2);
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\average_FRAP_mutant.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\boot6mutresults.mat')
bootavg_FRAP = mean(boot6FRAP,1);
resnormMut6 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICMut6 = 2*12+ 32*log((resnormMut6/(32-12))^2);
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\AsiA\Analysis\FRAP_analysis_181121\average_FRAP_all_cells.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\boot6asiaresults.mat')
bootavg_FRAP = mean(boot6FRAP,1);
resnormAsia6 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICAsia6 = 2*12 + 45*log((resnormAsia6/(45-12))^2);
load('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_M9\total_FRAP_all_cells_FH.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\boot6M9results.mat')
bootavg_FRAP = mean(boot6FRAP,1);
allFRAP = zeros(46,601);
for i = 1:46
    allFRAP(i,:) = total_FRAP_all_cells(i).ft2;
end
avg_FRAP = mean(allFRAP,1);
resnormM96 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICM96 = 2*12 + 46*log((resnormM96/(46-12))^2);
%model2
load('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM-Rif\Analysis\FRAP_analysis_181119\average_FRAP_all_cells.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\bootModel2rif.mat')
FRAP = horzcat(bootmodel2.FRAP);bootavg_FRAP = mean(vertcat(FRAP.FRAP),1);
resnormRif2 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICRif2= 2*6 + 45*log((resnormRif2/(45-6))^2);
load("D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_EZRDM\new_data_collection_method\Analysis\updated_analysis\average_FRAP_all_cells.mat")
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\bootModel2WT.mat')
FRAP = horzcat(bootmodel2.FRAP);bootavg_FRAP = mean(vertcat(FRAP.FRAP),1);
resnormWT2 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICWT2= 2*6 + 45*log((resnormWT2/(45-6))^2);
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\average_FRAP_Chlor.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\bootModel2Chlor.mat')
FRAP = horzcat(bootmodel2.FRAP);bootavg_FRAP = mean(vertcat(FRAP.FRAP),1);
resnormChlor2 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICChlor2= 2*6 + 45*log((resnormChlor2/(45-6))^2);
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\average_FRAP_mutant.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\bootModel2mut.mat')
FRAP = horzcat(bootmodel2.FRAP);bootavg_FRAP = mean(vertcat(FRAP.FRAP),1);
resnormMut2 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICMut2= 2*6 + 32*log((resnormMut2/(32-6))^2);
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Manuscripts\RNAP\AsiA\Analysis\FRAP_analysis_181121\average_FRAP_all_cells.mat')
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\bootModel2Asia.mat')
FRAP = horzcat(bootmodel2.FRAP);bootavg_FRAP = mean(vertcat(FRAP.FRAP),1);
resnormAsia2 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICAsia2= 2*6 + 45*log((resnormAsia2/(45-6))^2);
load('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\RNAP_GFP_M9\total_FRAP_all_cells_FH.mat')
allFRAP = zeros(46,601);
for i = 1:46
    allFRAP(i,:) = total_FRAP_all_cells(i).ft2;
end
avg_FRAP = mean(allFRAP,1);load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\bootModel2M9.mat')
FRAP = horzcat(bootmodel2.FRAP);bootavg_FRAP = mean(vertcat(FRAP.FRAP),1);
resnormM92 = (full(bootavg_FRAP)-full(avg_FRAP))*(full(bootavg_FRAP)-full(avg_FRAP))';
AICM92= 2*6 + 46*log((resnormM92/(46-6))^2);

%% histogram of kinetic values
figure
subplot(3,2,1)
[y,x] = hist(model2_allparams(1,:),0:0.05:2);
bar(x,y)
title('k1 distribution')
ylabel('counts')

subplot(3,2,2)
[y,x] = hist(model2_allparams(2,:), 0:0.05:2);
bar(x,y)
title('k2 distribution')

subplot(3,2,3)
[y,x] = hist(model2_allparams(3,:), 0:0.05:2);
bar(x,y)
title('k3 distribution')
ylabel('counts')

subplot(3,2,4)
[y,x] = hist(model2_allparams(4,:), 0:0.05:2);
bar(x,y)
title('k4 distribution')

subplot(3,2,5)
[y, x] = hist(model2_allparams(5,:), 0:0.05:2);
bar(x,y)
title('k5 distribution')
xlabel('rate')
ylabel('counts')

subplot(3,2,6)
[y,x] = hist(model2_allparams(6,:), 0:0.05:2);
bar(x,y)
title('k6 distribution')
xlabel('rate')
ylabel('counts')

%%
figure
plot(model2_allparams(1,:))
hold on
plot(model2_allparams(2,:), '-r')
plot(model2_allparams(3,:), '-g')
plot(model2_allparams(4,:), '-c')
plot(model2_allparams(5,:), '-m')
plot(model2_allparams(6,:), '-k')
xlabel('fitting iteration #', 'FontSize', 14)
ylabel('params', 'FontSize', 14)
title('Model 2 various output models from doing LSQNONLIN fitting', 'FontSize', 14)
legend('k1', 'k2', 'k3', 'k4', 'k5', 'k6')
set(gca, 'FontSize', 14)
%%
figure
hold on
for idx = 1:length(model2)
    plot(time, model2(idx).FRAP, '-', 'Color', [rand rand rand], 'LineWidth', 2)
end
plot(time, avg_FRAP, '-b', 'LineWidth', 3)
xlabel('time', 'FontSize', 14)
ylabel('FRAP', 'FontSize', 14)
title('Model 2 FRAP simulated curves vs. real data', 'FontSize', 14)
set(gca, 'FontSize', 14)
axis([0 601 0 1.5])

%%

figure
hold on
plot(time, avg_FRAP, '-b')
plot(time, model1.FRAP, '-r')
xlabel('time (s)', 'FontSize', 14)
ylabel('fraction recovery', 'FontSize', 14)
title({'rpoC-GFP + rifampicin, model 2', ...
    sprintf('k1 = %.4f, k2 = %.4f, Pfree = %.2f, Pbound = %.2f', k1, k2, k1/(k1+k2), k2/(k1+k2))}, 'FontSize', 14)
axis([0 601 0 1])
set(gca, 'FontSize', 14)

%%
meank1 = mean(model2_allparams(1,:));
stdk1 = std(model2_allparams(1,:));
meank3 = mean(model2_allparams(3,:));
stdk3 = std(model2_allparams(3,:));
meank4 = mean(model2_allparams(4,:));
stdk4 = std(model2_allparams(4,:));
meank5 = mean(model2_allparams(5,:));
stdk5 = std(model2_allparams(5,:));
meank6 = mean(model2_allparams(6,:));
stdk6 = std(model2_allparams(6,:));

figure
hold on
bar(1, meank1)
errorbar(1, meank1, stdk1, 'ok')
bar(3, meank3)
errorbar(3, meank3, stdk3, 'ok')
bar(4, meank4)
errorbar(4, meank4, stdk4, 'ok')
bar(5, meank5)
errorbar(5, meank5, stdk5, 'ok')
bar(6, meank6)
errorbar(6, meank6, stdk6, 'ok')
axis([0.5 6.5 -0.02 1.5])
ylabel('average value', 'FontSize', 14)
title('Best fit values for Model 2', 'FontSize', 14)
set(gca, 'FontSize', 14, 'XTick', 1:6, 'XTickLabel', {'k1', 'k2', 'k3', 'k4', 'k5', 'k6'})
%%


%% histogram of kinetic values
model6_allparams = horzcat(model6.params);
figure
subplot(3,2,1)
[y,x] = hist(model6_allparams(1,:),0.79:0.0005:0.80);
bar(x,y)
title('k1 distribution')
ylabel('counts')

subplot(3,2,2)
[y,x] = hist(model6_allparams(2,:), 0:0.05:2);
bar(x,y)
title('k2 distribution')

subplot(3,2,3)
[y,x] = hist(model6_allparams(3,:), 0:0.05:2);
bar(x,y)
title('k3 distribution')
ylabel('counts')

subplot(3,2,4)
[y,x] = hist(model6_allparams(4,:), 0:0.05:2);
bar(x,y)
title('k4 distribution')

subplot(3,2,5)
[y, x] = hist(model6_allparams(5,:), 0:0.05:2);
bar(x,y)
title('k5 distribution')
xlabel('rate')
ylabel('counts')

subplot(3,2,6)
[y,x] = hist(model6_allparams(6,:), 0:0.0005:1.5);
bar(x,y)
title('k6 distribution')
xlabel('rate')
ylabel('counts')


%%
figure
plot(model6_allparams(1,:))
hold on
plot(model6_allparams(2,:), '-r')
plot(model6_allparams(3,:), '-g')
plot(model6_allparams(4,:), '-c')
plot(model6_allparams(5,:), '-m')
plot(model6_allparams(6,:), '-k')
xlabel('fitting iteration #', 'FontSize', 14)
ylabel('params', 'FontSize', 14)
title('Model 6 various output models from doing LSQNONLIN fitting', 'FontSize', 14)
legend('k1', 'k2', 'k3', 'k4', 'k5', 'k6')
set(gca, 'FontSize', 14)
%%
figure
hold on
for idx = 1:length(model6b)
    plot(time, model6b(idx).FRAP, '-', 'Color', [rand rand rand], 'LineWidth', 2)
end
plot(time, avg_FRAP, '-b', 'LineWidth', 1)
xlabel('time', 'FontSize', 14)
ylabel('FRAP', 'FontSize', 14)
title('Model 6 FRAP simulated curves vs. real data', 'FontSize', 14)
set(gca, 'FontSize', 14)
axis([0 601 0 1])

%%
figure
for idx = 1:length(model6)
    hold on
    plot(time, avg_FRAP, '-b', 'LineWidth', 2)
    plot(time, model6(idx).FRAP, '-r')
    title(sprintf('model number %i', idx))
    waitforbuttonpress
    hold off
end

%%
%%
%%
%% Model 6b -- assume k3 = 0, k6 = 0, DON'T know k2, float all other parameters, input Pclust
clear
load('average_FRAP_all_cells.mat');
load('FRAP_fit.mat');
f = waitbar(0,'Initiliazing...','Name','LSQNONLIN takes forever',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

known = [NaN NaN 0 NaN NaN 0];

values = coeffvalues(FRAP_fit);
Pclust = values(1);

N = 100;
idx = 0;
nrounds = 0;
while idx < N
    nrounds = nrounds + 1;
    try
        test = frap_lsm2(time, average_FRAP, std_FRAP, 0.25, known, Pclust);
        residual_threshold = 2;
        sqsum_threshold = 25;
        g = test.resnorm;
        passes_resids = sum(find(test.residual>2));
        endFRAP = mean(test.FRAP(end-40:end));
        if passes_resids == 0 && g <= sqsum_threshold && endFRAP <= 1
            idx = idx + 1;
            model6b(idx) = test;
        end
    catch
    end
    waitbar(idx/N,f,sprintf('On round %i of %i', idx, N))
end

delete(f);

model6b_allparams = horzcat(model6b.params);

model6b_meank1 = mean(model6b_allparams(1,:));
model6b_meank2 = mean(model6b_allparams(2,:));
model6b_meank3 = mean(model6b_allparams(3,:));
model6b_meank4 = mean(model6b_allparams(4,:));
model6b_meank5 = mean(model6b_allparams(5,:));
model6b_meank6 = mean(model6b_allparams(6,:));

save('model6b.mat', 'model6b', 'model6b_allparams')

%% histogram of kinetic values
figure
subplot(3,2,1)
[y,x] = hist(model6b_allparams(1,:),0:0.05:2);
bar(x,y)
title('k1 distribution')
ylabel('counts')

subplot(3,2,2)
[y,x] = hist(model6b_allparams(2,:), 0:0.05:2);
bar(x,y)
title('k2 distribution')

subplot(3,2,3)
[y,x] = hist(model6b_allparams(3,:), 0:0.05:2);
bar(x,y)
title('k3 distribution')
ylabel('counts')

subplot(3,2,4)
[y,x] = hist(model6b_allparams(4,:), 0:0.05:2);
bar(x,y)
title('k4 distribution')

subplot(3,2,5)
[y, x] = hist(model6b_allparams(5,:), 0:0.05:2);
bar(x,y)
title('k5 distribution')
xlabel('rate')
ylabel('counts')

subplot(3,2,6)
[y,x] = hist(model6b_allparams(6,:), 0:0.05:2);
bar(x,y)
title('k6 distribution')
xlabel('rate')
ylabel('counts')

%%
figure
plot(model6b_allparams(1,:))
hold on
plot(model6b_allparams(2,:), '-r')
plot(model6b_allparams(3,:), '-g')
plot(model6b_allparams(4,:), '-c')
plot(model6b_allparams(5,:), '-m')
plot(model6b_allparams(6,:), '-k')
xlabel('fitting iteration #', 'FontSize', 14)
ylabel('params', 'FontSize', 14)
title('Model 6b various output models from doing LSQNONLIN fitting', 'FontSize', 14)
legend('k1', 'k2', 'k3', 'k4', 'k5', 'k6')
set(gca, 'FontSize', 14)
%%
figure
hold on
for idx = 1:length(model6b)
    plot(time, model6b(idx).FRAP, '-', 'Color', [rand rand rand], 'LineWidth', 2)
end
plot(time, avg_FRAP, '-b', 'LineWidth', 1)
xlabel('time', 'FontSize', 14)
ylabel('FRAP', 'FontSize', 14)
title('Model 6b FRAP simulated curves vs. real data', 'FontSize', 14)
set(gca, 'FontSize', 14)
axis([0 601 0 1.5])


%%
%%
%%
%% Model 7 -- float all the parameters????
clear
load('average_FRAP_all_cells.mat');
load('FRAP_fit.mat');
% f = waitbar(0,'Initiliazing...','Name','LSQNONLIN takes forever',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

known = [NaN NaN 0 NaN NaN 0];

values = coeffvalues(FRAP_fit);
Pclust = values(1);

q = parallel.pool.DataQueue;
afterEach(q, @disp);

parfor ii = 1:100
	model7(ii) = frap_lsm2b(time, avg_FRAP, std_FRAP, 0.25, known, Pclust);
	% waitbar(ii/N,f,sprintf('On round %i of %i', idx, N))
    send(q, ii);
end
 
model7_allparams = horzcat(model7.params);

model7_meank1 = mean(model7_allparams(1,:));
model7_meank2 = mean(model7_allparams(2,:));
model7_meank3 = mean(model7_allparams(3,:));
model7_meank4 = mean(model7_allparams(4,:));
model7_meank5 = mean(model7_allparams(5,:));
model7_meank6 = mean(model7_allparams(6,:));

save('model7.mat', 'model7', 'model7_allparams')