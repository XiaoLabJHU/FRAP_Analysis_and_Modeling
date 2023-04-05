%%
clear
load('D:\Xiao Lab Dropbox\Lab Members\Harris_Fran\Imaging\FinalFRAPresults\average_FRAP_bootrif.mat');
addpath('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\MATLAB\')
addpath('D:\Xiao Lab Dropbox\Lab Members\Alumni\Bettridge_Kelsey\MATLAB\FRAP_simulation')

tic
%f = waitbar(0,'Initiliazing...','Name','LSQNONLIN takes forever',...
%    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
boot6Pclust = zeros(1,100);
boot6Pfree = zeros(1,100);
boot6Pelong = zeros(1,100);
boot6Params = zeros(100,6);
boot6FRAP = zeros(100,601);
boot6resmorm = zeros(1,100);
known = [NaN NaN 0 NaN NaN 0];
N = 100;

parfor ii = 1:100
    Pclust = bootPclust(ii);
    idx = 0;
    nrounds = 0;
    minresnorm=50;
    while idx < N
        nrounds = nrounds + 1;
        try
            test = frap_lsm2(time, avg_FRAP(ii,:), std_FRAP(ii,:), 0.25, known, Pclust);
            residual_threshold = 2;
            sqsum_threshold = 50;
            g = test.resnorm;
            passes_resids = sum(find(test.residual>2));
            endFRAP = mean(test.FRAP(end-40:end));
            if passes_resids == 0 && g <= sqsum_threshold && endFRAP <= 1
                idx = idx + 1;
                if test.resnorm<minresnorm
                    bestmodel6 = test;
                    minresnorm = test.resnorm;
                end
                %waitbar(idx/N,f,sprintf('On round %i of %i', idx, N))
            end
        catch
        end
    end
    boot6Pclust(ii) = bestmodel6.Pclust;
    boot6Pfree(ii) = bestmodel6.Pfree;
    boot6Pelong(ii) = bestmodel6.Pelong;
    boot6resnorm(ii) = minresnorm;
    boot6Params(ii,:) =bestmodel6.params;
    boot6FRAP(ii,:) = bestmodel6.FRAP;
end
%delete(f);
toc

%%
%model6_allparams = horzcat(model6.params);
figure
hold on
plot(time,avg_FRAP(1,:))
plot(time, model6(2).FRAP)