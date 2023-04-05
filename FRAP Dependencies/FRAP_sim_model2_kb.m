%% JX181116, FRAP, model 2
% model 2: there is interchange between R_elong and R_clust
% R_free --> R_elong, k1
% R_elong --> R_free, k2, not dissociation rate but transcription rate:
% elongation complex only returns to free RNAP by finishing transcription
% R_free --> R_clust, k3
% R_clust --> R_free, k4
% R_elong --> R_clust, k5
% R_clust --> R_elong, k6
% one can disbale one reaction by setting the corresponding k = 0;
%
% kb modification: can input known parameters and unknown parameters
%     put known in order of [k1 k2 k3 k4 k5 k6] where if they are unknown,
%     you put in "NaN" instead of the value
%     - the unknownparams is an nx1 array where n <=6

function FRAP = FRAP_sim_model2_kb(knownparams, unknownparams, p, n, Pclust)

% the following values were estimated from Gotta-JBact-1991
% k1 = 1/2; %(s-1) (initiation rate on rrn)
% k2 = 1/130; % (s-1, elongation rate: ~ 2.2 min for the entire rrn operon)
% 
% % Values of k3, k4, k5 and k6 are unknown and are rough estimates
% k3 = 1/20; %(s-1)
% k4 = 1/1000; %(s-1)
% 
% k5 = 1/20;
% k6 = 1/100; 

known_indices = find(~isnan(knownparams));
for idx = 1:length(known_indices)
    k(known_indices(idx)) = knownparams(known_indices(idx));
end

unknown_indices = find(isnan(knownparams));
for idx = 1:length(unknown_indices)
    k(unknown_indices(idx)) = unknownparams(idx);
end

k1 = k(1);
k2 = k(2);
k3 = k(3);
k4 = k(4);
k5 = k(5);
k6 = k(6);

% n = round(6/min(k(find(k)))); % total simulation time: 6 times of the slowest non zero kinetic constant to be considered complete
% n = 601;
 R_free = zeros(n, 1);
 R_elong = zeros(n, 1); 
 R_clust = zeros(n, 1); 
%%
% initial condition
% try to estimate based off the input parameters
N = 4000; 

if nargin== 4
    R_free(1)  = N; % molecules per cell
    R_elong(1)  = 0;
    R_clust(1)  = 0;
elseif nargin == 5
    R_clust(1) = N*Pclust;
    R_free(1) = N-N*Pclust;
    R_elong(1) = 0;
end



for t = 2 : 1: n
    R_free(t)  = R_free(t-1) + k2*R_elong(t-1) + k4*R_clust(t-1) - (k3 + k1)*R_free(t-1);
    R_elong(t) = R_elong(t-1) + k1*R_free(t-1) + k6*R_clust(t-1) - (k2 + k5)*R_elong(t-1);
    R_clust(t) = R_clust(t-1) + k3*R_free(t-1) + k5*R_elong(t-1) - (k4 + k6)*R_clust(t-1);
end

RNAP.free_total  = R_free/N;
RNAP.elong_total = R_elong/N;
RNAP.clust_total = R_clust/N;


% simulate FRAP
% the only exchange happens at the level of R_free
% R_free_ <--> R_free_f, k7 ~= 20 s-1, using a D at 0.4 um2/s (D =
% 0.88*r^2/(4*0.693/k), r is the radius of the photobleaching area)
% k7 = 20; therefore I assume rapid equilibrium between the two
% populations in area A1 and A2

% A1 is the photobleached area
% A2 is the non=bleached area
% p = A1/A2
% A1 = 1;
% A2 = 9; 
% p = A1/(A1+A2); % ratio of photobleached area to the whole cell area
% at t = 0;

% bleached molecules in area A1
R_free_b_A1(1)  = p*R_free(end);
R_elong_b_A1(1) = p*R_elong(end);
R_clust_b_A1(1) = p*R_clust(end);

% bleached molecules in area A2
R_free_b_A2(1)  = 0;
R_elong_b_A2(1)  = 0;
R_clust_b_A2(1)  = 0;

% fluorescent molecules in area A2
R_free_f_A2(1)  = (1-p)*R_free(end);
R_elong_f_A2(1) = (1-p)*R_elong(end);
R_clust_f_A2(1) = (1-p)*R_clust(end);

% fluorescent molecules in area A1
R_free_f_A1(1)  = 0;
R_elong_f_A1(1)  = 0;
R_clust_f_A1(1)  = 0;


% at time t: give more time to allow FRAP to complete
 for t = 2 : 1 : n  %2*n
%     t = 4
    R_free_b_A1(t) = p*(R_free_b_A1(t-1) + R_free_b_A2(t-1)); % assuming fast equilibrium between the two populations
    R_free_b_A2(t) = (1-p)*(R_free_b_A1(t-1) + R_free_b_A2(t-1));
    R_free_f_A1(t) = p*(R_free_f_A1(t-1) + R_free_f_A2(t-1));
    R_free_f_A2(t) = (1-p)*(R_free_f_A1(t-1) + R_free_f_A2(t-1));
    
%     R_free(t)  = R_free(t-1) + k2*R_elong(t-1) + k4*R_clust(t-1) - (k3 + k1)*R_free(t-1);
%     R_elong(t) = R_elong(t-1) + k1*R_free(t-1) + k6*R_clust(t-1) - (k2 + k5)*R_elong(t-1);
%     R_clust(t) = R_clust(t-1) + k3*R_free(t-1) + k5*R_elong(t-1) - (k4 + k6)*R_clust(t-1);
    
    
    R_free_b_A1(t)  = R_free_b_A1(t) + k2*R_elong_b_A1(t-1) + k4*R_clust_b_A1(t-1) - (k1 + k3)*R_free_b_A1(t);
    R_elong_b_A1(t) = R_elong_b_A1(t-1) + k1*R_free_b_A1(t) + k6*R_clust_b_A1(t-1) - (k2 + k5)*R_elong_b_A1(t-1);
    R_clust_b_A1(t) = R_clust_b_A1(t-1) + k3*R_free_b_A1(t) + k5*R_elong_b_A1(t-1) - (k4 + k6)*R_clust_b_A1(t-1);
    

    R_free_b_A2(t)  = R_free_b_A2(t) + k2*R_elong_b_A2(t-1) + k4*R_clust_b_A2(t-1) - (k1 + k3)*R_free_b_A2(t);
    R_elong_b_A2(t) = R_elong_b_A2(t-1) + k1*R_free_b_A2(t) + k6*R_clust_b_A2(t-1) - (k2 + k5)*R_elong_b_A2(t-1);
    R_clust_b_A2(t) = R_clust_b_A2(t-1) + k3*R_free_b_A2(t) + k5*R_elong_b_A2(t-1) - (k4 + k6)*R_clust_b_A2(t-1);


    R_free_f_A2(t)  = R_free_f_A2(t) + k2*R_elong_f_A2(t-1) + k4*R_clust_f_A2(t-1) - (k1 + k3)*R_free_f_A2(t);
    R_elong_f_A2(t) = R_elong_f_A2(t-1) + k1*R_free_f_A2(t) + k6*R_clust_f_A2(t-1) - (k2 + k5)*R_elong_f_A2(t-1);
    R_clust_f_A2(t) = R_clust_f_A2(t-1) + k3*R_free_f_A2(t) + k5*R_elong_f_A2(t-1) - (k4 + k6)*R_clust_f_A2(t-1);
    
    R_free_f_A1(t)  = R_free_f_A1(t) + k2*R_elong_f_A1(t-1) + k4*R_clust_f_A1(t-1) - (k1 + k3)*R_free_f_A1(t);
    R_elong_f_A1(t) = R_elong_f_A1(t-1) + k1*R_free_f_A1(t) + k6*R_clust_f_A1(t-1) - (k2 + k5)*R_elong_f_A1(t-1);
    R_clust_f_A1(t) = R_clust_f_A1(t-1) + k3*R_free_f_A1(t) + k5*R_elong_f_A1(t-1) - (k4 + k6)*R_clust_f_A1(t-1);
  
 end

 % assign each time trace to the output structure
 
% RNAP.free_b_A1  = R_free_b_A1/N;
% RNAP.free_b_A2  = R_free_b_A2/N;
% RNAP.free_f_A1  = R_free_f_A1/N;
% RNAP.free_f_A2  = R_free_f_A2/N;
% 
% RNAP.clust_b_A1  = R_clust_b_A1/N;
% RNAP.clust_b_A2  = R_clust_b_A2/N;
% RNAP.clust_f_A1  = R_clust_f_A1/N;
% RNAP.clust_f_A2  = R_clust_f_A2/N;
% 
% RNAP.elong_b_A1  = R_elong_b_A1/N;
% RNAP.elong_b_A2  = R_elong_b_A2/N;
% RNAP.elong_f_A1  = R_elong_f_A1/N;
% RNAP.elong_f_A2  = R_elong_f_A2/N;

% calculate FRAP fraction according to experiments

I_bf_A1 = p*(1-p)*N; % expected recovery of fluorescence in number of molecules, adjusted for photobleached molecules
R_f_A1 = (R_free_f_A1 + R_elong_f_A1 + R_clust_f_A1)/I_bf_A1; % fraction of recovery at time t

FRAP = R_f_A1;

% plotting

% figure
% 
% subplot (2, 3, 1)
% plot (RNAP.free_total, 'k', 'linewidth', 3)
% hold on
% plot(RNAP.elong_total, 'r', 'linewidth', 3)
% plot(RNAP.clust_total, 'g', 'linewidth', 3)
% hold off
% ax = gca;
% ax.FontSize = 16;
% ax.LineWidth = 2;
% xlabel ('Time (s)', 'FontSize', 16)
% ylabel ('Fraction', 'FontSize', 16)
% lgd = legend ('RNAP free total', 'RNAP elong total', 'RNAP clust total');
% lgd.FontSize = 16;
% lgd.EdgeColor = [1 1 1];
% lgd.Location = 'best';
% 
% 
% subplot(2, 3, 2)
% plot(RNAP.FRAP, 'b', 'linewidth', 3)
% title ('Total Fluorescence Recovery in Bleached Area', 'FontSize', 16)
% ax = gca;
% ax.FontSize = 16;
% ax.LineWidth = 2;
% xlabel ('Time (s)', 'FontSize', 16)
% ylabel ('Fraction of Total Fluorescence Recovery', 'FontSize', 16)
% 
% subplot(2, 3, 3)
% plot (RNAP.free_f_A1, 'k', 'linewidth', 3)
% hold on
% plot(RNAP.elong_f_A1, 'r', 'linewidth', 3)
% plot(RNAP.clust_f_A1, 'g', 'linewidth', 3)
% hold off
% title ('Fraction of Fluorescent Molcules in Bleached Area', 'FontSize', 16)
% ax = gca;
% ax.FontSize = 16;
% ax.LineWidth = 2;
% xlabel ('Time (s)', 'FontSize', 16)
% ylabel ('Fraction', 'FontSize', 16)
% lgd = legend ('RNAP free', 'RNAP elong', 'RNAP clust');
% lgd.FontSize = 16;
% lgd.EdgeColor = [1 1 1];
% lgd.Location = 'best';
% 
% subplot(2, 3, 4)
% plot (RNAP.free_b_A1, 'k', 'linewidth', 3)
% hold on
% plot(RNAP.elong_b_A1, 'r', 'linewidth', 3)
% plot(RNAP.clust_b_A1, 'g', 'linewidth', 3)
% hold off
% title ('Fraction of Bleached Molcules in Bleached Area', 'FontSize', 16)
% ax = gca;
% ax.FontSize = 16;
% ax.LineWidth = 2;
% xlabel ('Time (s)', 'FontSize', 16)
% ylabel ('Fraction', 'FontSize', 16)
% lgd = legend ('RNAP free', 'RNAP elong', 'RNAP clust');
% lgd.FontSize = 16;
% lgd.EdgeColor = [1 1 1];
% lgd.Location = 'best';
% 
% subplot(2, 3, 5)
% plot (RNAP.free_f_A2, 'k', 'linewidth', 3)
% hold on
% plot(RNAP.elong_f_A2, 'r', 'linewidth', 3)
% plot(RNAP.clust_f_A2, 'g', 'linewidth', 3)
% hold off
% title ('Fraction of Fluorescent Molcules in non-Bleached Area', 'FontSize', 16)
% ax = gca;
% ax.FontSize = 16;
% ax.LineWidth = 2;
% xlabel ('Time (s)', 'FontSize', 16)
% ylabel ('Fraction', 'FontSize', 16)
% lgd = legend ('RNAP free', 'RNAP elong', 'RNAP clust');
% lgd.FontSize = 16;
% lgd.EdgeColor = [1 1 1];
% lgd.Location = 'best';
% 
% subplot(2, 3, 6)
% plot (RNAP.free_b_A2, 'k', 'linewidth', 3)
% hold on
% plot(RNAP.elong_b_A2, 'r', 'linewidth', 3)
% plot(RNAP.clust_b_A2, 'g', 'linewidth', 3)
% hold off
% title ('Fraction of Bleached Molcules in non-Bleached Area', 'FontSize', 16)
% ax = gca;
% ax.FontSize = 16;
% ax.LineWidth = 2;
% xlabel ('Time (s)', 'FontSize', 16)
% ylabel ('Fraction', 'FontSize', 16)
% lgd = legend ('RNAP free', 'RNAP elong', 'RNAP clust');
% lgd.FontSize = 16;
% lgd.EdgeColor = [1 1 1];
% lgd.Location = 'best';

end