% KB 2018-11-16
%
% Attempt to solve for kinetic parameters from average FRAP
%
% basically: LSM of simulated FRAP to real FRAP
% std_FRAP = optional

% input known params with all other values in array as NaN

function output = FRAP_lsm(time, avg_FRAP, std_FRAP, p, known, Pclust)

w = zeros(length(std_FRAP),1);
for idt = 1:length(std_FRAP)
    w(idt) = 1/std_FRAP(idt);
end
w = w';

if nargin == 5
    Pclust = 0;
end

know = known;
n_unknown = 6 - length(know(~isnan(know)));
guess = rand(n_unknown, 1);

resid1 = @(dontknow) sum(w.*(FRAP_sim_model2_kb(know, dontknow, p, length(time), Pclust)-avg_FRAP).^2);

[Coeff, g, exitflag] = fminsearch(resid1, guess,...
    optimset('MaxIter',10000,'TolX',5e-7,'TolFun',5e-7,'MaxFunEvals',100000));

know_indices = find(~isnan(know));
dk_indices = find(isnan(know));

params = zeros(6,1);

for idx = 1:length(know_indices)
    params(know_indices(idx)) = know(know_indices(idx));
end

for idx = 1:length(dk_indices)
    params(dk_indices(idx)) = Coeff(idx);
end

FRAP_out = FRAP_sim_model2_kb_out(params, p, length(time), Pclust);

output.params = params;
output.FRAP = FRAP_out;

end