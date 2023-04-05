
%%
%Jie Xiao, 150512

RingI_bB = mean(FRAP_data.RingI(2:20))-mean(FRAP_data.BackI(40:end))*FRAP_data.AreaRing(1)/FRAP_data.AreaBack(1); 
% photobleaching area intensity before photobleaching (background
% substracted)

CellI_bB = mean(FRAP_data.CellI(2:20))-mean(FRAP_data.BackI(40:end))*FRAP_data.AreaCell(1)/FRAP_data.AreaBack(1);
% total cell area intensity before photobleaching (background
% substracted)

CellI_aB = FRAP_data.CellI(21:end)-mean(FRAP_data.BackI(40:end))*FRAP_data.AreaCell(1)/FRAP_data.AreaBack(1);
% total intensity of the cell after photobleaching (background substracted)

% try to normalize so that it starts at 0; scale from pre-bleach frames
% aka Fnew(t) = (F(t)-F(0))/(Finf-F(0))



% METHOD 1 FOR DETERMINING D FROM FRAP DATA (Elowitz, 1999)
% substract the background
Cell_FL = double(FRAP_data.CellFL(:,:,21:end))-mean(FRAP_data.BackI(40:end))/FRAP_data.AreaBack;

Cell_FL_prof = [];
FL_prof_norm = [];
I = [];
A = [];
xx = [];
yy = [];

% Find and average the intensity along the x-axis for each time frame
for i=1:size(Cell_FL,3)
    for j=1:size(Cell_FL,2)
        Cell_FL_prof(i,j) = mean(Cell_FL(:,j,i));
    end
end
% Now normalize the fluorescence in relation to the rest of the stack
for i=1:size(Cell_FL_prof,1)
    FL_prof_norm(i,:) = Cell_FL_prof(i,:)/sum(Cell_FL_prof(i,:));
end

% Define the length of the cell
umperpix = 0.1; % convert pix to um
L = size(Cell_FL,2)*umperpix;

% Define parameters for calculating Fourier amplitude of nth mode An(t)
 n = 1;
 qn = n*pi/(size(FRAP_data.CellFL,2)*umperpix);

% find Fourier amplitude, proportional to exp(-qn^2*D*t)
for j=1:size(FL_prof_norm,1)
    for i=1:size(FL_prof_norm,2)
        I(j,i) = cos(qn*i*umperpix)*FL_prof_norm(j,i); % calculate the function
    end
    A(j) = 2/L * trapz(I(j,:)); % numerical integration
end

yy = A(1:30); % why only fit to first 30 pts?
xx = [0:1:29];

Fourier_fit_function = fittype(sprintf('a*exp(-x*b*%04d^2)+c',qn),'independent',{'x'},'coefficients',{'a','b','c'});
[Fourier_fit, Fourier_gof, Fourier_output] = fit(xx(:),yy(:),Fourier_fit_function,'StartPoint',[0 0.5 0],'Lower',[-Inf 0 0],'Upper',[Inf 5 Inf]);

% get the values of the coefficients
fit_values = coeffvalues(Fourier_fit);
% 'b' is the diffusion coefficient
Elowitz_D = fit_values(2);
Elowitz_ampl = fit_values(1);

RingI_maxR = CellI_aB*RingI_bB/CellI_bB; 
% maximal level of fluorescence recovery based on the ratio of photobleach
% area intensity and total cell instensity before photobleaching;note it is
% not adjusted by the ratio of areas

y_ring = FRAP_data.RingI(21:end)-mean(FRAP_data.BackI(40:end))*FRAP_data.AreaRing(1)/FRAP_data.AreaBack(1);
for i=1:length(RingI_maxR)
    y_ringNor(i) = y_ring(i)/RingI_maxR(i);
    y_ringNor2(i) = (y_ringNor(i)-y_ringNor(1))/(1-y_ringNor(1));
end
y_cell = FRAP_data.CellI(21:end)-mean(FRAP_data.BackI(40:end))*FRAP_data.AreaCell(1)/FRAP_data.AreaBack(1);
y_ringExp = y_cell * RingI_bB/CellI_bB;
x = 1*[1:length(y_ring)]';

% METHOD 2 FOR DETERMINING D FROM FRAP DATA (Bratton, 2011)

f0 = y_ring(1)/y_cell(1);
fINF = RingI_bB/CellI_bB;
for i=1:length(y_ring)
    f(i) = y_ring(i)/y_cell(i);
end
ft = (f-f0)/(fINF-f0);
xxx = [0:29];
yyy = ft(1:30);

ddt = f(2)-f(1);
ff0 = f(1)-ddt;
ft2 = (f-ff0)/(fINF-ff0);

Brutton_fittype = 'a*exp(-x/tau)+c';
try
    [Brutton_fit, Brutton_gof, Brutton_output] = fit(xxx(:),yyy(:),Brutton_fittype,'StartPoint',[0 0.5 1],'Lower',[-Inf 0 -Inf],'Upper',[1 Inf Inf]);
    Brutton_values = coeffvalues(Brutton_fit);
    Brutton_D = Brutton_values(2)/(qn^2);
    Brutton_mobile_fraction = -(Brutton_values(1));
    Brutton_immobile_fraction = 1-Brutton_mobile_fraction;
catch
    Brutton_fit = [];
    Brutton_gof = [];
    Brutton_output = [];
    Brutton_values = [];
    Brutton_D = [];
    Brutton_mobile_fraction = [];
    Brutton_immobile_fraction = [];
end
