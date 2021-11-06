%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vanilla-IMR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vanilla-IMR (Inertial microcavitation rheometry) is to solve the 
% viscoelastic properties of surrounding materials, 
%
% Author: Jin Yang; Postdoc @UW-Madison '19
% Contact: jyang526@wisc.edu
% Modified based on older version of IMR codes written by:
% Author: Jon Estrada @Brown Univ. PhD '17 &
% Author: Carlos Barajas @Univ.Michigan-Ann Arbor. BS '16
%
% =========================================================================
% Reference:
% [1] Estrada, J. B., Barajas, C., Henann, D. L., Johnsen, E., & Franck, C. (2018). 
%     High strain-rate soft material characterization via inertial cavitation. 
%     Journal of the Mechanics and Physics of Solids, 112, 291-317.
% [2] Barajas, C., & Johnsen, E. (2017). The effects of heat and mass diffusion
%     on freely oscillating bubbles in a viscoelastic, tissue-like medium. 
%     The Journal of the Acoustical Society of America, 141(2), 908-918.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

%% Section 1. Load R-t curve
fprintf('---------- Section 1 Load R-t curve ---------- \n');
fprintf('Choose method to load IMR R-t curve: \n');
fprintf('(1)Load R(t) curve by post-processing IMR experiment images \n');
fprintf('(2)Load R(t) curve by generating numerical simulations \n');
prompt = 'Input here: '; IMR_img_method = input(prompt);

switch IMR_img_method
    case 1
        %% Section 1.1 Get RofT curve by loading IMR images
        % ====== Postprocessing exp images ======
        fprintf('---------- Section 1.1 Load R-t curve from experiment images ---------- \n');
        videoFolderName = ['./data/IMR_JE_StiffPA_exp2/']; cd(videoFolderName);
        videoName = 'PAstiff_xy002.mp4'; LoadVideo; % Load video 
        simNo = 1; TimeRes = 0; expts = [2]; expt=expts;
        folderNamePrefix = ['./data/IMR_JE_PA/',num2str(simNo),'/'];
        fileNamePrefix = ['PA',num2str(simNo),'_RofTdata_',num2str(TimeRes)];

        runCalcRofT; % Need experience and still art of work; 
        cd('../../'); % Back to main directory
        fprintf('---------- Section 1.1 Done ---------- \n\n');

    case 2
        %% Section 1.2 Get RofT curve by numerical simulations
        fprintf('---------- Section 1.2 Numerical simulation for R-t curve ---------- \n');
        % ====== Time duration ======
        tspan = 2.25e-4; 
        % ====== Material model ======
        % G = 2.97e3; mu = 0.101; model = 'neoHook'; simNo = 4;
        % G = 7.69e3; mu = 0.05; model = 'neoHook'; simNo = 4; % Stiff PA in Estrada. et al JMPS paper
        % G = 2.12e3; mu = 0.118; model = 'neoHook'; simNo = 4; % Soft PA in Estrada. et al JMPS paper
        G = 2.97e3; mu = 0.04; alpha = 0.6; lambda_nu = 1; model = 'fungnlvis'; simNo = 7; 
        
        % ====== Providing initial conditions ======
        %(1)Given bubble initial radius and equilibirum radius
        ICReqOrPinit = 0; R0 = 4.085968222180783e2*1e-6; 
        Req0 = 54.818020051508434*1e-6; %R0*ones(1,1)./[10]; 
        expNo = length(Req0);
        %(2)Given bubble initial radius and initial inside pressure
        %ICReqOrPinit = 1; % R0 = 225e-5; P_guessList = [200]; expNo = length(P_guessList);
        % ====== Other default settings ======   
        NT = 500; % Inside bubble spatial grid #; 
        NMt = 10; % Outside bubble spatial grid #;
        IMRsolver_RelTolX = 1e-7; % Matlab ode23tb solver relative tolerance;
        TimeRes = 10; % "TimeRes" is the sampling rate you want to set over 270,000/s.
        % E.g., "TimeRes = 10" means we sample R-t data points in the sampling rate of 10*270,000/s.
        % % All other physical constant parameters, see IMRcall_parameters.m and runIMR_num_simulation.m.
        
        % ====== Solver ======   
        runIMR_num_simulation; % IMR numerical simulation
        fprintf('---------- Section 1.2 Done ---------- \n\n');
        
  otherwise
end

%% Section 1.3 Fit RofT curve using spline interpolation
fprintf('---------- Section 1.3 R-t curve fitting ----------\n');
% ====== Please the same parameters in Section 1.1 or 1.2 ======
% ====== Curve fittings ======
tspan = 2.25e-4;  % Time span 
expts = [1];  % E.g.: expts = [1:1:N], where N is total exp #.
% folderNamePrefix = ['./data/numsim/',num2str(simNo),'/'];
% fileNamePrefix = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes)];
simNo = 4; % simNo is the folder# for different material model/different experimental materials
% E.g., "simNo = 4" is the Neo-Hookean Kelvin-Voigt model/or #4 tested material
TimeRes = 10; % "TimeRes" is the sampling rate you want to set over 270,000/s.
% E.g., "TimeRes = 10" means we sample R-t data points in the sampling rate of 10*270,000/s.
% E.g., "TimeRes = 0" if using experimenal images, to use the original time data points

% folderNamePrefix = ['./data/IMR_JE_PA/',num2str(simNo),'/'];
% fileNamePrefix = ['PA',num2str(simNo),'_RofTdata_',num2str(TimeRes)];
folderNamePrefix = ['./data/numsim/',num2str(simNo),'/'];
fileNamePrefix = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes)];
prompt = 'Do you want to fit ROfT curve? 0-Yes; 1-No. Input here: '; FitROfTOrNot = input(prompt);
if FitROfTOrNot<1, FitROfTCurve; end
fprintf('---------- Section 1.3 Done ----------\n\n');


%% Sections 2-3 IMR solver to fit surrounding material properties
fprintf('---------- Sections 2-3 Fit material property using vanilla IMR ---------- \n');
% ====== Please the same parameters in Section 1 ======
% ====== Define the number of peaks in fitted R-t curve ======
fprintf('How many peaks of R-t curve to be used for fitting? \n');
prompt = 'Input here: '; PickNo = input(prompt);

% ====== Assign experiment index ======
tspan = 2.25e-4;
expts = [1]; % E.g.: expts = [1:1:N], where N is total exp #.
TimeRes = 10; % "TimeRes" is the sampling rate you want to set over 270,000/s.
% E.g., "TimeRes = 10" means we sample R-t data points in the sampling rate of 10*270,000/s.
% E.g., "TimeRes = 0" if using experimenal images, to use the original time data points
simNo = 4; % simNo is the folder# for different material model/different experimental materials
% E.g., "simNo = 4" is the Neo-Hookean Kelvin-Voigt model/or #4 tested material
expNo = length(expts);  model='neoHook'; %
% folderNamePrefix = ['./data/IMR_JE_PA/',num2str(simNo),'/'];
% fileNamePrefix = ['PA',num2str(simNo),'_RofTdata_',num2str(TimeRes)];
folderNamePrefix = ['./data/numsim/',num2str(simNo),'/'];
fileNamePrefix = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes)];
% if exist('folderNamePrefix','var')==0, folderNamePrefix = ['./data/Natassa/']; end
% if exist('fileNamePrefix','var')==0, fileNamePrefix = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes)]; end
        
% ====== Choose IMR solver method ======
fprintf('---------- Choose IMR solver method: ---------- \n');
fprintf('(1) Section 2.1: IMR solver with LSQ fitting \n');
fprintf('(2) Section 2.2: IMR solver with Nelder-Mead fitting \n');
fprintf('(3) Section 3.1: Decoupled IMR solver with LSQ fitting in time segments \n');
fprintf('(4) Section 3.2: Decoupled IMR solver with Nelder-Mead fitting in time segments \n');
prompt = 'Input here: '; IMRsolver_method = input(prompt);

switch IMRsolver_method
    
    case 1
        %% Section 2.1 IMR solver with LSQ fitting
        fprintf('---------- Section 2.1 Fit material property using vanilla IMR LSQ fitting ---------- \n');
        % ====== Define LSQ fitting range ======
        G_ooms = log10(2000):0.5:log10(20000); % Shear modulus
        mu_ooms = -2.0:0.5:-1.2; % Viscosity
        alpha_ooms = -Inf; lambda_nu_ooms = -Inf; G1_ooms = inf; % By default other parameters for Neo-Hookean KV material
        % ====== Solver ======
        runIMR; % Compute LSQ error matrix
        % Fitted material property output is in variable: "matPropVarList"
        fprintf('---------- Section 2.1 Done ---------- \n\n');
        
    case 2
        %% Section 2.2 IMR solver with Nelder-Mead fitting
        fprintf('---------- Section 2.2 Fit material property using vanilla IMR Nelder-Mead fitting ---------- \n');
        % ====== Initial value for fitting material property ======
        GInit = 2e3; muInit = 0.05; matPropVar0 = [log10(GInit), log10(muInit)]; % Store in matPropVar0 vector
        % ====== Upper and lower bounds for material property ======
        lb = log10([100, 5e-3]); ub = log10([1.5e5, 0.2]);
        % ====== Solver ======
        runIMR_NelderMead;
        fprintf('---------- Section 2.2 Done ---------- \n\n');

    case 3
        %% Section 3.1 Decoupled IMR solver with LSQ fitting in time segmentss
        fprintf('---------- Section 3.1 Fit material property using decoupled IMR w/ LSQ fitting ---------- \n');
        % ====== Step 1: Solve {P,C,T} decoupled system ====== 
        % tspan = 3.0e-7;
        fprintf('Do you want to solve first sub-system? 0-No; 1-Yes, accurate; 2-Yes, fast. \n');
        prompt = 'Input here: '; whetherToSolveP = input(prompt);
        DefineLocsMaxPOrNot = 0; DefineTimetolListOrNot = 0; PlotDefineTimetolListOrNot = 0;
        if whetherToSolveP > 0, PRange = [175,200]; runIMR2; end
        
        % ====== Step 2: Solve {R} decoupled system ======
        G_ooms = log10(2000):0.02:log10(12000); % Shear modulus
        mu_ooms = -2.0:0.02:-0.3; % Viscosity
        G1_ooms = Inf;
        % Solver 
        whetherToSolveP = 0; DefineLocsMaxPOrNot = 0; DefineTimetolListOrNot = 1; PlotDefineTimetolListOrNot = 0;
        if PlotDefineTimetolListOrNot == 1, matPropVar_best = [log10(275), log10(0.003)]; end
        runIMR2;
        fprintf('---------- Section 3.1 Done ---------- \n\n');

    case 4
        %% Section 3.2 Decoupled IMR solver with Nelder-Mead fitting in time segments
        % Only correct for cold surrounding hydrogel assumption
        fprintf('---------- Section 3.2 Fit material property using decoupled IMR w/ Nelder-Mead method ---------- \n');
        % ====== Step 1: Solve {P,C,T} decoupled system ====== 
        fprintf('Do you want to solve first sub-system? 0-No; 1-Yes, accurate; 2-Yes, fast. \n');
        prompt = 'Input here: '; whetherToSolveP = input(prompt);
        DefineLocsMaxPOrNot = 0; DefineTimetolListOrNot = 0; PlotDefineTimetolListOrNot = 0;
        if whetherToSolveP > 0, PRange = [80,200]; runIMR2; end
        
        % ====== Step 2: Solve {R} decoupled system ======
        % Initial value for fitting material property
        GInit = 2e3; muInit = 0.006; matPropVar0 = [log10(GInit), log10(muInit)]; % Store in matPropVar0 vector
        % Upper and lower bounds for material property
        lb = log10([1000, 0.001]); ub = log10([4000, 0.1]);
        % Solver 
        whetherToSolveP = 0; DefineLocsMaxPOrNot = 1; DefineTimetolListOrNot = 1; PlotDefineTimetolListOrNot = 0;
        if PlotDefineTimetolListOrNot == 1, matPropVar_best = [log10(3948), log10(0.021146)]; end
        NelderMead_TolX = 1e-6; runIMR2_NelderMead;
        fprintf('---------- Section 3.2 Done ---------- \n\n');
 
    otherwise
        
end

%% Finally output
matPropVarList % 10^(matPropVarList) is the final result.






