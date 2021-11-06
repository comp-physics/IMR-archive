%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMR numerical simulation
% with setting material model and bubble collapse type
% Author: Jin Yang, jyang526@wisc.edu
% Date: 2019.01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all; close all;
warning('off','all')
 
%% Input information
% ====== Set time duration for the simulation (s) ======
%tspan = 2.25E-4; %1.3E-4;

% ====== Set material model ======
% ------ numsim1: G=0, linear viscosity ------
% G = 0; mu = 0.01; alpha = 0; lambda_nu = 0; G1 = inf; simNo = 1;  
% ------ numsim2: G=0, first-order nonlinear viscosity ------
% G = 0; mu = 0.001; alpha = 0; lambda_nu = 0.1; G1 = inf; simNo = 2;
% ------ numsim3: linear Kevin-Voigt model ------
% ------ numsim4: Neo-Hookean Kevin-Voigt model ------
% G = 2970; mu = 0.01; alpha = 0; lambda_nu = 0; G1 = inf; simNo = 4; model = 'neoHook';
% ------ numsim5: Fung Kevin-Voigt model ------
% G = 461; mu = 0.01; alpha = 10; lambda_nu = 0; G1 = inf; simNo = 5; model = 'fung';   
% ------ numsim6: Neo-Hookean Kevin-Voigt w/ nonlinear viscosity ------
% G = 2.97e3; mu = 0.01; alpha = 0; lambda_nu = 1; G1 = inf; simNo = 6;   
% ------ numsim7: Fung Kevin-Voigt w/ nonlinear viscosity ------
% G = 2.97e3; mu = 0.01; alpha = 0.1; lambda_nu = 1; G1 = inf; simNo = 7; 
if exist('G','var') == 0, G = 0; end
if exist('mu','var') == 0, mu = 0; end
if exist('alpha','var') == 0, alpha = 0; end
if exist('lambda_nu','var') == 0, lambda_nu = 0; end
if exist('G1','var') == 0, G1 = inf; end

matprop = struct('G',G,'mu',mu,'alpha',alpha,'lambda_nu',lambda_nu,'G1',G1);

% ====== Other physical parameters ======
% Load parameters: 
Pmt = IMRcall_parameters(1,1,1,1); % The first input as "1" doesn't matter.
P_inf = Pmt(19); T_inf = Pmt(20); rho = Pmt(23); Uc = Pmt(24);
% % Comment codes below for consistency.
% P_inf = 101325; % (Pa) Atmospheric Pressure
% rho = 998.2; % Some paper use 1060; % JY!!! % (Kg/m^3) Material Density
% T_inf = 298.15; ST = 0.056; % JY!!! % 0.056; % (N/m) Liquid Surface Tension 
% Uc = sqrt(P_inf/rho); % Characteristic velocity
Tgrad = Pmt(25); Cgrad = Pmt(26); Tmgrad = Pmt(27); ST = Pmt(28);
% Tgrad = 1; %1 Temperature transfer bt bubble and material
% Cgrad = 1; %1 Vapor-non-condensible gas diffusion
% Tmgrad = 0; %0 Off means cold liquid assumption


%% Initial partial pressure of the non-condensible gas

% Assign R0 and Req0
% R0 = 225e-6; 
if ICReqOrPinit == 0
    %Req0 = R0*ones(1,1)./[10:-1:7]; % % Req0 = 22.5e-6;  R0 = 10*Req0;
    P_guessList = (P_inf+2*ST./Req0-Pvsat(T_inf)).*((Req0/R0).^3);
else
% % Alterinatively, we can assign initial values of bubble inside pressure
% P_guessList = [125,150,175,200,225,250,275,300,325,350]; 
end

for expt = 1:length(P_guessList)
    
    %general folder % cd('/Users/yangjin/Documents/MATLAB/Franck/IMR-master');
    %cd('E:\Jin\Franck\IMR-master');
    
    P_guess = P_guessList(expt);
    
    %%
    % NT = 500; % Mesh points in bubble, resolution should be >=500
    NTM = 10; % Mesh points in the material, should be 10
    Pext_type = 'IC'; %'IC' for Flynn, 'ga' for gaussian bubble growth

    if strcmp(Pext_type,'IC')
        Pext_Amp_Freq = [P_guess; 0];%[226; 0]; % Tune first number to match equi radii
    elseif strcmp(Pext_type,'ga')
        Pext_Amp_Freq = [P_guess; dt; tw;];
    end

    disptime = 0; % 1 = Displays time to complete simulation
    Dim = 1;  % 1 = displays results in dimensional form
    comp = 1; % 0 uses Rayleigh-Plesset, 1 uses Keller-Miksis

    %%
    if strcmp(Pext_type,'ga')
        Rmax = R0;
        R0 = eqR;
    end
     
    %Variables:
    %t2 = simulation time
    %R2 = simulation radius
    %R2dot = velocity of bubble wall
    %P = bubble pressure
    %S = stress integral
    %T = temperature inside bubble
    %C = relative concentration of vapor
    %Tm = temperature of wall
    
    REqNew = Req0(expt)/R0; REqIter = 0;
    % while (abs(REqNew-REqOld)/REqOld > 1e-3) && (REqIter < 100)
        REqOld = REqNew; REqIter = REqIter+1; disptime = 0;
        % ====== IMR solver part ======
        [t2, R2, R2dot, P, S, T, C, Tm, tdel, Tdel, Cdel] = funIMRsolver( ...
            model, matprop, tspan, R0, NT, NTM, Pext_type, Pext_Amp_Freq, ...
            disptime, Tgrad, Tmgrad, Cgrad, Dim, comp, REqOld,IMRsolver_RelTolX);
       
        % REqNew = mean(R2(end-20:end))/R0;
        % =============================
    % end
    % toc

    %% Plot
    figure; plot(t2,R2); set(gca,'fontsize',18); axis([0,tspan,0,R0]);
    xlabel('Time(s)','interpreter','latex'); ylabel('Simulated bubble radius $R$(m)','interpreter','latex');
    a=gca; a.TickLabelInterpreter = 'latex';
    % figure; plot(t2,P); axis([0,tspan,0,1e6])

    %% Fit R2-t2 data and generate data
    % cd(['/Users/yangjin/Documents/MATLAB/Franck/IMR-master/data/numsim/',num2str(simNo),'/',num2str(expt),'/']);
    tempDirectory = ['./data/numsim/',num2str(simNo),'/',num2str(expt),'/'];
    mkdir(tempDirectory); addpath( (tempDirectory)); cd(tempDirectory);
    
    Rfit_pp = interp1(t2',R2','pchip','pp');
    
    deltat = 1/270000/TimeRes; t = 0:deltat:tspan; Rnew = ppval(Rfit_pp,t);
               filename = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes)]; save([filename '.mat'], 't', 'Rnew');
    
%     TimeRes=0; t = t2'; Rnew = R2(1:1:length(R2))'; 
%                filename = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes)]; save([filename '.mat'], 't', 'Rnew');
%     TimeRes=1; deltat = 1/270000/TimeRes; t = 0:deltat:tspan;  Rnew = ppval(Rfit_pp,t);
%                filename = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes)]; save([filename '.mat'], 't', 'Rnew');
%     TimeRes=2; deltat = 1/270000/TimeRes; t = 0:deltat:tspan; Rnew = ppval(Rfit_pp,t);
%                filename = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes)]; save([filename '.mat'], 't', 'Rnew');
%     TimeRes=3; deltat = 1/270000/TimeRes; t = 0:deltat:tspan; Rnew = ppval(Rfit_pp,t);
%                filename = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes)]; save([filename '.mat'], 't', 'Rnew');           
%     TimeRes=5; deltat = 1/270000/TimeRes; t = 0:deltat:tspan; Rnew = ppval(Rfit_pp,t);
%                filename = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes)]; save([filename '.mat'], 't', 'Rnew'); 
%     TimeRes=10; deltat = 1/270000/TimeRes; t = 0:deltat:tspan; Rnew = ppval(Rfit_pp,t);
%                filename = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes)]; save([filename '.mat'], 't', 'Rnew');
       
    % t = t2'; Rnew = R2(1:1:length(R2))'; filename = ['numsim',num2str(simNo),'_RofTdata_100']; save([filename '.mat'], 't', 'Rnew');
    % t = t2'; t=t(1:2:length(t)); Rnew = R2(1:2:length(R2))'; filename = ['numsim',num2str(simNo),'_RofTdata_50']; save([filename '.mat'], 't', 'Rnew');
    % t = t2'; t=t(1:10:length(t)); Rnew = R2(1:10:length(R2))'; filename = ['numsim',num2str(simNo),'_RofTdata_10']; save([filename '.mat'], 't', 'Rnew');
    % t = t2'; t=t(1:20:length(t)); Rnew = R2(1:20:length(R2))'; filename = ['numsim',num2str(simNo),'_RofTdata_5']; save([filename '.mat'], 't', 'Rnew');

    % Run FitRtCurve
    if simNo == 4
        file_name = ['numsim',num2str(simNo),'_Exact_G',num2str(G),'_mu',num2str(mu),'_NT',num2str(NT),'_logRelTol',num2str(-log10(IMRsolver_RelTolX))];
    elseif simNo == 5 
        file_name = ['numsim',num2str(simNo),'_Exact_G',num2str(G),'_alpha',num2str(alpha),'_mu',num2str(mu),'_NT',num2str(NT),'_logRelTol',num2str(-log10(IMRsolver_RelTolX))];
    elseif simNo == 7
        file_name = ['numsim',num2str(simNo),'_Exact_G',num2str(G),'_alpha',num2str(alpha),'_mu',num2str(mu),'_lambda-nu',num2str(lambda_nu),'_NT',num2str(NT),'_logRelTol',num2str(-log10(IMRsolver_RelTolX))];
    end
        
    save([file_name '.mat'],'t2','R2','R2dot','C','T','P','matprop');

    cd('../../../../');
     
end


