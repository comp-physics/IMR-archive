%% Runfile for inertial microcavitation rheometry

% Authors:
% Jon Estrada
% jonathan_estrada@alumni.brown.edu
% Brown Solid Mechanics, PhD '17
% Carlos Barajas
% carlobar@umich.edu
% Umich Mechanical Engineering BS '16

% clear; close all;
warning('off','all')

%% Set up the file paths

%For CCV, file path is different
%general folder
% cd('/Users/yangjin/Documents/MATLAB/Franck/IMR-master');
% cd('E:\Jin\Franck\IMR-master');
% addpath(genpath('./data/')); % addpath(genpath('/gpfs/data/cf5/jbestrad/FL_Cav_Data/'));

%data folder
% fp = ['/Users/yangjin/Documents/MATLAB/Franck/IMR-master/data/numsim/',num2str(simNo),'/'];
% fp = './data/11kPa_PA/'; % fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/160420/11kPa_PA/';
% fp = './data/1pt3kPa_PA/';
% fp = './data/collagen/1/'

% simNo = 4; SampleNo = 100; GaussNo = 5;
% fp = ['./data/numsim/',num2str(simNo),'/'];

%fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/160420/water/';
%fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/160511/collagen/1/';
%fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/170403/Collagen/';
%fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/170411/Collagen/20170302p1/';

%For Local comp usage (i.e. single runs)
%general folder
%addpath(genpath('V:\data\jbestrad\FL_Cav_Data\'));
%data folder
%fp = 'V:\data\jbestrad\FL_Cav_Data\160420\water\';
%fp = 'V:\data\jbestrad\FL_Cav_Data\170411\Collagen\20170302p1/';
%fp = 'C:\Users\Jon\Brown Google Drive\RESEARCH DATA\FL Cav Data\160420\11kPa_PA\';

%Load the file RofTdata.mat, which contains vars Rnew and t
%Rnew has size [num_expts num_video_frames]
% load([fp 'RofTdata.mat']);
% tempfilename  =  ['numsim',num2str(simNo),'_RofTdata_100.mat'];
% load([fp tempfilename]); Rnew = Rnew*1e6;
 

%There are 4 choices of model:
%linkv (linear Kelvin-Voigt)
%neoHook (neo-Hookean Kelvin-Voigt)
%sls (Standard Linear Solid)
%nhzen (neo-Hookean Standard Solid)
% model = 'neoHook'; % JY!!!
% savename = '190220_sweep_JY_';
%'170727_softsweep_coarse';%'170413_collagenKVcoarse';

%Parallel pool creation
% curCluster = parcluster('local');
% curCluster.NumWorkers = 8;
% saveProfile(curCluster);
% pool = parpool(8);

%Set time duration for the simulation (s)

%Set the index choice of experiments
%expts = [12 14:19]; %water
%expts = [2,3,4,5,8,10,14,15,16,18,20,23,24]; %collagen
timestep = 1e-7; 
% PickNo = 7; % GaussNo = 5; SampleNo = 1; % numsim4

%Set the range of G and mu (most conveniently done as powers of 10)
% temp_G_ooms_List = 2970; G_ooms = log10(temp_G_ooms_List); % G_ooms = log10(7.69e3); % JY!!! % 1:0.2:5;%3.0:0.1:4.0; %soft PA %3.65:0.05:4.15 stiff PA
% temp_mu_ooms_List = [0.33, 0.101, 0.2, 0.01, 0.001, 0.1 ]; mu_ooms = log10(temp_mu_ooms_List); % 
% mu_ooms = log10(0.01); % JY!!! 0.0533; %-1.4:0.05:-0.9;%[-inf -4:0.25:-1.25 -1.05];%[-1.65:0.05:-0.9]%-2.25:0.25:-0.5%[-inf -1.65:0.05:-0.9];
% temp_alpha_ooms_List = [0.001,0.01,0.10,0.33, 0]; alpha_ooms = log10(temp_alpha_ooms_List);
% alpha_ooms = -Inf;
% temp_lambda_nu_oooms_List =[0, 3.3e-3, 1e-2, 3.3e-2, 1e-1, 1e0, 1e1, 1e2]; lambda_nu_ooms = log10(temp_lambda_nu_oooms_List);
% lambda_nu_ooms = -Inf;
% G1_ooms = inf;
%Note, mu/G1 should go to 0 in the limit, so G1 for K-V models should be infinite


%% Run IMR
matPropVarList = cell(length(expts),1);

for tempi = 1:length(expts)
    expt = expts(tempi);
    % mkdir([fp num2str(expt)]);
    % cd([fp num2str(expt)]);
    
    % JY!!!
    fp = ['./data/numsim/',num2str(simNo),'/',num2str(expt),'/'];
    tempfilename = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes),'.mat'];
    load([fp tempfilename]); Rnew = Rnew*1e6;
    
    tempfilename = ['Rfit_numsim',num2str(simNo),'_exp',num2str(expt),'_',num2str(TimeRes),'.mat']; 
    load([fp tempfilename]);  
    
    % Curve Rfit-t2 is the one we do fittings.
    
    %% Set time duration for the fittings (s)
    %tspan = 2.25E-4; %1.3E-4;
    %allRmax = max(Rnew,[],2);
    
    % mkdir([fp num2str(expt)]);
    %cd([fp]);
    [pksUp_P, tspanLocsMinR] = findpeaks( -(Rfit), t2, 'MinPeakDist', 1e-6 );
    [pksMaxR, tspanLocsMaxR] = findpeaks( Rfit, t2, 'MinPeakDist', 1e-6 );
    pksMaxR = [Rfit(1), pksMaxR]; tspanLocsMaxR = [t2(1), tspanLocsMaxR];
    
    %%
    % soln_mx = cell(length(mu_ooms),length(alpha_ooms),length(lambda_nu_ooms));
    % JY!!! % soln_mx = cell(length(G_ooms),length(mu_ooms),length(G1_ooms));
    % P_inf = 101325; % (Pa) Atmospheric Pressure
    % rho = 998.2; % JY!!! % rho = 998.2; % (Kg/m^3) Material Density
    % Uc = sqrt(P_inf/rho); % Characteristic velocity
    
    % Load parameters: 
    Pmt = IMRcall_parameters(1,1,1,1); % The first input as "1" doesn't matter.
    P_inf = Pmt(19); T_inf = Pmt(20); rho = Pmt(23); Uc = Pmt(24); ST = Pmt(28);
    
    Tgrad = Pmt(25); Cgrad = Pmt(26); Tmgrad = Pmt(27);
       
    %%
    R0 = Rfit(1); %[R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t); %need to get the inital radius from a fit (R)
    eqR = median(Rnew(expt, end-20:end)) *10^-6; % Solve for equilibirium R_inf
    R_eq = eqR/R0;
                  
    %% Bisection method to get initial partial pressure of the non-condensible gas
    P_guess = (P_inf+2*ST./eqR-Pvsat(T_inf)).*((eqR/R0).^3);
   
%     tic
%     % Initial guess (purely empirical) and range (e.g. ~226 for
%     % 11kPa, 92 for 1.3kPa, 20.6 for water)
%     % JY!!! which for collagen?
%     P_guess = 200;
%     P = [8, 400];
%     
%     % temp G, mu, G1
%     G = GInit; mu = muInit; G1 = Inf;
%     
%     R_eq_lims = [IMRCalc_Req(R0,Tgrad,Cgrad,P(1),G,G1,mu), IMRCalc_Req(R0,Tgrad,Cgrad,P(2),G,G1,mu)];
%     R_eqf = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%     error = 1; BisectionIter = 1;
%     while abs(error) > 0.00001 && BisectionIter<1e6
%         if R_eq > R_eqf
%             P(1) = P_guess;
%             R_eq_lims(1) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%         else
%             P(2) = P_guess;
%             R_eq_lims(2) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%         end
%         P_guess = mean(P);
%         R_eqf  = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%         error = abs(R_eq-R_eqf);
%         BisectionIter = BisectionIter+1;
%     end
%     %toc
%     if abs(error) > 1e-5 || P_guess < 8 || P_guess > 400
%         disp(['Fail to find initial partial pressure of gas in ', num2str(toc),' seconds.']);
%     else
%         disp(['P_guess = ',num2str(P_guess)]);
%         disp(['Find initial partial pressure of gas in ', num2str(toc),' seconds.']);
%     end
                
    %%
    NT = 500; % Mesh points in bubble, resolution should be >=500
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
    
    %% Variables:
    %t2_num = simulation time
    %R2_num = simulation radius
    %R2dot = velocity of bubble wall
    %P = bubble pressure
    %S = stress integral
    %T = temperature inside bubble
    %C = relative concentration of vapor
    %Tm = temperature of wall
    
    %% Solver
    
    %alpha_ooms = -Inf; lambda_nu_ooms = Inf; G1_ooms = Inf; 
    % Initial guess
    %matPropVarList0 = [2, -2];
    
    % Upper and lower bound
    %lb = log10([500, 5e-3]);
    %ub = log10([1.5e5, 0.1]);
    
    % options = optimset('Display','off' );
    options = optimset('Display','iter','PlotFcns','optimplotfval','TolX',1e-4);
    
%     close all;  
%     [LSQErr,LSQErrNo] = funIMRLSQErr_Seg(model,matPropVar0,tspanLocsMinR(1:PickNo),pksMaxR(1:PickNo),tspanLocsMaxR(1:PickNo),NT/10,NTM,...
%     Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp,eqR,t2,Rfit);
    [row1,col1] = find(t2<tspanLocsMinR(PickNo)); 
    t2_exp = t2(1:max(col1));
    Rfit_exp = Rfit(1:max(col1));

    NelderMead_SegOrNot = 0;
    
    tic;
    if NelderMead_SegOrNot == 1
        [matPropVar_best,fval,exitflag,output] = fminsearchbnd( @(matPropVar) ...
           funIMRLSQErr_Seg(model,matPropVar,tspanLocsMinR(1:PickNo),pksMaxR(1:PickNo),tspanLocsMaxR(1:PickNo),NT,NTM,...
           Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp,eqR,t2_exp,Rfit_exp), ...
           matPropVar0, lb,ub,options);
    else
        [matPropVar_best,fval,exitflag,output] = fminsearchbnd( @(matPropVar) ...
           funIMRLSQErr(model,matPropVar,tspanLocsMinR(PickNo),pksMaxR(1),tspanLocsMaxR(1),NT,NTM,...
           Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp,eqR,t2_exp,Rfit_exp), ...
           matPropVar0, lb,ub,options);
    end
    toc
  
    disp(['expt = ',num2str(expt),'; G_best = ',num2str(10^matPropVar_best(1)),'; mu_best = ', num2str(10^matPropVar_best(2))]);
            
    
    matPropVarList{expt} = matPropVar_best;
    
    
end

