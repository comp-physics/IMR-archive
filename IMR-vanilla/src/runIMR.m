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
%cd('E:\Jin\Franck\IMR-master');
%addpath(genpath('./data/')); % addpath(genpath('/gpfs/data/cf5/jbestrad/FL_Cav_Data/'));

%data folder
% fp = ['/Users/yangjin/Documents/MATLAB/Franck/IMR-master/data/numsim/',num2str(simNo),'/'];
% fp = './data/11kPa_PA/'; % fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/160420/11kPa_PA/';
% fp = './data/1pt3kPa_PA/';
% fp = './data/collagen/1/'

%simNo = 4; SampleNo = 100; GaussNo = 5;
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
%model = 'neoHook'; % JY!!!
%savename = '190220_sweep_JY_';
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
%expts = 1; timestep = 1e-7; PickNo = 7; % GaussNo = 5; SampleNo = 1; % numsim4

%% Set the range of G and mu (most conveniently done as powers of 10)
%temp_G_ooms_List = []; G_ooms = log10(temp_G_ooms_List);  
%G_ooms = 3.3:0.02:3.7; %3.0:0.1:4.0; %soft PA %3.65:0.05:4.15 stiff PA
%temp_mu_ooms_List = []; mu_ooms = log10(temp_mu_ooms_List);  
%mu_ooms = -2.2:0.02:-1.8; % JY!!! 0.0533; %-1.4:0.05:-0.9;%[-inf -4:0.25:-1.25 -1.05];%[-1.65:0.05:-0.9]%-2.25:0.25:-0.5%[-inf -1.65:0.05:-0.9];
%temp_alpha_ooms_List = [0.001,0.01,0.10,0.33,0]; alpha_ooms = log10(temp_alpha_ooms_List);
%alpha_ooms = -Inf;
%temp_lambda_nu_oooms_List =[0, 3.3e-3, 1e-2, 3.3e-2, 1e-1, 1e0, 1e1, 1e2]; lambda_nu_ooms = log10(temp_lambda_nu_oooms_List);
%lambda_nu_ooms = -Inf;
%G1_ooms = inf;
%Note, mu/G1 should go to 0 in the limit, so G1 for K-V models should be infinite


%% Run IMR
for expt = expts
    
    close all;
    
    % mkdir([fp num2str(expt)]);
    % cd([fp num2str(expt)]);
    
    % JY!!!
    %fp = ['E:\Jin\Franck\IMR-master\data\numsim\',num2str(simNo),'\',num2str(expt),'\'];
    %fp = ['./data/numsim/',num2str(simNo),'/',num2str(expt),'/'];
    fp = [folderNamePrefix,num2str(expt),'/'];
    %tempfilename = ['numsim',num2str(simNo),'_RofTdata_',num2str(TimeRes),'.mat'];
    tempfilename  =  [fileNamePrefix,'.mat'];
    load([fp tempfilename]); Rnew = Rnew*1e6;
    
    tempfilename = ['Rfit_numsim',num2str(simNo),'_exp',num2str(expt),'_',num2str(TimeRes),'.mat']; 
    load([fp tempfilename]);  
    
    % Curve Rfit-t2 is the one we do fittings.
    
    %% Set time duration for the simulation (s)
    %tspan = 2.25E-4; %1.3E-4;
    %allRmax = max(Rnew,[],2);
    
    % mkdir([fp num2str(expt)]);
    %cd([fp]);
    [pksUp_P, tspanLocsMinR] = findpeaks( -(Rfit), t2, 'MinPeakDist', 1e-6 );
    [pksMaxR, tspanLocsMaxR] = findpeaks( Rfit, t2, 'MinPeakDist', 1e-6 );
    pksMaxR = [Rfit(1), pksMaxR]; tspanLocsMaxR = [t2(1), tspanLocsMaxR];
    
    %%
    % JY!!!%soln_mx = cell(length(mu_ooms),length(alpha_ooms),length(lambda_nu_ooms));
    soln_mx = cell(length(G_ooms),length(mu_ooms),length(G1_ooms));
    LSQErr = cell(length(G_ooms),length(mu_ooms),length(G1_ooms)); LSQErrNo = LSQErr;
    % P_inf = 101325; % (Pa) Atmospheric Pressure
    % rho = 998.2; % JY!!! % rho = 998.2; % (Kg/m^3) Material Density
    % Uc = sqrt(P_inf/rho); % Characteristic velocity
    % 
    % Tgrad = 1; %1 Temperature transfer bt bubble and material
    % Cgrad = 1; %1;  %1 Vapor-non-condensible gas diffusion
    % Tmgrad = 0; %0 Off means cold liquid assumption
    
    % Load parameters: 
    Pmt = IMRcall_parameters(1,1,1,1); % The first input as "1" doesn't matter.
    P_inf = Pmt(19); T_inf = Pmt(20); rho = Pmt(23); Uc = Pmt(24); ST = Pmt(28);
    
    Tgrad = Pmt(25); Cgrad = Pmt(26); Tmgrad = Pmt(27);
       
    %%
    for k = 1:length(G1_ooms)
        for j= 1:length(mu_ooms)
            for i= 1:length(G_ooms)
                
                %[i,j,k]
                
                soln_mx{i,j,k} = struct('G',10^G_ooms(i),'mu',10^mu_ooms(j),'G1',10^G1_ooms(k),...
                    'tcs_star',[],'Rratios',[],'tmaxs_star',[],'t2_num',[],'R2_num',[],'U',[],'P',[],'J',[],'T',[],'C',[],...
                    'Cdel',[],'tdel',[],'Tdel',[]);
                
                G = soln_mx{i,j,k}.G;
                mu = soln_mx{i,j,k}.mu;
                G1 = soln_mx{i,j,k}.G1;
                alpha = 0; %soln_mx{i,j,k}.alpha;
                lambda_nu = 0; %soln_mx{i,j,k}.lambda_nu;
                
                matprop = struct('G',G,'mu',mu,'alpha',alpha,'lambda_nu',lambda_nu,'G1',G1);
                 
                %%
%                 R0 = Rfit(1); %[R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t); %need to get the inital radius from a fit (R)
%                 eqR = % median(Rnew(expt, end-20:end)) *10^-6; % Solve for equilibirium R_inf
%                 R_eq = eqR/R0;
 R0=Rfit(1); t0=0; 
        %[R0,t0] = calcR0(Rnew(1,:)*1E-6,t); %need to get the inital radius from a fit (R)
        eqR = mean(Rnew(1, end-20:end )) *10^-6; % Solve for equilibirium R_inf
        R_eq = eqR/R0;
                
                %% Bisection method to get initial partial pressure of the non-condensible gas
                P_guess = (P_inf+2*ST./eqR-Pvsat(T_inf)).*((eqR/R0).^3);
%                 
%                 tic
%                 % Initial guess (purely empirical) and range (e.g. ~226 for
%                 % 11kPa, 92 for 1.3kPa, 20.6 for water)
%                 % JY!!! which for collagen?
%                 P_guess = 200;
%                 P = [8, 400];
%                 
%                 R_eq_lims = [IMRCalc_Req(R0,Tgrad,Cgrad,P(1),G,G1,mu), IMRCalc_Req(R0,Tgrad,Cgrad,P(2),G,G1,mu)];
%                 R_eqf = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%                 error = 1; BisectionIter = 1;
%                 while abs(error) > 0.00001 && BisectionIter<1e6
%                     if R_eq > R_eqf
%                         P(1) = P_guess;
%                         R_eq_lims(1) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%                     else
%                         P(2) = P_guess;
%                         R_eq_lims(2) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%                     end
%                     P_guess = mean(P);
%                     R_eqf  = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
%                     error = abs(R_eq-R_eqf);
%                     BisectionIter = BisectionIter+1;
%                 end
%                 toc
%                 if abs(error) < 1e-5
%                     disp(['Find initial partial pressure of gas in ', num2str(toc),' seconds.']);
%                 else
%                     disp(['Fail to find initial partial pressure of gas in ', num2str(toc),' seconds.']);
%                 end
                
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
                %tic;
                IMRsolver_RelTolX = 1e-7;
                
                [row1,col1] = find(t2<tspanLocsMinR(PickNo));
                t2_exp = t2(1:max(col1));
                Rfit_exp = Rfit(1:max(col1));
                
                [t2_num, R2_num, R2dot, P, S, T, C, Tm, tdel, Tdel, Cdel] = funIMRsolver ...
                   (model, matprop, tspanLocsMinR(PickNo), R0, NT, NTM, ...
                   Pext_type, Pext_Amp_Freq, disptime, Tgrad, Tmgrad, Cgrad, Dim, comp, R_eq, IMRsolver_RelTolX);
                
                Rfit_pp = interp1(t2_num,R2_num,'pchip','pp'); Rfit_num = ppval(Rfit_pp,t2_exp);
                 
               
                %toc;
               
                %file_name = ['simnum',num2str(simNo),'_NHKV_G',num2str(G),'_mu',num2str(k),'_a',num2str(j),'_ln',num2str(i)];
     
                %cd('/Users/yangjin/Documents/MATLAB/Franck/IMR-master');
                %close all; plot_Fung_JY; 
                %figure(2); 
                % % fig_name = ['fig_exp',num2str(expt),'_Fung_mu',num2str(k),'_a',num2str(j),'_ln',num2str(i),'_pkloc.fig'];
                %fig_name = ['fig_simnum',num2str(simNo),'_NHKV_G',num2str(G),'_mu',num2str(k),'_a',num2str(j),'_ln',num2str(i),'_pkloc.fig'];
                %savefig(fig_name);
                 
                %%
                Rdiff = diff(R2_num);
                inc_idx = find(Rdiff>0);
                diffzeros_idx = find(Rdiff(1:(end-1)).*Rdiff(2:(end))<0)+1;
                
                localmax_idx = [1; diffzeros_idx(2:2:end)];
                localmin_idx = diffzeros_idx(1:2:end);
                Rratios = R2_num(localmax_idx)/R2_num(1);
                tcs = t2_num(localmin_idx);
                tcs_star = tcs*Uc/R0;
                tmaxs_star = t2_num(localmax_idx)*Uc/R0;
                %tpeaksall_star = [0; reshape([tcs_star,tmaxs_star(2:end)],size(diffzeros_idx))];
                
                soln_mx{i,j,k}.tcs_star = tcs_star;
                soln_mx{i,j,k}.Rratios = Rratios;
                soln_mx{i,j,k}.tmaxs_star = tmaxs_star;
                soln_mx{i,j,k}.t2_num = t2_num;
                soln_mx{i,j,k}.R2_num = R2_num;
                
                % soln_mx{i,j,k}.U = U;
                soln_mx{i,j,k}.P = P;
                % soln_mx{i,j,k}.J = J;
                soln_mx{i,j,k}.T = T;
                soln_mx{i,j,k}.C = C;
                soln_mx{i,j,k}.Cdel = Cdel;
                soln_mx{i,j,k}.tdel = tdel;
                soln_mx{i,j,k}.Tdel = Tdel;
                  
                LSQErr{i,j,k} = norm( Rfit_exp-Rfit_num ,2 );
                LSQErrNo{i,j,k} = length(Rfit_exp);
                
                % figure,plot(t2_exp,Rfit_exp)
                % hold on; plot(t2_exp,Rfit_num)
                 
            end
             
             
            %If you want to save intermediate files during runs, uncomment this line
            % save([savename 'G_' num2str(G_ooms(i)) '.mat'],'soln_mx','G_ooms','mu_ooms','Rnew','t');
        end
    end
    
    %cd(['/Users/yangjin/Documents/MATLAB/Franck/IMR-master/data/numsim/',num2str(simNo),'/']);
     
    % save([file_name '.mat'],'t2_num','R2_num','C','T','P','tdel','Tdel','Cdel');
    % save([savename '.mat'],'soln_mx','G_ooms','mu_ooms','Rnew','t');
    % cd(fp)
    
    runMinLSQErr; % Find minimum LSQ error
    
    matPropVar_best = [G_Peak,mu_Peak];
    
    disp(['expt = ',num2str(expt),'; G_best = ',num2str(10^matPropVar_best(1)),'; mu_best = ', num2str(10^matPropVar_best(2))]);
    
    matPropVarList{expt} = matPropVar_best;

end

%cd('/Users/yangjin/Documents/MATLAB/Franck/IMR-master');

 

%% Old codes
% sz = size(soln_mx);
% sG1 = sz(1);
% try
%     sM = sz(2);
% catch
%     sM = 1;
% end
% try
%     sG = sz(3);
% catch
%     sG = 1;
% end
%%
%
% LSQ = cell(sG1,sM,sG);
% LSQminidx = [inf 0 0 0];
% try
%     for i=1:sG1
%        %figure(i)
%         for j=1:sM
%             for k=1:sG
%
%             %figure(10000*i+100*j+k)
%             %figure(10000*i+100+k)
%              %figure(505)
%
%                 [R2max idx] = max(soln_mx{i,j,k}.R2_num);
%
%                 hold on;
%
%                 t2_num = soln_mx{i,j,k}.t2_num-soln_mx{i,j,k}.t2_num(idx);
%                 R2_num = soln_mx{i,j,k}.R2_num;
%                 %plot(t2_num,R2_num,'Color',[(i-1)/sG1(1) (k-1)/sG(1) (j-1)/sM(1)]);
%
%                 [R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t);
%                 t2_numexp=t-t0;
%                 R2exp=Rnew(expt,:)*1E-6;
%                 %plot(t2_numexp,R2exp, ' *');
%
%                 distRmat = bsxfun(@minus,R2_num,R2exp);
%                 disttmat = bsxfun(@minus,t2_num,t2_numexp);
%                 distmat = sqrt(disttmat.^2+distRmat.^2);
%                 weight = ones(1,101);
%                 %weight([1:10 18:end]) = 0; %for fitting the first peak
%                 %weight([1:10 19:26]) = 0; %for ignoring the second peak
%                 weight([1:9 70:end]) = 0;
%
%
%                 LSQ{i,j,k} = nansum(weight.*min(distmat,[],1));
%
%                 if LSQ{i,j,k}<LSQminidx(1)
%                     LSQminidx = [LSQ{i,j,k} i j k];
%                 end
%             end
% %             for n=k:-1:1
% %                 [maxRNHZ(n) idx(n)] = max(soln_mx{i,j,n}.R2_num);
% %             end
%
%         end
%
%     end
% catch
% end
%
% i=LSQminidx(2); j=LSQminidx(3); k=LSQminidx(4);
% log10([soln_mx{i,j,k}.G soln_mx{i,j,k}.mu soln_mx{i,j,k}.G1])
% [k j i]
% %%
% %k=7; j=8; i=16;
% %figure(9999); hold on;
% R2_num = soln_mx{i,j,k}.R2_num;
% [R2_nummax idx] = max(soln_mx{i,j,k}.R2_num);
% t2_num = soln_mx{i,j,k}.t2_num-soln_mx{i,j,k}.t2_num(idx);
% %plot(t2_num,R2_num);%,'Color','red');
% %plot(t2_num/max(R2_num)*10,R2_num/max(R2_num)./(1-((1-0.1301).*(exp(-t2_num/max(R2_num)*10/1.15))+0.1301)));%,'Color','red');
% %hold on; scatter(t2_numexp,R2_numexp,16,[0 0 0],'filled');
%
% [R2_nummax idx] = max(soln_mx{i,j,k}.R2_num);
%
% %save('dataprocoutputs.mat','t2_num' , 'R2_num','U','P','Z', 'T','C', 'Tm','tdel','Tdel','Cdel');
% end