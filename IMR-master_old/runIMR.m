%% Runfile for inertial microcavitation rheometry

% Authors:
% Jon Estrada
% jonathan_estrada@alumni.brown.edu
% Brown Solid Mechanics, PhD '17
% Carlos Barajas
% carlobar@umich.edu
% Umich Mechanical Engineering BS '16

clear; close all;
warning('off','all')

%% Set up the file paths


%For CCV, file path is different
%general folder
addpath(genpath('/gpfs/data/cf5/jbestrad/FL_Cav_Data/'));
%data folder
fp = '/gpfs/data/cf5/jbestrad/FL_Cav_Data/160420/11kPa_PA/';

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
load([fp 'RofTdata.mat']);

%There are 4 choices of model:
%linkv (linear Kelvin-Voigt)
%neoHook (neo-Hookean Kelvin-Voigt)
%sls (Standard Linear Solid)
%nhzen (neo-Hookean Standard Solid)
model = 'neoHook';
savename = '170821_sweep';
%'170727_softsweep_coarse';%'170413_collagenKVcoarse';

%Parallel pool creation
% curCluster = parcluster('local');
% curCluster.NumWorkers = 8;
% saveProfile(curCluster);
% pool = parpool(8);

%Set time duration for the simulation (s)
tspan = 2.25E-4;%1.3E-4;
allRmax = max(Rnew,[],2);

%Set the index choice of experiments
%expts = [12 14:19]; %water
%expts = [2,3,4,5,8,10,14,15,16,18,20,23,24]; %collagen
expts = 6;

%Set the range of G and mu (most conveniently done as powers of 10)
G_ooms = 1:0.2:5;%3.0:0.1:4.0; %soft PA %3.65:0.05:4.15 stiff PA
mu_ooms = -1.4:0.05:-0.9;%[-inf -4:0.25:-1.25 -1.05];%[-1.65:0.05:-0.9]%-2.25:0.25:-0.5%[-inf -1.65:0.05:-0.9];
G1_ooms = inf;
%Note, mu/G1 should go to 0 in the limit, so G1 for K-V models should be infinite

%% Run IMR
for expt = expts
    mkdir([fp num2str(expt)]);
    cd([fp num2str(expt)]);
    
    %%
    soln_mx = cell(length(G_ooms),length(mu_ooms),length(G1_ooms));
    P_inf = 101325; % (Pa) Atmospheric Pressure
    rho = 998.2; % (Kg/m^3) Material Density
    Uc = sqrt(P_inf/rho); % Characteristic velocity
    
    Tgrad = 1; %1 Temperature transfer bt bubble and material
    Cgrad = 1;%1;  %1 Vapor-non-condensible gas diffusion
    Tmgrad = 0; %0 Off means cold liquid assumption
       
    %%
    for k = 1:length(G1_ooms)
        for j= 1:length(mu_ooms)
            tic
            parfor i=1:length(G_ooms)
                soln_mx{i,j,k} = struct('G',10^G_ooms(i),'mu',10^mu_ooms(j),'G1',10^G1_ooms(k),...
                    'tcs_star',[],'Rratios',[],'tmaxs_star',[],'t2',[],'R2',[],'U',[],'P',[],'J',[],'T',[],'C',[],...
                    'Cdel',[],'tdel',[],'Tdel',[]);
                
                G = soln_mx{i,j,k}.G;
                mu = soln_mx{i,j,k}.mu;
                G1 = soln_mx{i,j,k}.G1;
                
                [R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t); %need to get the inital radius from a fit (R)
                eqR = median(Rnew(expt,62:end))*10^-6;
                R_eq = eqR/R0;
                
                % Bisection method to get initial partial pressure of the non-condensible gas
                
                % Initial guess (purely empirical) and range (e.g. ~226 for
                % 11kPa, 92 for 1.3kPa, 20.6 for water)
                P_guess = 226;
                P = [8, 400];
                
                R_eq_lims = [IMRCalc_Req(R0,Tgrad,Cgrad,P(1),G,G1,mu), IMRCalc_Req(R0,Tgrad,Cgrad,P(2),G,G1,mu)];
                R_eqf  = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
                error = 1;
                while abs(error) > 0.00001
                    if R_eq > R_eqf
                        P(1) = P_guess;
                        R_eq_lims(1) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
                    else
                        P(2) = P_guess;
                        R_eq_lims(2) = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
                    end
                    P_guess = mean(P);
                    R_eqf  = IMRCalc_Req(R0,Tgrad,Cgrad,P_guess,G,G1,mu);
                    error = abs(R_eq-R_eqf);
                end
                
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
                tic
                
                %Variables:
                %t2 = simulation time
                %R2 = simulation radius
                %R2dot = velocity of bubble wall
                %P = bubble pressure
                %S = stress integral
                %T = temperature inside bubble
                %C = relative concentration of vapor
                %Tm = temperature of wall
                [t2, R2, R2dot, P, S, T, C, Tm, tdel, Tdel, Cdel] = IMRsolver...
                    (model, G, G1, mu, tspan, R0, NT, NTM, ...
                    Pext_type, Pext_Amp_Freq , disptime, Tgrad, Tmgrad, Cgrad, Dim, comp);
                toc
                
                %%
                Rdiff = diff(R2);
                inc_idx = find(Rdiff>0);
                diffzeros_idx = find(Rdiff(1:(end-1)).*Rdiff(2:(end))<0)+1;
                
                localmax_idx = [1; diffzeros_idx(2:2:end)];
                localmin_idx = diffzeros_idx(1:2:end);
                Rratios = R2(localmax_idx)/R2(1);
                tcs = t2(localmin_idx);
                tcs_star = tcs*Uc/R0;
                tmaxs_star = t2(localmax_idx)*Uc/R0;
                %tpeaksall_star = [0; reshape([tcs_star,tmaxs_star(2:end)],size(diffzeros_idx))];
                
                soln_mx{i,j,k}.tcs_star = tcs_star;
                soln_mx{i,j,k}.Rratios = Rratios;
                soln_mx{i,j,k}.tmaxs_star = tmaxs_star;
                soln_mx{i,j,k}.t2 = t2;
                soln_mx{i,j,k}.R2 = R2;
                
                soln_mx{i,j,k}.U = U;
                soln_mx{i,j,k}.P = P;
                soln_mx{i,j,k}.J = J;
                %soln_mx{i,j,k}.T = T;
                %soln_mx{i,j,k}.C = C;
                soln_mx{i,j,k}.Cdel = Cdel;
                soln_mx{i,j,k}.tdel = tdel;
                soln_mx{i,j,k}.Tdel = Tdel;
                [i j k]
            end
            toc
            %If you want to save intermediate files during runs, uncomment this line
            %save([savename 'G_' num2str(G_ooms(i)) '.mat'],'soln_mx','G_ooms','mu_ooms','Rnew','t');
        end
    end
    save([savename '.mat'],'soln_mx','G_ooms','mu_ooms','Rnew','t');
    cd(fp)
end

%%
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
%                 [R2max idx] = max(soln_mx{i,j,k}.R2);
%
%                 hold on;
%
%                 t2 = soln_mx{i,j,k}.t2-soln_mx{i,j,k}.t2(idx);
%                 R2 = soln_mx{i,j,k}.R2;
%                 %plot(t2,R2,'Color',[(i-1)/sG1(1) (k-1)/sG(1) (j-1)/sM(1)]);
%
%                 [R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t);
%                 t2exp=t-t0;
%                 R2exp=Rnew(expt,:)*1E-6;
%                 %plot(t2exp,R2exp, ' *');
%
%                 distRmat = bsxfun(@minus,R2,R2exp);
%                 disttmat = bsxfun(@minus,t2,t2exp);
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
% %                 [maxRNHZ(n) idx(n)] = max(soln_mx{i,j,n}.R2);
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
% R2 = soln_mx{i,j,k}.R2;
% [R2max idx] = max(soln_mx{i,j,k}.R2);
% t2 = soln_mx{i,j,k}.t2-soln_mx{i,j,k}.t2(idx);
% %plot(t2,R2);%,'Color','red');
% %plot(t2/max(R2)*10,R2/max(R2)./(1-((1-0.1301).*(exp(-t2/max(R2)*10/1.15))+0.1301)));%,'Color','red');
% %hold on; scatter(t2exp,R2exp,16,[0 0 0],'filled');
%
% [R2max idx] = max(soln_mx{i,j,k}.R2);
%
% %save('dataprocoutputs.mat','t2' , 'R2','U','P','Z', 'T','C', 'Tm','tdel','Tdel','Cdel');
% end