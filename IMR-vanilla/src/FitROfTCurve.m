%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting R-t curve to get DRDt and D^2R/Dt^2
% with experimental data or numerical simulation data
% Author: Jin Yang, jyang526@wisc.edu
% Date: 2019.01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%% Fitting R to get DR and DDR
%close all; warning('off','all');

% add data folder
%cd('E:\Jin\Franck\IMR-master'); % cd('/Users/yangjin/Documents/MATLAB/Franck/IMR-master');
%addpath(genpath('.\data\')); %data folder % fp = './data/11kPa_PA/'; 
 
% set numsim number
% expts = 1:1:expNo;
  
for tempi = 1:length(expts)

    expt=expts(tempi);   
    %fp = ['E:\Jin\Franck\IMR-master\data\numsim\',num2str(simNo),'\',num2str(exptNo),'\']; 
    % folderNamePrefix = ['./data/numsim/',num2str(simNo),'/'];
    fp = [folderNamePrefix,num2str(expt),'/'];

    %Load the file RofTdata.mat, which contains vars Rnew and t
    %Rnew has size [num_expts num_video_frames]
    cd([fp]);
    % fileNamePrefix = ['numsim',num2str(simNo)];
    load( [fileNamePrefix,'.mat'] ); t=t(1:1:length(t));  
 
    %Set time duration for the simulation (s)
    %tspan = 2.25E-4; %1.3E-4;
    allRmax = max(Rnew,[],2);

    % ********* For exp images ***********
    % expts = 1; 
    % tList = t; RList = Rnew(expts,:); 
    % expt = expts;
    % [R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t); %need to get the inital radius from a fit (R)
    % eqR = median(Rnew(expt,62:end))*10^-6; % Solve for equilibirium R_inf
    % R_eq = eqR/R0;
    % figure; plot(t,Rnew(expts,:)*1e-6,'rs');
    % ********* For numsim ************
    tList = t; RList = Rnew;
    R0 = Rnew(1); t0 = 0;
    eqR = mean(Rnew(end-20:end) ) ; % Solve for equilibirium R_inf
    R_eq = eqR/R0;

    plot(t,Rnew ,'rs');

 
    %% Polynomial fittings
    % [DR,DDR] = calcDR(Rnew(expt,:)*1E-6,t);

    Rtemp = RList ;
    d1R = (Rtemp(3:end)-Rtemp(1:end-2))/(2*(t(5)-t(4)));
    d2R = diff(diff(Rtemp))/((t(5)-t(4))^2);
    %figure; plot(t(2:end-1),d1R' ); title('d1R');
    %figure; plot(t(2:end-1),d2R); title('d2R');
    % % Find all the DDR==0 points
    %hold on; plot(t,0*t);

    deltat = 1e-7; t2 = 0:deltat:tspan;
    Rfit_pp = interp1(t(1:end)',Rtemp(1:end)','pchip','pp');
    Rfit_pp_d1 = fnder(Rfit_pp,1);
    % Rfit_pp_d2 = fnder(Rfit_pp,2);

    % Valuate Rfit
    Rfit = ppval(Rfit_pp,t2); d1Rfit = ppval(Rfit_pp_d1,t2); 

    Rdotfit_pp = interp1(t2',d1Rfit','pchip','pp');
    Rfit_pp_d2 = fnder(Rdotfit_pp,1);
    d2Rfit = ppval(Rfit_pp_d2,t2);

    close all;
    figure; plot(t,Rnew,'rs'); %hold on; plot(t,Rnew ,'rs'); 
    set(gca,'fontsize',18); axis([0,tspan,0,R0]); 
    xlabel('Time($s$)','interpreter','latex'); ylabel('Sampled $R$','interpreter','latex');
    
    figure; plot(t2,Rfit); hold on; plot(t,Rnew ,'rs'); set(gca,'fontsize',18);  axis([0,tspan,0,R0]);
    xlabel('Time($s$)','interpreter','latex'); ylabel('Fitted bubble radius $R$','interpreter','latex');
    a=gca; a.TickLabelInterpreter = 'latex';
    set(gca,'fontsize',18); a=gca; a.TickLabelInterpreter = 'latex';
    figure; plot(t2,d1Rfit); hold on; plot( t(2:end-1),d1R','rs' );  axis([0,tspan,min(d1Rfit),max(d1Rfit)]);
    xlabel('Time($s$)','interpreter','latex'); ylabel('Fitted ${\partial}{R}/{\partial}{t}$','interpreter','latex');
    set(gca,'fontsize',18); a=gca; a.TickLabelInterpreter = 'latex';
    figure; plot(t2,d2Rfit); hold on; plot( t(2:end-1),d2R','rs' ); axis([0,tspan,min(d2Rfit),max(d2Rfit)]);
    xlabel('Time($s$)','interpreter','latex'); ylabel('Fitted ${\partial^2}{R}/{\partial}{t^2}$','interpreter','latex');
    set(gca,'fontsize',18); a=gca; a.TickLabelInterpreter = 'latex';

    t2IndexStart = min(find(t2>t0)); t2IndexEnd = length(t2);
    t2Start = t2(t2IndexStart)-t2(1); t2End = t2(t2IndexEnd)-t2(1);
    t2 = t2(t2IndexStart:t2IndexEnd)-t2(1);
    Rfit = Rfit(t2IndexStart:t2IndexEnd);
    d1Rfit = d1Rfit(t2IndexStart:t2IndexEnd);
    d2Rfit = d2Rfit(t2IndexStart:t2IndexEnd);

    % cd('/Users/yangjin/Documents/MATLAB/Franck/IMR-master/data/numsim/4/');
    % save Rfit_numsim5_100.mat t2 Rfit d1Rfit d2Rfit Rfit_pp Rfit_pp_d1 Rfit_pp_d2 t2Start t2End
    % filename=['Rfit_stiff_exp',num2str(expts),'.mat'];
    % filename=['Rfit_water_exp',num2str(expts),'.mat'];
    % filename=['Rfit_collagen1_exp',num2str(expts),'.mat'];
    filename=['Rfit_numsim',num2str(simNo),'_exp',num2str(expt),'_',num2str(TimeRes),'.mat'];
    save(filename,'t2','Rfit', 'd1Rfit', 'd2Rfit', 'Rfit_pp', 'Rfit_pp_d1', 'Rfit_pp_d2', 't2Start', 't2End');

    cd('../../../../');
    % cd('/Users/yangjin/Documents/MATLAB/Franck/IMR_code');
    %cd('E:\Jin\Franck\IMR-master');
 

end

% ====== Try to regularized d1R and d2R ========
% ====== But almost no improvements ===========
% FDMat = zeros(length(d1Rfit),length(d1Rfit));
% for tempi = 2:length(d1Rfit)-1
% FDMat(tempi,tempi-1) = -0.5; FDMat(tempi,tempi+1) = -0.5;
% end
% FDMat(1,1)=-1; FDMat(1,2)=1; FDMat(length(d1Rfit),length(d1Rfit)-1)=-1; FDMat(length(d1Rfit),length(d1Rfit))=1;
% d1Rfit_rg = (eye(length(d1Rfit))+1e-2*FDMat'*FDMat)\d1Rfit';
% figure; plot(t2,d1Rfit_rg); hold on; plot( t(2:end-1),d1R','rs' ); title('d1R');

% a1 = 1/12; a2 = -2/3; a4 = 2/3; a5 = -1/12;
% d1Rfit = ( a1*Rfit(1:end-4)+a2*Rfit(2:end-3)+a4*Rfit(4:end-1)+a5*Rfit(5:end) )/(deltat);
% figure; plot(t2(3:end-2),d1Rfit); hold on;
% plot( t(2:end-1),d1R','rs' ); title('d1R');

% ====== Using FD to compute derivatives ========
% a1 = 1/90; a2 = -3/20; a3 = 3/2; a4 = -49/18; a5 = 3/2; a6 = -3/20; a7 = 1/90;
% d2Rfit = (a1*Rfit(1:end-6) + a2*Rfit(2:end-5) + a3*Rfit(3:end-4) + a4*Rfit(4:end-3) + ...
%     a5*Rfit(5:end-2) + a6*Rfit(6:end-1) + a7*Rfit(7:end)) / (deltat^2);
% figure; plot(t2(4:end-3),d2Rfit); hold on;
% plot( t(2:end-1),d2R','rs' ); title('d2R')









