
clear all; close all;

%stiff
%fp = 'V:\data\jbestrad\FL_Cav_Data\160420\11kPa_PA\';
%fn = '170724_stiffsweep.mat';
% fn = '1101coarse_11kPa_PA_NHKV.mat';

%soft
fp = 'V:\data\jbestrad\FL_Cav_Data\160420\1pt3kPa_PA\';
fn = '170104mps_soft_sweep.mat';

%water
%fp = 'V:\data\jbestrad\FL_Cav_Data\160420\water\';
%fn = '170814_water.mat';%'170104mps_water_sweep.mat';
%fn = '170104mps_stiff_SNS_G_4_mu_-1.05.mat';
 
% fp = 'V:\data\jbestrad\FL_Cav_Data\160420\11kPa_PA\';
% fn = '1117Gfit_11kPa_PA_G.mat';

% %collagen
% fp = 'V:\data\jbestrad\FL_Cav_Data\160511\collagen\1\';
% %fn = '170206mps_collagen_sweep.mat';
% fn = '170411_collagenKVall.mat';

%collagen
%fp = 'V:\data\jbestrad\FL_Cav_Data\170403\Collagen\';
%fn = '170403_collagenSNS0404_G_5_mu_-1.3.mat';
%fn = '170408_collagenKVall2_.mat';


%close all;
expts=[1:5,7,9:20];
%expts = [12,14:19]; %chosen 7 water tests
%expts = [2,3,5:7]%[4,5,8,10,14,15,16,18,20,23,24];%15;%1:24%[1:7 9:20];
%expts = 6;
%for chsnmu = 1:6%1:

for expt=expts
    cd([fp num2str(expt)]);
    load(fn);
    
    sz = size(soln_mx);
    sG = sz(1);
    try
        sM = sz(2);
    catch
        sM = 1;
    end
    try
        sG1 = sz(3);
    catch
        sG1 = 1;
    end
    
    
    LSQ = cell(sG,sM,sG1);
    LSQminidx = [inf 0 0 0];
    
    for i=1:(sG)
        for j=1:sM %1:sM
            for k=1:(sG1)
                if ~isempty(soln_mx{i,j,k})
                t2 = interp(soln_mx{i,j,k}.t2,10);
                R2 = interp(soln_mx{i,j,k}.R2,10);
                
                warning('off','all')
                [R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t);
                t2exp=t-t0;
                R2exp=Rnew(expt,:)*1E-6;
                
                distRmat = bsxfun(@minus,R2,R2exp);
                disttmat = bsxfun(@minus,t2,t2exp);
                distmat = sqrt(disttmat.^2+distRmat.^2);
                weight = ones(1,101);
                
                
                
                w_idxs = find(([0 diff(R2exp)].*[diff(R2exp) 0])<0)-1;
                
                %%for fitting just the collapse
                %weight([1:(w_idxs(1)) (w_idxs(2)):end]) = 0;
                
                %for fitting first 3 peaks
                weight([1:(w_idxs(1)) (w_idxs(6)):end]) = 0;
                
                LSQ{i,j,k} = nansum(weight.*min(distmat,[],1));
                
                if LSQ{i,j,k}<LSQminidx(1)
                    LSQminidx = [LSQ{i,j,k} i j k];
                end
                
                figure(100*i);
                t2 = interp(soln_mx{i,j,k}.t2,10);
                R2 = interp(soln_mx{i,j,k}.R2,10);
                plot(t2,R2,'-','Color',[k/11 1-i/sG 1-j/sM],'LineWidth',1);
                hold on;
                plot(t2exp,R2exp,' s','Color','blue');
                set(gca,'XLim',[0 2.25E-4]);
                end
            end

        end
    end
    plot(t2exp,R2exp,' s','Color','blue');
     
    newLSQ = cell2mat(LSQ);
    [~, ix] = min(mean(newLSQ,3));
    %G_ooms(ix);
    
    i=LSQminidx(2); j=LSQminidx(3); k=LSQminidx(4);
    log10([soln_mx{i,j,k}.G soln_mx{i,j,k}.mu soln_mx{i,j,k}.G1]);
    [i j k];
    G_all(expt) = soln_mx{i,j,k}.G;
    mu_all(expt) = soln_mx{i,j,k}.mu;
    LSQ_all(expt) = LSQ{i,j,k};
    maxU_all(expt) = max(abs(soln_mx{i,j,k}.U));
    T_all(expt,:) = soln_mx{i,j,k}.T;
    
      figure(999);
      hold on;
      plot(t2exp/max(R2)*10.0751,R2exp/max(R2),' .','Color','blue');
      
      t2 = interp(soln_mx{i,j,k}.t2,10);
      R2 = interp(soln_mx{i,j,k}.R2,10);
      plot(t2,R2,'-','Color','red');
      plot(t2/max(R2)*10.0751,R2/max(R2),'-','Color','red');
t_all{expt} = t2; R_all{expt} = R2; allLSQ{expt} = cell2mat(LSQ);
%close all;
end

temp = zeros(1,41);
figure(3)
for p=expts
temp = temp+allLSQ{p};
plot(allLSQ{p}); hold on;
end
plot(mu_ooms,temp);

sumLSQ(chsnmu) = sum(LSQ_all);

%end
[minsumLSQ, idx] = min(sumLSQ);


mu_all = mu_all(expts);
LSQ_all = LSQ_all(expts);
G_all = G_all(expts);



%LSQ_all = LSQ_all.^2;
[mean(G_all) std(G_all)]

WmeanG = sum(G_all./(LSQ_all))/sum(1./LSQ_all);
WstdG = sqrt(sum((G_all-WmeanG).^2./LSQ_all)/(sum(1./(LSQ_all))-sum(1./(LSQ_all.^2))./sum(1./(LSQ_all))));

[WmeanG WstdG]

[mean(mu_all) std(mu_all)]

Wmeanmu = sum(mu_all./(LSQ_all))/sum(1./LSQ_all);
Wstdmu = sqrt(sum((mu_all-Wmeanmu).^2./LSQ_all)/(sum(1./(LSQ_all))-sum(1./(LSQ_all.^2))./sum(1./(LSQ_all))));

[Wmeanmu Wstdmu]

c = 1484;
maxU_all/c