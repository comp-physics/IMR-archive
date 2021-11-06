clear all; close all;
 
fp = 'V:\data\jbestrad\FL_Cav_Data\160420\11kPa_PA\';
fn1 = 'MakeFigure09a.mat';
fn2 = 'MakeFigure09b.mat';
expts = 6;

wht = [1 1 1];
grn = [0 68 27]/255;
prp = [64 0 75]/255;
org = [241 90 41]/255;

for expt=expts
    cd([fp num2str(expt)]);
    load(fn1,'mu_ooms');
    load(fn2,'G_ooms');
    
    sz = [length(G_ooms) length(mu_ooms)];
    sG = sz(1);
    sM = sz(2);
    
    dgrn = (wht-grn)/(sG-1);
    dprp = (wht-prp)/(sM-1);    
    
    fig = figure(1);
    fig.Units = 'inches';
    fig.Position = [7 2 7 8];
    
%% Plot for 9b is varying G
load(fn2);
    for i=1:(sG)
        for j=1 %1:sM
                if ~isempty(soln_mx{i,j})
                t2 = interp(soln_mx{i,j}.t2,10);
                R2 = interp(soln_mx{i,j}.R2,10);
                
                warning('off','all')
                [R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t);
                t2exp=t-t0;
                R2exp=Rnew(expt,:)*1E-6;
                
                distRmat = bsxfun(@minus,R2,R2exp);
                disttmat = bsxfun(@minus,t2,t2exp);
                distmat = sqrt(disttmat.^2+distRmat.^2);


                subplot(2,1,2);
                t2 = interp(soln_mx{i,j}.t2,10);
                R2 = interp(soln_mx{i,j}.R2,10);
                plot(t2,R2,'-','Color',grn+(sG-i)*dgrn,'LineWidth',1);
                hold on;
                if i == sG
                plot(t2exp,R2exp,' s','Color',org);
                end
                set(gca,'XLim',[0 2.25E-4],'YLim',[0 350E-6]);
                end
           
        end
    end
  
    load(fn1);
    for i=1
        for j=1:sM %1:sM

                if ~isempty(soln_mx{i,j})
                t2 = interp(soln_mx{i,j}.t2,10);
                R2 = interp(soln_mx{i,j}.R2,10);
                
                warning('off','all')
                [R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t);
                t2exp=t-t0;
                R2exp=Rnew(expt,:)*1E-6;
                
                distRmat = bsxfun(@minus,R2,R2exp);
                disttmat = bsxfun(@minus,t2,t2exp);
                distmat = sqrt(disttmat.^2+distRmat.^2);
                

                subplot(2,1,1);
                t2 = interp(soln_mx{i,j}.t2,10);
                R2 = interp(soln_mx{i,j}.R2,10);
                plot(t2,R2,'-','Color',prp+(sM-j)*dprp,'LineWidth',1);
                hold on;
                if j == sM
                plot(t2exp,R2exp,' s','Color',org);
                end
                set(gca,'XLim',[0 2.25E-4],'YLim',[0 350E-6]);
                end
           
        end
    end
    
    
    
end