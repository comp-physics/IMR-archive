clear; clc; %close all;

%fp = 'E:\Jon\Google Drive\RESEARCH DATA\FL Cav Data\160420\water\';
%fp = 'C:\Users\Jon\Brown Google Drive\RESEARCH DATA\FL Cav Data\160420\water\';
%load('Rtest.mat');
%load([fp 'RofTdata.mat']);

figure(38102); hold on;
for n=4
    expt = n;
    
    %%
    %G = num2str(log10(soln_mx{1,1,1}.G));
    %mu = num2str(log10(soln_mx{1,j,1}.mu))
    
    G = num2str(3.8);
    mu = num2str(-1.5);
    G1 = num2str(8.5);
    %load(['1021_11kPa_PA_NHSNS_G_' G '_mu_' mu '.mat']);
    load(['1026FB_11kPa_PA_NHSLS_G_' G '_G1_' G1 '.mat']);
    %%
    for k=length(G1_ooms):-1:1
        %figure(k);
        plot(soln_mx{1,1,k}.t2,soln_mx{1,1,k}.R2,'Color',[0 1-k/length(G1_ooms) 0]);
        hold on;
        %a{k} = convhull(soln_mx{1,1,k}.t2,soln_mx{1,1,k}.R2);
        %plot(soln_mx{1,1,k}.t2(a{k}),soln_mx{1,1,k}.R2(a{k}));
        
        Rp=diff(soln_mx{1,1,k}.R2);
        Rp1=Rp(2:end);
        Rp2=Rp(1:end-1);
        a{k,2}=[find(Rp1.*Rp2 < 0 & Rp1-Rp2 < 0)+1];
        a{k,1}=find(Rp1.*Rp2 < 0 & Rp1-Rp2 > 0)+1;
        
        t_minpeaks{k} = soln_mx{1,1,k}.t2(a{k,1}); R_minpeaks{k} = soln_mx{1,1,k}.R2(a{k,1});
        t_maxpeaks{k} = soln_mx{1,1,k}.t2(a{k,2}); R_maxpeaks{k} = soln_mx{1,1,k}.R2(a{k,2});
        
        %plot(t_minpeaks{k},R_minpeaks{k});
        %plot(t_maxpeaks{k},R_maxpeaks{k});
    end
    
    %%
    %figure(10);
    G1 = num2str(2.5);
    load(['1026UB_11kPa_PA_NHSLS_G_' G '_G1_' G1 '.mat']);
    k=1;
    plot(soln_mx{1,1,k}.t2,soln_mx{1,1,k}.R2,'Color',[1 0 0 ]); hold on;
    
    
    G1 = num2str(12.5);
    load(['1026LB_11kPa_PA_NHSLS_G_' G '_G1_' G1 '.mat']);
    k=4;
    plot(soln_mx{1,1,k}.t2,soln_mx{1,1,k}.R2,'Color',[0 0 1 ]); hold on;
    
    
    
    %for k=1:length(G1_ooms)
    %plot(soln_mx{1,j,1}.t2,soln_mx{1,j,1}.R2,'Color',[0 0 j/length(mu_ooms)]); hold on;
    %   plot(soln_mx{1,1,k}.t2,soln_mx{1,1,k}.R2,'Color',[1 0 0 ]); hold on;
    %end
    %  sz2 = length(mu_ooms);
    % j=0;
    % sz2 = 15;
    
    
    [Rmax, tmax] =  calcR0(Rnew(expt,:)*1E-6,t); %need to get the inital radius from a fit (R)
    eqR = median(Rnew(expt,62:end))*10^-6;
    R_eq = eqR/Rmax;
    
    %plot(t-tmax,Rnew(4,:)*1E-6, ' *','Color','blue'); hold on;
    
end

%t_maxpeaksfirst = t_maxpeaks{1};
%R_maxpeaksfirst = R_maxpeaks{1};

plot(soln_mx{1,1,1}.t2,soln_mx{1,1,1}.R2,'Color',[1 0 0]);

%t_maxpeakstemp = t_maxpeaks{46};
%R_maxpeakstemp = R_maxpeaks{46};
figure(4); hold on;

for c = 1:length(R_maxpeaks)
    plot(t_maxpeaks{c}(:),R_maxpeaks{c}(:),' o')
    
    R_maxpeaks_vec = R_maxpeaks{c};
    R_maxpeaks_pct{c} = (R_maxpeaks_vec - R_maxpeaks_vec(end))/(R_maxpeaks_vec(1) - R_maxpeaks_vec(end));
end


%%
%close all;


%clear all;
expt=10;

p_mapdec = [hex2dec('40') 0 hex2dec('4b');...
    hex2dec('76') hex2dec('2a') hex2dec('83');...
    hex2dec('99') hex2dec('70') hex2dec('ab');...
    hex2dec('c2') hex2dec('a5') hex2dec('cf');...
    hex2dec('e7') hex2dec('d4') hex2dec('e8');...
    hex2dec('e7') hex2dec('e7') hex2dec('e7');]/255;

g_mapdec = [0 hex2dec('44') hex2dec('1b');...
    hex2dec('1b') hex2dec('78') hex2dec('37');...
    hex2dec('5a') hex2dec('ae') hex2dec('61');...
    hex2dec('a6') hex2dec('db') hex2dec('a0');...
    hex2dec('d9') hex2dec('f0') hex2dec('d3');...
    hex2dec('e7') hex2dec('e7') hex2dec('e7');]/255;

fullmap = zeros(6,6,3);
fullmap(:,:,1) = [hex2dec('e7') hex2dec('e7') hex2dec('c2') hex2dec('99') hex2dec('76') hex2dec('40');...
    hex2dec('d9') hex2dec('d4') hex2dec('b9') hex2dec('91') hex2dec('7e') hex2dec('6d');...
    hex2dec('a6') hex2dec('bc') hex2dec('a2') hex2dec('7d') hex2dec('54') hex2dec('4c');...
    hex2dec('5a') hex2dec('9e') hex2dec('87') hex2dec('62') hex2dec('40') hex2dec('27');...
    hex2dec('1b') hex2dec('70') hex2dec('5c') hex2dec('40') hex2dec('22') hex2dec('0f');...
    hex2dec('00') hex2dec('45') hex2dec('35') hex2dec('20') hex2dec('0b') hex2dec('00');]/255;

fullmap(:,:,2) = [hex2dec('e7') hex2dec('d4') hex2dec('a5') hex2dec('70') hex2dec('2a') hex2dec('00');...
    hex2dec('f0') hex2dec('e3') hex2dec('c6') hex2dec('91') hex2dec('65') hex2dec('44');...
    hex2dec('db') hex2dec('d7') hex2dec('c2') hex2dec('93') hex2dec('5b') hex2dec('37');...
    hex2dec('aa') hex2dec('ce') hex2dec('b1') hex2dec('7f') hex2dec('4e') hex2dec('34');...
    hex2dec('78') hex2dec('ba') hex2dec('90') hex2dec('6f') hex2dec('48') hex2dec('2b');...
    hex2dec('44') hex2dec('9e') hex2dec('8b') hex2dec('5d') hex2dec('3e') hex2dec('1b');]/255;

fullmap(:,:,3) = [hex2dec('e7') hex2dec('e8') hex2dec('cf') hex2dec('ab') hex2dec('83') hex2dec('4b');...
    hex2dec('d3') hex2dec('ec') hex2dec('e1') hex2dec('cb') hex2dec('b0') hex2dec('96');...
    hex2dec('a0') hex2dec('e0') hex2dec('d5') hex2dec('bf') hex2dec('a4') hex2dec('8a');...
    hex2dec('61') hex2dec('c2') hex2dec('c3') hex2dec('a7') hex2dec('92') hex2dec('78');...
    hex2dec('37') hex2dec('a1') hex2dec('af') hex2dec('98') hex2dec('7d') hex2dec('64');...
    hex2dec('1b') hex2dec('88') hex2dec('93') hex2dec('7c') hex2dec('62') hex2dec('47');]/255;

x1 = 8; x2 = 12;
[xq,yq] = meshgrid(linspace(1,6,x1),linspace(1,6,x2));

figure(30);
fullmapmx = cat(3,interp2(fullmap(:,:,1),xq,yq),interp2(fullmap(:,:,2),xq,yq),interp2(fullmap(:,:,3),xq,yq));
image(fullmapmx);
axis image;
%
% hsb = zeros(6,6,3);
% hsb(1:6,1,1) = [291.2 297 281.43 281.69 291.24 291.2];
% hsb(1:6,1,2) = [0 8.62 20.29 34.5 67.94 100];
% hsb(1:6,1,3) = [90.59 90.98 81.18 67.06 51.37 29.41];
%
% hsb(1,1:6,1) = [291.2 107.59 113.9 125 138.06 143.82];
% hsb(1,1:6,2) = [0 12.08 27 48.28 77.5 100];
% hsb(1,1:6,3) = [90.59 94.12 85.88 68.24 47.06 26.67];
%
% for m = 2:6
%     for n = 2:6
%         for p=1:3
%         hsb(m,n,p) = (hsb(1,n,p)+hsb(m,1,p))/2;
%         end
%     end
% end
% amat = meshgrid(linspace(0,6,6),linspace(0,6,6))/6;
% for n=1:3
% b_mapdec(:,:,n) = bsxfun(@times,amat,p_mapdec(:,n))+ bsxfun(@times,1-amat,p_mapdec(:,n))
% end
%close all; 

figure(8887);

for i=1:length(G_ooms)%[1 7:16]
    [xqg,yqg] = meshgrid(1:3,linspace(1,length(g_mapdec),length(G_ooms)));
    gmapnew = interp2(g_mapdec,xqg,yqg);
    
    [xqp,yqp] = meshgrid(1:3,linspace(1,length(p_mapdec),length(mu_ooms)));
    pmapnew = interp2(p_mapdec,xqp,yqp);
    
    
    %figure(100+i);
    for j=1:length(mu_ooms)
        plot(soln_mx{i,j,1}.t2,soln_mx{i,j,1}.R2,'Color',pmapnew(length(mu_ooms)-(j-1),:));
        %plot(soln_mx{i,j,1}.t2,soln_mx{i,j,1}.R2,'Color',pmapnew(length(mu_ooms)-(j-1),:));
        hold on;
    end
    
    [Rmax, tmax] =  calcR0(Rnew(expt,:)*1E-6,t); %need to get the inital radius from a fit (R)
    eqR = median(Rnew(expt,62:end))*10^-6;
    R_eq = eqR/Rmax;
    
    plot(t-tmax,Rnew(expt,:)*1E-6, ' s','Color','blue'); hold on;
    set(gca,'XLim',[0 2.25E-4]);
    
end


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

for i=1:sG
    figure(i)
    for j=1:sM
        for k=1:sG1
            
            %figure(10000*i+100*j+k)
            %figure(10000*i+100+k)
            %figure(505)
            
            %[R2max idx] = max(soln_mx{i,j,k}.R2);
            
            hold on;
            
            t2 = interp(soln_mx{i,j,k}.t2,10);
            R2 = interp(soln_mx{i,j,k}.R2,10);
            %plot(t2,R2,'Color',[(i-1)/sG(1) (k-1)/sG1(1) (j-1)/sM(1)]);
            plot(t2,R2,'Color',[0 (i-1)/sG(1) 0]);
            
            [R0,t0] =  calcR0(Rnew(expt,:)*1E-6,t);
            t2exp=t-t0;
            R2exp=Rnew(expt,:)*1E-6;
            plot(t2exp,R2exp, ' *');
            
            distRmat = bsxfun(@minus,R2,R2exp);
            disttmat = bsxfun(@minus,t2,t2exp);
            distmat = sqrt(disttmat.^2+distRmat.^2);
            weight = ones(1,101);
            
            %for fitting just the collapse
            
            %[~,w_idx0] = max(R2exp);
            w_idxs = find(([0 diff(R2exp)].*[diff(R2exp) 0])<0)-1;
            weight([1:(w_idxs(1)) (w_idxs(2)-1):end]) = 0; %1.3kPa
            %weight([1:9 18:end]) = 0; %11kPa
            
            %weight([1:9 10:17 70:end]) = 0;
            
            %for fitting first 3 peaks
            %weight([1:9 31:69 70:end]) = 0; %11 kPa
            %weight([1:11 36:69 70:end]) = 0; %1.3 kPa
            
            
            %for fitting the whole curve
            %weight([1:9 70:end]) = 0; %11 kPa
            
            LSQ{i,j,k} = nansum(weight.*min(distmat,[],1));
            
            if LSQ{i,j,k}<LSQminidx(1)
                LSQminidx = [LSQ{i,j,k} i j k];
            end
        end
        %             for n=k:-1:1
        %                 [maxRNHZ(n) idx(n)] = max(soln_mx{i,j,n}.R2);
        %             end
        
    end
    
end

newLSQ = cell2mat(LSQ);
[~, ix] = min(mean(newLSQ,2));
G_ooms(ix)

figure(120);

i=LSQminidx(2); j=LSQminidx(3); k=LSQminidx(4);
log10([soln_mx{i,j,k}.G soln_mx{i,j,k}.mu soln_mx{i,j,k}.G1])
[i j k]
plot(soln_mx{i,j,k}.t2,soln_mx{i,j,k}.R2,'Color','blue');
hold on;

[Rmax, tmax] =  calcR0(Rnew(expt,:)*1E-6,t); %need to get the inital radius from a fit (R)
eqR = median(Rnew(expt,62:end))*10^-6;
R_eq = eqR/Rmax;

plot(t-tmax,Rnew(expt,:)*1E-6, ' s','Color','blue','Linewidth',1);
set(gca,'XLim',[0 2.25E-4]);
%%


figure(4909);
R = soln_mx{i,j,k}.R2;
U = soln_mx{i,j,k}.U;
%Req = R_eq*max(R);
Req = median(R((end-100):end));
S_rrE = 2*soln_mx{i,j,k}.G/3*((Req./R).^4-(R./Req).^2);
S_rrV = -4*U*soln_mx{i,j,k}.mu./R;
S_rr = S_rrE+S_rrV;
plot(soln_mx{i,j,k}.t2,log10(abs(S_rr)),'Color','blue');
title('Stress vs. Time')
hold on;

figure(4908);
plot(soln_mx{i,j,k}.t2,log10(abs(S_rrE)),'Color','green');
hold on;
plot(soln_mx{i,j,k}.t2,log10(abs(S_rrV)),'Color','magenta');

figure(4900);
plot(soln_mx{i,j,k}.t2,S_rr,'Color','blue');
title('Stress vs. Time')
hold on;

figure(4902);
S_rrpos = S_rr(S_rr>0); t2pos = soln_mx{i,j,k}.t2(S_rr>0);
S_rrneg = S_rr(S_rr<0); t2neg = soln_mx{i,j,k}.t2(S_rr<0);

plot(t2neg,log10(abs(S_rrneg)),' .','Color','red');
hold on;
plot(t2pos,log10(abs(S_rrpos)),' .','Color','blue');
title('Stress vs. Time')
b3 = bsxfun(@times,ones(length(S_rr),3),[44, 123, 182]/255); b3 = bsxfun(@times,b3,S_rr<0);
r3 = bsxfun(@times,ones(length(S_rr),3),[215, 25, 28]/255); r3 = bsxfun(@times,r3,S_rr>0);
rgb = b3+r3;
gsc = double(S_rr>0);

E_rr = -2*log(R/Req);

figure(5902);
plot(soln_mx{i,j,k}.t2,E_rr,'Color','blue');
title('Strain vs. Time');
hold on;

figure(5903);
cline(soln_mx{i,j,k}.t2,log10(S_rr),[],gsc); 
nmap = [44, 123, 182; 215, 25, 28]/255;
colormap(nmap);
hold on;

%hold on;
figure(49001);
%plot(soln_mx{i,j,k}.t2,soln_mx{i,j,k}.U./R);
E_rrdot = -2*soln_mx{i,j,k}.U./R;
cline(soln_mx{i,j,k}.t2,log10(abs(E_rrdot)),[],double(E_rrdot>0)); colormap(nmap);
hold on;
set(gca,'XLim',[0 1.25E-4]);

r0_Rfrag = 1.201*Rmax;
r_Rfrag = (Rmax^3-Req^3+r0_Rfrag.^3).^(1/3);
E_rrdot_Rfrag = -2*soln_mx{i,j,k}.U.*R.^2./(r_Rfrag.^3);
cline(soln_mx{i,j,k}.t2,log10(abs(E_rrdot_Rfrag)),[],double(E_rrdot_Rfrag>0)); colormap(nmap);



figure(5914);
r0 = linspace(Req,10*Rmax,10000);
r = (Rmax^3-Req^3+r0.^3).^(1/3);
plot(log10(r/Req),log10(-2*log(r./r0)),'Color','red');
hold on;
r_ = (min(R)^3-Req^3+r0.^3).^(1/3);
plot(log10(r_/Req),log10(-2*log(r_./r0)),'Color','blue');
title('Strain vs. Position');

figure(4911);
%Req = R_eq*max(R);
S_rrMax = 2*soln_mx{i,j,k}.G/3*((r0./r).^4-(r./r0).^2);
plot(log10(r/Req),log10(abs(S_rrMax)),'Color','red');
title('Stress vs. Position')
hold on;
S_rrMin = 2*soln_mx{i,j,k}.G/3*((r0./r_).^4-(r_./r0).^2);
plot(log10(r_/Req),log10(abs(S_rrMin)),'Color','blue');
title('Stress vs. Position')

figure(294);
plot(log10(r/Req),log10(r0/Req),'Color','red');
hold on;
plot(log10(r_/Req),log10(r0/Req),'Color','blue');


