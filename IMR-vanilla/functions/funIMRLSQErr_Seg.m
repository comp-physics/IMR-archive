function [LSQErr,LSQErrNo] = funIMRLSQErr_Seg(model,matPropVarList,tspanLocsMinR,pksMaxR,tspanLocsMaxR,NT,NTM,...
    Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp,eqR,t2,Rfit_exp)
% ================================================
% Variable List 
% ================================================
% PickNo: useful for segments computation
% model:  material model
% matPropVarList: material properties variable list in the order of 
%                 {G_ooms,mu_ooms,alpha_ooms,lambda_nu_ooms,G1_ooms}
% ================================================

DefineTimetolListOrNot = 1; PlotDefineTimetolListOrNot = 0; PlotExact = 0;
Rfit_pp = interp1(t2,Rfit_exp,'pchip','pp');

try G = 10^matPropVarList(1); catch G = 0; end
try mu = 10^matPropVarList(2); catch mu = 0; end
try alpha = 10^matPropVarList(3); catch alpha = 0; end
try lambda_nu = 10^matPropVarList(4); catch lambda_nu = 0; end
try G1 = 10^matPropVarList(5); catch G1 = Inf; end

matprop = struct('G',G,'mu',mu,'alpha',alpha,'lambda_nu',lambda_nu,'G1',G1);

% Initialize
t3_num = []; R3_num = []; R3dot = []; LSQErr = 0; LSQErrNo = 0;
tempi = 1; timetol = 1e-6; % default value %timetol = timetolList(expt,1);

%% First collapse
[t3_num1, R3_num1, R3dot,P3,S3,T3,C3,Tm3,tdel,Tdel,Cdel] = IMRsolver( ...
    model, matprop, tspanLocsMaxR(tempi), pksMaxR(tempi), NT, NTM, Pext_type, Pext_Amp_Freq, ...
    disptime, Tgrad, Tmgrad, Cgrad, Dim, comp, eqR/pksMaxR(tempi)); %eqR should be ratio REq/R0.
% Find results of t3_num1 considering some timetol
[rowtemp,coltemp] = find(t3_num1<(tspanLocsMinR(tempi)-timetol)); 
t3_num1 = t3_num1(1:max(ceil(rowtemp))); R3_num1 = R3_num1(1:max(ceil(rowtemp)));
[t3_num1,ia] = unique(t3_num1); R3_num1 = R3_num1(ia);
t3_num = [t3_num;t3_num1]; R3_num = [R3_num;R3_num1]; 

if DefineTimetolListOrNot==1, 
    disp('hi');hold on; plot(t3_num1 ,R3_num1,'b-'); title(['G=',num2str(G),'; \mu=',num2str(mu),';']); end
if PlotDefineTimetolListOrNot==1,
    if PlotExact == 1
        hold on; plot(t3_num1,R3_num1,'b-'); title(['G=',num2str(G),'; \mu=',num2str(mu),';']);
    else
        hold on; plot(t3_num1,R3_num1,'k-.'); title(['G=',num2str(G),'; \mu=',num2str(mu),';']);
    end
end
 
timestep = t2(2)-t2(1);
tNew1 = tspanLocsMaxR(tempi):timestep:min(tspanLocsMinR(tempi)-timetol); Rfit_exp1 = ppval(Rfit_pp,tNew1);
Rfit_pp1 = interp1(t3_num1,R3_num1,'pchip','pp'); Rfit_num1 = ppval(Rfit_pp1,tNew1);
LSQErr = LSQErr + sum((Rfit_exp1-Rfit_num1).^2);
LSQErrNo = LSQErrNo + length(Rfit_exp1);

PickNo = length(tspanLocsMinR);
for tempi = 1:PickNo-1
   % Use default timetol %timetol = timetolList(expt,2*tempi+1); 
   
   %% From each Rmax peak, going to time + direction
   [t3_num1, R3_num1, R3dot,P3,S3,T3,C3,Tm3,tdel,Tdel,Cdel] = IMRsolver( ...
       model, matprop, tspanLocsMaxR(tempi+2)-tspanLocsMaxR(tempi+1), pksMaxR(tempi+1), ...
       NT, NTM, Pext_type, Pext_Amp_Freq, disptime, Tgrad, Tmgrad, Cgrad, Dim, comp, eqR/pksMaxR(tempi+1));
   % Find results of t3_num1 considering some timetol
   [rowtemp,coltemp] = find(t3_num1<(tspanLocsMinR(tempi+1)-timetol));
   t3_num1 = t3_num1(1:max(ceil(rowtemp))); R3_num1 = R3_num1(1:max(ceil(rowtemp)));
   [t3_num1,ia] = unique(t3_num1); R3_num1 = R3_num1(ia);
   t3_num = [t3_num;t3_num1]; R3_num = [R3_num;R3_num1];
   
   if DefineTimetolListOrNot==1, hold on; plot(t3_num1,R3_num1,'b-'); end
   if PlotDefineTimetolListOrNot==1,
       if PlotExact == 1, hold on; plot(t3_num1,R3_num1,'b-');
       else  hold on; plot(t3_num1,R3_num1,'k-.'); end
   end

   timestep = t2(2)-t2(1);
   tNew1 = tspanLocsMaxR(tempi+1):timestep:min(tspanLocsMinR(tempi+1)-timetol); Rfit_exp1 = ppval(Rfit_pp,tNew1);
   Rfit_pp1 = interp1(t3_num1,R3_num1,'pchip','pp'); Rfit_num1 = ppval(Rfit_pp1,tNew1);
   LSQErr = LSQErr + sum((Rfit_exp1-Rfit_num1).^2);
   LSQErrNo = LSQErrNo + length(Rfit_exp1);
   
   %% From each Rmax peak, going to time - direction
   [t3_num2, R3_num2, R3dot,P3,S3,T3,C3,Tm3,tdel,Tdel,Cdel] = IMRsolver( ...
        model, matprop, tspanLocsMaxR(tempi)-tspanLocsMaxR(tempi+1), pksMaxR(tempi+1), ...
        NT, NTM, Pext_type, Pext_Amp_Freq, disptime, Tgrad, Tmgrad, Cgrad, Dim, comp, eqR/pksMaxR(tempi+1));
    t3_num2 = t3_num2 + tspanLocsMaxR(tempi+1);
    % Find results of t3_num1 considering some timetol
    [rowtemp,coltemp] = find(t3_num2>(tspanLocsMinR(tempi)+timetol)); 
    t3_num2 = t3_num2(1:max(ceil(rowtemp)));R3_num2 = R3_num2(1:max(ceil(rowtemp)));
    [t3_num2,ia] = unique(t3_num2); R3_num2=R3_num2(ia);
    t3_num=[t3_num;t3_num2]; R3_num=[R3_num;R3_num2];  
    
    if DefineTimetolListOrNot==1, hold on; plot(t3_num2,R3_num2,'b-'); end
    if PlotDefineTimetolListOrNot==1,
        if PlotExact == 1, hold on; plot(t3_num2,R3_num2,'b-');
        else  hold on; plot(t3_num2,R3_num2,'k-.'); end
    end
    
    timestep = t2(2)-t2(1);
    tNew2 = tspanLocsMaxR(tempi+1):-timestep:(tspanLocsMinR(tempi)+timetol); Rfit_exp2 = ppval(Rfit_pp,tNew2);
    Rfit_pp2 = interp1(t3_num2,R3_num2,'pchip','pp'); Rfit_num2 = ppval(Rfit_pp2,tNew2);
    LSQErr = LSQErr + sum((Rfit_exp2-Rfit_num2).^2);
    LSQErrNo = LSQErrNo + length(Rfit_exp2);
   
   
end


 


