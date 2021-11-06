function [LSQErr,LSQErrNo] = funIMR2LSQErr(model,matPropVarList,tspanLocsMaxP,pksMaxR,tspanLocsMaxR,NT,NTM,...
    Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp,Pfit_pp,Pfit_pp_d1,eqR,timestep,timetolList,Rfit_pp)
% ================================================
% Variable List 
% ================================================
% PickNo: useful for segments computation
% model:  material model
% matPropVarList: material properties variable list in the order of 
%                 {G_ooms,mu_ooms,alpha_ooms,lambda_nu_ooms,G1_ooms}
% ================================================

try G = 10^matPropVarList(1); catch G = 0; end
try mu = 10^matPropVarList(2); catch mu = 0; end
try alpha = 10^matPropVarList(3); catch alpha = 0; end
try lambda_nu = 10^matPropVarList(4); catch lambda_nu = 0; end
try G1 = 10^matPropVarList(5); catch G1 = Inf; end

matprop = struct('G',G,'mu',mu,'alpha',alpha,'lambda_nu',lambda_nu,'G1',G1);
%G = 10^NH_KV_varList(1); mu = 10^NH_KV_varList(2);
%G1 = 10^G1_ooms; alpha = 10^alpha_ooms; lambda_nu = 10^lambda_nu_ooms;

% Initialize
expt = 1; % in each local function 
t3_num = []; R3_num = []; R3dot = []; LSQErr = 0; LSQErrNo = 0;
tempi = 1; timetol = timetolList(expt,1);

%% First collapse
[t3_num1, R3_num1, R3dot1, tdel, Tdel, Cdel] = funIMRsolver2R_Seg( ...
    model, matprop, tspanLocsMaxR(tempi), tspanLocsMaxR(tempi+1), pksMaxR(tempi), ...
    NT, NTM, Pext_type, Pext_Amp_Freq, disptime, Tgrad, Tmgrad, Cgrad, ...
    Dim, comp, Pfit_pp, Pfit_pp_d1, eqR);
%[rowtemp,coltemp]=find(R3_num1<1e-5); t3_num1 = t3_num1(1:min(floor(rowtemp)));R3_num1 = R3_num1(1:min(floor(rowtemp)));
[rowtemp,coltemp]=find(t3_num1<(tspanLocsMaxP(tempi)-timetol)); 
t3_num1 = t3_num1(1:max(ceil(rowtemp))); R3_num1 = R3_num1(1:max(ceil(rowtemp)));
[t3_num1,ia] = unique(t3_num1); R3_num1=R3_num1(ia);
t3_num=[t3_num;t3_num1]; R3_num=[R3_num;R3_num1]; %R3dot = [R3dot;R3dot1];

% if DefineTimetolListOrNot==1, hold on; plot(t3_num1 ,R3_num1,'b-'); title(['G=',num2str(G),'; \mu=',num2str(mu),';']);end
% if PlotDefineTimetolListOrNot==1,
%     if PlotExact == 1
%         hold on; plot(t3_num1,R3_num1,'b-'); title(['G=',num2str(G),'; \mu=',num2str(mu),';']);
%     else
%         hold on; plot(t3_num1,R3_num1,'k-.'); title(['G=',num2str(G),'; \mu=',num2str(mu),';']);
%     end
% end

tNew1 = tspanLocsMaxR(tempi):timestep:min(tspanLocsMaxP(tempi)-timetol); Rfit_exp1 = ppval(Rfit_pp,tNew1);
Rfit_pp1 = interp1(t3_num1,R3_num1,'pchip','pp'); Rfit_num1 = ppval(Rfit_pp1,tNew1);
LSQErr = LSQErr + sum((Rfit_exp1-Rfit_num1).^2);
LSQErrNo = LSQErrNo + length(Rfit_exp1);

PickNo = length(tspanLocsMaxP);
for tempi = 1:PickNo-1
    
    timetol = timetolList(expt,2*tempi+1);
    
    %% From each Rmax peak, going to time + direction
    [t3_num1, R3_num1, R3dot1, tdel, Tdel, Cdel] = funIMRsolver2R_Seg( ...
        model, matprop, tspanLocsMaxR(tempi+1), tspanLocsMaxR(tempi+2), pksMaxR(tempi+1), NT, NTM, ...
        Pext_type, Pext_Amp_Freq, disptime, Tgrad, Tmgrad, Cgrad, Dim, comp, Pfit_pp, Pfit_pp_d1, eqR);
    % Find results of t3_num1 considering some timetol
    [rowtemp,coltemp] = find(t3_num1<(tspanLocsMaxP(tempi+1)-timetol)); 
    t3_num1 = t3_num1(1:max(ceil(rowtemp))); R3_num1 = R3_num1(1:max(ceil(rowtemp)));
    [t3_num1,ia] = unique(t3_num1); R3_num1 = R3_num1(ia);
    t3_num = [t3_num;t3_num1]; R3_num = [R3_num;R3_num1]; %R3dot = [R3dot;R3dot1];
%     if DefineTimetolListOrNot==1, hold on; plot(t3_num1,R3_num1,'b-'); end
%     if PlotDefineTimetolListOrNot==1,
%         if PlotExact == 1, hold on; plot(t3_num1,R3_num1,'b-');
%         else  hold on; plot(t3_num1,R3_num1,'k-.'); end
%     end
    
    tNew1 = tspanLocsMaxR(tempi+1):timestep:min(tspanLocsMaxP(tempi+1)-timetol); Rfit_exp1 = ppval(Rfit_pp,tNew1);
    Rfit_pp1 = interp1(t3_num1,R3_num1,'pchip','pp'); Rfit_num1 = ppval(Rfit_pp1,tNew1);
    LSQErr = LSQErr + sum((Rfit_exp1-Rfit_num1).^2);
    LSQErrNo = LSQErrNo + length(Rfit_exp1);
    
    %% From each Rmax peak, going to time - direction
    timetol=timetolList(expt,2*tempi);
    [t3_num2, R3_num2, R3dot2, tdel, Tdel, Cdel] = funIMRsolver2R_Seg( ...
        model, matprop, tspanLocsMaxR(tempi+1) , tspanLocsMaxR(tempi) , pksMaxR(tempi+1), NT, NTM, ...
        Pext_type, Pext_Amp_Freq, disptime, Tgrad, Tmgrad, Cgrad, Dim, comp, Pfit_pp, Pfit_pp_d1, eqR);
    % Find results of t3_num1 considering some timetol
    [rowtemp,coltemp]=find(t3_num2>(tspanLocsMaxP(tempi)+timetol)); 
    t3_num2 = t3_num2(1:max(ceil(rowtemp)));R3_num2 = R3_num2(1:max(ceil(rowtemp)));
    [t3_num2,ia] = unique(t3_num2); R3_num2=R3_num2(ia);
    t3_num=[t3_num;t3_num2]; R3_num=[R3_num;R3_num2]; %R3dot = [R3dot;R3dot2];
%     if DefineTimetolListOrNot==1, hold on; plot(t3_num2,R3_num2,'b-'); end
%     if PlotDefineTimetolListOrNot==1,
%         if PlotExact == 1, hold on; plot(t3_num2,R3_num2,'b-');
%         else  hold on; plot(t3_num2,R3_num2,'k-.'); end
%     end

    tNew2 = tspanLocsMaxR(tempi+1):-timestep:(tspanLocsMaxP(tempi)+timetol); Rfit_exp2 = ppval(Rfit_pp,tNew2);
    Rfit_pp2 = interp1(t3_num2,R3_num2,'pchip','pp'); Rfit_num2 = ppval(Rfit_pp2,tNew2);
    LSQErr = LSQErr + sum((Rfit_exp2-Rfit_num2).^2);
    LSQErrNo = LSQErrNo + length(Rfit_exp2);
   

end

 % end of function
