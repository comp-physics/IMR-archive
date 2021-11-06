function [LSQErr,LSQErrNo] = funIMRLSQErr(model,matPropVarList,tspan_NMS,pksMaxR,tspanLocsMaxR,NT,NTM,...
    Pext_type,Pext_Amp_Freq,disptime,Tgrad,Tmgrad,Cgrad,Dim,comp,eqR,t2,Rfit_exp)
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

% Initialize
LSQErr = 0; LSQErrNo = 0;
% tempi = 1; timetol = timetolList(1,1);

% [t3_num1, R3_num1, R3dot1, tdel, Tdel, Cdel] = IMRsolver2R_S ...
%     (model, G, G1, mu, alpha, lambda_nu, peakLoc_exp(tempi) , ...
%     peakLoc_exp(tempi+1) , peakVal_exp(tempi), NT/10, NTM, ...
%     Pext_type, Pext_Amp_Freq, disptime, Tgrad, Tmgrad, Cgrad, ...
%     Dim, comp, Pfit_pp, Pfit_pp_d1, eqR);
[t3_num, R3_num, R3dot, P3,S3,T3,C3,Tm3,tdel,Tdel,Cdel] = funIMRsolver( ...
    model, matprop, tspan_NMS, pksMaxR, NT, NTM, Pext_type, Pext_Amp_Freq, ...
    disptime, Tgrad, Tmgrad, Cgrad, Dim, comp, eqR/pksMaxR, 1e-7); %eqR should be ratio REq/R0.
 

Rfit_pp = interp1(t3_num,R3_num,'pchip','pp'); Rfit_num1 = ppval(Rfit_pp,t2);
LSQErr = LSQErr + sum((Rfit_exp - Rfit_num1).^2);
LSQErrNo = LSQErrNo + length(Rfit_exp);


