
close all; 
RstList1 = 0:0.01:3;
RstList = 10.^(RstList1);
Eterm_Fung = 0*RstList; 
Eterm_NH = 0*RstList;
Eterm_YC = 0*RstList;
Eterm_Lin = 0*RstList;

alpha=0.01;
 for tempi = 1:length(RstList)
     Rst = RstList(tempi);
     Eterm_Fung(tempi) = (1-3*alpha)*(5/4-1/4/Rst^4-1/Rst) + 2*alpha*(-27/40-1/8/Rst^8-1/5/Rst^5-1/Rst^2+2*Rst);
     Eterm_NH(tempi) = (5/4-1/4/Rst^4-1/Rst);
     Eterm_YC(tempi) = 2/3*(1-1/Rst^3);
     Eterm_Lin(tempi) = 1-1/Rst^2;
 end
 
figure; loglog(RstList,Eterm_Fung);
hold on; loglog(RstList,Eterm_NH);
hold on; loglog(RstList,Eterm_YC);
hold on; loglog(RstList,Eterm_Lin);
set(gca,'fontsize',20);
xlabel('R/REq')
ylabel('S/(2G)')
title('Dimensionless elastic term','fontweight','normal')
lgd = legend('Fung model','Neo-Hookean','Yang & Church','Linear elastic');
% axis([10^-3,100,0,1e20])