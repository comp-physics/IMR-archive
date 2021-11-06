  
if strcmp(model,'fung')==1
    G_ooms = alpha_ooms;
end

LSQErrMatrix = zeros(length(G_ooms),length(mu_ooms)); LSQErrNoMatrix = LSQErrMatrix;
for tempi = 1:length(G_ooms)
    for tempj = 1:length(mu_ooms)
        LSQErrMatrix(tempi,tempj) = LSQErr{tempi,tempj};
        LSQErrNoMatrix(tempi,tempj) = LSQErrNo{tempi,tempj};
    end
end

        
% LSQErrMatrix = zeros(length(G_ooms),1);
% for tempi = 1:length(G_ooms)
%     LSQErrMatrix(tempi) = soln_mx{tempi,1,5,1}.LSQErr;
% end
% figure; plot(G_ooms,LSQErrMatrix); set(gca,'yscale','log');
        
% LSQErrMatrix = zeros(length(mu_ooms),1);
% for tempi = 1:length(mu_ooms)
%     LSQErrMatrix(tempi) = soln_mx{1,tempi,5,1}.LSQErr;
% end
% figure; plot(mu_ooms,LSQErrMatrix); set(gca,'yscale','log');
        


[mu_oomsGrid,G_oomsGrid] = meshgrid(mu_ooms,G_ooms);
figure; surf(G_oomsGrid,mu_oomsGrid,log10(LSQErrMatrix));
axis([ (G_ooms(1)), (G_ooms(end)),(mu_ooms(1)), (mu_ooms(end)),-9, -7]); caxis([-9, -7]);
 axis auto; caxis auto;
colormap gray; set(gca,'fontsize',20); xlabel('log(G)');ylabel('log(\mu)');zlabel('LSQ error');


[mu_Peak,G_Peak,LSQErrMin] = findpeak((-log10(LSQErrMatrix)),1);
G_oomsfit_pp = interp1(1:length(G_ooms),G_ooms,'pchip','pp');
G_Peak = ppval(G_oomsfit_pp,G_Peak); 10^G_Peak
mu_oomsfit_pp = interp1(1:length(mu_ooms),mu_ooms,'pchip','pp');
mu_Peak = ppval(mu_oomsfit_pp,mu_Peak); 10^mu_Peak
% LSQErrMin
        
hold on; p1=plot3(G_Peak,mu_Peak,-LSQErrMin+0.00001,'ro','MarkerFaceColor','r');
LSQErrMinExact  = interp2(mu_oomsGrid,G_oomsGrid,-log10(LSQErrMatrix),log10(0.01), log10(2.97e3) );
% hold on; p2=plot3(log10(2.97e3),log10(0.01),-LSQErrMinExact+0.02,'bo','MarkerFaceColor','b');
% lgd = legend([p1,p2],'Numerical','Exact');

axis tight; view([-1, -2.5, 3.5])
%axis([ (G_ooms(1)), (G_ooms(end)),(mu_ooms(1)), (mu_ooms(end)),-LSQErrMin-0.5, -LSQErrMin+2]); caxis([-LSQErrMin-0.5, -LSQErrMin+2]);
% axis([3.4,  3.82,(mu_ooms(1)), (mu_ooms(end)),-LSQErrMin-0.5, -LSQErrMin+2]); caxis([-LSQErrMin-0.5, -LSQErrMin+2]);
        
        
            % save(['LSQErrMatrix_stiff_expt',num2str(expt),'_PeakAll3_',num2str(PickNo),'.mat'],'LSQErrMatrix','LSQErrNoMatrix','G_ooms','mu_ooms');
            % save(['LSQErrMatrix_numsim',num2str(simNo),'_PeakAll.mat'],'LSQErrMatrix');
        
            % ****** Generate G_accept and mu_accept LSQMatrix mask ******
        %      % LSQErrMatrix
        %     [ G_accept, mu_accept ] = find( -log10(LSQErrMatrix)-LSQErrMin > - 0.5);
        %
        %     LSQErr_acceptMask = zeros(length(G_ooms),length(mu_ooms));
        %     for tempi = 1:length(mu_accept)
        %         LSQErr_acceptMask(G_accept(tempi),mu_accept(tempi)) = 1;
        %     end
        %     figure; surf(G_oomsGrid,mu_oomsGrid,LSQErr_acceptMask );
        %     set(gca,'fontsize',20); xlabel('log(G)');ylabel('log(\mu)');
        %     axis tight;    view(2); colormap gray; %  view([-1,-2.5, 3.5]);
        %
        %     LSQErr_acceptMaskGlobal = LSQErr_acceptMaskGlobal+LSQErr_acceptMask.*log10(LSQErrMatrix);
        %
        %     end
        %
        %     figure; surf(G_oomsGrid,mu_oomsGrid,LSQErr_acceptMaskGlobal);
        %     set(gca,'fontsize',20); xlabel('log(G)');ylabel('log(\mu)');
        %     axis tight;    view(2); colormap gray; %  view([-1,-2.5, 3.5]);
        
        
        %
%             save(['LSQErrMatrix_numsim',num2str(simNo),'_NH_expt',num2str(expt),'_PeakAll_',num2str(PickNo),'_100.mat'],'LSQErrMatrix','LSQErrNoMatrix','alpha_ooms','mu_ooms');
        
        % LSQErrMatrix Hessian matrix analysis
        
%             if mu_Peak<-2.9999
%                 mu_Peak=mu_Peak+0.1;
%             end
%             Glin = G_ooms(1):G_ooms_step:G_ooms(end); mulin=mu_ooms(1):mu_ooms_step:mu_ooms(end);
%             [muNodes,GNodes]=meshgrid(mulin,Glin);
%             [LSQErrMatrix2] = griddata(G_oomsGrid(:) ,mu_oomsGrid(:) ,log10(LSQErrMatrix(:)) ,GNodes,muNodes,'cubic');
%             [LSQErrMatrix2] = gridfit(G_oomsGrid(:) ,mu_oomsGrid(:) ,log10(LSQErrMatrix(:)),Glin,mulin,'smoothness',10); LSQErrMatrix2=LSQErrMatrix2';
%             figure; surf(GNodes ,muNodes ,LSQErrMatrix2);
%             hold on; p1=plot3(G_Peak,mu_Peak,-LSQErrMin+0.0001,'ro','MarkerFaceColor','r');
%              LSQErrMinExact  = interp2(mu_oomsGrid,G_oomsGrid,-log10(LSQErrMatrix),log10(0.01), log10(2.97e3) );
%             axis tight; view([-1,-2.5, 3.5])
%             axis([ (G_ooms(1)), (G_ooms(end)),(mu_ooms(1)), (mu_ooms(end)),-LSQErrMin-0.01, -LSQErrMin+1]); caxis([-LSQErrMin-0.2, -LSQErrMin+1]);
%         
%             temp = ((GNodes-G_Peak).^2+(muNodes-mu_Peak).^2) ;
%             [temprow,tempcol]=find(temp==min(min( temp )));
%             tempDist=6;
%             DfAxis = [temprow-tempDist,temprow+tempDist,tempcol-tempDist,tempcol+tempDist];
%             [DLSQErrDx,DLSQErrDy,DDLSQErrDxx,DDLSQErrDyy,DDLSQErrDxy] = funGradient(DfAxis,LSQErrMatrix2);
%         
%         	Hessian11=median(reshape(DDLSQErrDxx(tempDist+1,tempDist+1),1,1))/(0.02^2);
%             Hessian22=median(reshape(DDLSQErrDyy(tempDist+1,tempDist+1),1,1))/(0.02^2);
%             Hessian12=median(reshape(DDLSQErrDxy(tempDist+1,tempDist+1),1,1))/(0.02^2);
%         
%             invHessian = inv([Hessian11,Hessian12;Hessian12,Hessian22]);
%             stdHessianG = sqrt(invHessian(1,1));
%             stdHessianmu = sqrt(invHessian(2,2));
%         
%             disp([num2str(10^G_Peak),'; ',num2str(10^mu_Peak),'; ',num2str(G_Peak),'; ',num2str(stdHessianG),'; ',num2str(mu_Peak),'; ',...
%                 num2str(stdHessianmu),'; ']);
        