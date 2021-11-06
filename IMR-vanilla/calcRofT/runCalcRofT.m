% clear all; close all;

mag = 20; 
%files = dir('PA1000_*.tif');
files = dir('Img_*.tif');
mkdir('CircOver');
um2px = 1.38; %calibration from image, in um2px

% Rnew = zeros(length(files),1);
% CircleFitPar = zeros(length(files),2); 

%% code by JYangm

im = cell(length(files),1);
for i = 1:length(files)
    im{i} = files(i).name;
end

Rnew = zeros(length(files),1);
for tempk = [1:length(files)]
    
    tempk
    
    SatisfyOrNot = 1;
    while SatisfyOrNot > 0
        
        close all;
        
        TiffRef = imgaussfilt(  double(rgb2gray(imread(im{1},1)))  );
        TiffRef = imgaussfilt( TiffRef ); TiffRef = imgaussfilt( TiffRef );
        TiffCurrent = imgaussfilt(  double(rgb2gray(imread(im{tempk},1)))  );
        TiffCurrent = imgaussfilt( TiffCurrent ); TiffCurrent = imgaussfilt( TiffCurrent );
        
        figure, imshow(uint8(TiffRef));
        figure, imshow(uint8(TiffCurrent));
        TiffDiff = TiffCurrent ;
        
        TiffDiffNorm = 255/(max(TiffDiff(:))-min(TiffDiff(:)))*(TiffDiff-min(TiffDiff(:)));
        figure, imshow(uint8(TiffDiffNorm));
        
        TiffDiffGauss = imgaussfilt(TiffDiffNorm,1);
        TiffDiffGauss = imgaussfilt(TiffDiffGauss,1);
        TiffDiffGauss = imgaussfilt(TiffDiffGauss,1);
        figure, imshow(uint8(TiffDiffGauss));
        
       % figure,  imshow(uint8(imread(im{tempk},1))); pause;
        
        %TiffDiff2 = imadjust(uint8(TiffDiffGauss)); figure, imshow(TiffDiff2);
        %[centers, radii, metric] = imfindcircles((TiffDiff2),[20,60]);
        
        % Threshold the Image
        % BW1 = edge(uint8(TiffDiffGauss),'Prewitt',0.025); figure, imshow(BW1);
        
        % [centers,radii] = imfindcircles(uint8(TiffDiffGauss),[320 340]/2 ,'ObjectPolarity','dark', ...
        %     'Sensitivity',0.8)
        
        prompt='BigOrSmall bubble: Small-0 or Big-1 ? '; BigOrSmall = input(prompt);
        
        close all; figure, imshow(uint8(TiffDiffGauss));
        
        ParameterTuneOrNot = 0; removeobjradius = 10;
        if BigOrSmall == 0, removeobjradius = 0; end
        gridx = [1,size(TiffDiffGauss,1)]; gridy = [1,size(TiffDiffGauss,2)];
        
%         if BigOrSmall == 0,
%             % Choose ROI domain
%             [gridy(1),gridx(1)] = ginput(1);
%             [gridy(2),gridx(2)] = ginput(1);
%             gridx = round(gridx); gridy = round(gridy);
%             % Adjust ROI domain
%             if gridx(1)<1, gridx(1)=1; end
%             if gridy(1)<1, gridy(1)=1; end
%             if gridx(2)>size(TiffDiffGauss,1), gridx(2)=size(TiffDiffGauss,1);end
%             if gridy(2)>size(TiffDiffGauss,2), gridy(2)=size(TiffDiffGauss,2);end
%         end
         
        while ParameterTuneOrNot < 1
            
            bw=imbinarize(uint8(TiffDiffGauss(gridx(1):gridx(2),gridy(1):gridy(2))),'adaptive','ForegroundPolarity','dark','Sensitivity',0.55);
            figure,imshow(bw)
            %[centers, radii, metric] = imfindcircles(uint8(bw),[20,60]);
            
            % remove all object containing fewer than 30 pixels
            bw = bwareaopen(bw,90); imshow(bw)
            
            if removeobjradius>1
                % fill a gap in the pen's cap
                se = strel('disk',removeobjradius);
                bw = imclose(bw,se); imshow(bw)
                
                % fill any holes, so that regionprops can be used to estimate
                % the area enclosed by each of the boundaries
                bw = imfill(bw,'holes'); imshow(bw)
            end
            
            prompt='ParameterTuneOrNot: Yes-0 or No-1 ? '; ParameterTuneOrNot = input(prompt);
            if ParameterTuneOrNot < 1
                prompt='What is removeobjradius? Input integer here: '; removeobjradius = floor(input(prompt));
            end
            
        end
        
        %% For large bubbles out of view
        if BigOrSmall == 1 % Big bubble
            
            [B,L] = bwboundaries(bw);
            
            figure,
            % Display the label matrix and draw each boundary
            imshow(label2rgb(L, @jet, [.5 .5 .5]))
            hold on
            % for k = 1:length(B)
            %   boundary = B{k};
            %   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
            % end
            XY=[];
            if length(B) ==3
                
                for k = [1,3]
                    boundary = B{k};
                    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
                    XY=[XY;[boundary(:,1)+gridx(1)-1,boundary(:,2)+gridy(1)-1]];
                end
            elseif length(B) ==2
                
                for k = [1,2]
                    boundary = B{k};
                    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
                    XY=[XY;[boundary(:,1)+gridx(1)-1,boundary(:,2)+gridy(1)-1]];
                end
                
            elseif length(B) ==1
                
                for k = [1]
                    boundary = B{k};
                    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
                    XY=[XY;[boundary(:,1)+gridx(1)-1,boundary(:,2)+gridy(1)-1]];
                end
                
            else
                
                for k = [1, length(B)]
                    boundary = B{k};
                    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
                    XY=[XY;[boundary(:,1)+gridx(1)-1,boundary(:,2)+gridy(1)-1]];
                end
            end
            
            try [row1,~] = find(XY(:,1)<gridx(1)+10); catch row1 = []; end
            try [row2,~] = find(XY(:,1)>gridx(2)-10); catch row2 = []; end
            try [row3,~] = find(XY(:,2)<gridy(1)+10); catch row3 = []; end
            try [row4,~] = find(XY(:,2)>gridy(2)-10); catch row4 = []; end
            row = unique([row1;row2;row3;row4]);
            XY(row,:) = [];
            
            
        %% For small bubbles totally inside view
        else
            bws = conv2(double(bw),ones(5,5),'same');
            figure, surf(bws)
            bws2 = bws<5;
            figure, imshow(bws2);
            
            se = strel('disk',10);
            bws2 = imclose(bws2,se); imshow(bws2)
            
            props = regionprops(bws2,'Area','PixelIdxList','MajorAxisLength','MinorAxisLength');
            [~,indexOfMax] = max([props.Area]);
            approximateRadius =  props(indexOfMax).MajorAxisLength/2;
            
            largestBlobIndexes  = props(indexOfMax).PixelIdxList;
            bws2 = false(size(bws));
            bws2(largestBlobIndexes) = 1;
            bws2 = imfill(bws2,'holes');
            figure;imshow(bws2);title('Leaving only largest blob and filling holes');
            edgebws2 = edge(bws2);
            figure;imshow(edgebws2);title('Edge detection');
            
            [row,col] = find(edgebws2==1);
            
            XY = [row+gridx(1)-1,col+gridy(1)-1];
            
        end
        
        try Par=CircleFitByTaubin(XY); catch Par = [0,0,0]; end
        
       
        
        Rnew(tempk) = Par(3);
        CircleFitPar(tempk,:) = [Par(1),Par(2)];
        
        %figure, imshow(uint8(TiffCurrent)); 
        figure, imshow(uint8(imread(im{tempk},1))); 
        hold on;
        viscircles([Par(2),Par(1)],Par(3),'Color','b');
        
        prompt='SatisfyOrNot for this frame: Yes-0 or No-1 ? '; SatisfyOrNot = input(prompt);
        
        pause(1);
        
    end
    
    
     
end

deltat = 1/270000; deltax = 2.76e-6;
%Rnew = Rnew*deltax;

%t2=1:1:length(files);
t=1:1:101;  % Rnew=Rnew';
figure,  plot(t, Rnew); hold on; plot(t, Rnew,'bd')
set(gca,'fontsize',18)
a=gca; a.TickLabelInterpreter = 'latex';
xlabel('Time(s)','interpreter','latex'); ylabel('Bubble radius $R$($\mu$m)','interpreter','latex');

[RmaxAll,RmaxTimeLoc] = max(Rnew);
try
    [fitobj] = fit(t(RmaxTimeLoc-2:RmaxTimeLoc+2),Rnew(RmaxTimeLoc-2:RmaxTimeLoc+2),'poly2');
    p = coeffvalues(fitobj); 
    RmaxTimeLoc = -p(2)/2/p(1);
    RmaxAll = (4*p(1)*p(3)-p(2)^2)/(4*p(1));
catch
    % Don't change RmaxAll and RmaxTimeLoc;
end

KeepIndex = find(t>RmaxTimeLoc);
t = [RmaxTimeLoc, t(KeepIndex)] ;
Rnew = [RmaxAll, Rnew(KeepIndex)];

t=t*deltat;
t=t-t(1);

fp = [folderNamePrefix,num2str(expt),'/'];

    %Load the file RofTdata.mat, which contains vars Rnew and t
    %Rnew has size [num_expts num_video_frames]
    cd([fp]);
    
    filename = [fileNamePrefix,'.mat'];
    save(filename,'t','Rnew');
    cd('../../../../');
    
    pause;

%% Make movie
% files = dir('Img_*.tif');
% im = cell(length(files),1);
% for i = 1:length(files)
%     im{i} = files(i).name;
% end
% 
% v = VideoWriter('RofT.mp4');
% open(v); 
% figure, 
% for tempk = 4:length(files)
%     clf
%     imshow(uint8(imread(im{tempk},1))); hold on;
%         viscircles([CircleFitPar(tempk,2),CircleFitPar(tempk,1)],Rnew(tempk),'Color','b');
%         axis([1,512,1,128]);
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%         tempk
%        
% end
% close(v);


%%
% for i=1  %:length(files)
% close all;
% [R,tmax,centroid] = calcRofT(files(i).name); 
% 
% t = (1:tmax)*1/270000; %filming speed
% %R = R*1.24; %1.24 um/px
% figure;
% plot(t,R*um2px*20/mag,' o'); %1.38 um/px for 20x, 2.76 um/px for 10x
% xlabel('Time [\mus]');
% ylabel('Bubble radius [\mum]');
% %save([files(i).name(1:end-4) '.mat'],'R','t','centroid');
% %save([files(i).name(1:end-4) 'nc_nan.mat'],'R','t','centroid');
% Rnew(i,:) = R*um2px*20/mag;
% movefile('*CircOver*',[pwd '\CircOver']);
% 
% end
% 
% %Rnew is in microns, t is in seconds
% save('RofTdata.mat','Rnew','t');
% 
% %% plot the graphs
% close all;
% matfiles = dir('*.mat');
% load('RofTdata.mat');
% figure;
% 
% for i = 1:length(files)
%    %load(matfiles(i).name)
%    %figure;
%    hold on;
%    %plot(t,Rnew(i,:)/max(Rnew(i,5:end)),' o'); set(gca,'YLim',[0 1]);
%    plot(t,Rnew(i,:),' o'); %set(gca,'YLim',[0 1]);
%    xlabel('Time [\mus]');
%    ylabel('Bubble radius [\mum]');   
% end