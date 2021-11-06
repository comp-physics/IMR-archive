clear all; close all;

mag = 20; 

%files = dir('PA1000_*.tif');
files = dir('*.tif');
mkdir('CircOver');
um2px = 1.38; %calibration from image, in um2px


Rnew = zeros(length(files),101);
for i=1:length(files)
close all;
[R,tmax,centroid] = calcRofT(files(i).name); 

t = (1:tmax)*1/270000; %filming speed
%R = R*1.24; %1.24 um/px
figure;
plot(t,R*um2px*20/mag,' o'); %1.38 um/px for 20x, 2.76 um/px for 10x
xlabel('Time [\mus]');
ylabel('Bubble radius [\mum]');
%save([files(i).name(1:end-4) '.mat'],'R','t','centroid');
%save([files(i).name(1:end-4) 'nc_nan.mat'],'R','t','centroid');
Rnew(i,:) = R*um2px*20/mag;
movefile('*CircOver*',[pwd '\CircOver']);

end

%Rnew is in microns, t is in seconds
save('RofTdata.mat','Rnew','t');

%% plot the graphs
close all;
matfiles = dir('*.mat');
load('RofTdata.mat');
figure;

for i = 1:length(files)
   %load(matfiles(i).name)
   %figure;
   hold on;
   %plot(t,Rnew(i,:)/max(Rnew(i,5:end)),' o'); set(gca,'YLim',[0 1]);
   plot(t,Rnew(i,:),' o'); %set(gca,'YLim',[0 1]);
   xlabel('Time [\mus]');
   ylabel('Bubble radius [\mum]');   
end