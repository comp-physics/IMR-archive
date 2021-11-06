function [R, tmax, centroid] = calcRofT(im)

%load tiff into a gray stack
%name = '512x128_270k_14';
%imext = '.tif';
%im = [name '.' imext];

T = Tiff(im,'r');
%I0 = read(T);
tiffinfo = imfinfo(im);
tiffstack = double(rgb2gray(imread(im,1)));
differential_images = tiffstack - double(rgb2gray(imread(im,1)));
for i = 2:1:size(tiffinfo,1)
    temptiff = double(rgb2gray(imread(im,i)));
    tiffstack = cat(3,tiffstack,temptiff);
    differential_images = cat(3,differential_images,temptiff-double(rgb2gray(imread(im,1))));
end
tmax = size(tiffstack,3);

%% initialize tiff stack for overlaid circle pictures 
f = figure; imshow(tiffstack(:,:,1),'Border','tight'); axis image; colormap gray; caxis([0 255]);
imwrite(frame2im(getframe(f)),[im(1:(end-4)) 'CircOver.tif']);
f = figure; imshow(squeeze(tiffstack(:,:,2)),'Border','tight'); axis image; colormap gray; caxis([0 255]);
imwrite(frame2im(getframe(f)),[im(1:(end-4)) 'CircOver.tif'],'WriteMode','append');
%% Fit circles on each timestep
R = zeros(tmax,1);
centroid = zeros(tmax,2);

for a=3:tmax %first image is ref, second is laser flash, third onward are bubbles
%threshold to get bubble dimensions
newim = abs(differential_images(:,:,a));
threshold = 0.40; %10/max(newim(:)); %set by looking at the histogram and finding the peak width
%threshold adjusted to 60% of the original difference, 20160111 
newim = ceil(newim/max(newim(:))-threshold);

CC = bwconncomp(newim);
CClabel = bwlabel(newim);
numPixels = cellfun(@numel,CC.PixelIdxList);
%if the second-largest convex hull is larger than half the largest convex
%hull, include it also
[biggestSize,idx] = max(numPixels);
secondBiggestSize = max(numPixels(numPixels<biggestSize));

BW2 = zeros(size(newim));
BW2(CC.PixelIdxList{idx}) = 1;
if secondBiggestSize*2>biggestSize
    BW2(CC.PixelIdxList{numPixels==secondBiggestSize}) = 1;    
end
% 
% figure;
% %subplot(4,1,1);
% imshow(uint8(tiffstack(:,:,a)),'Border','tight'); colormap gray; axis image;
% %subplot(4,1,2);
% imshow(uint8(abs(differential_images(:,:,a))),'Border','tight'); colormap gray; axis image;
% 
% imshow(newim,'Border','tight'); colormap gray; axis image;
% %subplot(4,1,3);
% imshow(BW2,'Border','tight'); colormap gray; axis image;
% %subplot(4,1,4);
% imshow(BW2.*tiffstack(:,:,a)); colormap gray; axis image;

stats = regionprops('table', BW2, 'Centroid', 'MajorAxisLength','MinorAxisLength','Orientation','ConvexHull');
circpts = stats.ConvexHull{1}; %takes convex hull data 
edgeidx = find((circpts(:,2)==0.5) | (circpts(:,2)==128.5) | (circpts(:,1)==0.5) | (circpts(:,1)==512.5)); %locating the edges on the images which shouldn't go into the circle fit
circpts(edgeidx,:) = []; %removing the non-circle image edges
fitcirc = CircleFitByTaubin(circpts); %use Taubin circle fit method to get the circle fit for the bubble
%figure; plot(circpts(:,1),circpts(:,2),' o'); axis image;

R(a) = fitcirc(3);
centroid(a,:) = fitcirc(1:2);

f = figure;
imshow(squeeze(tiffstack(:,:,a)),'Border','tight'); axis image; axis tight; colormap gray;  caxis([0 255]);
over_circ = viscircles(centroid(a,:),R(a),'LineStyle','--','EnhanceVisibility',false)
axis([1,512,1,128]);
imwrite(frame2im(getframe(f)),[im(1:(end-4)) 'CircOver.tif'],'WriteMode','append');
close all;
end

%movefile('*CircOver*','E:\Jon\Google Drive\Cavitation Data\20160107\PA001\1000\tiffs\circmovies')
 %centdiff = bsxfun(@minus,centroid,centroid(3,:));
 %norms = sqrt(sum(abs(centdiff).^2,2));
 %R(norms > 30) = NaN;

end
