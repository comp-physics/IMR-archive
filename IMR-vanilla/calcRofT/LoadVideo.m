% Load images from video
% videoName = 'PAstiff_xy002.mp4';
v = VideoReader(videoName);

k = 1;
while hasFrame(v)
    mov(k).cdata = readFrame(v);
    
    
    %figure, imshow(mov(k).cdata);
    if k<10
        imgname = ['Img_000',num2str(k),'.tif'];
    elseif k<100
        imgname = ['Img_00',num2str(k),'.tif'];
    elseif k<1000
        imgname = ['Img_0',num2str(k),'.tif'];
    else
        imgname = ['Img_',num2str(k),'.tif'];
    end
    imwrite(mov(k).cdata,imgname);
    
    k = k + 1;
    
end

