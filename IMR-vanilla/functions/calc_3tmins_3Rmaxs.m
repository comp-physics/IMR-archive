function [tmins,Rmaxs,tpeak0] = calc_3tmins_3Rmaxs(R,t,p0_locidx,p1_locidx,p2_locidx,p3_locidx)
%[first 3 collapse times, first 3 max radii, time of first peak] = calc_3_tmins_3Rmaxs(R, t, index locations of first 4 individual curves)
dRnew = diff(R);
check = [dRnew nan].*[nan dRnew]; 
extrema = find(check<0);
mins = extrema(2:2:end);
tmins = t(mins); tmins = tmins(1:3);

tf = linspace(t(1),t(end),10000);
A0 = polyfit(t(p0_locidx),R(p0_locidx),4);
cfitx0 = linspace(t(p0_locidx(1)-2),t(p0_locidx(end)+2),500);
cfit0 = A0(1)*cfitx0.^4+A0(2)*cfitx0.^3+A0(3)*cfitx0.^2+A0(4)*cfitx0+A0(5);
cfit0_full = A0(1)*tf.^4+A0(2)*tf.^3+A0(3)*tf.^2+A0(4)*tf+A0(5);
[~,peakidx] = max(cfit0_full);
tpeak0 = tf(peakidx);
%plot(cfitx0,cfit0,exptstar(idx,:),expRstar(idx,:));

A1 = polyfit(t(p1_locidx),R(p1_locidx),4);
cfitx1 = linspace(t(p1_locidx(1)-2),t(p1_locidx(end)+2),500);
cfit1 = A1(1)*cfitx1.^4+A1(2)*cfitx1.^3+A1(3)*cfitx1.^2+A1(4)*cfitx1+A1(5);
cfit1_full = A1(1)*tf.^4+A1(2)*tf.^3+A1(3)*tf.^2+A1(4)*tf+A1(5);
%plot(cfitx0,cfit0,cfitx1,cfit1,exptstar(idx,:),expRstar(idx,:));

A2 = polyfit(t(p2_locidx),R(p2_locidx),2);
cfitx2 = linspace(t(p2_locidx(1)-2),t(p2_locidx(end)+2),500);
cfit2 = A2(1)*cfitx2.^2+A2(2)*cfitx2+A2(3);
cfit2_full = A2(1)*tf.^2+A2(2)*tf+A2(3);
plot(cfitx0,cfit0,cfitx1,cfit1,cfitx2,cfit2,t,R);

A3 = polyfit(t(p3_locidx),R(p3_locidx),2);
cfitx3 = linspace(t(p3_locidx(1)-2),t(p3_locidx(end)+2),500);
cfit3 = A3(1)*cfitx3.^2+A3(2)*cfitx3+A3(3);
cfit3_full = A3(1)*tf.^2+A3(2)*tf+A3(3);
plot(cfitx0,cfit0,cfitx1,cfit1,cfitx2,cfit2,t,R);

[~, loc1] = min(abs(tf-tmins(1)));
[~, loc2] = min(abs(tf-tmins(2)));
[~, loc3] = min(abs(tf-tmins(3)));
crop = 50;
idx1 = (loc1-crop):(loc1+crop);
idx2 = (loc2-crop):(loc2+crop);
idx3 = (loc3-crop):(loc3+crop);

[~, tmin1_idx] = min(abs(cfit0_full(idx1)-cfit1_full(idx1)));
[~, tmin2_idx] = min(abs(cfit1_full(idx2)-cfit2_full(idx2)));
[~, tmin3_idx] = min(abs(cfit2_full(idx3)-cfit3_full(idx3)));

tmins = [tf(tmin1_idx+loc1-(crop+1)) tf(tmin2_idx+loc2-(crop+1)) tf(tmin3_idx+loc3-(crop+1))];
Rmaxs = [max(cfit0) max(cfit1) max(cfit2)];
end