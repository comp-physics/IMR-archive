%normalize Rnew
clear all;

load('RofTdata.mat');

sizeRnew = size(Rnew);

s = 12;

normRnew = Rnew;
Rmax = max(Rnew,[],2);
Rmaxstats = [mean(Rmax),std(Rmax),min(Rmax),max(Rmax)];
normRnew = bsxfun(@rdivide,Rnew,Rmax);
shift_normTnew = bsxfun(@rdivide,repmat(t,sizeRnew(1),1),Rmax*(10^-6))*sqrt(101325/1060);
shift_normTnew = bsxfun(@minus,shift_normTnew,shift_normTnew(:,s));
[loc_1s_j, loc_1s_i] = find(normRnew'==1);
shiftj = s-loc_1s_j;

for i=1:sizeRnew(1)
    shift_normRnew(i,:) = circshift(normRnew(i,:),shiftj(i),2);
    if shiftj(i) > 0
        shift_normRnew(i,1:shiftj(i)) = nan(1,length(1:shiftj(i)));
    end
    if shiftj(i) < 0
        shift_normRnew(i,(sizeRnew(2)+shiftj(i)+1):sizeRnew(2)) = nan(1,length((sizeRnew(2)+shiftj(i)+1):sizeRnew(2)));
    end
end


shift_normRnew_stats = [nanmean(shift_normRnew,1); nanstd(shift_normRnew,1);];

figure(1) 
errorbar(t,shift_normRnew_stats(1,:),shift_normRnew_stats(2,:),' o');
hold on;
figure(2) 
plot(t,shift_normRnew_stats(1,:),' s','MarkerSize',12);
hold on;

figure(3)
for i=1:sizeRnew(1)
plot(shift_normTnew(i,:),shift_normRnew(i,:))
hold on;
end

save('RofT_shifted.mat','t','shift_normRnew','shift_normTnew','Rmax','Rmaxstats','shift_normRnew_stats');