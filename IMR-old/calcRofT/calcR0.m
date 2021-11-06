function [R0, t0] = calcR0(R,t)
d2R = diff(diff(R));
maxidx = min(find([-inf -inf d2R(3:end)]>0));
minidx = max(find(isnan(R)))+1;
if isempty(minidx)
    minidx = min(find(R~=0));
end
idxs = minidx:maxidx;

A0 = polyfit(t(idxs),R(idxs),4);
cfitx0 = linspace(t(minidx),t(maxidx),5000);
cfit0 = A0(1)*cfitx0.^4+A0(2)*cfitx0.^3+A0(3)*cfitx0.^2+A0(4)*cfitx0+A0(5);

[R0 idx]=max(cfit0);
t0=cfitx0(idx); 
end