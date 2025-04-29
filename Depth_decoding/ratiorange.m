function[limit_floor,limit_ceil]=ratiorange(aa)
bb=diff(aa);
n=ceil(size(bb,1)/1000); 
mean = ones(1,n)./n;  
y = conv(bb,mean,'same');  
%% find 
compy=find(y>std(y)*2);
compy_dix=find(diff(compy)>(n*10));
%%
limit_floor=aa(compy(compy_dix));
limit_ceil=aa(compy(compy_dix+1));
end
