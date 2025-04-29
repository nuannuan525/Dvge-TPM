function [G_fxfy] = FourierTransform2D(x0,y0,g0,fx,fy)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~(size(fx,1)==1 || size(fx,2)==1)
    G_fxfy = zeros(size(fx));
end
if ~(size(x0,1)==1 || size(x0,2)==1)
    Dx0 =  x0(1,2)-x0(1,1);
    Dy0 =  y0(2,1)-y0(1,1);
end
 
for ii = 1:size(fx,1)
    for jj = 1:size(fx,2)
        temp_F = exp(-1i*2*pi*(fx(ii,jj)*x0+ fy(ii,jj)*y0));
%         figure;imagesc(angle(temp_F));
        temp_integ = g0.*temp_F;        
        G_fxfy(ii,jj) = sum(temp_integ(:));
    end
end
G_fxfy = G_fxfy*Dx0*Dy0;


end

