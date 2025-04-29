function [psf_ratio,kernel_1,kernel_2,distance] = psf_parameters(psf_1,psf_2)

n_pc=size(psf_2,3);
select_pc=[1:1:n_pc];
n_pc=size(select_pc,2);
psf_1=psf_1(:,:,select_pc);
psf_2=psf_2(:,:,select_pc);

ctr_x=round(size(psf_1,1)/2);
kernel_1= zeros(size(psf_1));
kernel_2= zeros(size(psf_1));
for k = 1:size(psf_1, 3)
    [max_value_1, max_index_1] = max(psf_1(:,:,k), [], 'all', 'linear');
    [x1, y1] = ind2sub(size(psf_1(:,:,k)), max_index_1);
    kernel_1(ctr_x, y1, k) = 1;
    [max_value_2, max_index_2] = max(psf_2(:,:,k), [], 'all', 'linear');
    [x2, y2] = ind2sub(size(psf_2(:,:,k)), max_index_2);   
    kernel_2(ctr_x, y2, k) = 1;
    distance(k)=sqrt((x1-x2)^2+(y1-y2)^2);
end
   
    psf_ratio=sqrt(squeeze(sum(sum(psf_2,1),2)))./(sqrt(squeeze(sum(sum(psf_1,1),2)))+sqrt(squeeze(sum(sum(psf_2,1),2)))); 

end

