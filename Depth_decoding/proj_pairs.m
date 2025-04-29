function [select_pairs_all,select_pairs_unique,NeurTT,s1,psf_ratio_new,distance] = proj_pairs(NeurTT_L,NeurTT_R,c_L,c_R,Im_L,Im_R,psf_11,psf_22)
%This function is used to find the pairs of neuron from images by Beam 1 and Beam 2
% All of these neurons distributions ans Ga traces are getting from CNMF code
% the paris are found by the corrlection of Ga traces and the y positons distance 

% The inputs to this function are
%     NeurTT_L:  the Ga traces of all neurons by Beam 1  
%     NeurTT_R:  the Ga traces of all neurons by Beam 2
%     c_L     :  the distributions of all neurons by Beam 1  
%     c_R     :  the distributions of all neurons by Beam 2 
%     distance:  the distance of beam 1 and beam 2 along z
%     psf_ratio: the ratio [sqrt(Beam_2)/(sqrt(Beam_1)+sqrt(Beam_2))]
%     along z

% Used parameters

%     R_squared: Parameters used to evaluate correlation

%
% The Outputs of this function are
%      final_combinations : selcet neurons
%      NeurTT            : Ga traces;

%
% 20240124 - xiaona
%% psf
[psf_ratio,kernel_1,kernel_2,distance]=psf_parameters(psf_11,psf_22); % srqt(Im_R)./(sqrt(Im_L)+sqrt(Im_R))
n_pc=size(psf_ratio,1);
z=[1:1:size(distance,2)];
p = polyfit(z, distance, 1);
distance=polyval(p, z);

min_dis=min(distance);
max_dis=max(distance);

%% find proj pairs by separate distances

pos_A=[];pos_B=[];
for ii=1:size(c_L,3)
    mask_region_A=c_L(:,:,ii)>0;
    stats = regionprops(mask_region_A, 'Centroid', 'Area','MajorAxisLength', 'MinorAxisLength');
    pos=stats.Centroid;
    pos_A=[pos_A;pos];
end
for ii=1:size(c_R,3)
    mask_region_B=c_R(:,:,ii)>0;
    stats = regionprops(mask_region_B, 'Centroid', 'Area','MajorAxisLength', 'MinorAxisLength');
    pos=stats.Centroid;
    pos_B=[pos_B;pos];
end


min_dis=min(distance);
max_dis=max(distance);
combinations=[];
for hh = 1:size(pos_A, 1)
    for gg=1:size(pos_B, 1)
        x_diff = pos_B(gg, 1) - pos_A(hh, 1); 
        y_diff = abs(pos_B(gg, 2) - pos_A(hh, 2));

        if y_diff <= 4 && x_diff >= 0 && x_diff <= 25         
           dis_pais=sqrt((x_diff)^2+(y_diff)^2); 
           combinations=[combinations;hh,gg,dis_pais];
        end
    end
end

%% check all possibel proj pairs according to the correlation of fluorescence time traces

    for hh=1:size(combinations,1)
        tt1 = sqrt(NeurTT_L(combinations(hh,1),:));
        tt2 = sqrt(NeurTT_R(combinations(hh,2),:));
        p = polyfit(tt1, tt2, 1); % Perform a linear fit
        B_fit = polyval(p, tt1); % Calculate the fitted values
        B_mean = mean(tt2);
        SS_tot = sum((tt2 - B_mean).^2);
        SS_res = sum((tt2 - B_fit).^2);
        R_squared = 1 - (SS_res / SS_tot); %  R-squared value represent the accuracy of the mode
        combinations(hh,4)=R_squared;
    end

    sig_ROIs=combinations(:,4)>0.9;  % select nueron according to R-squared value
    select_combinations = combinations(sig_ROIs,:);                                       % Isolate significant ROIs
    [counts, values] = hist(select_combinations(:,1), unique(select_combinations(:,1)));  
    pairs_1=[];
    for ii=1:size(values,1)
        if counts(ii)==1
            pairs=select_combinations(select_combinations(:,1)==values(ii),:);
            pairs_1=[pairs_1;pairs];
        end 
    end
%% adjust the upper and lower limits of intensity ratio

    ratio_stack=[];
    for gg=1:size(pairs_1,1)
        dis_pairs=pairs_1(gg,3);
        mask_A=mat2gray(c_L(:,:,pairs_1(gg,1)))>0;
        mask_B=mat2gray(c_R(:,:,pairs_1(gg,2)))>0;    
        Flu_A=squeeze(sum(sum(Im_L.*mask_A,1),2));
        [~,loc_max]=max(Flu_A);
        img_L=Im_L(:,:,loc_max).*(mask_A>0);
        if loc_max<2
            loc_max=2;
        end
        img_R=mean(Im_R(:,:,loc_max-1:loc_max+1),3).*(mask_B>0);
        index_dis=round(min(find(abs(dis_pairs-distance) == min(abs(dis_pairs-distance))))); 

        y1 = conv2(img_L, kernel_2(:,:,index_dis),'same'); % corrcoef
        y2 = conv2(img_R, kernel_1(:,:,index_dis),'same'); % corrcoef
        ratio=sqrt(y2)./(sqrt(y1)+sqrt(y2));

        rr= reshape(ratio,1,[]);
        rr(isinf(rr))=[];rr(rr==0)=[];rr(rr==1)=[];rr(isnan(rr))=[];
        ratio_pairs_all(gg,1)=most_Element(rr);
        ratio_stack=cat(3,ratio_stack,ratio);
    end

    ratio_stack(ratio_stack==0)=nan;
    aa=reshape(ratio_stack,[],1);
    aa(isnan(aa))=[];
    aa1=sort(aa);
    [limit_floor_1,limit_ceil_1]=ratiorange(aa1);

    index_min=round(min(find(abs(min(pairs_1(:,3))-distance) == min(abs(min(pairs_1(:,3))-distance))))); % depth get from distance
    index_max=round(min(find(abs(max(pairs_1(:,3))-distance) == min(abs(max(pairs_1(:,3))-distance))))); % depth get from distance

    limit_new=interp1([index_min,index_max],[limit_floor_1,limit_ceil_1],(1:n_pc),'spline');
    psf_ratio_Max=max(limit_new);
    psf_ratio_Min=min(limit_new);
    psf_ratio_new = ((psf_ratio - min(psf_ratio)) * (psf_ratio_Max - psf_ratio_Min)) / (max(psf_ratio) - min(psf_ratio)) + psf_ratio_Min;

%%  select true proj pairs by Ga traces 
select_pairs_all = select_proj_pairs(combinations,kernel_1,kernel_2,distance,psf_ratio_new,c_L,c_R,Im_L,Im_R,NeurTT_L,NeurTT_R);
[~, unique_indices] = unique(select_pairs_all(:, 1), 'first');
select_pairs_unique = select_pairs_all(unique_indices, :);

[~, unique_indices] = unique(select_pairs_unique(:, 1), 'first');
select_pairs_unique = select_pairs_unique(unique_indices, :);

NeurTT=[];
for jj=1:size(select_pairs_unique,1)
    if select_pairs_unique(jj,5)==0
        time_trace=(NeurTT_L(select_pairs_unique(jj,1),:)+NeurTT_R(select_pairs_unique(jj,2),:))./2;
    elseif select_pairs_unique(jj,5)==1
        time_trace=NeurTT_L(select_pairs_unique(jj,1),:);
    else
        time_trace=NeurTT_R(select_pairs_unique(jj,2),:);
    end
    NeurTT=[NeurTT;time_trace]; 
end
 
%% get depth and Ga traces of neuron pairs
for gg=1:size(select_pairs_unique,1)
    dis_pairs=select_pairs_unique(gg,3);
    index_dis=round(min(find(abs(dis_pairs-distance) == min(abs(dis_pairs-distance))))); % depth get from distance
    
    mask_A=mat2gray(c_L(:,:,select_pairs_unique(gg,1)))>0.1;
    mask_B=mat2gray(c_R(:,:,select_pairs_unique(gg,2)))>0.1;
    Flu_A=squeeze(sum(sum(Im_L.*mask_A,1),2));
    [~,loc_max]=max(Flu_A);
    img_L=Im_L(:,:,loc_max).*(mask_A>0);
    img_R=mean(Im_R(:,:,loc_max-1:loc_max+1),3).*(mask_B>0);

    y1 = conv2(img_L, kernel_2(:,:,index_dis),'same'); % corrcoef
    y2 = conv2(img_R, kernel_1(:,:,index_dis),'same'); % corrcoef
    ratio=sqrt(y2)./(sqrt(y1)+sqrt(y2));
    
    rr= reshape(ratio,1,[]);
    rr(isinf(rr))=[];rr(rr==0)=[];rr(rr==1)=[];rr(isnan(rr))=[];
    index_ratio=round(mean(find(abs(most_Element(rr)-psf_ratio_new) == min(abs(most_Element(rr)-psf_ratio_new)))));% depth get from ratio
    
    if abs(index_dis-index_ratio)<3
        index=mean([index_dis,index_ratio]);
    else   
        img_L=mean(Im_L.*mask_A,3);
        img_R=mean(Im_R.*mask_B,3);
        area_neuron=min(sum(sum((img_L>0),1),2),sum(sum((img_R>0),1),2));
        area_overlap=0;
        max_iterations=3;
        iterations = 0;
        min_index=min(index_dis,index_ratio);
        max_index=max(index_dis,index_ratio);
        index=0;
        while area_overlap<=0.9*area_neuron && iterations < max_iterations
            for hh=min_index:max_index
                y1 = conv2(img_L, kernel_2(:,:,hh),'same'); % corrcoef
                y2 = conv2(img_R, kernel_1(:,:,hh),'same'); % corrcoef
                if sum(sum((y1>0).*(y2>0),1),2)> area_overlap
                    index=hh;
                    area_overlap=sum(sum((y1>0).*(y2>0),1),2);
                end
            end
            iterations=iterations+1;
            min_index=max(1,min_index-2);
            max_index=min(max_index+2,n_pc);
        end
    end
    select_pairs_unique(gg,6)=index;
end
s1=[];
end

