function [select_pairs] = select_proj_pairs(combinations,kernel_1,kernel_2,distance,psf_ratio_new,c_L,c_R,Im_L,Im_R,NeurTT_L,NeurTT_R)

    normal_pairs=[];
    combinations(:,5)=0; 

    sig_ROIs=combinations(:,4)>0.75;  % select nueron
    select_combinations = combinations(sig_ROIs,:);    % Isolate significant ROIs
    [counts, values] = hist(select_combinations(:,1), unique(select_combinations(:,1)));
    
    for ii=1:size(values,1)
        Num_neuron=values(ii);
        pairs=select_combinations(select_combinations(:,1)==Num_neuron,:);
        if counts(ii)>1
            for hh=1:size(pairs,1)
                mask_L=mat2gray(c_L(:,:,pairs(hh,1)))>0;
                mask_R=mat2gray(c_R(:,:,pairs(hh,2)))>0;   
                dis_pairs=pairs(hh,3);
                Flu_L=squeeze(sum(sum(Im_L.*mask_L,1),2));
                [~,loc_max]=max(Flu_L);
                img_L=Im_L(:,:,loc_max).*(mask_L>0);
                img_R=mean(Im_R(:,:,loc_max-1:loc_max+1),3).*(mask_R>0);
                index_dis=round(min(find(abs(dis_pairs-distance) == min(abs(dis_pairs-distance))))); % depth get from distance
                y1 = conv2(img_L, kernel_2(:,:,index_dis),'same'); % corrcoef
                y2 = conv2(img_R, kernel_1(:,:,index_dis),'same'); % corrcoef
                ratio=sqrt(y1)./(sqrt(y1)+sqrt(y2));
                rr= reshape(ratio,1,[]);
                rr(isinf(rr))=[];rr(rr==0)=[];rr(rr==1)=[];rr(isnan(rr))=[];
                index_ratio=round(mean(find(abs(most_Element(rr)-psf_ratio_new) == min(abs(most_Element(rr)-psf_ratio_new)))));
                pairs(hh,5)=abs(index_dis-index_ratio);
            end

            [~,pc_choose]=min(pairs(:,5));
            pairs_t=pairs(pc_choose,:);
            pairs_t(:,5)=0;
            normal_pairs=[normal_pairs;pairs_t];

        else
            normal_pairs=[normal_pairs;pairs];
        end 
    end


Fres_combinations= combinations(~sig_ROIs,:); % Fres ROIs
retain_pairs = setdiff(combinations(:,1), select_combinations(:,1));
rows_to_clear = ~ismember(Fres_combinations(:, 1), retain_pairs);
Fres_combinations(rows_to_clear,:)=[];
retain_pairs = setdiff(combinations(:,2), select_combinations(:,2));
rows_to_clear = ~ismember(Fres_combinations(:, 2), retain_pairs);
Fres_combinations(rows_to_clear,:)=[];

special_pairs=[];
for hh=1:size(Fres_combinations,1)
    pairs_1=Fres_combinations(Fres_combinations(:,1)==Fres_combinations(hh,1),:);
    pairs_2=Fres_combinations(Fres_combinations(:,2)==Fres_combinations(hh,2),:);
    pairs=[pairs_1;pairs_2];
    pairs = unique(pairs, 'rows', 'stable'); 
    if size(pairs,1)==1 % one to one
        mask_L=mat2gray(c_L(:,:,pairs(1,1)))>0.1;
        mask_R=mat2gray(c_R(:,:,pairs(1,2)))>0.1;  

        tt1 =squeeze(sum(sum(Im_L.*mask_L,1),2));
        tt2 =squeeze(sum(sum(Im_R.*mask_R,1),2));
        p = polyfit(tt1, tt2, 1);  
        B_fit = polyval(p, tt1);
        B_mean = mean(tt2);
        SS_tot = sum((tt2 - B_mean).^2);
        SS_res = sum((tt2 - B_fit).^2);
        R_squared = 1 - (SS_res / SS_tot);
        if R_squared>0.7
            pairs(:,4)=R_squared;
            normal_pairs=[normal_pairs;pairs];
        end
    elseif pairs(1,1)==pairs(2,1)  % more proj in im_L to one proj in im_R
        [R_squared_max,pc_choose]=max(pairs(:,4));
        pairs_t=pairs(pc_choose,:);        
        tt1 = sqrt(NeurTT_L(pairs(1,1),:));
        tt2 = sqrt(NeurTT_R(pairs(1,2),:)+NeurTT_R(pairs(2,2),:)); % with overlapping proj
        p = polyfit(tt1, tt2, 1);
        B_fit = polyval(p, tt1);   
        B_mean = mean(tt2);
        SS_tot = sum((tt2 - B_mean).^2);
        SS_res = sum((tt2 - B_fit).^2);
        R_squared = 1 - (SS_res / SS_tot);
        if  R_squared >max(R_squared_max,0.7) 
            pairs(:,4)=R_squared;
            pairs(:,5)=1; 
            special_pairs=[special_pairs;pairs];
        else
            mask_L=mat2gray(c_L(:,:,pairs_t(1,1)))>0.1; % without overlapping proj
            mask_R=mat2gray(c_R(:,:,pairs_t(1,2)))>0.1;  

            tt1 =squeeze(sum(sum(Im_L.*mask_L,1),2));
            tt2 =squeeze(sum(sum(Im_R.*mask_R,1),2));
            p = polyfit(tt1, tt2, 1); 
            B_fit = polyval(p, tt1);   
            B_mean = mean(tt2);
            SS_tot = sum((tt2 - B_mean).^2);
            SS_res = sum((tt2 - B_fit).^2);
            R_squared = 1 - (SS_res / SS_tot);
            if R_squared>0.7
                pairs(:,4)=R_squared;
                normal_pairs=[normal_pairs;pairs_t];
            end
        end
                    
     elseif pairs(1,2)==pairs(2,2) % one proj in im_L to more proj in im_R
            [R_squared_max,pc_choose]=max(pairs(:,4));
            pairs_t=pairs(pc_choose,:);        
            tt1 = sqrt(NeurTT_L(pairs(1,1),:)+NeurTT_L(pairs(2,1),:));% with overlapping proj
            tt2 = sqrt(NeurTT_R(pairs(1,2),:));
            p = polyfit(tt1, tt2, 1);
            B_fit = polyval(p, tt1);
            B_mean = mean(tt2);
            SS_tot = sum((tt2 - B_mean).^2);
            SS_res = sum((tt2 - B_fit).^2);
            R_squared = 1 - (SS_res / SS_tot);
        if  R_squared >max(R_squared_max,0.7)
            pairs(:,4)=R_squared;
            pairs(:,5)=2; 
            special_pairs=[special_pairs;pairs];
        else
            mask_L=mat2gray(c_L(:,:,pairs_t(1,1)))>0.1; % without overlapping proj
            mask_R=mat2gray(c_R(:,:,pairs_t(1,2)))>0.1;  

            tt1 =squeeze(sum(sum(Im_L.*mask_L,1),2));
            tt2 =squeeze(sum(sum(Im_R.*mask_R,1),2));
            p = polyfit(tt1, tt2, 1); 
            B_fit = polyval(p, tt1);  
            B_mean = mean(tt2);
            SS_tot = sum((tt2 - B_mean).^2);
            SS_res = sum((tt2 - B_fit).^2);
            R_squared = 1 - (SS_res / SS_tot);
            if R_squared>0.7
                pairs_t(:,4)=R_squared;
                normal_pairs=[normal_pairs;pairs_t];
            end
        end
    end
end
   
special_pairs=unique(special_pairs, 'rows', 'stable'); 
normal_pairs = unique(normal_pairs, 'rows', 'stable'); 
select_pairs=[normal_pairs;special_pairs];

end

