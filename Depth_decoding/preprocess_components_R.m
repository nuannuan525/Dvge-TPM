function [cell_1,NeurTT_1] = preprocess_components_R(mask_cell,C_or,percent_1,percent_2,Df)

    cell=[];NeurTT=[];

    [nx,ny,nz]=size(mask_cell);
    
    for ii=1:nz
        img=mask_cell(:,:,ii);
        threshold=prctile(nonzeros(img), percent_1); 
        img(img<threshold)=0;
        labeled = bwlabel(img);
        for jj=1:max(max(labeled))
            Binary=zeros(nx,ny);
            Binary(labeled==jj)=1;
            conn_comp = bwconncomp(Binary); 
            props = regionprops(conn_comp, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
            ellipticity = props.MajorAxisLength/ props.MinorAxisLength;

              if sum(Binary(:))>50 && ellipticity<5
                img_new=img.*Binary;
                cell=cat(3,cell,img_new);
                NeurTT=[NeurTT;C_or(ii,:)./Df(ii)];
              end
        end
    end
    

    cell_1=[];    
    NeurTT_1=[];
    
    for hh=1:size(cell,3)
        img=cell(:,:,hh);
        value= prctile(nonzeros(img), percent_2); 
        s1=img>value;
        CC = bwconncomp(s1);
        stats = regionprops(CC, 'Centroid');
        
        if length(stats)==1
            cell_1=cat(3,cell_1,cell(:,:,hh));
            NeurTT_1=[NeurTT_1;NeurTT(hh,:)];
            
        elseif length(stats)==2
            centroid1 = stats(1).Centroid;
            centroid2 = stats(2).Centroid;
            slope = (centroid2(1) - centroid1(1)) / (centroid2(2) - centroid1(2)); 
            if slope == 0
                xBoundary = (centroid1(1) + centroid2(1)) / 2; 
                img(round(xBoundary))=0; 
                labeled = bwlabel(img);
                img_1=zeros(nx,ny);
                img_1(labeled==1)=img(labeled==1);
                img_2=zeros(nx,ny);
                img_2(labeled==2)=img(labeled==2);    
            else
                perpendicularSlope = -1 / slope;
                midPoint = (centroid1 + centroid2) / 2; 
                [X, Y] = meshgrid(1:nx, 1:ny); 
                YBoundary = perpendicularSlope * (X - midPoint(1)) + midPoint(2);
                region1Mask = Y < YBoundary; 
                region2Mask = Y > YBoundary;
                img_1=img.*region1Mask;
                img_2=img.*region2Mask;
            end
                cell_1=cat(3,cell_1,img_1);
                cell_1=cat(3,cell_1,img_2);
                NeurTT_1=[NeurTT_1;NeurTT(hh,:);NeurTT(hh,:)];
                
            end
       
    end
            
% figure(999);imagesc(sum(cell_1,3));axis image;
%%    %% check cell
select_neuron=[1,3,37];
    for gg=select_neuron
        img=cell_1(:,:,gg);
        value= prctile(nonzeros(img), 90); 
        s1=img>value;
        CC = bwconncomp(s1);
        stats = regionprops(CC, 'Centroid');
        
        if length(stats)==1
            cell_1=cat(3,cell_1,cell(:,:,hh));
            NeurTT_1=[NeurTT_1;NeurTT(hh,:)];
            
        elseif length(stats)==2
            centroid1 = stats(1).Centroid;
            centroid2 = stats(2).Centroid;
            slope = (centroid2(1) - centroid1(1)) / (centroid2(2) - centroid1(2)); 
            if slope == 0
                xBoundary = (centroid1(1) + centroid2(1)) / 2; 
                img(round(xBoundary))=0; 
                labeled = bwlabel(img);
                img_1=zeros(nx,ny);
                img_1(labeled==1)=img(labeled==1);
                img_2=zeros(nx,ny);
                img_2(labeled==2)=img(labeled==2);  
            elseif slope>10
                img(:,round((centroid1(1) + centroid2(1)) / 2))=0;
                labeled = bwlabel(img);
                img_1=zeros(nx,ny);
                img_1(labeled==1)=img(labeled==1);
                img_2=zeros(nx,ny);
                img_2(labeled==2)=img(labeled==2);   
            else
                perpendicularSlope = -1 / slope;
                midPoint = (centroid1 + centroid2) / 2; 
                [X, Y] = meshgrid(1:nx, 1:ny); 
                YBoundary = perpendicularSlope * (X - midPoint(1)) + midPoint(2);
                region1Mask = Y < YBoundary;
                region2Mask = Y > YBoundary;
                img_1=img.*region1Mask;
                img_2=img.*region2Mask;                              
            end
            cell_1=cat(3,cell_1,img_1);
            cell_1=cat(3,cell_1,img_2);
            NeurTT_1=[NeurTT_1;NeurTT(hh,:);NeurTT(hh,:)];
        end

    end
    cell_1(:,:,select_neuron)=[];
    NeurTT_1(select_neuron,:)=[];

end

