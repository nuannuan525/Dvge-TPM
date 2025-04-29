classdef FourierOptics
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        ObjPl ; 
        Img_plane; 
        lambda;
        n; % water,1.33 the refractive index of immersion material
        Obj_f;  % the focal length of Obj
        NA;  % include the numerical aperture of objective : NA;
        
    end
    
    methods
        function obj = FourierOptics()  
            radius_exp_obj=8.28e3/2; 
            N_EXP_obj=648;
            obj.ObjPl.xx = linspace(-radius_exp_obj,radius_exp_obj,N_EXP_obj);  % unit micron pupil size
            obj.ObjPl.yy = obj.ObjPl.xx;
            [obj.ObjPl.X0,obj.ObjPl.Y0] = meshgrid(obj.ObjPl.xx,obj.ObjPl.yy);
            obj.ObjPl.Pupil = sqrt(obj.ObjPl.X0.^2 + obj.ObjPl.Y0.^2);
            obj.ObjPl.Pupil = obj.ObjPl.Pupil<(1.*max(obj.ObjPl.xx)); %logical

            [obj.ObjPl.Theta0,obj.ObjPl.R0] = cart2pol(obj.ObjPl.X0,obj.ObjPl.Y0);
            obj.ObjPl.R0 = mat2gray(obj.ObjPl.R0,[0 max(obj.ObjPl.xx)]); 
            obj.ObjPl.phase = zeros(size(obj.ObjPl.X0));

            obj.ObjPl.nmcoeff(1,:) = [0  1  1  2  2  2  3  3  3  3];
            obj.ObjPl.nmcoeff(2,:) = [0 -1  1 -2  0  2 -3 -1  1  3];
            obj.ObjPl.nmcoeff(3,:) = [0  0  0  2  0  0  0  0  0  0];
            for ii = 1:size(obj.ObjPl.nmcoeff,2)
                obj.ObjPl.phase(obj.ObjPl.Pupil) =obj.ObjPl.phase(obj.ObjPl.Pupil) + ...
                    obj.ObjPl.nmcoeff(3,ii)*zernfun(obj.ObjPl.nmcoeff(1,ii),obj.ObjPl.nmcoeff(2,ii),obj.ObjPl.R0(obj.ObjPl.Pupil),obj.ObjPl.Theta0(obj.ObjPl.Pupil)); % phase
            end
            obj.ObjPl.LightField = obj.ObjPl.Pupil.*exp(1i*obj.ObjPl.phase);

            %% basic parameters
            obj.Obj_f = 8e3;  % 
            obj.lambda = 0.92; % unit micron
            obj.NA=1.1;
            obj.n=1.33;
            %% sigle image plane
            obj.Img_plane.xx = 0;
            obj.Img_plane.yy = 0;
            [obj.Img_plane.X0,obj.Img_plane.Y0] = meshgrid(obj.Img_plane.xx,obj.Img_plane.yy);
            eps = 1e-6;
            obj.Img_plane.zlist = -100:1:0 + eps; % unit micros
            obj.Img_plane.DefocusF = reshape(1./(1./obj.Obj_f-1./(obj.Img_plane.zlist+obj.Obj_f))+obj.Obj_f,1,1,length(obj.Img_plane.zlist)); % 
            obj.Img_plane.DefocusFactor = exp(-1i*pi*obj.lambda/obj.n./obj.Img_plane.DefocusF.*(obj.ObjPl.X0.^2+obj.ObjPl.Y0.^2)); %
            obj.Img_plane.LightField_zstack = zeros(size(obj.Img_plane.X0,1),size(obj.Img_plane.X0,2),length(obj.Img_plane.zlist)); % 3D
        end
        
        function [obj,outInt]= getLightIntensity(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.Img_plane.LightField_zstack = zeros(size(obj.Img_plane.X0,1),size(obj.Img_plane.X0,2),length(obj.Img_plane.zlist));
            for ii = 1:length(obj.Img_plane.zlist)
                obj.Img_plane.LightField_zstack(:,:,ii) = 1/(obj.lambda/obj.n*obj.Obj_f)*FourierTransform2D(...
                obj.ObjPl.X0,obj.ObjPl.Y0,...
                obj.ObjPl.LightField.*obj.Img_plane.DefocusFactor(:,:,ii),...
                obj.Img_plane.X0/(obj.lambda/obj.n*obj.Obj_f),obj.Img_plane.Y0/(obj.lambda/obj.n*obj.Obj_f));      
            end
            outInt =abs(obj.Img_plane.LightField_zstack).^2; % intensity
            
        end
        function [obj,outField]= getLightField(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.Img_plane.LightField_zstack = zeros(size(obj.Img_plane.X0,1),size(obj.Img_plane.X0,2),length(obj.Img_plane.zlist));
            for ii = 1:length(obj.Img_plane.zlist)
                obj.Img_plane.LightField_zstack(:,:,ii) = 1/(obj.lambda/obj.n*obj.Obj_f)*FourierTransform2D(...
                obj.ObjPl.X0,obj.ObjPl.Y0,...
                obj.ObjPl.LightField.*obj.Img_plane.DefocusFactor(:,:,ii),...
                obj.Img_plane.X0/(obj.lambda/obj.n*obj.Obj_f),obj.Img_plane.Y0/(obj.lambda/obj.n*obj.Obj_f));      
            end
            outField =obj.Img_plane.LightField_zstack;  % complex amplitude
            
        end
        
        function  obj = setZernike_nmcoeff(obj,n_m_coeff)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if size(n_m_coeff,1)==3
            
                obj.ObjPl.nmcoeff = n_m_coeff;
                obj.ObjPl.phase = zeros(size(obj.ObjPl.X0));
                for ii = 1:size(obj.ObjPl.nmcoeff,2)
                    obj.ObjPl.phase(obj.ObjPl.Pupil) =obj.ObjPl.phase(obj.ObjPl.Pupil) + ...
                        obj.ObjPl.nmcoeff(3,ii)*zernfun(obj.ObjPl.nmcoeff(1,ii),obj.ObjPl.nmcoeff(2,ii),obj.ObjPl.R0(obj.ObjPl.Pupil),obj.ObjPl.Theta0(obj.ObjPl.Pupil));
                end
                obj.ObjPl.LightField = obj.ObjPl.Pupil.*exp(1i*obj.ObjPl.phase); 
                
%                 obj.lambda = 0.92;
            else
                error('wrong input: n_m_coeff must be 3xn 2D matrix');
            end
        end    % function   setZernike_nmcoeff


        
        function drwaing(obj,Intstack)
            if ~isreal(Intstack)
                Intstack = abs(Intstack).^2;
            end
            figure;
            imagemin = min(Intstack(:));
            imagemax = max(Intstack(:));
            for ii = 1:length(obj.Img_plane.zlist)
                subplot(1,length(obj.Img_plane.zlist),ii)
                imagesc(obj.Img_plane.xx,obj.Img_plane.yy, Intstack(:,:,ii),[imagemin,imagemax]);
                axis square;    
                title(['z = ',num2str(obj.Img_plane.zlist(ii)),' um']);
            end
        end
        

    end
end











