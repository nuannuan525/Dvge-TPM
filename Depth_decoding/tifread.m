function [ image1 ] =tifread( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


    aa=imfinfo(input_args);
    widthpix=aa(1).Width;
    heightpix=aa(1).Height;
    depthpix=size(aa,1);
    image1=zeros(heightpix,widthpix,depthpix);
    %% ---  read image
    for ii=1:depthpix
        image1(:,:,ii)=(imread(input_args,ii));
    end
end

