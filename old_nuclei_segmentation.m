%code taken from https://blogs.mathworks.com/steve/2006/06/02/cell-segmentation
%NOT USED!

function [boundaries] = cell_segmentation(im_path)
	I = imread(im_path);
	I = rgb2gray(I);
	I = imcomplement(I);
	%figure; imshow(I)

	%increase contrast
	I = adapthisteq(I, 'ClipLimit', 1, 'Distribution', 'exponential');

	%create a disk shaped morphological object with a radius of 5 pixels
	d = strel('disk',5);
	%morpholically open the image with the disk, basically just makes the cells
	% more prominant and the background less so 
	I_blur = imopen(I, d);
	
	%convert the image to binary
	bin = im2bw(I_blur);

	%remove any noise pixels that may still exist (connectivity < 10)
	bin = bwareaopen(bin, 10);

	%compute the boundaries of the cells
	boundaries = bwperim(bin);

%original code by:
% _Steve Eddins_
% _Copyright 2006 The MathWorks, Inc._