function [lumen_boundaries] = lumen_segmentation2(slide_im, nuclei_boundaries)
	hsv = rgb2hsv(slide_im);
	%hsv thresholding to get lumen 
	lumen_hsv = hsv(:,:,2);
	lumen_hsv(lumen_hsv > 0.12) = 1;
	lumen_hsv(lumen_hsv <= 0.12) = 0;
	%convert to binary
	bin = imcomplement(im2bw(lumen_hsv, 0.1));
	%remove noise 
	bin = bwareaopen(bin, 50);
	%fill in nucleir
	bin = imfill(bin, 'holes');
	%erode then dilate the boundaries to smooth
	bin = bwmorph(bin, 'open', 5);
	%segment into superpixels 
	[sp,num_sp] = superpixels(uint8(bin*255), 2000);
	scored_superpixels = zeros(size(sp));
	%loop over all the superpixels
	for i = 1:num_sp
		%for each sp, get the fraction of white pixels
		pixels = bin(sp==i);
		ratio = size(pixels(pixels==1),1) / size(pixels,1);
		scored_superpixels(sp==i) = ratio;
	end
	%threshold the scores to create a binary image of the candidate lumen
	lumen = im2bw(scored_superpixels, 0.9);
	%perform spatial analysis to determine whether candidate lumen regions
	%reside inside epithelial regions (IE, surrounded by nuclei)
	lumen = bwareaopen(lumen, 800);
	nuclei_centers = regionprops(nuclei_boundaries, 'centroid');
	nuclei_centers = cat(1, nuclei_centers.Centroid);
	%note that the indices of these 2 CC's, even though one is for the
	%booundaries and one is for the whole lumen, are the same. This was
	%tested, although it doesn't seem clear why this holds always
	CC_boundary = bwconncomp(bwperim(lumen));
	CC = bwconncomp(lumen);
	%loop over lumen candidates boundaries and determine the number of nearby nuclei
	%t will store the coordinates of the lumen centers and their
	%corresponding number of nearby nuclei
	%col1 = x, col2 = y, col3 = number of nearby nuclei
	distance = 15;		%200
	nearby_nuclei_threshold = 200;	%70
	for i = 1:CC_boundary.NumObjects
		surrounding_nuclei = 0;
		try
			if size(CC.PixelIdxList{i},1) > 175000
				lumen(CC.PixelIdxList{i}) = 0;
				continue;
			end
		%corner case where CC and CC_boundary sizes don't line up. 
		catch e
			continue;
		end
		for j = 1:size(CC_boundary.PixelIdxList{i},1)
			[bound_r, bound_c] = ind2sub(size(lumen),CC_boundary.PixelIdxList{i}(j));
			for k = 1:size(nuclei_centers, 1)
				if pdist([bound_r bound_c ; nuclei_centers(k,:)]) < distance
					surrounding_nuclei = surrounding_nuclei + 1;
				end
			end
		end
		%{
		a = zeros(size(lumen));
		a(CC.PixelIdxList{i}) = 1;
		figure; imshow(a);
		disp(strcat("Surrounding nuclei: ", string(surrounding_nuclei)));
		%disp(string(size(CC.PixelIdxList{i}, 1)));
		%}
		if surrounding_nuclei < nearby_nuclei_threshold
			lumen(CC.PixelIdxList{i}) = 0;
		end
	end
	lumen_boundaries = bwperim(lumen);
	figure; imshow(imoverlay(slide_im, lumen_boundaries, 'r'));
end