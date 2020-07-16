function [lumen_boundaries] = lumen_segmentation(slide_im, nuclei_boundaries)
	hsv = rgb2hsv(slide_im);
	%hsv thresholding to get lumen 
	lumen_hsv = hsv(:,:,2);
	lumen_hsv(lumen_hsv > 0.12) = 1;
	lumen_hsv(lumen_hsv <= 0.12) = 0;
	%convert to binary
	bin = imcomplement(im2bw(lumen_hsv, 0.1));
	%remove noise 
	bin = bwareaopen(bin, 50);
	%fill in nuclei
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
	lumen_centers = regionprops(lumen, 'centroid');
	lumen_centers = cat(1, lumen_centers.Centroid);
	nuclei_centers = regionprops(nuclei_boundaries, 'centroid');
	nuclei_centers = cat(1, nuclei_centers.Centroid);
	CC = bwconncomp(lumen);
	%loop over lumen candidates and determine the number of nearby nuclei
	%t will store the coordinates of the lumen centers and their
	%corresponding number of nearby nuclei
	%col1 = x, col2 = y, col3 = number of nearby nuclei
	t = [];
	distance = 300;		%200
	nearby_nuclei = 70;	%70
	for i = 1:size(lumen_centers, 1)
		surrounding_nuclei = 0;
		for j = 1:size(nuclei_centers, 1)
			if pdist([lumen_centers(i,:);nuclei_centers(j,:)]) < distance
				surrounding_nuclei = surrounding_nuclei + 1;
			end
		end
		t = [t;lumen_centers(i,1) lumen_centers(i,2) surrounding_nuclei];
	end
	%only consider lumen centers with >nearby_nuclei nearby nuclei
	t = t(t(:,3)>nearby_nuclei,1:2);
	t = uint16(t);
	%remove all candidate lumen who do not have enough nearby nuclei
	new_lumen = zeros(size(lumen));
	for i = 1:size(t, 1)
		%note that row/column are swapped here for some reason. when trying
		%to visualize, plotting the points in t as row/column seems
		%correct, as it plots the center points of the lumen. but when
		%indexing the lumen variable, it needs to be column/row. no clue
		%why
		%a better solution to this would be to link CC centers with the
		%index of that CC in the struct. 
		lin = sub2ind(size(lumen),t(i,2),t(i,1));
		for j = 1:CC.NumObjects
			if ismember(lin,CC.PixelIdxList{j})
				new_lumen(CC.PixelIdxList{j}) = 1;
			end
		end
	end
	lumen=new_lumen;
	lumen_boundaries = bwperim(lumen);
	%figure; imshow(imoverlay(slide_im, bwmorph(lumen_boundaries,'thicken',3), 'r'));
end