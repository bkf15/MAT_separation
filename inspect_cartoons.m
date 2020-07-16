%code to inspect akash cartoons
%	LABELS - 
%			adh08-1b_seg19 = FEA
%			adh19-1a_seg19 = Columnar
%			adh19-1a_seg27 = ADH
%			adh23-1a_seg76 = Normal
cartoons = {"adh08-1b_seg19_cartoon.jpg", "adh19-1a_seg19_cartoon.jpg", "adh19-1a_seg27_cartoon.jpg", "adh23-1a_seg76_cartoon.jpg"};
cartoon_labels = {"FEA", "Columnar", "ADH", "Normal"};
num_nuclei = 40;
t = table('Size', [4,3], 'VariableTypes', {'double','double','double'}, 'VariableNames', {'Ribbon', 'Taper', 'Separation'});

for i = 1:size(cartoons, 2)
	slide_im = imcomplement(im2bw(imread(char(pwd + "/akash_segmentation/" + cartoons{i})),0.5));
	slide_im = thin_skeleton(slide_im);
	CC = bwconncomp(slide_im);
	for j = 1:CC.NumObjects
		if size(CC.PixelIdxList{j}, 1) > 200
			slide_im(CC.PixelIdxList{j}) = 0;
		end
	end
	[skeleton] = generate_skeletons_from_img(slide_im, 'invert', true);
	skeleton = thin_skeleton(skeleton);
	CC = bwconncomp(skeleton);
	skeleton = remove_non_nucleus_skeleton_points(skeleton, slide_im, CC);
	%figure; imagesc(skeleton);
	dists = bwdist(slide_im);
	CC = bwconncomp(skeleton);
	sk_scores = compute_symmetry_for_sk_pts(dists, skeleton, 4, 0.1, CC);
	CC = bwconncomp(slide_im);
	b_scores = compute_boundary_symmetry(sk_scores(:,:,1), sk_scores(:,:,2), sk_scores(:,:,3), skeleton, slide_im, CC);
	avg_ribbon_score = find_high_symmetry_nuclei(b_scores(:,:,1), [],[],[]);
	avg_ribbon_score = avg_ribbon_score(:,1);
	avg_ribbon_score = avg_ribbon_score(1:num_nuclei);
	avg_ribbon_score = mean(avg_ribbon_score);
	
	avg_taper_score = find_high_symmetry_nuclei(b_scores(:,:,2), [],[],[]);
	avg_taper_score = avg_taper_score(:,1);
	avg_taper_score = avg_taper_score(1:num_nuclei);
	avg_taper_score = mean(avg_taper_score);
	
	avg_sep_score = find_high_symmetry_nuclei(b_scores(:,:,3), [],[],[]);
	avg_sep_score = avg_sep_score(:,1);
	avg_sep_score = avg_sep_score(1:num_nuclei);
	avg_sep_score = mean(avg_sep_score);
	
	t(i,:) = {avg_ribbon_score avg_taper_score avg_sep_score};
end


txt = table2cell(t);
txt = cellfun(@num2str, txt, 'UniformOutput', false); 
%X = [1 2 3 1 2 3 1 2 3 1 2 3];
%Y = [1 1 1 2 2 2 3 3 3 4 4 4];
X = [1 1 1 1 2 2 2 2 3 3 3 3];
Y = [1 2 3 4 1 2 3 4 1 2 3 4];
figure; imagesc(table2array(t));
text(X, Y, txt, 'HorizontalAlignment', 'Center'); colorbar;
title(sprintf('Cartoon Scores \n Ribbon           Taper            Separation'));
ylabel('Normal             ADH             Columnar           FEA');

