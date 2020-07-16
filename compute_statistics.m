%computes average symmetry values for the top num_nuclei most symmetric nuclei in
%the image specified by im_path. nuclei_scores is 1x3:.
%	nuclei_scores(1, 1) = avg ribbon symmetry
%	nuclei_scores(1, 2) = avg taper symmetry
%	nuclei_scores(1, 3) = avg separation symmetry
function [nuclei_scores, lumen_scores] = compute_statistics(im_path, num_nuclei)
	im_path = char(im_path);
	im_path = char(pwd + "/all_slides/" + im_path(1:strfind(im_path,'_')-1) + "/" + im_path);
	[~, nuclei_boundary_scores, lumen_boundary_scores] = compute_skeleton(im_path);
	try
		sorted_ribbon = find_high_symmetry_nuclei(nuclei_boundary_scores(:,:,1),[],[],[]);
		sorted_ribbon = sorted_ribbon(:, 1);
		sorted_ribbon = sorted_ribbon(1:num_nuclei);
		%sorted_average_ribbon = sorted_average_ribbon(1:end);
		sorted_taper = find_high_symmetry_nuclei(nuclei_boundary_scores(:,:,2),[],[],[]);
		sorted_taper = sorted_taper(:, 1);
		sorted_taper = sorted_taper(1:num_nuclei);
		%sorted_average_taper = sorted_average_taper(1:end);
		sorted_separation = find_high_symmetry_nuclei(nuclei_boundary_scores(:,:,3),[],[],[]);
		sorted_separation = sorted_separation(:, 1);
		sorted_separation = sorted_separation(1:num_nuclei);
		%sorted_average_separation = sorted_average_separation(1:end);
		nuclei_scores = [sum(sorted_ribbon)/size(sorted_ribbon, 1) sum(sorted_taper)/size(sorted_taper, 1) ...
			sum(sorted_separation)/size(sorted_separation, 1)];
		lumen_scores = [compute_lumen_scores(lumen_boundary_scores(:,:,1)), ...
			compute_lumen_scores(lumen_boundary_scores(:,:,2)), ...
			compute_lumen_scores(lumen_boundary_scores(:,:,3))];
	catch e
		disp("Error :(");
	end
end