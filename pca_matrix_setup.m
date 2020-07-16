%function to construct the pca_matrix. 
%%pca_matrix will be used as input to the pca function. it will be
%(3xnum_nuclei) x num_slides. Columns correspond to data samples, IE
%slides. row 1 will be the score of the top ribbon symmetric nuclei from
%that image, row 2 will be the top taper score, row 3 top separation, then
%it goes into 2nd and 3rd top scores etc. 



%DEPRECATED !!!!!!!!!!!!!!!!!!!!!!
function [pca_matrix] = pca_matrix_setup(ribbon, taper, separation, num_nuclei)
	pca_matrix = zeros(3*num_nuclei, 1);
	ind = 1;
	for i = 1:num_nuclei
		pca_matrix(ind) = ribbon(i);
		pca_matrix(ind+1) = taper(i);
		pca_matrix(ind+2) = separation(i);
		ind = ind+3;
	end
end