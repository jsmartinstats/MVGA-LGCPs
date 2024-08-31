function script_MVanisotLGCPsim(integerlabel)
% This is a script for generating multivariate, geometric anisotropic log-Gaussian 
% Cox point patterns, with Matern second-order dependence structures. The
% anisotropy, mean and covariance parameters are specified below, and the simulated
% data is printed to a binary file in the form of a cell. In the current setup, 
% 10 independent MVGALGCPs are simulated each time the script is executed. 
% 
% INPUT
% integerlabel     	a nonnegative, integer-valued index, allowing multiple executions
%					of the script with distinct output files
%
% OUTPUT FILES
% seed_datagen_%i.bin			a binary file containing the seed used for the rng
% trueparams_MVGALGCP_%i.bin	a binary file containing the true mean, covariance and 
% 								and anisotropy parameters for the generating point 
%								point process model
%								*** only generated if integerlabel = 1 ***
% data_MVGALGCP_%i.bin			a binary file containing 10 independent MVGALGCPs.
% 								Each point pattern consists of 4 columns of data when
%								properly restructured: the first column contains the 
%								component labels; the second is a column of NaNs; the
%								third and fourth are the x- and y-coordinates, respectively.
% 
% last modified by jsmartin.stats@gmail.com in Feb 2018
%%
	rng('shuffle');
	s = rng;
	fid = fopen(sprintf('seed_datagen_%i.bin',integerlabel),'a');
	fwrite(fid,s.Seed,'double');
	fclose(fid);
%%
 	truemu = [5.75,5.625];
	truealpha = [0.09,0.1;0.1,0.12];
	truenu	  = [0.5,0.5;0.5,0.5];
 	truezeta  = [0.2,0.35;0.35,0.2];
	truetheta = [pi*0.2,pi*0.3;pi*0.3,pi*0.4];
	truesig11UB = (pi/truezeta(1,1))*(4*truenu(1,1)./(truealpha(1,1)^2))^(-truenu(1,2)) * gamma(truenu(1,1));
	truesig22UB = (pi/truezeta(2,2))*(4*truenu(2,2)./(truealpha(2,2)^2))^(-truenu(1,2)) * gamma(truenu(2,2));
	truesig12UB = (pi/truezeta(1,2))*(4*truenu(1,2)./(truealpha(1,2)^2))^(-truenu(1,2)) * gamma(0.5*(truenu(1,1)+truenu(2,2))+1)*gamma(truenu(1,2))/gamma(truenu(1,2)+1);
	
    truesigma = [truesig11UB,truesig12UB;truesig12UB,truesig22UB] .* [1,0.775;0.775,1.125*truesig11UB./truesig22UB];
                                                               % matrix multiplier must be nonnegdef
	truesigma = 2*truesigma./(truesigma(1,1));
    expabund = exp(truemu+diag(truesigma)'/2);
    
%  To check the validity of the multivariate geometric anisotropic Matern covariance structure,
%  using the conditions in the accompanying manuscript, uncomment the following line:
%     MVGALGCP_condchecker({truealpha,truenu,truesigma,truetheta,truezeta}); 
%%    
	if integerlabel==1
		fid = fopen(sprintf('trueparams_MVGALGCP_%i.bin',integerlabel),'a');
		fwrite(fid,[truemu',truealpha,truenu,truesigma,truetheta,truezeta],'double');
		fclose(fid);
	end
	data = MVGALGCPsimulate(10,[0,1,0,1],{expabund,truealpha,truenu,truesigma,truetheta,truezeta});
	
	for i=1:numel(data)
		fid = fopen(sprintf('data_MVGALGCP_%i.bin',(integerlabel-1)*10+i),'a');
		fwrite(fid,data{i},'double');
		fclose(fid);
	end

end


