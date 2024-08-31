function script_MVGALGCPfit(integerlabel)
% This is a script for fitting multivariate, geometric anisotropic log-Gaussian 
% Cox point processes, with Matern second-order dependence structures. The
% anisotropy, mean and covariance parameters are estimated and printed to a binary 
% file in the form of a cell. In the current setup, 100 datasets have processes
% fitted to them each time the script is executed; these datasets are labelled
% integerlabel+0:99.
% 
% INPUT
% integerlabel     	a nonnegative, integer-valued index, allowing multiple executions
%					of the script with distinct input/output
%
% OUTPUT FILES
% seed%i.bin			a binary file containing the seed used for the rng
% allestmu_%i.bin		a binary file containing 100 matrices of mu estimates
% allesttheta_%i.bin	a binary file containing 100 matrices of theta estimates
% allestzeta_%i.bin		a binary file containing 100 matrices of zeta estimates
% allestalpha_%i.bin	a binary file containing 100 matrices of alpha estimates
% allestnu_%i.bin		a binary file containing 100 matrices of nu estimates
% allestsigma_%i.bin	a binary file containing 100 matrices of sigma estimates
%
% last modified by jsmartin.stats@gmail.com in March 2018
%%
	for i=1:integerlabel
		rng('shuffle');
	end

rng('shuffle');
s = rng;
fid = fopen(sprintf('seed%i.bin',integerlabel),'a');
fwrite(fid,s.Seed,'double');
fclose(fid);

allesttheta = nan(2,2,100);
allestzeta = nan(2,2,100);

allestmu = nan(2,100);
allestalpha = nan(2,2,100);
allestnu = nan(2,2,100);
allestsigma = nan(2,2,100);

	for index=1:100
		% simulate multivariate data
%		script_MVanisotLGCPsim(index);
		datasetindex = (integerlabel-1)*100+index;
		try
			% read in multivariate anisotropic data
			fid = fopen(sprintf('data_MVGALGCP_%i.bin',datasetindex),'r');
			data = reshape(fread(fid,'double'),[],4);
			fclose(fid);
			
			fprintf('estimating anisotropy parameters...\n');
			[~,esttheta,estzeta] = evalc('PP_GAest(data,[0,1,0,1]);');

			allesttheta(:,:,index) = esttheta;
			allestzeta(:,:,index) = estzeta;

			fprintf('estimating Matern parameters...\n');
			% estimate using MPLE
			[~,estmu,~,estalpha,estnu,estsigma] = evalc('PP_Maternest_validMV(data,[0,1,0,1],esttheta,estzeta,[0,10.0],[0,5.0],[0,50.0]);');
			allestmu(:,index) = diag(estmu);
			allestalpha(:,:,index) = estalpha;
			allestnu(:,:,index) = estnu;
			allestsigma(:,:,index) = estsigma;
		catch ME
			warning(strcat(sprintf('dataset %i: ',datasetindex),ME.message));
			continue;
		end
	end

% print out estimates and true parameters
	fid = fopen(sprintf('allestmu_%i.bin',integerlabel),'a');
	fwrite(fid,allestmu,'double');
	fclose(fid);
	fid = fopen(sprintf('allestalpha_%i.bin',integerlabel),'a');
	fwrite(fid,allestalpha,'double');
	fclose(fid);
	fid = fopen(sprintf('allestnu_%i.bin',integerlabel),'a');
	fwrite(fid,allestnu,'double');
	fclose(fid);
	fid = fopen(sprintf('allestsigma_%i.bin',integerlabel),'a');
	fwrite(fid,allestsigma,'double');
	fclose(fid);
	fid = fopen(sprintf('allesttheta_%i.bin',integerlabel),'a');
	fwrite(fid,allesttheta,'double');
	fclose(fid);
	fid = fopen(sprintf('allestzeta_%i.bin',integerlabel),'a');
	fwrite(fid,allestzeta,'double');
	fclose(fid);
    
end