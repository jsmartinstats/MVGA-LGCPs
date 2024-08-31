function varargout = MinContrastSPP_coarse(data,W,ra,rb,rstep,expon,alpha,nu,sigma,varargin)
% This function performs parameter estimation for a Matern LGCP, via a coarse 
% minimum contrast procedure (the grid of candidate parameter values is coarse).

% INPUT: 
% data      spatial point pattern data
% W 		the observation window for the data, in the form
% 				[lower x bound, upper x bound, lower y bound, upper y bound]
% ra        the smallest distance at which the empirical and theoretical
%           moments are contrasted
% rb        the largest distance at which the empirical and theoretical
%           moments are contrasted
% rstep     the size of the intervals in the sequence of distances at which
%           the empirical and theoretical moments are contrasted
% expon     exponent for 
% alpha 	gives either:
% 				the range of estimation [alpha_1,alpha_2] in which to seek estimates; or
%				the fixed value of alpha if it is not to be estimated.
% nu		similar to alpha.
% sigma		similar to alpha.
% varargin  potentially contains an angle (varargin{1}) and ratio (varargin{2})
%           of anisotropy, which we may wish to use to transform the data
%           before finding the pcf/covariance
 
% OUTPUT:
% varargout contains the parameters that have been estimated

% last modified by jsmartin.stats@gmail.com in Dec 2017

%%
	rvec = ra:rstep:rb;

	if numel(alpha)>1
% 		alpha = linspace(alpha(1),alpha(2),100);
		alpha = alpha(1):0.2:alpha(2);
	end
	if numel(nu)>1
% 		nu = linspace(nu(1),nu(2),10);
% 		nu = nu(1):0.5:nu(2);
        nu = [0.05,0.5,5]; % use discrete candidate smoothnesses; these differ from Moller & Toftaker (2014) as we use a different parameterisation
        nu = nu(nu>=nu(1));   % ...as long as they are not below our input lower bound
        nu = nu(nu<=nu(2));   % ...and as long as they are not above our input upper bound
	end
	if numel(sigma)>1
% 		sigma = linspace(sigma(1),sigma(2),10);
		sigma = sigma(1):1.0:sigma(2);
	end
		
    [rgrid,alphagrid,nugrid,sigmagrid] = ndgrid(rvec,alpha,nu,sigma);
    numgridels = numel(rgrid);
	griddims = size(rgrid);
	
    if nargin>8
        empir_g = Gsthat_isot_genW(data,W,rvec,'box',0.01,varargin{1},varargin{2});
        theor_cov = covfun_geoanisot(rgrid(:),zeros(numgridels,1),'Matern_HW94',sigmagrid(:),alphagrid(:),nugrid(:),varargin{1},1,varargin{2});
        % can't use covfun_geoanisot_quick(. . .) here
    else
        empir_g = Gsthat_isot_genW(data,W,rvec,'box',0.01);
        theor_cov = covfun_geoanisot(rgrid(:),zeros(numgridels,1),'Matern_HW94',sigmagrid(:),alphagrid(:),nugrid(:),0,1,1);
    	% can't use covfun_geoanisot_quick(. . .) here
    end
    if numel(empir_g)>numel(rvec)
        empir_g = empir_g(2,1,:);
    end
    empir_g = reshape(repmat(empir_g(:),1,numgridels./numel(rvec)),griddims);
    empir_cov = log(empir_g);
	theor_cov = reshape(theor_cov,size(rgrid));
    
    integrand = (empir_cov.^expon - theor_cov.^expon).^2;
    minimand = sum(rstep.*(integrand(1:end-1,:,:,:)+integrand(2:end,:,:,:))./2);
    linind = find(minimand==min(minimand(:)));
    
    alphagrid = squeeze(alphagrid(1,:,:,:));
    nugrid = squeeze(nugrid(1,:,:,:));
    sigmagrid = squeeze(sigmagrid(1,:,:,:));

    alphahat = alphagrid(linind(1));
    nuhat = nugrid(linind(1));
    sigsqhat = sigmagrid(linind(1));
    
    varargout = {alphahat,nuhat,sigsqhat};
end