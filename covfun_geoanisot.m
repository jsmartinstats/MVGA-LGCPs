function covmatoutput = covfun_geoanisot(u1,u2,covtype,sig,alp,gam,the,ome,zet)
% This function evaluates one of three geometric anisotropic covariance functions 
% over the spatial lags specified by the inputs u1,u2. This is achieved by evaluating
% the corresponding isotropic covariance function on a geometrically deformed space. 
% The three isotropic covariance functions available are:
% 		- Power Exponential, with power, scale and shape (exponent) parameters
%		- Matern (parameterised as in Stein, 1999), with power, scale and shape parameters
%		- Matern (parameterised as in Handcock & Wallis 1994), with  power, scale and shape
%				parameters
% 
% The dimension of each covariance matrix is determined by the dimensions of the 
% parameter inputs; no checks are carried out to ensure consistency of the dimensions of 
% the inputs. The covariance matrices are concatenated and returned in an n*d-by-d array, 
% where n is the number of spatial lags.
% 
% 
% INPUT
% u1,u2     	coordinates for the spatial lags at which the covariance matrices are to be evaluated
% covtype		string determining the covariance model to use. 
% 				One of: 'powexp', 'Matern_Stein99' or 'Matern_HW94' 
% sig			power parameter matrix (d-by-d)
% alp			scale parameter matrix (d-by-d)
% gam			shape (exponent) parameter matrix (d-by-d)
% the			angle of anistotropy parameter matrix (d-by-d)
% ome			scale of anisotropy parameter matrix (d-by-d)
%				(should be set to ones(d) to avoid identifiability issues with alp)
% zet			ratio of anisotropy parameter matrix (d-by-d)
%
% OUTPUT
% covmatoutput	n*d-by-d output array containing n concatenated d-by-d covariance matrices
% 
% 
% last modified by jsmartin.stats@gmail.com in Nov 2017
%%
    s2 = sin(the).^2;
    c2 = cos(the).^2;
    sc = sin(the).*cos(the);
    invzetsq = 1./(zet.^2);
    rmat = (1./ome).*sqrt( (c2+invzetsq.*s2).*u1.^2 + 2.*(1-invzetsq).*sc.*u1.*u2 + (s2+invzetsq.*c2).*u2.^2 );
    if strcmp(covtype,'powexp')
        covmatoutput = sig .* exp(-(alp.*rmat).^gam);
    elseif strcmp(covtype,'Matern_Stein99')
        phi = sig;
        nu = gam;
        if any(rmat(:)==0)
            zinds = rmat==0;
            nzinds = rmat~=0;
            covmatoutput = nan(size(rmat));
            covmatoutput(nzinds) = (sqrt(pi)*phi(nzinds) .*((alp(nzinds).*rmat(nzinds)).^nu(nzinds)).*arrayfun(@(x,y) besselk(x,y),nu(nzinds),alp(nzinds).*rmat(nzinds)))./(2.^(nu(nzinds)-1).*gamma(nu(nzinds)+0.5).*(alp(nzinds).^(2*nu(nzinds))));
            covmatoutput(zinds) = (sqrt(pi).*phi(zinds).*gamma(nu(zinds)))./(alp(zinds).^(2.*nu(zinds)).*gamma(nu(zinds)+0.5));
        else
            covmatoutput = (sqrt(pi)*phi' .*((alp.*rmat).^nu).*arrayfun(@(x,y) besselk(x,y),nu,alp.*rmat))./(2.^(nu-1).*gamma(nu+0.5).*(alp.^(2*nu)));
        end
    elseif strcmp(covtype,'Matern_HW94')
        phi = sig;
        nu = gam;
        if any(rmat(:)==0)
            zinds = rmat==0;
            nzinds = rmat~=0;
            covmatoutput = nan(size(rmat));
            covmatoutput(nzinds) = (phi(nzinds).*(2.*sqrt(nu(nzinds)).*rmat(nzinds)./alp(nzinds)).^nu(nzinds) .* arrayfun(@(x,y) besselk(x,y),nu(nzinds),2.*sqrt(nu(nzinds)).*rmat(nzinds)./alp(nzinds)))./(2.^(nu(nzinds)-1).*gamma(nu(nzinds)));
            covmatoutput(zinds) = phi(zinds);
        else
            covmatoutput = (phi.*(2.*sqrt(nu).*rmat./alp).^nu .* arrayfun(@(x,y) besselk(x,y),nu,2.*sqrt(nu).*rmat./alp))./(2.^(nu-1).*gamma(nu));
        end
    else
        error('invalid covariance structure type - choose "powexp", "Matern_Stein99" or "Matern_HW94"');
    end
end