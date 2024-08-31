function [estmu,estlam_pqfit,estalpha,estnu,estsigma] = PP_Maternest_validMV(data,W,thetamat,zetamat,alpha,nu,sigma)
% This function estimates the Matern parameters for the covariance structure underlying the 
% potentially multivariate LGCP provided in data. The LGCP is assumed to be geometric 
% anisotropic, with GA parameters provided in thetamat and zetamat.

% We transform the LGCP into an isotropic setting and use approximate maximum likelihood to
% obtain estimates of the Matern parameters. After this transformation, the (presumably
% rectangular) observation window W will become a parallelogram. We use new results to 
% calculate the required approximate log-likelihood.


% INPUT: 
% data      potentially multivariate GA LGCP data
% W 		the observation window for the data, in the form
% 				[lower x bound, upper x bound, lower y bound, upper y bound]

% thetamat  a matrix of parameters specifying the angles of anisotropy for each of 
%           the marginal and bivariate dependence structures
% zetamat   a matrix of parameters specifying the ratios of anisotropy for each of 
%           the marginal and bivariate dependence structures
% alpha 	gives either:
% 				the range of estimation [alpha_1,alpha_2] in which to seek estimates; or
%				the fixed value of alpha if it is not to be estimated.
% nu		similar to alpha.
% sigma		similar to alpha.
 
% OUTPUT
% 
% estalpha is a matrix of estimated Matern scale parameters
% estnu is a matrix of estimated Matern smoothness (shape) parameters
% estsigma is a matrix of estimated Matern power parameters
% 
% last modified by jsmartin.stats@gmail.com in May 2018

%% Set options for optimisation procedure, fmincon
        options = optimoptions('fmincon');
        options = optimoptions(options,...
            'Algorithm','interior-point',...
            'HonorBounds',true,...
            'SubproblemAlgorithm','factorization',...
            'FiniteDifferenceType','central',...
            'OptimalityTolerance',1e-6... % this is the default, and is sufficient; tried 1e-10 also
            ); %#ok<NASGU>     
%%
P = max(data(:,1));

estmu = nan(P,P);
estlam_pqfit = nan(P,P);
estalpha = 999*ones(P,P);
estnu = nan(P,P);
estsigma = 999*ones(P,P);

if or(alpha(end)>999,sigma(end)>999)
	error('upper bound for alpha or sigma is too high\n')
end

% Separate the x- and y-coordinates of the observation window and find the window area, before isotropization.
Wx = [W(1),W(2),W(2),W(1)];
Wy = [W(3),W(3),W(4),W(4)];
WA_orig = (max(Wx)-min(Wx))*(max(Wy)-min(Wy)); % area of the original observation window (input), as long as it is rectangular and parallel to the x-axis


% first, redefine rvec, for use in coarse minimum contrast procedure. This should be high-resolution,
% as it is being used to evaluate the trapezium rule, i.e. to approximate an integral. The 'coarse' 
% aspect of the min contrast procedure refers to the grid of parameters that we are optimising over.
    rvecmincon = 0.01:0.01:0.25;

    if numel(nu)>1
        inputnu = nu;
        nu = [0.05,0.5,5]; % use discrete candidate smoothnesses; these differ from Moller & Toftaker (2014) as we use a different parameterisation
        nu = nu(nu>=inputnu(1));   % ...as long as they are not below our input lower bound
        nu = nu(nu<=inputnu(2));   % ...and as long as they are not above our input upper bound
    end
    % calculate lower and upper bounds for the scale and power parameters
    lb = [alpha(1),sigma(1)]; 	  %#ok<NASGU> % if the input value for a parameter is a scalar, then the corresponding 
    ub = [alpha(end),sigma(end)]; %#ok<NASGU> % value for ub and lb will match, i.e. the parameter will be fixed.

	% Estimate Matern parameters for marginal covariance structures
	for p=1:P
		% First, find suitable initialisation points for the approx MLE procedure, using coarse min contrast
		[initalpha,~,initsigma] = MinContrastSPP_coarse(data(data(:,1)==p,:),W,rvecmincon(1),rvecmincon(end),rvecmincon(2)-rvecmincon(1),1,alpha,[nu(1),nu(end)],sigma,thetamat(p,p),zetamat(p,p)); 
		count=0; % if the nonlinear optimisation degenerates, we re-run the parameter estimation (with randomly perturbed initialisations), up to a maximum of 10 iterations
		while and(or(any(estalpha(p,p)>0.95*alpha(end)),any(estsigma(p,p)>0.95*sigma(end))),count<10)
			count=count+1;
			if count==1
				% first time round, just use the min-contrast initialisation
				pertinitalpha = initalpha;
				pertinitsigma = initsigma;
			else
				% otherwise, add random perturbation to initial values
				pertinitalpha = initalpha + 0.025*randn;
				% make sure that random perturbation doesn't take initial values outside of valid interval
				while or(pertinitalpha<alpha(1),pertinitalpha>alpha(end))
					pertinitalpha = initalpha + 0.025*randn;
				end
				pertinitsigma = initsigma + 0.25*randn;
				while or(pertinitsigma<sigma(1),pertinitsigma>sigma(end))
					pertinitsigma = initsigma + 0.25*randn;
				end
			end
			% estimate the Matern parameters through numerically optimising the
			% Palm log-likelihood, using fmincon. The mean of the underlying 
			% GRF (parameterised by mu) can be estimated directly as the MPLE
			% for the intensity is available in closed form...if we first 
			% calculate the inner region and the difference process for the 
			% Palm log-likelihood:
				R=0.1;
				PLLcell = PLLinputcell(data(data(:,1)==p,:),W,R,thetamat(p,p),zetamat(p,p));
				
			mles=nan(numel(nu),2);
			fvals = nan(numel(nu),1);
			for nuindex = 1:numel(nu)
				candnu = nu(nuindex);
	%             % using initialising values of the Matern parameters,
	%             % calculate the corresponding MPLE for mu; use this as the
	%             % initial value for mu
				[~,mles(nuindex,:),fvals(nuindex,1)] = evalc('fmincon(@(params) -palmloglik_genW(PLLcell,W,R,''MPLE'',params(1),candnu,params(2),thetamat(p,p),zetamat(p,p)),[pertinitalpha,pertinitsigma],[],[],[],[],lb,ub,[],options);');
			end
			[~,mlenuind] = min(fvals);
			estalpha(p,p) = mles(mlenuind,1);
			estnu(p,p) = nu(mlenuind);
			estsigma(p,p) = mles(mlenuind,2);
		end
		fprintf('%i,',count);
		% use the estimated values of the Matern params, along with the 
		% difference pattern in PLLcell, to calculate the
		% MPLE of mu
		numinnerfocalpts = size(PLLcell{2},1);
		diffPPr = PLLcell{3};
		estmu(p,p) = log(MPLElam(numel(diffPPr),numinnerfocalpts,R,estalpha(p,p),estnu(p,p),estsigma(p,p))) - estsigma(p,p)/2; % % MPLE of lambda, see Dvorak & Prokesova (2012)                    
	end

	if P==2
		% estimate off-diagonal entries - only good for P=2
		p = 2; q = 1;
		% For the fitted model to be valid, we require \nu_pq to be bounded below by the mean of the
		% corresponding marginal parameters (Condition 1 in the manuscript). We therefore discard all
		% candidate values of nu_pq that contravene this requirement.
		nupqLB = (estnu(p,p)+estnu(q,q))/2;% lower bound for valid specifications of nu_pq
		validnupqvec = nu(nu>=nupqLB);
			
		% we also require the inverse correlation length (4*nu/alpha^2) to satisfy a similar constraint;
		% this is detailed in Remark 1 of the manuscript.
		iclpp = 4*estnu(p,p)/estalpha(p,p)^2; % marginal inverse correlation length for the process p
		iclqq = 4*estnu(q,q)/estalpha(q,q)^2; % marginal inverse correlation length for the process q
		iclpqLB = (iclpp + iclqq)./2; % lower bound for the inverse cross-correlation length between processes p & q
		% for fixed nu, the lower bound on the icl translates to an upper bound on alpha_pq
		alphapqUBvec = sqrt(4.*validnupqvec./iclpqLB); % upper bound for alpha_pq
		
		% we also need the fitted value of sigma_pq to satisfy a constraint (provided in Remark 2 of the
		% manuscript). This constraint is in the form of an upper bound for sigma_pq, and this UB is 
		% dependent on alpha_pq. Since alpha_pq and sigma_pq are being estimated together via fmincon, 
		% we can achieve this through specification of a nonlinear constraint in fmincon.
		[initalpha,~,initsigma] = MinContrastSPP_coarse(data(ismember(data(:,1),[p,q]),:),W,rvecmincon(1),rvecmincon(end),rvecmincon(2)-rvecmincon(1),1,alpha,[nu(1),nu(end)],sigma,thetamat(p,q),zetamat(p,q)); 
		count=0; % if the nonlinear optimisation degenerates, we re-run the parameter estimation (with randomly perturbed initialisations), up to a maximum of 10 iterations
		while and(or(any(estalpha(p,q)>0.95*alpha(end)),any(estsigma(p,q)>0.95*sigma(end))),count<10)
			count=count+1;
			if count==1
				% first time round, just use the min-contrast initialisation
				pertinitalpha = initalpha;
				pertinitsigma = initsigma;
			else
				% otherwise, add random perturbation to initial values
				pertinitalpha = initalpha + 0.025*randn;
				% make sure that random perturbation doesn't take initial values outside of valid interval
				while or(pertinitalpha<alpha(1),pertinitalpha>alpha(end))
					pertinitalpha = initalpha + 0.025*randn;
				end
				pertinitsigma = initsigma + 0.25*randn;
				while or(pertinitsigma<sigma(1),pertinitsigma>sigma(end))
					pertinitsigma = initsigma + 0.25*randn;
				end
			end
				
			lamp_MPLE = exp(estmu(p,p) + estsigma(p,p)/2);
			lamq_MPLE = exp(estmu(q,q) + estsigma(q,q)/2);
			
			mles=nan(numel(validnupqvec),2);
			fvals = nan(numel(validnupqvec),1);
			for nuindex = 1:numel(validnupqvec)
				candnu = validnupqvec(nuindex);
				
				% The nonlinear constraint for fmincon will depend on additional parameters W_pp and W_qq, which we 
				% calculate here - these ensure that the marginal power parameters relate to the marginal
				% scale and smoothness parameters in the appropriate way. W_pp and W_qq are calculated outside of 
				% the fmincon implementation for reasons of computational efficiency. They are calculated inside the
				% nu loop as Del_nu is dependent on the current candidate value of nu.
				Delnu = candnu-(estnu(p,p)+estnu(q,q))/2; % value of Delta_nu corresponding to estimated marginal values of nu and candidate cross-covariance nu
				V = nan(P,1);
				for p=1:P
					V(p) = sqrt( (estsigma(p,p)*zetamat(p,p)/pi) * (4*estnu(p,p)/(estalpha(p,p)^2))^(Delnu+estnu(p,p)) * 1/gamma(estnu(p,p)) );
				end
				
				lb = [alpha(1),sigma(1)]; %#ok<NASGU> % lower bounds on lambda_p and lambda_q are both 10; lower bounds on alpha and sigma are inherited from inputs
				ub = [alphapqUBvec(nuindex),sigma(end)]; %#ok<NASGU> % renew upper bound for alpha; restriction on sigmapq is left to the nonlinear constraint in fmincon
														  % upper bounds on lambda_p and lambda_q are both 1000;
				% check that the initialising values of alpha and sigma satisfy their respective constraints
				if initalpha>alphapqUBvec(nuindex)
					initalpha = alphapqUBvec(nuindex);
				end
				initiclpq = 4*candnu/(initalpha^2);
				initsigUB = (V(p)*V(q)*pi/zetamat(p,q)) * initiclpq.^(-candnu) * gamma((estnu(p,p)+estnu(q,q))/2 + 1)*gamma(candnu)/gamma(candnu + 1);
				if initsigma>initsigUB
					initsigma = initsigUB;
				end
				[~,mles(nuindex,:),fvals(nuindex,1)] = evalc('fmincon(@(params) -palmloglik_genW_bv_sym(data(ismember(data(:,1),[p,q]),:),W,R,lamp_MPLE,lamq_MPLE,params(1),candnu,params(2),thetamat(p,q),zetamat(p,q)),[pertinitalpha,pertinitsigma],[],[],[],[],lb,ub,@(params) mycon(params,candnu,estnu(p,p),estnu(q,q),zetamat(p,q),V(p),V(q)),options);');
			end
			[~,mlenuind] = min(fvals);
			estalpha(p,q) = mles(mlenuind,1);
			estnu(p,q) = validnupqvec(mlenuind);
			estsigma(p,q) = mles(mlenuind,2);
		end
		fprintf('%i...',count);
	end
end
%%    

function [c,ceq] = mycon(candparams,nupq,nupp,nuqq,zetapq,Vp,Vq) 
% candparams = [alpha_pq,sigma_pq]
% nupq is the current value of nu_pq being used
% fmincon performs constrained minimisation, using this nonlinear constraint
% needs revising for P\ge3 (Delta_nu will be different)
iclpq = 4*nupq/(candparams(1)^2); % inverse correlation length (p,q)
sigmapqUB = (Vp*Vq*pi/zetapq) * iclpq.^(-nupq) * gamma((nupp+nuqq)/2 + 1)*gamma(nupq)/gamma(nupq + 1);
c = candparams(2) - sigmapqUB;
ceq = [];
end

function mple = MPLElam(ndiffr,n_inner,R,alp,nu,sig)
% this calculates the maximum Palm-likelihood estimator for the intensity
% of the point pattern with ndiffr points in its first-order difference 
% and with n_inner points in the inner region, defined by buffer radius R.
    rvec = (0:0.01:1).*R;
    pcfatrvec_pq = exp(covfun_geoanisot_quick(rvec',zeros(numel(rvec),1),'Matern_HW94',sig,alp,nu,0,1,1));
    RipKatR_p = trapz(rvec,pcfatrvec_pq.*2.*pi.*rvec');
    mple = ndiffr/(RipKatR_p*n_inner);
end