function pll = palmloglik_genW_bv_sym(data,W,R,lam_p,lam_q,alpha,nu,sigma,varargin)
% This function calculates the Palm log-cross-likelihood, cf the Palm log-
% likelihood detailed by Dvorak & Prokesova (2012) (using the border edge 
% correction method)

% Note that, if estimating the cross-likelihood, the parameter mu will pertain to
% the second process, i.e. process q when using the palm intensity lam_{0,pq}

% INPUT
% data      the original bivariate Matern LGCP
% W         contains the limits of the observation window, which is assumed to be rectangular and parallel to the x-axis:
%               W = [lower x limit, upper x limit, lower y limit, upper y limit]
% R         the max separation distance at which we consider point pairs
% lam_p     the intensity of the first component of the bivariate SPP
% lam_q     the intensity of the second component of the bivariate SPP
% alpha 	the value of alpha at which to evaluate the Palm log-likelihood
% nu        the value of nu at which to evaluate the Palm log-likelihood
% sigma 	the value of sigma at which to evaluate the Palm log-likelihood
% varargin  may contain angle and ratio variables, which are assumed to describe the anisotropy of the data:
%               varargin{1} - the angle of anisotropy (in radians, in [0,pi])
%               varargin{2} - the ratio of anisotropy (in (0,1))

% OUTPUT:
% pll    contains the estimated Palm log likelihood

% last modified by jsmartin.stats@gmail.com in Oct 2017

%% ~~~~~~~~~~~~~~~~~~~~

% Separate the x- and y-coordinates.
Wx = [W(1),W(2),W(2),W(1)];
Wy = [W(3),W(3),W(4),W(4)];

% find the number of distinct types of point in the multi-type SPP
P = numel(unique(data(:,1)));
if P>2
    error('this function cannot deal with more than two types of point in the observed PP');
end
    


% If anisotropy parameters are given, isotropise the data and the observation window, then 
if nargin>8
    %%

    theta = varargin{1};
    zeta = varargin{2};

    theta_deg = rad2deg(theta); % use degrees to get exact values of sin & cos when mod(theta,pi/2)=0
    U = [cosd(theta_deg),-sind(theta_deg);sind(theta_deg),cosd(theta_deg)]; % use degrees to get exact values of sin & cos when mod(theta,pi/2)=0
    
    data(:,3:4) = data(:,3:4)*U*diag([1,1./zeta]); % is this right?
 
    W = [Wx',Wy']*U*diag([1,1./zeta]);

    Wx = W(:,1)';
    Wy = W(:,2)';
else
    theta = 0;
    zeta = 1;
end

% define the inner region of the window (for edge-correction purposes)
[buffx,buffy]=bufferm([Wx,Wx(1)],[Wy,Wy(1)],R,'in');
% this returns points on both the observation window and the inner region
% use points in [buffx,buffy] that are on the inner region boundary only;
% these are separated from the points on the observation window by [Nan,Nan]
nanindx = find(isnan(buffx),1);
nanindy = find(isnan(buffy),1);
if nanindx==nanindy
    innerregion = [buffx(1:nanindx-1),buffy(1:nanindy-1)];
else
    error('error in finding the innerregion');
end

% establish which datapoints are inside the inner region
innerdata = data(inpolygon(data(:,3),data(:,4),innerregion(:,1),innerregion(:,2)),:);
%%

% define the point pairs y-x, with x, a point in process p, located in the inner region,
% and y, a point in process q, located in W. NB: the notation here is slightly different from
% above - for the univariate case, x and y were somewhat interchangeable.
    diffPP1 = repmat(data(data(:,1)==2,3)',sum(innerdata(:,1)==1),1)-repmat(innerdata(innerdata(:,1)==1,3),1,sum(data(:,1)==2));
    diffPP2 = repmat(data(data(:,1)==2,4)',sum(innerdata(:,1)==1),1)-repmat(innerdata(innerdata(:,1)==1,4),1,sum(data(:,1)==2));
    diffPPr = hypot(diffPP1(:),diffPP2(:));
% restrict attention to those point pair vectors with magnitude in (0,R)
    restrictedinds = diffPPr>0 & diffPPr<R;
%     diffPP = [diffPP1(restrictedinds),diffPP2(restrictedinds)];
    diffPPr_pq = diffPPr(restrictedinds);
%     ndiffPPr_pq = numel(diffPPr_pq);
% repeat these steps for the opposite difference vectors
    diffPP1 = repmat(data(data(:,1)==1,3)',sum(innerdata(:,1)==2),1)-repmat(innerdata(innerdata(:,1)==2,3),1,sum(data(:,1)==1));
    diffPP2 = repmat(data(data(:,1)==1,4)',sum(innerdata(:,1)==2),1)-repmat(innerdata(innerdata(:,1)==2,4),1,sum(data(:,1)==1));
    diffPPr = hypot(diffPP1(:),diffPP2(:));
    restrictedinds = diffPPr>0 & diffPPr<R;
%     diffPP = [diffPP1(restrictedinds),diffPP2(restrictedinds)];
    diffPPr_qp = diffPPr(restrictedinds);
%     ndiffPPr_qp = numel(diffPPr_qp);

diffPPr = sort([diffPPr_pq;diffPPr_qp]); % the magnitude of all point pairs in both difference patterns
ndiffPPr = numel(diffPPr);
numinnerpts_p = sum(innerdata(:,1)==1); % number of focal points (denoted here as x) of type p in the inner region
numinnerpts_q = sum(innerdata(:,1)==2); % number of focal points (denoted here as x) of type q in the inner region

% approximately evaluate K_pq(R) [=K_qp(R)], i.e. Ripley's bivariate K function at R
% this is approximate, as it uses the trapezium rule to evaluate the integral
rvec = (0:0.01:1).*R;
pcfatrvec_pq = exp(covfun_geoanisot_quick(rvec',zeros(numel(rvec),1),'Matern_HW94',sigma,alpha,nu,0,1,1));
RipKatR_pq = trapz(rvec,pcfatrvec_pq.*2.*pi.*rvec');

% evaluate the palm intensities at each point pair in the union of the difference patterns defined above;
% the bivariate Palm intensity, lambda_{0,pq}(r) = lambda_q g_{pq}(r)
% since the pcf is symmetric, lambda_{0,pq}(r)+lambda_{0,qp}(r) = (lam_p+lam_q)g_{pq}(r)
% we calculate here log(lambda_{0,pq}(r)+lambda_{0,qp}(r)) for r given by each point pair in the difference
% pattern union
logpalmintensities = log(lam_p + lam_q) + covfun_geoanisot_quick(diffPPr,zeros(ndiffPPr,1),'Matern_HW94',sigma,alpha,nu,0,1,1);

% hold on;
% plot(diffPPr,logpalmintensity,'.');
% plot(rvec,mu+sigma/2 + pcfatrvec,'r');
% hold off;
% evaluate the palm log likelihood at each of the above point pairs
pll = sum(logpalmintensities) - (numinnerpts_p.*lam_q + numinnerpts_q.*lam_p).*RipKatR_pq;
end