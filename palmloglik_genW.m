function pll = palmloglik_genW(data,W,R,mu,alpha,nu,sigma,varargin)
% This function calculates the Palm log-likelihood, as first defined by
% Tanaka et al (2008), but as detailed by Dvorak & Prokesova (2012) (esp
% with respect to the border edge correction method)

% Note that, if estimating the cross-likelihood, the parameter mu will pertain to
% the second process, i.e. process q when using the palm intensity lam_{0,pq}

% INPUT
% data      either a matrix containing only the original Matern LGCP,
%           or a cell containing three elements:
%               data{1} - the original Matern LGCP
%               data{2} - a matrix containing the subset of the data that 
%                         lies within the inner region (calculated using
%                         the same R as is used in this program.
%               data{3} - a vector containing the absolute values of the
%                         distances between all point pairs (i.e. |x-y| and
%                         |y-x|)
%           NB: if a cell is provided, it is assumed to contain isotropic,
%           or already-isotropised data; if it is not a cell, then the data
%           is isotropised using the given anisotropy parameters
% W         contains the limits of the observation window, which is assumed to be rectangular and parallel to the x-axis:
%               W = [lower x limit, upper x limit, lower y limit, upper y limit]
% R         the max separation distance at which we consider point pairs
% mu        either: the string 'MPLE', in which case the maximum
%           Palm log-likelihood estimator is used; or a value of at which
%           to evaluate the Palm log-likelihood
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
WA = (max(Wx)-min(Wx))*(max(Wy)-min(Wy)); % area of the original observation window (input), as long as it is rectangular and parallel to the x-axis

% If data is a cell, extract the various elements
if iscell(data)
    if numel(varargin)>0
        warning('Anisotropy parameters are unused: input cell for PLL is assumed to contain isotropic data');
    end
    datamat = data{1};
    innerdata = data{2};
    diffPPr = data{3};
    data = datamat;
    % ensure that diffPPr is sorted, and contains only values in (0,R)
    diffPPr = sort(diffPPr);
    restrictedinds = diffPPr>0 & diffPPr<R;
    diffPPr = diffPPr(restrictedinds);
    % set flag to avoid calculating innerdata and diffPPr below
    calculatediffPPr = false;
else
    calculatediffPPr = true;
end

numdata=size(data,1);
% find the number of distinct types of point in the multi-type SPP
P = numel(unique(data(:,1)));
if P>1
    error('this function cannot deal with more than one type of point in the observed PP');
end


if calculatediffPPr
    if numel(varargin)>0
    % If anisotropy parameters are given, isotropise the data and the observation window, 
        theta = varargin{1};
        zeta = varargin{2};

        theta_deg = rad2deg(theta); % use degrees to get exact values of sin & cos when mod(theta,pi/2)=0
        U = [cosd(theta_deg),-sind(theta_deg);sind(theta_deg),cosd(theta_deg)]; % use degrees to get exact values of sin & cos when mod(theta,pi/2)=0

        data(:,3:4) = data(:,3:4)*U*diag([1,1./zeta]); % is this right?

        W = [Wx',Wy']*U*diag([1,1./zeta]);

        Wx = W(:,1)';
        Wy = W(:,2)';
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
    %     [inpts,onpts] = inpolygon(buffx,buffy,Wx,Wy);
    %     innerregion = [buffx(inpts&~onpts),buffy(inpts&~onpts)];
    end

    innerWA = polyarea(innerregion(:,1),innerregion(:,2)); % area of the inner region
    % establish which datapoints are inside the inner region
    innerdata = data(inpolygon(data(:,3),data(:,4),innerregion(:,1),innerregion(:,2)),:);
    numinnerfocalpts = size(innerdata,1);

    if P==1
        % define the point pairs x-y, with x in the inner region and y in W
    %     diffPP1 = repmat(innerdata(:,3),1,numdata)-repmat(data(:,3)',numinnerdata,1);
        diffPP1 = repmat(data(:,3)',numinnerfocalpts,1)-repmat(innerdata(:,3),1,numdata);
    %     diffPP2 = repmat(innerdata(:,4),1,numdata)-repmat(data(:,4)',numinnerdata,1);
        diffPP2 = repmat(data(:,4)',numinnerfocalpts,1)-repmat(innerdata(:,4),1,numdata);
        diffPPr = hypot(diffPP1(:),diffPP2(:));
    %     numinnerfocalpts = numinnerdata;% number of focal points (denoted here as y) in the inner region
    else % i.e. if P==2
        % define the point pairs y-x, with x, a point in process p, located in the inner region,
        % and y, a point in process q, located in W. NB: the notation here is slightly different from
        % above - for the univariate case, x and y were somewhat interchangeable.
        diffPP1 = repmat(data(data(:,1)==2,3)',sum(innerdata(:,1)==1),1)-repmat(innerdata(innerdata(:,1)==1,3),1,sum(data(:,1)==2));
        diffPP2 = repmat(data(data(:,1)==2,4)',sum(innerdata(:,1)==1),1)-repmat(innerdata(innerdata(:,1)==1,4),1,sum(data(:,1)==2));
        diffPPr = hypot(diffPP1(:),diffPP2(:));
        numinnerfocalpts = sum(innerdata(:,1)==1); % number of focal points (denoted here as x) in the inner region
    end
    % restrict attention to those point pair vectors with magnitude in (0,R)
    restrictedinds = diffPPr>0 & diffPPr<R;
    diffPP = [diffPP1(restrictedinds),diffPP2(restrictedinds)];
    diffPPr = diffPPr(restrictedinds);
    diffPPr = sort(diffPPr);
    ndiffPPr = numel(diffPPr);
else
    numinnerfocalpts = size(innerdata,1);
    ndiffPPr = numel(diffPPr);
end
%%
% approximately evaluate K(R), i.e. Ripley's K function at R, for the
% (isotropised) LGCP
% this is approximate, as it uses the trapezium rule to evaluate the integral
rvec = (0:0.01:1).*R;
pcfatrvec_pq = exp(covfun_geoanisot_quick(rvec',zeros(numel(rvec),1),'Matern_HW94',sigma,alpha,nu,0,1,1));
RipKatR_pq = trapz(rvec,pcfatrvec_pq.*2.*pi.*rvec');

if strcmp(mu,'MPLE')
    % If the input value mu is the string 'MPLE', then we use the value of
    % lambda that maximises the Palm log-likelihood for the given values of
    % alpha, nu and sigma. For details, see e.g. Dvorak & Prokesova (2012)
    lambda = ndiffPPr/(RipKatR_pq*numinnerfocalpts); % MPLE of lambda, see Dvorak & Prokesova (2012)
else
    lambda = exp(mu + sigma/2); % this is the LGCP intensity corresponding to the input values mu and sigma
end

% evaluate the palm intensity at each point pair; this is equal to the
% product of the intensity of the (isotropised) LGCP and its pcf
logpalmintensity = log(lambda) + covfun_geoanisot_quick(diffPPr,zeros(ndiffPPr,1),'Matern_HW94',sigma,alpha,nu,0,1,1);

% evaluate the palm log likelihood at each of the above point pairs
pll = sum(logpalmintensity) - numinnerfocalpts.*lambda.*RipKatR_pq;
end