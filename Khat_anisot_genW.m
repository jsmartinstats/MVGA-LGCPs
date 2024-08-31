function varargout = Khat_anisot_genW(SPP,W,rvec,phivec,bw_phi,varargin)

% Returns an estimated matrix of sector-K functions (these are anisotropic
% analogues of Ripley's K function) for the stationary multivariate spatial
% point pattern in SPP. If varargin contains anisotropy parameters, SPP is 
% transformed according to these before proceeding.

% This function allows for data that is originally geometric anisotropic, 
% and post-transformation, is able to deal with the resulting 
% parallelogram-shaped observation window.

% INPUT
% SPP       contains the spatial point pattern of interest. This is an
%           n-by-4 matrix: the first column contains the variable labels;
%           the second is a column of nan values and the third and fourth
%           columns contain the x- and y-coordinates for the points.
% W         contains the limits of the observation window:
%               W = [lower x limit, upper x limit, lower y limit, upper y limit]
% rvec      contains the vector of distances at which to perform estimation
% phivec    contains the vector of angles at which to center the sectors
%           of interest
% bw_phi    bandwidth parameter for the sectors of interest; this is the 
%           angle subtended by each sector
% vargin    potentially contains an angle of anisotropy and ratio of
%           anisotropy, which we may wish to use to transform the data
%           before finding the pcf

% OUTPUT
% vargout   contains an estimated matrix of sector-K functions
%

% last modified by jsmartin.stats@gmail.com in Oct 2017
%%
% Order the data, first according to their x-coordinate, then according to their process labels (this is required)
data = sortrows(SPP,3);
data = sortrows(data,1);

% should this ordering happen after isotropisation?

% separate the x- and y-coordinates (needed for Ripley-based edge-correction)
Wx = [W(1),W(2),W(2),W(1)];
Wy = [W(3),W(3),W(4),W(4)];

% If anisotropy parameters are given, isotropise the data and the observation window
    origdata = data;
    origWx = Wx;
    origWy = Wy;
if nargin>5
    %%
    theta = varargin{1};
    zeta = varargin{2};
    theta_deg = rad2deg(theta); % use degrees to get exact values of sin & cos when mod(theta,pi/2)=0
    U = [cosd(theta_deg),-sind(theta_deg);sind(theta_deg),cosd(theta_deg)]; % use degrees to get exact values of sin & cos when mod(theta,pi/2)=0
    data(:,3:4) = data(:,3:4)*U*diag([1,1./zeta]);
    dummy = [Wx',Wy']*U*diag([1,1./zeta]);
    Wx = dummy(:,1)';
    Wy = dummy(:,2)';
end

% find the number of distinct types of point in the multi-type SPP, as well
% as the number of values of r and phi that we wish to estimate the
% sector-K functions
P = numel(unique(data(:,1)));
numrs = length(rvec);
numphis = length(phivec);

% create variables to contain the isotropised set covariance and our estimators
sectorK_classicint = zeros(P,P,numrs,numphis);

% The estimator used here is periodic with period pi, and it involves averaging two functions evaluated at phi and phi+pi; these functions cannot be
% evaluated at phi>=2pi. We therefore need to ensure that phivec is in the interval [0,pi):
phivec = mod(phivec,2*pi); 
%     this line ensures that phivec is in [0,2pi)...not [0,pi).
%     that's not really a problem though - the anisotropic xpcf is not pi-periodic for p\ne q

for rind=1:numrs
    fprintf('%i...',rind);
    for phiind=1:numphis
        r = rvec(rind);
        phi = phivec(phiind);
        
        n = histc(data(:,1),unique(SPP(:,1)));
        sumn = sum(n);

        dists = pdist(data(:,3:4));
        distsarray = squareform(dists);

        xdistsarray = repmat(data(:,3)',sumn,1) - repmat(data(:,3),1,sumn); % targetpoints-focalpoints
        ydistsarray = repmat(data(:,4)',sumn,1) - repmat(data(:,4),1,sumn); % targetpoints-focalpoints

        % calculate the angle between the vector connecting each pair and the abcissa axis - this is measured in radians, and is calculated by
        % arctan as being in the interval (-pi,pi); we modify the output such that it is instead in the interval [0,2pi)

        anglearray = atan2(ydistsarray,xdistsarray);
        anglearray(anglearray<0) = anglearray(anglearray<0) + 2*pi;
            % note: lower triangle of anglearray will only contain elements in [pi/2,3pi/2)...
            %       and the upper triangle will contain only elements in [0,pi/2) and [3pi/2,2pi)...
            %       due to the fact that the points are ordered by their x-coordinate;
        anglearray(logical(eye(sumn))) = nan;
            % we ensure that the principal diagonal is nan, and not 0; if it is zero, then this will be interpreted below as corresponding to a 
            % pair of points separated by a vector at angle 0.

        %% estimate numerator for each element of K - using similar numerator as used by Moller and Toftaker (2014, Appendix C2)
            % note that we take the edge correction factor into the numerator such that it is itself an estimator for the second-order intensity
            % function.

            % NOTE: we use kernel functions below to ascertain whether the relative bearings between two points are close
            % enough to phi. We must take care when evaluating kernel functions on the relative bearings of two points: if the
            % relative bearing between two points is 2pi-a and we wish to judge the proximity of this to phi=b, for small enough a and b, then
            % the two quantities will be close, but will not show up if we just consider the difference (or squared difference) between the two
            % quantities; we therefore evaluate the angular kernel under two scenarios.
        % use a box kernel to define the sector with half-angle bw_phi
        kernelarray = double(abs(distsarray)<=r & (abs(anglearray-phi)<=0.5*bw_phi | abs(anglearray-phi)>=2*pi-0.5*bw_phi));
        % repeat for mod(phi+pi,2*pi): we take the mean below for marginal processes
        % note: we take the modulus here to ensure that, for A at bearing phi\in[pi,2pi) from B, the bearing of B from A lies in [0,pi) 
        kernelarray2 = double(abs(distsarray)<=r & (abs(anglearray-mod(phi+pi,2*pi))<=0.5*bw_phi | abs(anglearray-mod(phi+pi,2*pi))>=2*pi-0.5*bw_phi));
            
        % for p=q, the geometric anisotropic pcf is periodic with period pi; exploit this by using the relative position of point B to point A
        % as well as the relative position of point A to point B.
        % we take care of this here, and not at the end of the function; it would only be possible to do so at the end of the function if both
        % phi and phi+pi were given as inputs;
        % note: first, we create a conditioning matrix to ensure that we are only considering this modification where p=q
        margindarray = []; % initialise an array to indicate the entries pertaining to the marginal dependencies
        for p=1:P        % this will have a block diagonal structure, and we construct it iteratively
            margindarray = blkdiag(margindarray,true(n(p),n(p)));
        end
        margindarray=logical(margindarray);
        kernelarray(margindarray) = (kernelarray(margindarray) + kernelarray2(margindarray))./2;

        kernelarray(logical(eye(size(kernelarray))))=0;
        % force the leading diagonal to be logical(0); we do not consider the interaction of a point with itself

        % In the estimator of the integrated second order intensity function, ie the numerator of the (cross) K-function, we include an edge
        % correction factor, which is calculated as the scaled area of the overlap between the original observation window W and its translated 
        % self, where the translation is according to the vector separating the point pair of interest and the scaling is according to the given
        % ratio of anisotropy.

        Wxlen = max(origWx)-min(origWx); % note: we use the original window (ie pre-isotropisation) here
        Wylen = max(origWy)-min(origWy);
        xdistsarray2 = squareform(pdist(origdata(:,3))); % note: we use the original data (ie pre-isotropisation) here
        ydistsarray2 = squareform(pdist(origdata(:,4)));
        edgecorrection = (Wxlen-xdistsarray2).*(Wylen-ydistsarray2);
        if exist('zeta','var') % if we have isotropised the data
            edgecorrection = edgecorrection./zeta;
        end
        correctedscaledkernel = kernelarray./edgecorrection;

        % in order to find the sums of the various submatrices of numarray, we
        % create a summed area table (the matrix analogue of a cumulative sum)
        K_num= sat(correctedscaledkernel,n);

%%
        diagels = logical(eye(P,P));
        WA = Wxlen*Wylen;
        if exist('zeta','var') % if we have isotropised the data
            WA = WA/zeta;
        end
        n = histc(data(:,1),unique(SPP(:,1)));
        npmat = repmat(n,1,P);
        nqmat = repmat(n',P,1);
        classicintensity = npmat.*nqmat./(WA^2);
        classicintensity(diagels) = classicintensity(diagels).*(npmat(diagels)-1)./npmat(diagels);
        
        sectorK_classicint(:,:,rind,phiind) = K_num./classicintensity;
    end
end
fprintf('\n');

varargout{1} = sectorK_classicint;
end