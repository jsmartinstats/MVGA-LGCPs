function varargout = Gsthat_anisot_genW(SPP,W,rvec,phivec,kerneltype,bandwidths,varargin)

% Returns a matrix of anisotropic (cross-)pair correlation functions for the
% stationary multivariate spatial point pattern in SPP. If varargin contains 
% anisotropy parameters, SPP is transformed according to these before proceeding.

% Unlike Gsthat_anisot.m, this function allows for data that is originally 
% geometric anisotropic, and post-transformation, is able to deal with the 
% resulting parallelogram-shaped observation window.

% INPUT
% SPP       contains the spatial point pattern of interest. This is an
%           n-by-4 matrix: the first column contains the variable labels;
%           the second is a column of nan values and the third and fourth
%           columns contain the x- and y-coordinates for the points.
% W         contains the limits of the observation window:
%               W = [lower x limit, upper x limit, lower y limit, upper y limit]
% rvec      contains the vector of distances at which to estimate the (x)pcf
% phivec    contains the vector of angles at which to estimate the (x)pcf
% vargin    potentially contains an angle of anisotropy and ratio of
%           anisotropy, which we may wish to use to transform the data
%           before finding the pcf

% OUTPUT
% vargout   contains some subset of the following, in this order:
%     varargout{1} = Gsthats_classicint;
%     varargout{2} = Gsthats_adapint;
%

% last modified by jsmartin.stats@gmail.com in Nov 2017

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
if nargin>6
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

% find the number of distinct types of point in the multi-type SPP, as well as the number of values of r and phi that we wish to measure the
% anisotropic xpcf for.
P = numel(unique(data(:,1)));
numrs = length(rvec);
numphis = length(phivec);

% create variables to contain the isotropised set covariance and our estimators
isotsetcov = nan(numrs,numphis);
Gsthats_adapint = zeros(P,P,numrs,numphis);
Gsthats_classicint = zeros(P,P,numrs,numphis);

% The estimator used here is periodic with period pi, and it involves averaging two functions evaluated at phi and phi+pi; these functions cannot be
% evaluated at phi>=2pi. We therefore need to ensure that phivec is in the interval [0,pi):
phivec = mod(phivec,2*pi); 
%     this line ensures that phivec is in [0,2pi)...not [0,pi).
%     that's not really a problem though - the anisotropic xpcf is not pi-periodic for p\ne q

% Set bandwidths for distance- and angle-measuring kernel functions
bw_r = bandwidths(1);
bw_phi = bandwidths(2);
if bw_r==0 || bw_phi==0
    error('Cannot have zero bandwidth for the kernel function');
end

for rind=1:numrs
    fprintf('%i...',rind);
    for phiind=1:numphis
%         if mod(phiind,ceil(numphis/4))==0
%             fprintf('.');
%         end
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

        %% estimate numerator for each element of G - using similar numerator as used by Moller and Toftaker (2014, Appendix C2)
            % note that we take the edge correction factor into the numerator such that it is itself an estimator for the second-order intensity
            % function.

            % NOTE: we use kernel functions below to ascertain whether the separation distance or relative bearings between two points are close
            % enough to r and phi, resp. We must take care when evaluating kernel functions on the relative bearings of two points: if the
            % relative bearing between two points is 2pi-a and we wish to judge the proximity of this to phi=b, for small enough a and b, then
            % the two quantities will be close, but will not show up if we just consider the difference (or squared difference) between the two
            % quantities; we therefore evaluate the angular kernel under two scenarios.
        switch kerneltype
            case 'box'
                % use 1d uniform (i.e. box) kernels with standard deviation bw_r and bw_phi
                kernelarray = (1/(sqrt(12)*bw_r*sqrt(12)*bw_phi)).*double(abs(distsarray-r)<=0.5*sqrt(12)*bw_r & (abs(anglearray-phi)<=0.5*sqrt(12)*bw_phi | abs(anglearray-phi)>=2*pi-0.5*sqrt(12)*bw_phi));
                % repeat for mod(phi+pi,2*pi): we take the mean below for marginal processes
                % note: we take the modulus here to ensure that, for A at bearing phi\in[pi,2pi) from B, the bearing of B from A lies in [0,pi) 
                kernelarray2 = (1/(sqrt(12)*bw_r*sqrt(12)*bw_phi)).*double(abs(distsarray-r)<=0.5*sqrt(12)*bw_r & (abs(anglearray-mod(phi+pi,2*pi))<=0.5*sqrt(12)*bw_phi | abs(anglearray-mod(phi+pi,2*pi))>=2*pi-0.5*sqrt(12)*bw_phi));
            case 'Epan'
                % or use Epanechnikov kernels with std dev bw_r and bw_phi
                kernelarray_r = (0.75/(sqrt(5)*bw_r)).*(1-((distsarray-r)./(sqrt(5)*bw_r)).^2);
                kernelarray_phi_1 = (0.75/(sqrt(5)*bw_phi)).*(1-((anglearray-phi)./(sqrt(5)*bw_phi)).^2);
                kernelarray_phi_2 = (0.75/(sqrt(5)*bw_phi)).*(1-((2*pi-abs(anglearray-phi))./(sqrt(5)*bw_phi)).^2);
                kernelarray_r(kernelarray_r<0) = 0;
                kernelarray_phi_1(kernelarray_phi_1<0) = 0;
                kernelarray_phi_2(kernelarray_phi_2<0) = 0;
                kernelarray = kernelarray_r.*(kernelarray_phi_1 + kernelarray_phi_2);
                % repeat calculation of phi kernel for mod(phi+pi,2*pi): we take the mean below for marginal processes
                kernelarray2_phi_1 = (0.75/(sqrt(5)*bw_phi)).*(1-((anglearray-mod(phi+pi,2*pi))./(sqrt(5)*bw_phi)).^2);
                kernelarray2_phi_2 = (0.75/(sqrt(5)*bw_phi)).*(1-((2*pi-abs(anglearray-mod(phi+pi,2*pi)))./(sqrt(5)*bw_phi)).^2);
                kernelarray2_phi_1(kernelarray2_phi_1<0) = 0;
                kernelarray2_phi_2(kernelarray2_phi_2<0) = 0;
                kernelarray2 = kernelarray_r.*(kernelarray2_phi_1 + kernelarray2_phi_2);
            otherwise
                error('choose a valid kernel type/n');
        end

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

        scaledkernelarray = kernelarray./r;

        % In the estimator of the second order intensity function, ie the numerator of the (cross) pair intensity function, we include an edge
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
        correctedscaledkernel = scaledkernelarray./edgecorrection;

        % in order to find the sums of the various submatrices of numarray, we
        % create a summed area table (the matrix analogue of a cumulative sum)
        g_num= sat(correctedscaledkernel,n);

        %% estimation of the surface-adapted intensity (Illian et al, 2008, p194)
%         [ec1array,~] = edgecorr_disc(data,Wx,Wy,r); % does this have to be Wx_cw and Wy_cw?
%         boundoverlapcell = mat2cell(ec1array,n);
%         boundoverlappropnarray = cellfun(@(M1) sum(M1),boundoverlapcell);
% 
%         % ISOTROPISED SET COVARIANCE DOES NOT VARY WITH rind OR phiind...IT CAN BE CALCULATED BEFOREHAND
%         
%         % isotropised set covariance, calculated as in Illian et al (2008, p486)
%         a=Wylen;
%         b=Wxlen;
%         if r<=a
%             isotsetcov(rind,phiind) = a*b - 2*r*(a+b)/pi + (r^2)/pi;
%         else
%             isotsetcov(rind,phiind) = a*b*(2*asin(a/r) - a/b - 2*(r/a-sqrt((r/a)^2-1)))/pi;
%         end
% 
%         Gadaptedintensityarray = boundoverlappropnarray./isotsetcov(rind,phiind);
%         
%         Gsthats_adapint(:,:,rind,phiind) = g_num./(nansum(Gadaptedintensityarray,2)*nansum(Gadaptedintensityarray,2)');
%         
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
        
        Gsthats_classicint(:,:,rind,phiind) = g_num./classicintensity;
    end
end
fprintf('\n');

varargout{1} = Gsthats_classicint;
varargout{2} = 1;%Gsthats_adapint;
end