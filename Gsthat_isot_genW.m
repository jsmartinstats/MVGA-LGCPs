function varargout = Gsthat_isot_genW(SPP,W,rvec,kerneltype,bw_r,varargin)

% Returns a matrix of isotropic (cross-)pair correlation functions for the
% stationary multivariate spatial point pattern in SPP. If varargin contains 
% anisotropy parameters, SPP is assumed to be anisotropic and is therefore
% isotropized before proceeding. 

% Unlike Gsthat_isot.m, this function allows for data that is originally 
% geometric anisotropic, and post-transformation, is able to deal with the 
% resulting parallelogram-shaped observation window.

% INPUT
% SPP       contains the spatial point pattern of interest. This is an
%           n-by-4 matrix: the first column contains the variable labels;
%           the second is a column of nan values and the third and fourth
%           columns contain the coordinates for the points.
% W         contains the limits of the observation window:
%               W = [lower x limit, upper x limit, lower y limit, upper y limit]
% rvec      contains the vector of distances at which to estimate the (x)pcf
% kerneltype  is the type of kernel used for the (x)pcf estimation:
%               'Epan' results in the use of the Epanechnikov kernel;
%               'box' results in the use of the box kernel.
% bw_r      is a bandwidth parameter for the kernel function..
%           The kernel is scaled such that bw_r is the standard deviation; this
%           matches the convention in spatstat.
% varargin  potentially contains an angle of anisotropy and ratio of
%           anisotropy, which we may wish to use to transform the data
%           before finding the pcf

% OUTPUT
% vargout   contains some subset of the following, in this order:
%     varargout{1} = Gsthats_classicint;
%     varargout{2} = Gsthats_adapint;
%

% last modified by jsmartin.stats@gmail.com in Dec 2017

%%
% Order the data, first according to their process labels, then according to their x-coordinate (this is required)
data = sortrows(SPP,3);
data = sortrows(data,1);

Wx = [W(1),W(2),W(2),W(1)];
Wy = [W(3),W(3),W(4),W(4)];

if nargin>5
    %%
    theta = varargin{1};
    zeta = varargin{2};
    theta_deg = rad2deg(theta); % use degrees to get exact values of sin & cos when mod(theta,pi/2)=0
    U = [cosd(theta_deg),-sind(theta_deg);sind(theta_deg),cosd(theta_deg)]; % use degrees to get exact values of sin & cos when mod(theta,pi/2)=0
    origdata = data;
    origWx = Wx;
    origWy = Wy;
    data(:,3:4) = data(:,3:4)*U*diag([1,1./zeta]);
    W = [Wx',Wy']*U*diag([1,1./zeta]);
    Wx = W(:,1)';
    Wy = W(:,2)';
    if length(varargin)>2
        % an optional argument 'ignoresmallr' is included
        ignoresmallr = varargin{3};
    else
		ignoresmallr = 0; % default is False/0
	end
end

% find the number of distinct types of point in the multi-type SPP, as well as the number of distances at which we wish to measure the isotropic xpcf
P = numel(unique(SPP(:,1)));
numrs = length(rvec);

% create variables to contain the isotropised set covariance and our estimators
isotsetcov = nan(1,numrs);
Gsthats_adapint = zeros(P,P,numrs); % (x)pcf estimator using a G-adapted intensity estimator
Gsthats_classicint = zeros(P,P,numrs); % (x)pcf estimator using a classical intensity estimator

% We also need to ensure that the vertices of the observation window are listed clockwise...
centre = [mean(Wx),mean(Wy)];
Wleftinds = Wx<centre(1);
Wrightinds = Wx>=centre(1);
Wbelowinds = Wy<centre(2);
Waboveinds = Wy>=centre(2);
angles = nan(size(Wx));
angles(Wrightinds&Waboveinds) = atan((Wy(Wrightinds&Waboveinds)-centre(2))./(Wx(Wrightinds&Waboveinds)-centre(1)));
angles(Wleftinds&Waboveinds) = pi + atan((Wy(Wleftinds&Waboveinds)-centre(2))./(Wx(Wleftinds&Waboveinds)-centre(1)));
angles(Wleftinds&Wbelowinds) = pi + atan((Wy(Wleftinds&Wbelowinds)-centre(2))./(Wx(Wleftinds&Wbelowinds)-centre(1)));
angles(Wrightinds&Wbelowinds) = 2*pi + atan((Wy(Wrightinds&Wbelowinds)-centre(2))./(Wx(Wrightinds&Wbelowinds)-centre(1)));
[~,Worder] = sort(angles);
% ...and, to start with, we choose to make the first point in the list of vertices that which has the lowest y-coordinate (if there are two, then we
% choose the one that also has the lowest x-coordinate).
firstpoint = find(Wy==min(Wy));
if numel(firstpoint)>1
%     firstpoint = firstpoint(find(Wx(firstpoint)==min(Wx(firstpoint))));
    firstpoint = firstpoint(Wx(firstpoint)==min(Wx(firstpoint)));
end
Worder = circshift(Worder',-find(Worder==firstpoint)+1)';
Wx_cw = Wx(Worder);
Wy_cw = Wy(Worder);

% to set an arbitrary convention for the purposes of this program, we wish to have the data oriented such that the longest edge of the observation
% window is horizontal. This will be useful for our calculations in the case where W is a parallelogram. Therefore if, as currently
% labelled, the edge between vertices 1 and 2 is longer than the edge between vertices 1 and 4, we rotate the data and the observation window such
% that this is the horizontal base; otherwise, we rotate the other way such that 
% WHY ORIENTATE THE DATA SUCH THAT THE LONGEST EDGE OF THE WINDOW IS HORIZONTAL? We do this so that the edge correction is easier - the overlapping 
% area in translated windows can then be found using the length of the base and the height of the parallelogram, along with the x-dists and y-dists
% between points

if Wy_cw(1)~=Wy_cw(2) % if the edge between vertices 1 and 2 is not horizontal
    edgelengths = sqrt(diff(Wx_cw).^2 + diff(Wy_cw).^2);
    if edgelengths(1)>=edgelengths(2) % we compare edge12 with edge23 (which is the same as edge41)
        artrotang = atan((Wy_cw(2)-Wy_cw(1))./(Wx_cw(2)-Wx_cw(1))); % artificial angle of rotation
        artrotmat = [cos(artrotang),-sin(artrotang);sin(artrotang),cos(artrotang)];
        data(:,3:4) = data(:,3:4)*artrotmat;
        W_cw = ([Wx_cw',Wy_cw']*artrotmat);
        Wx_cw = W_cw(:,1)';
        Wy_cw = W_cw(:,2)';
    elseif edgelengths(1)<edgelengths(2)
        artrotang = atan((Wy_cw(4)-Wy_cw(1))./(Wx_cw(1)-Wx_cw(4))); % artificial angle of rotation
        artrotmat = [cos(artrotang),sin(artrotang);-sin(artrotang),cos(artrotang)];
        data(:,3:4) = data(:,3:4)*artrotmat;
        W_cw = ([Wx_cw',Wy_cw']*artrotmat);
        Wx_cw = circshift(W_cw(:,1),1)';
        Wy_cw = circshift(W_cw(:,2),1)';
    end
else
    artrotang = 0;
end

% establish the dimensions and shape of the observation window

    % find the longest edge - this is the major diagonal in a parallelogram.
    edgesanddiags = pdist([Wx_cw',Wy_cw']);
    majdiagind = edgesanddiags==max(edgesanddiags);
    % If majdiagind indicates two major diagonals...then W is a rectangle (or a square)
    if sum(majdiagind)>1
        Wshape='rectangle';
    else
        Wshape='parallelogram';
    end

if strcmp(Wshape,'rectangle')
    Wx1 = min(Wx_cw);
    Wxend = max(Wx_cw);
    Wxlen = Wxend-Wx1;
    Wylen = max(Wy_cw)-min(Wy_cw);
% establish the dimensions of the observation window that will be needed for e.g. calculating its area or the isotropised set covariance
    a = Wylen; % height of rectangle
    b = Wxlen; % width of rectangle
elseif strcmp(Wshape,'parallelogram')
    Wx1 = min(Wx_cw);
    Wxend = max(Wx_cw);
    Wxlen = Wxend-Wx1;
    Wylen = max(Wy_cw)-min(Wy_cw);
% establish the dimensions of the observation window that will be needed for e.g. calculating its area or the isotropised set covariance
    a = Wy_cw(4)-Wy_cw(1); % height of parallelogram
    b = Wx_cw(2)-Wx_cw(1); % length of base of parallelogram
    c = hypot(Wx_cw(4)-Wx_cw(1),Wy_cw(4)-Wy_cw(1));% length of diagonal edge of parallelogram
    edge12 = [Wx_cw(2)-Wx_cw(1),Wy_cw(2)-Wy_cw(1)];
    edge23 = [Wx_cw(3)-Wx_cw(2),Wy_cw(3)-Wy_cw(2)];
    extanglevertex2 = acos(dot(edge12,edge23)/(norm(edge12)*norm(edge23)));
    if extanglevertex2>pi/2
        extanglevertex2_acute = pi-extanglevertex2;
    else
        extanglevertex2_acute = extanglevertex2;
    end
    extanglevertex2_acute_deg = extanglevertex2_acute*180/pi;
end

if any(rvec>min(Wxlen,Wylen)/2 & ~ignoresmallr)
    warning('the pcf estimator is known to be unreliable for r greater than half the shortest side of W: choose smaller r!');
end

if bw_r==0
    error('Cannot have zero bandwidth for the kernel function');
end

for rind=1:numrs
    fprintf('%i...',rind);
    r = rvec(rind);
    n = histc(data(:,1),unique(SPP(:,1)));
    sumn = sum(n);

    dists = pdist(data(:,3:4));
    distsarray = squareform(dists);
    
    switch kerneltype
        case 'box'
            % use a 1d uniform (i.e. box) kernel with standard deviation bw
            kernelarray = (1/(sqrt(12)*bw_r)).*double(abs(r - distsarray)./(sqrt(12)*bw_r) <= 0.5);
        case 'Epan'
            % or use the Epanechnikov kernel with std dev bw
            kernelarray = (0.75/(sqrt(5)*bw_r)).*(1-((r - distsarray)./(sqrt(5)*bw_r)).^2); 
            kernelarray(kernelarray<0) = 0;
        otherwise
            error('choose a valid kernel type/n');
    end
    kernelarray(logical(eye(size(kernelarray))))=0;
%         % force the leading diagonal to be logical(0); we do not consider the interaction of a point with itself

    scaledkernelarray = kernelarray./(2*pi*r); % scale kernel (could also use distsarray in place of r here)

    % carry out edge correction via the translation method; divide the scaled kernel by the area of the overlap
    % created by translating the observation window by the separation vector
    xdistsarray = squareform(pdist(data(:,3)));
    ydistsarray = squareform(pdist(data(:,4)));
    if strcmp(Wshape,'rectangle')
        % When W is a rectangle, Wxlen and Wylen are the lengths of W's edges; 
        % the edge correction factor is the overlap of W with its translated self
        edgecorrection = (Wxlen-xdistsarray).*(Wylen-ydistsarray);
    elseif strcmp(Wshape,'parallelogram')
        % When W is a parallelogram formed by isotropizing the original data and observation window,
        % a is the height of the parallelogram, and b is the length of the base. 
        % The edge correction factor is again the overlap of W with its translated self, 
        % which we obtain straightforwardly using a and b
        edgecorrection = (b-xdistsarray).*(a-ydistsarray);
%         edgecorrection = edgecorrection./zeta;
    end
    correctedscaledkernel = scaledkernelarray./edgecorrection;
    
    % in order to find the sums of the various submatrices of the edge-corrected,
    % scaled kernel, we use a summed area table (the matrix analogue of a cumulative sum)
    g_num = sat(correctedscaledkernel,n); % this is the numerator of the estimator
                                               % (i.e. before standardising with the intensities)
            
%% estimation of the surface-adapted intensity (Illian et al, 2008, p194)
    if strcmp(Wshape,'rectangle')
        [ec1array,~] = edgecorr_disc(data,Wx_cw,Wy_cw,r);
    elseif strcmp(Wshape,'parallelogram')
        [ec1array,~] = edgecorr_disc_prlgm(data,Wx_cw,Wy_cw,r);
    end
    boundoverlapcell = mat2cell(ec1array,n);
    boundoverlappropnarray = cellfun(@(M1) sum(M1),boundoverlapcell);
    
    % isotropised set covariance
    if strcmp(Wshape,'rectangle')
        if r<=a % a and b should be calculated near the top of this function; these are the orthogonal dimensions of W
            isotsetcov(rind) = a*b - 2*r*(a+b)/pi + (r^2)/pi;
        else
            isotsetcov(rind) = a*b*(2*asin(a/r) - a/b - 2*(r/a-sqrt((r/a)^2-1)))/pi;
        end
    elseif strcmp(Wshape,'parallelogram')
        % if W is a parallelogram, there will be a piecewise equation for calculating the isotropised set covariance, however, assuming that the
        % parallelogram-shaped W is generated by rotating and stretching by some factor, it is perhaps more straightforward to calculate the
        % isotropised set covariance for the original rectangular window (ie with matching edges) and subsequently scale using the stretch factor
        % - this will be the anisotropy ratio for us.
        if nargin>5
            % If we're here, then we have access to the original rectilinear W, along with the required scaling factor zeta.
            a_new = min(max(origWx)-min(origWx),max(origWy)-min(origWy));
            b_new = max(max(origWx)-min(origWx),max(origWy)-min(origWy));
            if r<=a_new
                isotsetcov(rind) = a_new*b_new - 2*r*(a_new+b_new)/pi + (r^2)/pi;
            else
                isotsetcov(rind) = a_new*b_new*(2*asin(a_new/r) - a_new/b_new - 2*(r/a_new-sqrt((r/a_new)^2-1)))/pi;
            end
            isotsetcov(rind) = isotsetcov(rind)/zeta;
        else
            error('parallelogram-shaped W input into function...find a way of calculating the isotropised set covariance');
        end
    end
    Gadaptedintensityarray = boundoverlappropnarray./isotsetcov(rind);
    Gsthats_adapint(:,:,rind) = g_num./(nansum(Gadaptedintensityarray,2)*nansum(Gadaptedintensityarray,2)');
    
    diagels = logical(eye(P,P));
    WA = a*b;
    n = histc(data(:,1),unique(SPP(:,1)));
    npmat = repmat(n,1,P);
    nqmat = repmat(n',P,1);
    classicintensity = npmat.*nqmat./(WA^2);
    classicintensity(diagels) = classicintensity(diagels).*(npmat(diagels)-1)./npmat(diagels);
    
    Gsthats_classicint(:,:,rind) = g_num./classicintensity;
        
%         Ghats_kernestint(:,:,rind) = sum(g_kei_array,3);
end
fprintf('\n');

varargout{1} = Gsthats_classicint;
varargout{2} = Gsthats_adapint;
% varargout{3} = Ghats_kernestint;
% varargout{4} = [min(data(:,3)),max(data(:,3)),min(data(:,4)),max(data(:,4))];
end
