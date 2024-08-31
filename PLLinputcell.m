function PLLcell = PLLinputcell(data,W,R,varargin)
% This function calculates a cell of inputs for palmloglik_genW.m. This is
% useful for maximising the PLL, as it allows the user to calculate the
% difference point pattern outside of each evaluation of the PLL, reducing
% computational budget.
% 
% If anisotropy parameters are provided, the input data is assumed to be
% anisotropic, and is isotropised. The data in the output cell is therefore
% isotropic.
% 
% INPUTS
% data      A matrix, containing a (potentially geometric anisotropic) LGCP
% W         The original, rectangular observation window
% R         the max separation distance at which we consider point pairs
% varargin  may contain angle and ratio variables, which are assumed to 
%           describe the anisotropy of the data:
%             varargin{1} - the angle of anisotropy (in radians, in [0,pi])
%             varargin{2} - the ratio of anisotropy (in (0,1))
% 
% OUTPUTS
% PLLcell   A cell containing three elements:
%               PLLcell{1} - the original Matern LGCP
%               PLLcell{2} - a matrix containing the subset of the data
%                            that lies within the inner region (calculated
%                            using the provided R.
%               PLLcell{3} - a vector containing the absolute values of the
%                            distances between all point pairs (i.e. |x-y|
%                            and |y-x|) 

% last modified by jsmartin.stats@gmail.com in Oct 2017

%%

    % Separate the x- and y-coordinates of the observation window and find the window area, before isotropization.
    Wx = [W(1),W(2),W(2),W(1)];
    Wy = [W(3),W(3),W(4),W(4)];
    
    if size(varargin)>0
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
    
    numdata=size(data,1);
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
    numinnerfocalpts = size(innerdata,1);
    % define the point pairs x-y, with x in the inner region and y in W
    diffPP1 = repmat(data(:,3)',numinnerfocalpts,1)-repmat(innerdata(:,3),1,numdata);
    diffPP2 = repmat(data(:,4)',numinnerfocalpts,1)-repmat(innerdata(:,4),1,numdata);
    diffPPr = hypot(diffPP1(:),diffPP2(:));
    % restrict attention to those point pair vectors with magnitude in (0,R)
    restrictedinds = diffPPr>0 & diffPPr<R;
    diffPPr = diffPPr(restrictedinds);
    diffPPr = sort(diffPPr);
    PLLcell = {data,innerdata,diffPPr};
end