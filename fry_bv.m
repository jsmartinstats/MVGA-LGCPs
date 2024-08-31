function varargout = fry_bv(SPP,labels,varargin)
% This function outputs the "cross-Fry plot" for the spatial point pattern in SPP
% This is a plot of the points x2-x1 such that x1 is in the first component of 
% the bivariate point pattern "SPP", and x2 is in the second component
% 
% If instructed, the function will include the reflected points x1-x2

% varargin can contain a minimum and maximum separation distance

% INPUT
% SPP       contains the spatial point pattern of interest
%           This is an n-by-4 matrix containing the x- and y-coords in the third 
%               and fourth columns, resp., and the process label in the first column
% labels    the numeric labels of the components of SPP; this is given explicitly, in
%           order to efficiently deal separately with the asymmetry between p vs q and
%           q vs p.
% varargin  potentially contains:
%               a - the minimum separation distance for contributory point pairs
%               b - the maximum separation distance for contributory point pairs
%               reflect - a boolean, telling us whether to consider the reflected Fry points too
% 
% OUTPUT
% varargout contains an m-by-2 matrix with the x- and y-coords for the Fry points
% 
% last modified by jsmartin.stats@gmail.com in Feb 2017

SPP1 = SPP(SPP(:,1)==labels(1),3:4);
SPP2 = SPP(SPP(:,1)==labels(2),3:4);

    if nargin==5
        a = varargin{1};
        b = varargin{2};
        reflect = varargin{3};
    elseif nargin==4
        a = varargin{1};
        b = varargin{2};
        reflect = false
    elseif nargin==3
        a = 0;
        b = varargin{1};
        reflect = false;
    elseif nargin==2
        a = 0;
        b = realmax;
        reflect = false;
    end
        

    n1 = size(SPP1,1);
    n2 = size(SPP2,1);
    
    point1 = nan(n1,n2,2);
    point2 = nan(n1,n2,2);
    point1(:,:,1) = repmat(SPP1(:,1),1,n2);
    point1(:,:,2) = repmat(SPP1(:,2),1,n2);
    point2(:,:,1) = repmat(SPP2(:,1)',n1,1);
    point2(:,:,2) = repmat(SPP2(:,2)',n1,1);

    % % calculate the vector from point 1 to point 2
    xsepdist = point2(:,:,1)-point1(:,:,1);
    ysepdist = point2(:,:,2)-point1(:,:,2);
    
    % % restrict according to max and min separation distances
    restrictedinds = hypot(xsepdist,ysepdist)>=a & hypot(xsepdist,ysepdist)<=b;
    xsepdist = xsepdist(restrictedinds);
    ysepdist = ysepdist(restrictedinds);
    xsepdist = xsepdist(:);
    ysepdist = ysepdist(:);
    if reflect
        xsepdist = [xsepdist;-xsepdist];
        ysepdist = [ysepdist;-ysepdist];
    end

    varargout{1} = [xsepdist,ysepdist];
end