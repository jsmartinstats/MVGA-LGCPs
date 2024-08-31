function varargout = MVGALGCPsimulate(n,W,params)
% Function for simulating log-Gaussian Cox processes
%   -   atm, we use the Matern covariance function; to be extended

% INPUT
% n         the number of datasets to be simulated
% W         contains the limits of the simulation window:
%               W = [lower x limit, upper x limit, lower y limit, upper y limit]
% params    a cell containing the mean, covariance and anisotropy parameters
%               Mean parameters:
%                   params{1} contains the expected number of points
%               Covariance parameters:
%                   params{2} contains a matrix of scale parameters (alpha)
%                   params{3} contains a matrix of smoothness parameters (nu)
%                   params{4} contains a matrix of variance parameters (sigma)
%               Anisotropy parameters:
%                   params{5} contains a matrix of anisotropy angle parameters (theta)
%                   params{6} contains a matrix of anisotropy ratio parameters (zeta)


% OUTPUT
% 
%     varargout{1}  simdata,   a cell containing n simulated MVGA LGCPs with the specified parameters;
%     varargout{2}  will contain the number of iterations that the simulation was terminated for computational concerns.
%

% last modified by jsmartin.stats@gmail.com in Nov 2017
%%
    % extract the parameter matrices
    expabundance = params{1};
    alphamat = params{2};
    numat = params{3};
    sigmamat = params{4};
    betamat = ones(size(alphamat)); 
                 % the simulation procedure below requires a scale parameter for the geometric anisotropy;
                 % this is redundant given the scale parameter in the covariance, so is set to 1 wlog
    thetamat = params{5};
    zetamat  = params{6};
    
    % establish the dimension of the LGCPs to be created
    P = size(alphamat,1);
    
    % establish remaining simulation settings
    Mx = 1024; % the resolution of the GRFs to be simulated
    My = 1024; % the resolution of the GRFs to be simulated
    covtype = 'Matern_HW94';

    % carry out simulation
    muscalar = log(expabundance/((W(2)-W(1))*(W(4)-W(3))))-diag(sigmamat)'/2;

    simdata = cell(n,1);
    index = 0;
    nanindex = 0;
    fprintf('\nGenerating process');
    while index < n
        index = index+1;
        fprintf('...%i',index);

%         % simulate a (possibly anisotropic) Gaussian Random Field with the specified covariance structure
        [~,GRF] = evalc('GRFcreate(W,[Mx,My],covtype,alphamat,numat,sigmamat,betamat,thetamat,zetamat);');
        
        % create the intensity using an expected abundance and the simulated GRF
        muvec = muscalar.*ones(1,P);
        inhomlambda = exp(reshape(repmat(muvec,Mx*My,1),Mx,My,P) + GRF);

        Aforeachlatticecell = ((W(2)-W(1))*(W(4)-W(3)))/(Mx*My);
        nineachlatticecell = poissrnd(inhomlambda.*Aforeachlatticecell);

        if sum(nineachlatticecell(:))>2000 % for the sake of memory requirements/computational budget
            index = index-1;
            nanindex = nanindex+1;
        else
            ptsineachlatticecell = rand(sum(nineachlatticecell(:)),2);
            [~,lattcelllininds]=histc(1:sum(nineachlatticecell(:)),cumsum(nineachlatticecell(:))+1);
            lattcelllininds = lattcelllininds'+1;
            yind = rem(lattcelllininds,My);
            yind(yind==0)=My;

            xind = mod(ceil(lattcelllininds./My),Mx);
            xind(xind==0)=Mx;

            Coxpts = [W(1)+(W(2)-W(1)).*(xind+ptsineachlatticecell(:,1)-1)./Mx,W(3)+(W(4)-W(3)).*(yind+ptsineachlatticecell(:,2)-1)./My];
            totaln = squeeze(sum(sum(nineachlatticecell,1),2));
            sumtotaln = sum(totaln);
            specnos = nan(sumtotaln,1);
            for p=1:P
                specnos(sum(totaln(1:p-1))+1:sum(totaln(1:p))) = p;
            end
            simdata{index} = [specnos,nan(sumtotaln,1),Coxpts];
        end
    end
    
    varargout{1} = simdata;
    varargout{2} = nanindex;
end