function [MCaxisdiff,oppositefitcount] = ellipseMCCI(MClength,numpoints,rotmat,scalemat,sh)

% This function generates a Monte Carlo sample of the difference in
% semi-axis lengths for the estimated ellipse specified by Ah and sh, 
% where 
%       Ah = rotmat*diag(scalemat)*rotmat'

% We proceed by using the transformation matrix Ah and sh to sample
% MClength different sets of points (each containing numpoints points).
% Each point is a noisy realisation 
%           e+epsilon       epsilon~mvnorm(0,sh*eye(2))
% of the ellipse described by 
%           e'*Ah*e = 1

% INPUT
% MClength      is the number of MC samples to use
% numpoints     is the number of points in each MC-sampled noisy ellipse
% Ah            is a psd matrix that specifies an ellipse through the above
% sh            is the standard deviation for both components of the 
%               bivariate noise.
%
% OUTPUT
% MCaxisdiff	is the MC sample of axis length differences
% oppositefitcount counts the number of MC iterations for which the fitted ellipse
%				has an angle of rotation opposite to that implied by the input rotmat
%
% last modified by jsmartin.stats@gmail.com in Apr 2017
%%
    lmax = 1;
    MC_ls = 1;
    MCaxisdiff = nan(MClength,lmax);
    
    % we wish to generate MClength MC samples of length numpoints from the
    % quadratic form. Since all MC samples are identically distributed, we 
    % can create the entire sample of length MClength*numpoints at once, 
    % and split up the resulting sample afterwords
    
    % generate MClength*numpoints locations on the unit circle;
    sampledangles = 2*pi*rand(1,MClength*numpoints);
    sampledcircpts = [cos(sampledangles);sin(sampledangles)];
    % now transform these, such that they are random points on the l-th fitted ellipse
    sampledellpts = reshape(rotmat*scalemat*sampledcircpts,2,numpoints,MClength);

    % now add Gaussian errors to both coordinates, with shared sd given by the l-th component of estnoisesds
    noisyellpts = sampledellpts + reshape(mvnrnd(zeros(1,2),(sh^2)*eye(2),MClength*numpoints)',2,numpoints,MClength); %#ok<NASGU>
    
    % The order in which eig calculates the evals/evecs for each noisy MC ellipse is not necessarily the same; this 
    % would affect the estimated angle, and the sign of the difference between the estimated semi-axis lengths.
    % To guard against this, first calculate the angle of rotation for the ellipse represented by the input transformation.
    estangrot = atan(rotmat(2,1)/rotmat(1,1));  % returns a value in [-pi/2,pi/2]...this is all we care about; we don't care
                                                % about whether the rotation is by theta or theta+pi, as in both these scenarios, 
                                                % the semi-axes are the same.
    % the angle of rotation (restricted to [-pi/2,pi/2]) should be the same sign in each MC iteration, as it is here
    % count the number of times that it is not...
    oppositefitcount=0;
    
    % now fit ellipses for each MC sample:
    count=0;
    for MCiter=1:MClength
        [~,A_MC,dummy1,dummy2,projind] = evalc('ellipsefit_als(noisyellpts(:,:,MCiter));');
        while projind % if the ellipse is fitted by projecting onto a set of psd matrices, one of the evals will be 0; we don't want this, we want 
                      % our evals to be positive because we're looking for these to describe the axes of the fitted ellipse.
                      % ...so we redo this MC iteration
            noisyellpts(:,:,MCiter) = sampledellpts(:,:,MCiter) + reshape(mvnrnd(zeros(1,2),(sh^2)*eye(2),numpoints)',2,numpoints);
            [~,A_MC,dummy1,dummy2,projind] = evalc('ellipsefit_als(noisyellpts(:,:,MCiter));');
%             fprintf('redoing MC iteration %i...\n',MCiter);
            % If the number of times we have to redo the fitting is greater than 50% of the size of the MC sample, raise a warning
            count = count+1;
            if count==ceil(0.5*MClength)
                warning('We are having to re-fit ellipses in more than 50 percent of the MC iterations');
            end
        end
        % A = QLQ^T, with Q the evec matrix and L the eval matrix.
        [evecs_MC,evals_MC]=eig(A_MC);
        if det(evecs_MC)<0 % if the evecs matrix represents an improper rotation, make it proper
            evecs_MC = evecs_MC(:,[2,1]);
            evals_MC = diag(flip(diag(evals_MC)));
        end
        % this proper rotation matrix will have rotation angle given by...
        estangrot_MC = atan(evecs_MC(2,1)/evecs_MC(1,1)); % (...or estangrot_MC = atan(evecs_MC(2,1)/evecs_MC(1,1))+pi, if evecs(2,1)<0 ... but this is ok! )
        if estangrot_MC*estangrot<0
            oppositefitcount=oppositefitcount+1;
        end
        
        % estimate scale of ellipticity of the MC-sampled ellipse
        % NB: we do not need to use back rotation/pcf approach here as
        % we are not estimating the compression of the original PP,
        % just the ellipticity of the l-th contour.
        MCestscale = diag(sqrt(1./diag(evals_MC)));
        MCnormestscale = MCestscale./sqrt(det(MCestscale));
        MCaxisdiff(MCiter,1) = MCnormestscale(1)-MCnormestscale(4);
        
        % %
        % look at the angles of rotation in both the input ellipse and the current MC ellipse: if they are separated by \le pi/4, treat the two
        % ellipses as being of the same true orientation; otherwise, we consider the true angles to be phase-shifted (i.e. separated by pi/2), 
        % and so the difference in semi-axes will be the negative of what has been calculated
        if abs(estangrot)<=pi/4 % if the input ellipse has an angle of rotation with abs value closer to 0 than pi/2, 
            %.... we take the angle of rotation in the current MC ellipse to be estangrot          
            if abs(estangrot-estangrot_MC)>pi/4
                % true orientations are phase-shifted, so negate the calculated axis diff
                MCaxisdiff(MCiter,1) = -MCaxisdiff(MCiter,1);
            % else
                % true orientations are similar
            end
        else 
            % ...otherwise (i.e. if the input ellipse has a.o.r. with abs val closer to pi/2 than 0
            % ... we add pi to either the input ellipse's a.o.r, the current ellipse's a.o.r., or both, 
            % ... such that they are both positive (and in the interval [pi/4,3pi/4]
            if estangrot<0
                posestangrot = estangrot+pi; % positive estimated a.o.r.
            else
                posestangrot = estangrot;
            end
            if estangrot_MC<0
                posestangrot_MC = estangrot_MC+pi; % positive estimated a.o.r.
            else
                posestangrot_MC = estangrot_MC;
            end
            if abs(posestangrot-posestangrot_MC)>pi/4
                % true orientations are phase-shifted, so negate the calculated axis diff
                MCaxisdiff(MCiter,1) = -MCaxisdiff(MCiter,1);
            % else
                % true orientations are similar
            end
        end
    end
end
