function [esttheta,estzeta] = PP_GAest(data,W)
% This function estimates the angle and scale of anisotropy displayed by the potentially  
% multivariate SPP in data. The approach is described in the manuscript. The (symmetric) Fry 
% points are calculated for each dimension and each pair of dimensions. Contours are calculated, 
% where the l-th contour contains the l-th closest point in a number of constructed directions.
% For anisotropic processes, these contours will be elliptical, with assumed Gaussian noise;
% an ellipse is therefore fitted using adjusted least squares, and the deformation of this
% ellipse is calculated. 

% INPUT: 
% data      potentially multivariate spatial point pattern data
% W 		the observation window for the data, in the form
% 				[lower x bound, upper x bound, lower y bound, upper y bound]
% 
% OUTPUT
% 
% esttheta is a matrix of estimated angles of anisotropy
% estzeta is a matrix of estimated scales of anisotropy 
% 
% last modified by jsmartin.stats@gmail.com in May 2018
%%


minr = 0.0;
maxr = 0.25; % maximum separation distance considered when calculating Fry points

P = max(data(:,1));

esttheta = nan(P,P);
estzeta = nan(P,P);

for p=1:P
    for q=1:p

        if p==q
			numsec = floor(sum(data(:,1)==p)/6); % number of sectors/directions to test; this is a resolution parameter - the fitted ellipses each use numsec points
												% this is chosen using RSRS16's rule of thumb: numsec~abundance/6
            % compute Fry points
            frypts = fry(data(data(:,1)==p,3:4),minr,maxr);
        else
			% numsec = floor(min(sum(data(:,1)==p)/6,sum(data(:,1)==q)/6)); 
			WA = (W(2)-W(1))*(W(4)-W(3));
			lamp = sum(data(:,1)==p)/WA;
			lamq = sum(data(:,1)==q)/WA;
			numsec = floor((lamp*lamq*WA)/(6*(lamp+lamq)));
												% number of sectors/directions to test; this is a resolution parameter - the fitted ellipses each use numsec points
												% this is chosen using an extension of RSRS16's rule of thumb: we use numsec ~ lam1*lam2*|W|/6(lam1+lam2)
            % compute bivariate Fry points
            % use reflected Fry points as well - ie x2-x1 as well as x1-x2
            frypts = fry_bv(data(ismember(data(:,1),[p,q]),:),[p,q],minr,maxr,true);
        end
%%
        % represent Fry points in polar coords
        fryptsr = hypot(frypts(:,1),frypts(:,2));
        fryptsth = atan2(frypts(:,2),frypts(:,1));
            % convert from [-pi,pi] to [0,2pi]
            fryptsth(fryptsth<0) = fryptsth(fryptsth<0) + 2*pi;
        % assign each Fry point to a sector (with width 2pi/numsec)
%         [~,~,secind] = histcounts(fryptsth,0:2*pi/numsec:2*pi); % v2016 (I think) onwards
        [counts,secind] = histc(fryptsth,0:2*pi/numsec:2*pi);

        % create matrices by repeating col vector (think about what information we're wanting to keep)
        fryrmat = repmat(fryptsr,1,numsec);
        fryxymat = permute(reshape(repmat(frypts,1,numsec),numel(fryptsr),2,numsec),[1,3,2]);
        % NA-out values that dont belong in each col
        fryrmat(repmat(secind,1,numsec)~=repmat(1:numsec,numel(fryptsr),1)) = nan;
        fryxymat(reshape(repmat(secind,1,2*numsec),size(fryxymat,1),numsec,2)~=reshape(repmat([1:numsec,1:numsec],size(fryxymat,1),1),size(fryxymat,1),numsec,2))=nan;
        % sort each col independently, in ascending order of r
        [fryrmat,sortinds] = sort(fryrmat);
        for secind=1:numsec
            fryxymat(:,secind,:) = fryxymat(sortinds(:,secind),secind,:);
        end
        % trim matrix
        fryrmat = fryrmat(1:max(sum(~isnan(fryrmat))),:);
        fryxymat = fryxymat(1:max(sum(~isnan(fryrmat))),:,:);

        % fit ellipses to the first lmax pseudo-Fry point sets, using adjusted least squares with unknown noise variance
        lmax = min(counts(1:numsec));
        estnoisesds = nan(1,lmax);
        estrotmats = nan(2,2,lmax);
        estscalemats = nan(2,2,lmax);
        invalid_ls = [];
        
		% fit an ellipse to the lmax-th contour first; use this to fix the specification of the ellipse (i.e. st zeta\in(0,1])
		[~,A,~,estnoisesds(lmax),notpd] = evalc('ellipsefit_als(reshape(fryxymat(lmax,:,:),numsec,2)'');');
		[evecs,evals]=eig(A);
		% if the evec matrix corresponds to an improper rotation, find the corresponding 'proper' eigendecomposition 
		if det(evecs)<0
			evecs = evecs(:,[2,1]);
			evals = diag(flip(diag(evals)));
		end
		estrotmats(:,:,lmax) = evecs;
		estscalemats(:,:,lmax) = diag(sqrt(1./diag(evals)));
			% set the direction of rotation for this ellipse as a benchmark against which all ellipses are calibrated for this point pattern
			benchangle = mod(atan(estrotmats(2,1,lmax)/estrotmats(1,1,lmax))+pi/2+pi/2,pi)-pi/2; 
					% introduce phaseshift by setting thetahat_lmax=atan(R_21/R_11)+pi/2
					% use mod(thetahat_lmax+pi/2,pi)-pi/2 to ensure final benchangle is in [-pi/2,pi/2) (for consistency, as this is what is returned by atan)
			% check whether the implied zeta is in (0,1]...i.e. the fitted ellipse is compressed top-down before rotating
			impliedzeta = estscalemats(2,2,lmax)/estscalemats(1,1,lmax);
			if impliedzeta>1 % if not, switch the specification of the ellipse:
				% overwrite the proper rotation matrix using the benchmark angle
				estrotmats(:,:,lmax) = [cos(benchangle),-sin(benchangle);sin(benchangle),cos(benchangle)];
				% overwrite the scale matrix to its new form
				estscalemats(:,:,lmax) = diag(flip(diag(estscalemats(:,:,lmax))));
			end
		if notpd
			invalid_ls = [invalid_ls,lmax]; %#ok<AGROW>
		end
        
        for l=1:(lmax-1)
            % fit the ellipse corresponding to each l-contour
             [~,A,~,estnoisesds(l),notpd] = evalc('ellipsefit_als(reshape(fryxymat(l,:,:),numsec,2)'');');
            [evecs,evals]=eig(A);
            % if the evec matrix corresponds to an improper rotation, find the corresponding 'proper' eigendecomposition 
            if det(evecs)<0
                evecs = evecs(:,[2,1]);
                evals = diag(flip(diag(evals)));
            end
            estrotmats(:,:,l) = evecs;
            estscalemats(:,:,l) = diag(sqrt(1./diag(evals)));
                % check whether the angle of rotation for this fitted ellipse is in the same direction as the benchmark; i.e. ensure that the abs
                % diff between the angles is \le pi/4
                if abs(benchangle)<pi/4
                % if the abs value of the benchmark angle itself is less than pi/4, we can proceed straightforwardly:
                    if abs(benchangle-atan(estrotmats(2,1,l)/estrotmats(1,1,l))) > pi/4
                        % if the abs diff between angles is >pi/4, change the specification of the current fitted ellipse
                        % phase-shift the angle of rotation
                        shiftedangle = mod(atan(estrotmats(2,1,l)/estrotmats(1,1,l))+pi/2+pi/2,pi)-pi/2; % see above for explanation of syntax
                        % overwrite the proper rotation matrix using the shifted angle
                        estrotmats(:,:,l) = [cos(shiftedangle),-sin(shiftedangle);sin(shiftedangle),cos(shiftedangle)];
                        % overwrite the scale matrix to its new form
                        estscalemats(:,:,l) = diag(flip(diag(estscalemats(:,:,l))));
                    end
                else
                % if the abs value of the benchmark angle is closer to pi/2 than 0, we add pi to either the benchmark, the current angle, or both,
                % such that they are both positive (and in the interval [pi/4,3pi/4]
                    if benchangle<0
                        pba = benchangle+pi;% pba = 'positive benchangle'
                    else
                        pba = benchangle;
                    end
                    if atan(estrotmats(2,1,l)/estrotmats(1,1,l))<0
                        pca = atan(estrotmats(2,1,l)/estrotmats(1,1,l)) + pi; % pca = 'positive current angle'
                    else
                        pca = atan(estrotmats(2,1,l)/estrotmats(1,1,l));
                    end
                    if abs(pba-pca) > pi/4
                        % if the abs diff between angles is >pi/4, change the specification of the current fitted ellipse
                        % phase-shift the angle of rotation
                        shiftedangle = mod(atan(estrotmats(2,1,l)/estrotmats(1,1,l))+pi/2+pi/2,pi)-pi/2; % see above for explanation of syntax
                        % overwrite the proper rotation matrix using the shifted angle
                        estrotmats(:,:,l) = [cos(shiftedangle),-sin(shiftedangle);sin(shiftedangle),cos(shiftedangle)];
                        % overwrite the scale matrix to its new form
                        estscalemats(:,:,l) = diag(flip(diag(estscalemats(:,:,l))));
                    end
                end                    
            if notpd
                invalid_ls = [invalid_ls,l]; %#ok<AGROW>
            end
        end
        valid_ls = setdiff(1:lmax,invalid_ls);
% maybe further restrict fitted ellipses to those with low noise sd?
% or figure out what is causing non-psd transformation matrices to be
% fitted in ellipsefit_als
        
%%    
    % for the {10:10:lmax}th fitted ellipses, find the Monte Carlo confidence intervals for the semi-axis differences:
        MClength=200;
        numpoints=50;
    % create vector of original l-indices for which we find the semi-axis difference MC CIs 
        if lmax<10
            MC_ls = 1:lmax;
        elseif lmax<20
            MC_ls = [1,2:2:lmax];
        elseif lmax<50
            MC_ls = [1,5:5:lmax]; 
        elseif lmax<200
            MC_ls = [1,10:10:lmax]; 
        else
            MC_ls = [1,100:100:lmax];
        end
        validMC_ls = MC_ls(ismember(MC_ls,valid_ls)); % subset of MC l-indices for which the corresponding contour has a pd fitted ellipse 
        axisdiff = nan(MClength,numel(validMC_ls));
        ofc = nan(numel(1,validMC_ls));
        for l_ind=1:numel(validMC_ls)
            l=validMC_ls(l_ind);
            disp([l,lmax]);
            [axisdiff(:,l_ind),ofc(l_ind)] = ellipseMCCI(MClength,numpoints,estrotmats(:,:,l),estscalemats(:,:,l),estnoisesds(l));
        end
%%   
        axisdiff = sort(axisdiff);
        axisdiff95pcMCCI_LB = axisdiff(floor(0.025*MClength),:);
        axisdiff95pcMCCI_UB = axisdiff(ceil(0.975*MClength),:);

        % We expect some of the 95% MC confidence intervals to envelope zero: we expect the l-contours to be roughly isotropic for low l. We therefore
        % identify these isotropic contours using the MCCIs. When the first MCCI indicating isotropy is known, we identify the l-contour corresponding
        % to the first subsequent CI that doesn't envelope 0, and we use all l-contours with l greater than this value. This is based on TARs 
        % observation that, for clustered processes, the fitted ellipses overestimate the 'roundness' of the transformation at smaller scales. 
        MCCInotenv0 = find(axisdiff95pcMCCI_LB>0 | axisdiff95pcMCCI_UB<0);% find which MCCIs do not envelope zero, i.e. indicating anisotropy
        firstMCCIenv0 = find(axisdiff95pcMCCI_LB<0 & axisdiff95pcMCCI_UB>0,1);% find the first MCCI that envelopes zero, i.e. indicating isotropy
        if numel(firstMCCIenv0)>0
            firstsubseqMCCInotenv0 = MCCInotenv0(find(MCCInotenv0>firstMCCIenv0,1));% find the index of the first MCCI to indicate anisotropy after those indicating isotropy
            firstanisotropic_valid_MC_l_ind = firstsubseqMCCInotenv0;
        else
            firstanisotropic_valid_MC_l_ind = 1;
        end
        if numel(firstanisotropic_valid_MC_l_ind)==0
            % if the above hasn't worked, just use all l-contours with l greater than that corresponding to the first CI that doesn't envelope 0.
            firstanisotropic_valid_MC_l_ind = MCCInotenv0(1);
        end
        
        useful_valid_ls = valid_ls(valid_ls>=validMC_ls(firstanisotropic_valid_MC_l_ind));
        % Now, the 'useful, valid l-contours' are those that have pos-def fitted ellipses and display a significant degree of anisotropy. We use these to
        % find the overall rotation estimator.

        % We sample 10000 points on the surfaces of the ellipses fitted to the useful, valid contours, using the reciprocal of the estimated noise standard
        % deviations as sampling weights (...noting that each fitted ellipse has already been normalised to have unit volume).
        sampweights = 1./estnoisesds; 
        sampweights(useful_valid_ls) = sampweights(useful_valid_ls)./sum(sampweights(useful_valid_ls));
        sampweights(setdiff(1:numel(sampweights),useful_valid_ls)) = 0;

        sampsizes = round(10000*sampweights); % sumofsamplesizes = sum(sampsizes)
        allsampledcircpts = [];
        allsampledellpts = [];
        for l=useful_valid_ls
            sampledangles = 2*pi*rand(1,sampsizes(l));
            sampledcircpts = [cos(sampledangles);sin(sampledangles)];
            sampledellpts = estrotmats(:,:,l)*(estscalemats(:,:,l)./sqrt(det(estscalemats(:,:,l))))*sampledcircpts;
            allsampledcircpts = [allsampledcircpts,sampledcircpts]; %#ok<AGROW>
            allsampledellpts = [allsampledellpts,sampledellpts]; %#ok<AGROW>
        end

        % fit an ellipse to these 10000 sampled points
        [~,A,~,~,notpd] = evalc('ellipsefit_als(allsampledellpts)');
        % the above A may have been given by projecting onto the space of psd matrices, i.e. it may contain an eigenvalue that is zero. Since the
        % quadratic form of an ellipse requires A to be positive definite, this is therefore unsuitable for further use. Degeneracy may have occurred
        % due to different scales of anisotropy being present in the data. We therefore repeat the above calculations, reducing the number of 
        % l-contours being used to fit the final ellipse; we sequentially remove the lowest-indexed l-contours.
        while notpd
            useful_valid_ls = useful_valid_ls(2:end);
            sampweights = 1./estnoisesds; 
            sampweights(useful_valid_ls) = sampweights(useful_valid_ls)./sum(sampweights(useful_valid_ls));
            sampweights(setdiff(1:numel(sampweights),useful_valid_ls)) = 0;
            sampsizes = round(10000*sampweights); % sumofsamplesizes = sum(sampsizes)
            allsampledcircpts = [];
            allsampledellpts = [];
            for l=useful_valid_ls
                sampledangles = 2*pi*rand(1,sampsizes(l));
                sampledcircpts = [cos(sampledangles);sin(sampledangles)];
                sampledellpts = estrotmats(:,:,l)*(estscalemats(:,:,l)./sqrt(det(estscalemats(:,:,l))))*sampledcircpts;
                allsampledcircpts = [allsampledcircpts,sampledcircpts]; %#ok<AGROW>
                allsampledellpts = [allsampledellpts,sampledellpts]; %#ok<AGROW>
            end
            [~,A,~,~,notpd] = evalc('ellipsefit_als(allsampledellpts)');
        end
       
        [evecs,evals]=eig(A);
		% if the evec matrix corresponds to an improper rotation, find the corresponding 'proper' eigendecomposition 
		if det(evecs)<0
			evecs = evecs(:,[2,1]);
			evals = diag(flip(diag(evals)));
		end
        estrotmat_MC = evecs;
        estscalemat_MC = diag(sqrt(1./diag(abs(evals))));
% %
		% check whether the angle of rotation for this fitted ellipse is in the same direction as the benchmark; i.e. ensure that the abs
		% diff between the angles is \le pi/4
		if abs(benchangle)<pi/4
			% if the abs value of the benchmark angle itself is less than pi/4, we can proceed straightforwardly:
			if abs(benchangle-atan(estrotmat_MC(2,1)/estrotmat_MC(1,1))) > pi/4
				% if the abs diff between angles is >pi/4, change the specification of the current fitted ellipse
				% phase-shift the angle of rotation
				shiftedangle = mod(atan(estrotmat_MC(2,1)/estrotmat_MC(1,1))+pi/2+pi/2,pi)-pi/2; % see above for explanation of syntax
				% overwrite the proper rotation matrix using the shifted angle
				estrotmat_MC = [cos(shiftedangle),-sin(shiftedangle);sin(shiftedangle),cos(shiftedangle)];
				% overwrite the scale matrix to its new form
				estscalemat_MC = diag(flip(diag(estscalemat_MC)));
			end
		else
			% if the abs value of the benchmark angle is closer to pi/2 than 0, we add pi to either the benchmark, the current angle, or both,
			% such that they are both positive (and in the interval [pi/4,3pi/4]
			if benchangle<0
				pba = benchangle+pi;% pba = 'positive benchangle'
			else
				pba = benchangle;
			end
			if atan(estrotmat_MC(2,1)/estrotmat_MC(1,1))<0
				pca = atan(estrotmat_MC(2,1)/estrotmat_MC(1,1)) + pi; % pca = 'positive current angle'
			else
				pca = atan(estrotmat_MC(2,1)/estrotmat_MC(1,1));
			end
			if abs(pba-pca) > pi/4
				% if the abs diff between angles is >pi/4, change the specification of the current fitted ellipse
				% phase-shift the angle of rotation
				shiftedangle = mod(atan(estrotmat_MC(2,1)/estrotmat_MC(1,1))+pi/2+pi/2,pi)-pi/2; % see above for explanation of syntax
				% overwrite the proper rotation matrix using the shifted angle
				estrotmat_MC = [cos(shiftedangle),-sin(shiftedangle);sin(shiftedangle),cos(shiftedangle)];
				% overwrite the scale matrix to its new form
				estscalemat_MC = diag(flip(diag(estscalemat_MC)));
			end
		end                    

        sampledangles = 2*pi*rand(1,1000);
        sampledcircpts = [cos(sampledangles);sin(sampledangles)];
        sampledellpts = estrotmat_MC*(estscalemat_MC./sqrt(det(estscalemat_MC)))*sampledcircpts;

        estangle = atan(estrotmat_MC(2,1)/estrotmat_MC(1,1));        
        esttheta(p,q) = mod(estangle+pi/2,pi);

%%        
        kerntype = 'box'; % NB: diameter of box kernel is sqrt(12)*bw_r; bw_r is the std dev of the kernel
        bandwidths_zetaest = [0.05,pi/20];
        bandwidths_zetaest = [0.05,(pi/2)]./sqrt(12); % try bw_phi such that all points contribute to Gsthat_anisot_genW in one direction or t'other.
        bandwidths_zetaest = [0.05,(pi/8)]./sqrt(12); % No - earlier expts showed that lower bw_phi is better. Use bws such that half-angle is pi/16 (more extreme than RRSS17) and half-distance is 0.05 
        bandwidths_zetaest = [nan,(pi/8)]; % Using Khat_anisot_genW - the input required *is* the bandwidth, not the std dev.
    
        numcandzetas = 199;
        rdel = 0.01;
		rvec = rdel:rdel:(maxr-rdel);
        phivec = [0,pi/2];
    
        zetavec = linspace(0,2,numcandzetas+2); 

        isoKarray = nan(1,1,numel(rvec),numel(phivec),numel(zetavec));
        for zetaind=1:length(zetavec)
            zeta = zetavec(zetaind);
            if p==q
                [~,anisoK_ci] = evalc('Khat_anisot_genW(data(data(:,1)==p,:),W,rvec,phivec,bandwidths_zetaest(2),esttheta(p,q),zeta);');
                isoKarray(:,:,:,:,zetaind) = anisoK_ci;
            else
                [~,anisoK_ci] = evalc('Khat_anisot_genW(data(ismember(data(:,1),[p,q]),:),W,rvec,phivec,bandwidths_zetaest(2),esttheta(p,q),zeta);');
                isoKarray(:,:,:,:,zetaind) = (anisoK_ci(1,2,:,:) + anisoK_ci(2,1,:,:))./2;            
            end
        end
%%
		absdiffs = abs(isoKarray(:,:,:,1,:)-isoKarray(:,:,:,2,:));
		midpts = squeeze((absdiffs(:,:,1:end-1,:,:)+absdiffs(:,:,2:end,:,:))/2);
		estV = sum(midpts.*rdel,1);
		[~,minloc]=min(estV);
		estzeta(p,q) = zetavec(minloc);        
		
		% Here, we estimate zeta using an input rvec. Previously, we have tried estimating zeta for a range of rvecs and choosing the modal estimate; the
		% different rvecs tried had the same endpoint and increment, but started at different values. 
        
    end
end

end
    
