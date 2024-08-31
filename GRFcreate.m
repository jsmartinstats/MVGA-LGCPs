function varargout = GRFcreate(W,M,covtype,varargin)

% Returns a potentially multidimensional zero-mean Gaussian random field 
% with spatial autocovariance structure specified by covtype.

% Note that this procedure is approximate, as it uses circulant embedding
% % Is this true? It is based on Wood and Chan (1994), which gives an 
% % exact-in-principle procedure if certain conditions are satisfied...

% INPUT:
% W         A vector giving the limits for the observation window in each
%           dimension
% M         A vector giving the discretisation resolution for each dimension
% covtype   A string naming the type of spatial autocovariance structure to
%           use. 'Matern_Stein' is recommended; this refers to the Matern
%           covariance structure, with parameterisation as in Stein 1999, 
%           eq .
%
% varargin  Contains parameter matrices for the spatial covariance structure. 
%           The parameters on the diagonals determine the autocovariance
%           structure, those on the off-diagonals the cross-cov structure.
%           varargin also contains parameters for the geometric anisotropy.
%           The covariance parameters come first and depend on covtype:
%           covtype             | Parameters
%           -----------------------------------------------
%           Matern_HW94         |   alpha   (varargin{1})
%                               |   nu      (varargin{2})
%                               |   sigma   (varargin{3})
%           --------------------|---------------------------
%                               |
%           The anisotropy parameters come next:
%                                   betamat  (varargin{-nargin-2})
%                                   thetamat (varargin{-nargin-1})
%                                   zetamat  (varargin{-nargin})
% OUTPUT:
% varargout contains a potentially multidimensional Gaussian random field (GRF)

% last modified by jsmartin.stats@gmail.com in Nov 2017


% We use the n-dimensional fft below; the output of which may be required 
% to be real and/or (if complex, Hermitian) symmetric. We therefore make
% use of the following tolerance parameters:
    imagtol = 10^(-10);% 10^(-5); % tolerance for imaginary component of any fft output
    asymtol = 10^(-10);% 10^(-1); % tolerance for asymmetry of output of n-dim fft.
    zerotol = 10^(-12);% tolerance for 

% Extract the parameter vectors from varargin, depending on covtype

if or(strcmp(covtype,'Matern_HW94'),strcmp(covtype,'powexp'))
    try
        narginchk(6,9);
    catch ME
        newME = MException(ME.identifier,...
                           strcat(ME.message,'3 Covariance parameters required, 3 anisotropy parameters optional.'));
        throw(newME);
    end
    alphamat = varargin{1};
    numat = varargin{2};
    sigmasigmamat = varargin{3};
    if nargin>6
        try
            narginchk(9,9);
            betamat = varargin{4};
            thetamat = varargin{5};
            zetamat = varargin{6};
        catch
            fprintf('Only %g anisotropy parameters provided, isotropy assumed.\n',length(varargin)-3);
            betamat = ones(size(alphamat));
            thetamat = zeros(size(alphamat));
            zetamat = ones(size(alphamat));
        end
%        if the number of inputs is more than 6, then varargin contains 6 inputs,
%        the covariance parameters and anisotropy parameters (and we know
%        that the field is anisotropic).
        geoanisot = true;
    else
        geoanisot = false;
    end
end

% Find out the dimension of the GRF to be created
    P = size(alphamat,1);

% Find out the resolution of the discretisation
    Mx = M(1);
    My = M(2);

% We use the Wood & Chan (1994) procedure for generating a planar Gaussian
% random field. In the variables defined here, M and N refer, respectively,
% to n and m in that paper.
%%
    if geoanisot
        covfun = @(u1,u2) covfun_geoanisot_quick(u1,u2,covtype,sigmasigmamat,alphamat,numat,thetamat,betamat,zetamat);
    else
        covfun = @(r) covfun_isot(r,covtype,sigmasigmamat,alphamat,numat);
    end
    
% check if the covariance function is even in all coordinates:
    dummy = -0.5:0.05:0.5;
    lendum = length(dummy);
    dummyx = repmat(dummy,lendum,1);
    dummyy = repmat(dummy',1,lendum);
    if geoanisot
        dummycell = arrayfun(covfun,dummyx,dummyy,'UniformOutput',false);
    else
        dummydists = sqrt(dummyx.^2 + dummyy.^2);
        dummycell = arrayfun(covfun,dummydists,'UniformOutput',false);
    end
    posinds = 1:lendum;
    neginds = lendum:-1:1;
    eveninx=nan(lendum,1);
    eveniny=nan(lendum,1);
    for i=1:length(posinds)
        eveninx(i) = all(all(cell2mat(dummycell(:,posinds(i)))==cell2mat(dummycell(:,neginds(i)))));
        eveniny(i) = all(all(cell2mat(dummycell(posinds(i),:))==cell2mat(dummycell(neginds(i),:))));
    end

% The Wood & Chan (1992) method generates a GRF over a discrete
% lattice that is at least twice as large (in every dimension) as the
% one required; this is so that the spatial covariance structure can be
% imbued using circulant embedding.
% I define here the possible resolutions of the extended GRF; it may be
% the case that for a particular resolution we calculate a circulant 
% embedding matrix that is not nonnegative definite. In such a
% scenario, we would move to a higher resolution. 
% To increase the efficiency of the fft, we define the resolution, Nx,
% such that it is a power of 2 (if the covariance fun is even in x) or 
% a power of 3 (otherwise)...
    if all(eveninx)
        possNxs = 2.^(1:20);
    else
        possNxs = 3.^(1:20);
%         possNxs = 2.^(1:10);
    end
    if all(eveniny)
        possNys = 2.^(1:20);
    else
        possNys = 3.^(1:20);
%         possNys = 2.^(1:10);
    end
% ...and such that Nx>=2(Mx-1)
    possNxs = possNxs(possNxs>=2*(Mx-1)); 
%     possNxs = possNxs(possNxs<=2^11); % censor the possible resolution of N according to available computational resources
    possNxs = possNxs(possNxs<=3^7); % censor the possible resolution of N according to available computational resources
% similarly for Ny
    Ny = 2*(My-1);
    possNys = possNys(possNys>=Ny);
%     possNys = possNys(possNys<=2^11); % censor the possible resolution of N according to available computational resources
    possNys = possNys(possNys<=3^7); % censor the possible resolution of N according to available computational resources
    
    evalsofevals=-1;
    Nxind = 0;
    Nyind = 0;
    while any(evalsofevals(:)<0)
        if or(Nxind==numel(possNxs),Nyind==numel(possNys)) 
            %% added June 2017
            % This section should only be hit if the while condition has succeeded (i.e. if we haven't managed to 
            % obtain Nbar P-by-P eigenvalue matrices that are all nonnegative definite), for all attempted block
            % circulant embedding matrices within the specified size limits for [Nx,Ny].
            % If this is the case, we resort to approximate circulant embedding (see also Wood & Chan, 1994).
            trLam = sum(evals,3);% trace of Lambda (WC94 notn) is here a PxP matrix, where the (p,q)-element is 
                                 % the sum of evals(p,q,:); see above for a reminder of the structure of evals.
            trLamplus = sum(max(evals,0),3); % trace of Lambda_+ (WC94 notn) - this is the same as above, but 
                                             % where the negative eigenvalues are ignored.
            rho1 = trLam./trLamplus; % approx embedding matrix is now 
            evals = reshape(repmat(rho1.^2,1,Nbar),P,P,Nbar).*max(evals,0); % eigenvalues of the approximate 
                                 % block-circulant matrix. Since this is constructed from the valid (i.e. real,
                                 % symmetric) eigenvalue matrices calculated for the largest attempted block-
                                 % circulant embedding matrix, it will also be real and symmetric.

        else
            %% Find the eigenvalues for the defined block circulant embedding matrices.       
            % strategy for deciding which dimension to extend:
            %   - use pcolor to depict the matrices of the covariances in covs;
            %   - one of the directions will contain a discontinuity - choose
            %   this dimension to extend

            Nxind = Nxind+1;
            Nyind = Nyind+1;
        %     Nx = 2^(log2(Nx)+1)
            Nx=possNxs(Nxind);
        %     Ny = 2^(log2(Ny)+1);
            Ny=possNys(Nyind);

            Nbar = Nx*Ny;

            Ndaggx = floor(Nx/2)+1;
            Ndaggy = floor(Ny/2)+1;
            Ndaggbar = Ndaggx*Ndaggy;
            [pdaggy,pdaggx] = meshgrid(0:Ndaggy-1,0:Ndaggx-1);
                % Note order of outputs in line above; this is due to the bijection
                % (3.4) in Wood & Chan (1994), which defines the relation between the
                % row index s, of the Nbar-by-2 matrix p, and the multi-index j, which
                % is an element of the set I(N).
%             pdagg = [pdaggx(:),pdaggy(:)]; % here, pdagg is the matrix I(mdagg), as noted by Wood & Chan (1994)
%                                            % pdagg needs to be integer-valued, as it is used to index later on
            pdaggx=pdaggx(:);pdaggy=pdaggy(:);

            [covcol1distsy,covcol1distsx] = meshgrid(0:Ny-1,0:Nx-1);
            covcol1distsx=covcol1distsx(:);
            covcol1distsy=covcol1distsy(:);
                % Here, we have [covcol1distsx,covcol1distsy] = p, in the
                % notation of Wood & Chan (1994)
                % Note that the order of the outputs in the line above follows the same
                % convention as before.            

        % We use the 2d FFT to calculate the eigenvalues of our block
        % circulant embedding matrix C; the input for the 2d FFT is the
        % covariance function evaluated at each row of the matrix I(m) (in
        % notn of Wood & Chan, 1994) after a transformation has been applied.

        % First, we define I(m) as [px(:),py(:)] - this obeys the bijection
        % specified by Wood & Chan (1994, p414). Then, for the first column
        % of this (Nx*Ny)-by-2 matrix, each element has the following
        % transformation applied: for element h,
        %   if -Nx+1 <= h <  -Nx/2 , h := h+Nx
        %   if -Nx/2 <= h <=  Nx/2 , no transformation is applied
        %   if  Nx/2 <  h <=  Nx-1 , h := h-Nx
        % and similarly for the second column, with Ny replacing Nx
            inds = covcol1distsx>(Nx/2) & covcol1distsx<=(Nx-1);
            covcol1distsx(inds) = covcol1distsx(inds)-Nx;
            inds = covcol1distsy>(Ny/2) & covcol1distsy<=(Ny-1);
            covcol1distsy(inds) = covcol1distsy(inds)-Ny;
            inds = -covcol1distsx>(Nx/2) & -covcol1distsx<=(Nx-1);
            covcol1distsx(inds) = covcol1distsx(inds)+Nx;
            inds = -covcol1distsy>(Ny/2) & -covcol1distsy<=(Ny-1);
            covcol1distsy(inds) = covcol1distsy(inds)+Ny;

        % we also scale the covcol1dists according to the size of the
        % observation window
            covcol1distsx = covcol1distsx.*(W(2)-W(1));
            covcol1distsy = covcol1distsy.*(W(4)-W(3));

            if geoanisot
                covs = covfun_geoanisot_quick(covcol1distsx./Mx,covcol1distsy./My,covtype,sigmasigmamat,alphamat,numat,thetamat,betamat,zetamat);
            else
                disp('need to change covfun_isot so that it can take a col vector as distance input');
            end

        %%
            evals = fft2(permute(reshape(covs',P*P,Nx,Ny),[2,3,1]));
            clear covs; % the array 'covs' is huge - one of the biggest arrays created in this program.
                        % delete it to save memory

        % The above cycles through the elements of the PxP covariance 
        % matrix at each distance d(1,1),d(1,2),...,d(1,P), 
        % d(2,1),...d(2,P),...,...,d(P,1),...,d(P,P). Technically, this is
        % the wrong way; this is cycling through the transpose. Since the
        % covariance matrix is symmetric, however, the end result is the
        % same.
        % we can use fft2 here because the argument is a block
        % circulant matrix
        % NB: 'covs' is the first column of blocks in the block circulant
        % embedding matrix; evals gives the eigenvalues of the P^2 Nbar-by-Nbar
        % matrices, each of which is formed by taking a single element from
        % each block. 

            tic;
            inds = abs(imag(evals(:)))<imagtol;
            evals(inds) = real(evals(inds));

            evals = permute(reshape(evals(:),Nbar,P,P),[3,2,1]); 

            if sum(sum(sum(evals - permute(evals,[2,1,3]))))<asymtol
                evals = (evals + permute(evals,[2,1,3]))./2;
            end
            time2=toc;
        end
        % two possible approaches for checking that all eval matrices are
        % nnd: the first takes longer, but uses less memory; the second
        % (might be) quicker, but definitely takes more memory. 
        %% Option 1:
        i=0;
        while i<size(evals,3)
            i=i+1;
            % we want all slices of evals to be nnd; iterate through 3rd
            % dim of evals, calculating eigenvalues, and if we get a
            % negative eigenvalue, set evalsofevals to be negative and exit
            % this loop
            evalsofevalslice=eig(evals(:,:,i));
            if any(evalsofevalslice<-zerotol) % allow these eigenvalues to 
                                              % be negative, but v v small.
                evalsofevals=-1;
                i=size(evals,3); % forces us to exit the while loop
            else
                evalsofevals=1;
            end
        end
        %% Option 2:
%         for i=1:size(evals,3)
%             estruct(i).mat=evals(:,:,i); %#ok<AGROW>
%         end
%         [~,evalsofevals]=arrayfun(@(M1) eig(M1.mat), estruct,'UniformOutput',false);
%         evalsofevals=cell2mat(evalsofevals);
%         inds = imag(evalsofevals(:))<imagtol;
%         evalsofevals(inds) = real(evalsofevals(inds));
% 
%     % we wish the evals matrices to be nonnegative definite; after forcing nnd-ness above, we may
%     % still have that the evals matrices have evalsofevals<0, but very small. Allow this up to some tolerance
%         inds = abs(evalsofevals(:))<zerotol;
%         evalsofevals(inds) = 0;
    end
    fprintf('Nx chosen: %i',Nx);
%%
    

        s = -1;
        amat = nan(Nbar,P);
        count0=0;
        count1=0;
        count2=0;


        while s<Ndaggbar-1
            s = s+1;
%             j = pdagg(s+1,:);
            j = [pdaggx(s+1),pdaggy(s+1)];
            Aindicator = j>0 & j < [Nx,Ny]./2;
            xi = sum(Aindicator);
            if xi==0
                amat(j(1)+Nx*j(2)+1,:) = mvnrnd(zeros(1,P),evals(:,:,j(1)+Nx*j(2)+1))./sqrt(Nbar);
                count0 = count0+1;
            elseif xi==1
                if Aindicator(1)
                    j1 = [Nx-j(1),j(2)];
                elseif Aindicator(2)
                    j1 = [j(1),Ny-j(2)];
                else
                   disp('error in Aindicator #1');
                end
                j2 = j;
                ind1 = j1(1)+Nx*j1(2) + 1; % find rowindex that tells us which row of pdagg j1 corresponds to
                ind2 = j2(1)+Nx*j2(2) + 1;
                U = mvnrnd(zeros(1,P),evals(:,:,ind1));
                V = mvnrnd(zeros(1,P),evals(:,:,ind1));
                amat([ind1,ind2],:) = [U+1i*V;U-1i*V]./sqrt(2*Nbar);
                count1=count1+2;
            elseif xi==2
                j1 = [Nx-j(1),Ny-j(2)];
                j2 = j;
                ind1 = j1(1)+Nx*j1(2) + 1; % find rowindex that tells us which row of pdagg j1 corresponds to
                ind2 = j2(1)+Nx*j2(2) + 1;

                U = mvnrnd(zeros(1,P),evals(:,:,ind1));
                V = mvnrnd(zeros(1,P),evals(:,:,ind1));
                amat([ind1,ind2],:) = [U+1i*V;U-1i*V]./sqrt(2*Nbar);

                j1 = [Nx-j(1),j(2)];
                j2 = [j(1),Ny-j(2)];
                ind1 = j1(1)+Nx*j1(2) + 1; % find rowindex that tells us which row of pdagg j1 corresponds to
                ind2 = j2(1)+Nx*j2(2) + 1;
                try
                U = mvnrnd(zeros(1,P),evals(:,:,ind1));
                V = mvnrnd(zeros(1,P),evals(:,:,ind1));
                catch
                    warning('blah');
                end
                amat([ind1,ind2],:) = [U+1i*V;U-1i*V]./sqrt(2*Nbar);
                count2=count2+4;
            else
                disp('error in Aindicator');
                return;
            end

        end
        clear evals; % the array 'evals' is huge - one of the biggest arrays created in this program.
                     % delete it to save memory
%%
        GRF = fft2(reshape(amat,Nx,Ny,P)); % this GRF is the extendedGRF, which we reduce 3 lines below
        inds = abs(imag(GRF))<imagtol;
        GRF(inds) = real(GRF(inds));
        GRF = GRF(1:Mx,1:My,:);

        % Above follows Wood & Chan (1994) for creating a d-dimensional GRF
        % where, in an (N_1*N_2*...*N_d)-by-d matrix of indices for each
        % dimension, the first column cycles first, the second cycles
        % second, etc. It is convenient to represent a 2-dimensional GRF as
        % a matrix, however for a matrix the first dimension that cycles 
        % refers to the row numbers, ie on the vertical. We therefore
        % transpose our random field, which we are treating as a matrix, so
        % that it appears correctly in plots.
        GRF = permute(GRF,[2,1,3]);

        varargout{1} = GRF;