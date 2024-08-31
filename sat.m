function output = sat(inputmat,nvec)
% This function uses a summed area table (the matrix analogue of a cumulative sum)
% to sum the specified submatrices of a given matrix input. At the minute, this
% only allows for a symmetric submatrix structure.

% INPUT:    inputmat - the data matrix with submatrices to be summed
%           nvec     - the split of elements in each dimension of the full matrix
%                       i.e. the full matrix has dimensions sum(nvec) by sum(nvec)

% last modified by jsmartin.stats@gmail.com in June 2016

    inputmat(isnan(inputmat)) = 0;
    inputmat = cumsum(cumsum(inputmat),2);
    cumsumn = cumsum(nvec);
    if sum(cumsumn==0)>0
        cumsumn_nz=cumsumn(cumsumn>0);
        gnumsummand1 = inputmat(cumsumn_nz,cumsumn_nz);
        gnumsummand2 = zeros(size(gnumsummand1));
        gnumsummand3 = zeros(size(gnumsummand1));
        gnumsummand4 = zeros(size(gnumsummand1));        inputmat(isnan(inputmat)) = 0;

        gnumsummand2(2:end,2:end) = inputmat(cumsumn_nz(1:end-1),cumsumn_nz(1:end-1));
        gnumsummand3(2:end,:) = inputmat(cumsumn_nz(1:end-1),cumsumn_nz);
        gnumsummand4(:,2:end) = inputmat(cumsumn_nz,cumsumn_nz(1:end-1));
        gnum2 = gnumsummand1+gnumsummand2-gnumsummand3-gnumsummand4;
        gnum2(gnum2<=eps) = 0;
        gnum3 = zeros(numel(cumsumn),numel(cumsumn));
        gnum3(sum(cumsumn==0)+1:end,sum(cumsumn==0)+1:end) = gnum2;
    else
        gnumsummand1 = inputmat(cumsumn,cumsumn);
        gnumsummand2 = zeros(size(gnumsummand1));
        gnumsummand3 = zeros(size(gnumsummand1));
        gnumsummand4 = zeros(size(gnumsummand1));
        gnumsummand2(2:end,2:end) = inputmat(cumsumn(1:end-1),cumsumn(1:end-1));
        gnumsummand3(2:end,:) = inputmat(cumsumn(1:end-1),cumsumn);
        gnumsummand4(:,2:end) = inputmat(cumsumn,cumsumn(1:end-1));
        gnum2 = gnumsummand1+gnumsummand2-gnumsummand3-gnumsummand4;
        gnum2(gnum2<=eps) = 0;
        gnum3 = gnum2;
    end
    output = gnum3;
end