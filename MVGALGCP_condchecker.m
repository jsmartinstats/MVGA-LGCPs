function allcondssatisfied = MVGALGCP_condchecker(paramests)
% A function to check whether Conditions 1-4 in our manuscript are
% satisfied by the MVGALGCP specified by the parameters in params

% For the moment, this assumes the model to be bivariate

% last modified by jsmartin.stats@gmail.com in June 2017

%%
% extract the parameter estimates
alpha_est = paramests{1};
    alpha_est11 = alpha_est(1,1);
    alpha_est22 = alpha_est(2,2);
    alpha_est12 = alpha_est(2,1);
nu_est = paramests{2};
    nu_est11 = nu_est(1,1);
    nu_est22 = nu_est(2,2);
    nu_est12 = nu_est(2,1);
sigma_est = paramests{3};
    sigma_est11 = sigma_est(1,1);
    sigma_est22 = sigma_est(2,2);
    sigma_est12 = sigma_est(2,1);
theta_est = paramests{4};
    theta_est11 = theta_est(1,1);
    theta_est22 = theta_est(2,2);
    theta_est12 = theta_est(2,1);
zeta_est = paramests{5};
    zeta_est11 = zeta_est(1,1);
    zeta_est22 = zeta_est(2,2);
    zeta_est12 = zeta_est(2,1);

% Checking Condition 1:
    nu12_LB = (nu_est11+nu_est22)./2;
    cond1 = nu_est12>=nu12_LB;
    if cond1
        fprintf('Condition 1 satisfied\n');
    else
        fprintf('Condition 1 not satisfied\n');
    end
    
% Checking Remark 2: [Rmk2 true => Cond2 true]
    icl11 = 4.*nu_est11./(alpha_est11.^2);
    icl22 = 4.*nu_est22./(alpha_est22.^2);
    icl12_LB = (icl11+icl22)./2;
    icl12 = 4.*nu_est12./(alpha_est12.^2);
    rmk2 = icl12>=icl12_LB;
    if rmk2
        fprintf('Remark 2 satisfied\n');
    else
        fprintf('Remark 2 not satisfied\n');
    end
    
% Checking Condition 3
    Delta_nu = nu_est12-nu12_LB;
    cond3mat11 = (zeta_est11*sigma_est11/pi) * gamma(nu_est11+1)/(gamma(nu_est11+1)*gamma(nu_est11)) * icl11^(Delta_nu+nu_est11);
    cond3mat22 = (zeta_est22*sigma_est22/pi) * gamma(nu_est22+1)/(gamma(nu_est22+1)*gamma(nu_est22)) * icl22^(Delta_nu+nu_est22);
    cond3mat12 = (zeta_est12*sigma_est12/pi) * gamma(nu_est12+1)/(gamma(nu12_LB+1)*gamma(nu_est12)) * icl12^(Delta_nu+nu12_LB);
    cond3 = all(eig([cond3mat11,cond3mat12;cond3mat12,cond3mat22])>=0);
    if cond3
        fprintf('Condition 3 satisfied\n');
    else
        fprintf('Condition 3 not satisfied\n');
    end
    
% Checking Remark 4: [Rmk4 true => Cond4 true]
    th=theta_est11; z=zeta_est11;
    Sig11 = reshape([cos(th).^2 + z.*z.*sin(th).^2,...
                    (1-z.^2).*sin(th).*cos(th),...
                    (1-z.^2).*sin(th).*cos(th),...
                     z.*z.*cos(th).^2 + sin(th).^2 ...
                     ]',2,2,numel(th));
    th=theta_est22; z=zeta_est22;
    Sig22 = reshape([cos(th).^2 + z.*z.*sin(th).^2,...
                    (1-z.^2).*sin(th).*cos(th),...
                    (1-z.^2).*sin(th).*cos(th),...
                     z.*z.*cos(th).^2 + sin(th).^2 ...
                     ]',2,2,numel(th));
    th=theta_est12; z=zeta_est12;
    Sig12 = reshape([cos(th).^2 + z.*z.*sin(th).^2,...
                    (1-z.^2).*sin(th).*cos(th),...
                    (1-z.^2).*sin(th).*cos(th),...
                     z.*z.*cos(th).^2 + sin(th).^2 ...
                     ]',2,2,numel(th));
    Sig12_LB = (Sig11 + Sig22)/2 ;
    Sig12_LB(2,1)=nan;
    Sig12_LB(1,2)=nan;
    rmk4 = all([Sig12(1,1),Sig12(2,2)] >= [Sig12_LB(1,1),Sig12_LB(2,2)]);
    if rmk4
        fprintf('Remark 4 satisfied\n');
    else
        fprintf('Remark 4 not satisfied\n');
    end

allcondssatisfied = all([cond1,rmk2,cond3,rmk4]);    


end