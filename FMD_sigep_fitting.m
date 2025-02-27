function data_cell = FMD_sigep_fitting(Original_UTSS_Data)
%%% ====================================================================== %%%
% This function will...
% Take UTSS data (stress-strain) from MTC_Analysis_LiveFunc()
% Interpolate
% Stress zero the data
% Fit a polynomial to the desired strain_range
% Calculate secant and tangent modulus off polynomial fit at strain_range.
%
%   Plotting:
%       1. Original data (black circle)
%       2. Yields (green arrows)
%       3. Interpolated data (orange dash line)
%       4. Polynomial fit (black line)
%       5. Secant and Tangent Modulus (blue lines)
%
%   OPTIONS:
%       1. interpolation_step
%       2. interpolation_type
%       3. polynomial_order (1 = first order (y=mx+b))
%       4. strain_range = [0.300, 0.650] is 30-65% strain
%       5. tangent_point_strain = e.g. 0.06 (must be within strain_range)
%
%
% SAVE AS: WT_T2_data_cell = MTC_sigep_fitting(WT_T2_UTSS)
% 
%%%% Robert J. Wiener (c) Oct. 2021 %%%%
%=========================================================================%

%%%%%%%%%%%%% OPTIONS %%%%%%%%%%%%%
interpolation_steps = 0.01;
interpolation_type = 'linear' ;
polynomial_order = 2;  % (1 = first order (y=mx+b))
strain_range = [0.300, 2.50];  % set strain range (e.g. 150% strain_range=1.500) %note this has interpolation_steps
tangent_point_strain = 2.500; % typically max strain_range (strain_range(2)), but can be e.g. 0.05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Store original data for data_cell
% Shift interpolated data to stress starting at 0 for comparing power fits
Ringdat_SSsigep_UTSS = Original_UTSS_Data;
Ringdat_SSsigep_UTSS(:,2) = Original_UTSS_Data(:,2) - Original_UTSS_Data(1,2);


%%% INTERPOLATE stress-strain data from UTSS in 1% strain intervals
Ringdat_SSsigep_UTSS(end-1:end,:) = []; %interp error, end of broken tissue reads non-unique strain values
x = Ringdat_SSsigep_UTSS(:,1);
v = Ringdat_SSsigep_UTSS(:,2);
xq = [0 : interpolation_steps : Ringdat_SSsigep_UTSS(end,1)];
vq1 = interp1(x,v,xq, interpolation_type); %default='linear', also try 'spline' or 'cubic'
%use this to plot original and interpolated data
%figure; plot(x,v,'o',xq,vq1,':.');
xq = xq';
vq1 = vq1';
intpdat = horzcat(xq, vq1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FITTING THE STRESS_STRAIN DATA %%%
%Polynomial fit to calculated data ranges
%Fit curve to multiple data segments for stiffness interpretation
intpdat(:,1) = round(intpdat(:,1),3); %round intpdat(strain) to 3 decimals for integer
strain_range_index = [ find( intpdat(:,1) == strain_range(1) ), find( intpdat(:,1) == strain_range(2) ) ];

%Linear fit "polyfit(~,~,1)" means p(1) is slope/stiffness(m)
[p, S] = polyfit(intpdat(strain_range_index(1):strain_range_index(2) ,1), intpdat(strain_range_index(1):strain_range_index(2) ,2), polynomial_order); %1=first order (y=mx+b)
[y_est, delta] = polyval(p, intpdat(strain_range_index(1):strain_range_index(2), 1), S); %polyfit y_estimates & CI
%Calculare R2 error
SSR = sum( ( vq1(strain_range_index(1):strain_range_index(2)) - y_est ) .^2);
SST = sum( ( vq1(strain_range_index(1):strain_range_index(2)) - mean(vq1(strain_range_index(1):strain_range_index(2))) ) .^2);
R2 = 1 - (SSR/SST);

%find peaks, first local max = yield stress (Y_stress)
[pks, pks_loc] = findpeaks(Ringdat_SSsigep_UTSS(:,2),Ringdat_SSsigep_UTSS(:,1));
Y_stress = pks(1); %yield stress (kPa)
Y_strain = pks_loc(1); %yield strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATING SECANT MODULUS FROM STRAIN_RANGE %%%
% Slope of line defined by (y2-y1)/(x2-x1)
% Use this matrix for plotting
sec(1,1) = strain_range(1);
sec(2,1) = strain_range(2);
sec(1,2) = y_est(1);
sec(2,2) = y_est(end);

% Use this as secant slope for reporting modulus
secant_slope = ( y_est(end) - y_est(1) ) / ( strain_range(2) - strain_range(1) ) ;


%%% CALCULATING TANGENT MODULUS FROM STRAIN_RANGE %%%
% differnetiate polyfit equation (which is over strain range)
syms T(x)

% equation of polyfit over strain range
T(x) = (p(1))*x.^2 + (p(2))*x + (p(3));
% differentiate T wrt x (x is strain)
dT(x) = diff(T,x);  %result should be dT(x)= (p(1))*2*x + (p(2))

% find slope of tangent line @ tangent_point
tangent_slope = (p(1))*2*(tangent_point_strain) + (p(2));  %dT(tangent_point_strain) = (p(1))*2*(tangent_point_strain) + (p(2))
% find y-intercept of tangent line @ tangent_point
% plugging into y=mx+b
% note 'b' is based off p_polyfit function
tangent_intercept = -1* ( (tangent_slope*(tangent_point_strain)) - ((p(1)*(tangent_point_strain)^2 + p(2)*(tangent_point_strain) + p(3) ) ) ) ;

% Use this matrix for plotting
% plot tang over entire strain_range
tang(1,1) = strain_range(1);
tang(2,1) = strain_range(2);
tang(1,2) = tangent_slope*strain_range(1) + tangent_intercept; %y=mx+b
tang(2,2) = tangent_slope*strain_range(2) + tangent_intercept; %y=mx+b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%% PLOTTING %%%%
% Plot orginal, interpolated, curve fits + CI, and local peaks
figure; plot(Ringdat_SSsigep_UTSS(:,1), Ringdat_SSsigep_UTSS(:,2), 'ko','LineWidth',1.5 ); %original
hold on; plot( xq(strain_range_index(1):strain_range_index(2)), vq1(strain_range_index(1):strain_range_index(2)) ); %interpolated
hold on; plot( intpdat((strain_range_index(1):strain_range_index(2)), 1), intpdat((strain_range_index(1):strain_range_index(2)), 2), 'r--','LineWidth',4 ); %zero_interpolated
hold on; plot( intpdat(strain_range_index(1):strain_range_index(2) ,1), y_est, 'k','LineWidth',1.5 ); %fit
%hold on; plot( intpdat(strain_range_index(1):strain_range_index(2) ,1), y_est+2*delta, 'k--', intpdat(strain_range_index(1):strain_range_index(2) ,1), y_est-2*delta, 'k--' ); %95_CI
hold on; findpeaks(Ringdat_SSsigep_UTSS(:,2),Ringdat_SSsigep_UTSS(:,1));
hold on; plot(sec(:,1), sec(:,2), 'b-', 'LineWidth',1.5 ); %secant
hold on; plot(tang(:,1), tang(:,2), 'b-', 'LineWidth',1.5 ); %tangent
xlabel('Strain ε'); ylabel('Stress σ [kPa]');
%xlim([ ]); ylim([0, 10])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store Data
data_cell{1,1} = 'original data'; data_cell{2,1} = Original_UTSS_Data;
data_cell{1,2} = 'interpolated data (zero-shifted)'; data_cell{2,2} = intpdat;
data_cell{1,3} = strcat('p_polyfits',mat2str(strain_range)); data_cell{2,3} = p;
data_cell{1,4} = 'y_est'; data_cell{2,4} = y_est;
data_cell{1,5} = 'R2'; data_cell{2,5} = R2;
data_cell{1,6} = 'delta_CI'; data_cell{2,6} = delta;
data_cell{1,7} = strcat('secant_E',mat2str(strain_range)); data_cell{2,7} = secant_slope;
data_cell{1,8} = strcat('tangent_E[',num2str(tangent_point_strain) ,']'); data_cell{2,8} = tangent_slope;
data_cell{1,9} = 'yield stress'; data_cell{2,9} = Y_stress;
data_cell{1,10} = 'yield strain'; data_cell{2,10} = Y_strain;

disp([data_cell{2,3}(1); data_cell{2,3}(2); data_cell{2,5}; data_cell{2,7}; data_cell{2,8}])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  MODIFICATION IDEAS  %%%
% Try Lo/Hi poly fit regions
% Try basing strain_range off Yield_stress
% Try higher order polyfits for UTSS data
% Integrate data to find AOC
end

