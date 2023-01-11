function data_cell = MTC_sigep_fitting(Original_UTSS_Data)
%%% ====================================================================== %%%
% This function will...
% Take UTSS data (stress-strain) from MTC_Analysis_LiveFunc()
% Interpolate
% Stress zero the data
% Fit a polynomial to the desired strain_range
% Calculate secant and tangent modulus off polynomial fit at strain_range.
%
% NOTE: code is written for second order (quadratic) modeling
%
% REQUIRED ADD-ONS: <Symbolic Math Toolbox>
%
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
%       3. polynomial_order ( 2=second order: y=p(1)x^2+p(2)x+p(3) )
%       4. strain_range = [0.300, 0.650] is 30-65% strain
%       5. tangent_point_strain_1 = e.g. 0.06 (must be within strain_range)
%       6. tangent_point_strain_2 = second tangent point option(must be within strain_range)
%
%
% SAVE AS: WT_T2_data_cell = MTC_sigep_fitting(WT_T2_UTSS)
% 
%%%% Robert J. Wiener (c) Oct. 2021 %%%%
%=========================================================================%

%%%%%%%%%%%%% OPTIONS %%%%%%%%%%%%%
interpolation_steps = 0.001; %should be around 1/100 of motor_step_size
interpolation_type = 'linear' ;
polynomial_order = 2;  % ( 2=second order: y=p(1)x^2+p(2)x+p(3) )
strain_range = [0.30, 0.65];  % set strain range (e.g. 150% strain_range=1.500) %note this has interpolation_steps
tangent_point_strain_1 = strain_range(2); % typically max strain_range (strain_range(2)), but can be e.g. 0.05
tangent_point_strain_2 = median(strain_range); %try defining as middle of strain_range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Store original data for data_cell
Original_UTSS_Data = Original_UTSS_Data{2,7};
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FITTING THE STRESS_STRAIN DATA %%%
%Polynomial fit to calculated data ranges
%Fit curve to multiple data segments for stiffness interpretation
intpdat(:,1) = round(intpdat(:,1),3); %round intpdat(strain) to 3 decimals for integer
strain_range_index = [ find( intpdat(:,1) == strain_range(1) ), find( intpdat(:,1) == strain_range(2) ) ];

%Quadratic fit "polyfit(~,~,2)" means p(1) is p(1)*x^2
[p, S] = polyfit(intpdat(strain_range_index(1):strain_range_index(2) ,1), intpdat(strain_range_index(1):strain_range_index(2) ,2), polynomial_order); %2=second order: y=p(1)x^2+p(2)x+p(3)
[y_est, delta] = polyval(p, intpdat(strain_range_index(1):strain_range_index(2), 1), S); %polyfit y_estimates & CI
%Calculare R2 error
SSR = sum( ( vq1(strain_range_index(1):strain_range_index(2)) - y_est ) .^2);
SST = sum( ( vq1(strain_range_index(1):strain_range_index(2)) - mean(vq1(strain_range_index(1):strain_range_index(2))) ) .^2);
R2 = 1 - (SSR/SST);

%find peaks, first local max = yield stress (Y_stress)
[pks, pks_loc] = findpeaks(Ringdat_SSsigep_UTSS(:,2),Ringdat_SSsigep_UTSS(:,1));
Y_stress = pks(1); %yield stress (kPa)
Y_strain = pks_loc(1); %yield strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  CACULATE MODULI  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATE SECANT MODULUS FROM STRAIN_RANGE %%%
% Slope of line defined by (y2-y1)/(x2-x1)
% Use this matrix for plotting
sec(1,1) = intpdat(1,1); %initial strain (= 0) [sec x1]
sec(2,1) = strain_range(2); %end strain (at end of strain range) [sec x2]
sec(1,2) = intpdat(1,2); %initial stress = 0 [sec y1]
sec(2,2) = y_est(end); %end stress (at end of strain range) of fitted curve [sec y2]

% Use this as secant slope for reporting modulus
secant_slope = ( y_est(end) - intpdat(1,2) ) / ( strain_range(2) - intpdat(1,1) ) ;

%%%%%%%%%%%%%%%%
%%% CALCULATE TANGENT MODULUS FROM STRAIN_RANGE %%%
% differnetiate polyfit equation (which is over strain range)
syms Fstrain x

% equation of polyfit over strain range
Fstrain = @(x) (p(1))*x.^2 + (p(2))*x + (p(3));
% differentiate ds wrt x (x is strain)
ds(x) = diff(Fstrain,x);  %result should be ds(x)= (p(1))*2*x + (p(2))

% find slope of tangent line @ tangent_point
% manual power law diff of fit (change if higher order poly_fit)
% TANGENT 1 %
tangent_slope_1 = (p(1))*2*(tangent_point_strain_1) + (p(2));  %dT(tangent_point_strain_1) = (p(1))*2*(tangent_point_strain_1) + (p(2))
% find y-intercept of tangent line @ tangent_point
% plugging into y=mx+b
% note 'b' is based off p_polyfit function
tangent_intercept_1 = -1* ( (tangent_slope_1*(tangent_point_strain_1)) - ((p(1)*(tangent_point_strain_1)^2 + p(2)*(tangent_point_strain_1) + p(3) ) ) ) ;

% TANGENT 2 %
tangent_slope_2 = (p(1))*2*(tangent_point_strain_2) + (p(2));  %dT(tangent_point_strain_2) = (p(1))*2*(tangent_point_strain_2) + (p(2))
tangent_intercept_2 = -1* ( (tangent_slope_2*(tangent_point_strain_2)) - ((p(1)*(tangent_point_strain_2)^2 + p(2)*(tangent_point_strain_2) + p(3) ) ) ) ;


% Use this matrix for plotting
% plot tang over entire strain_range
% TANG 1 %
tang_1(1,1) = strain_range(1);
tang_1(2,1) = strain_range(2);
tang_1(1,2) = tangent_slope_1*strain_range(1) + tangent_intercept_1; %y=mx+b
tang_1(2,2) = tangent_slope_1*strain_range(2) + tangent_intercept_1; %y=mx+b

% TANG 2 %
tang_2(1,1) = strain_range(1);
tang_2(2,1) = strain_range(2);
tang_2(1,2) = tangent_slope_2*strain_range(1) + tangent_intercept_2; %y=mx+b
tang_2(2,2) = tangent_slope_2*strain_range(2) + tangent_intercept_2; %y=mx+b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATE STRAIN ENERGY FROM POLY_FIT ACROSS STRAIN_RANGE %%%

%syms Fstrain x % redundant this is initialized in TANGENT
%Fstrain = @(x) (p(1))*x.^2 + (p(2))*x + (p(3));

% Integrate poly_fit function, only over strain_range (offset to account for bending)
% calculate the numeric integration of input function ‘Fx’, which in turn signifies the area under a curve.
% A = integral (Fx, Xminimum, Xmaximum)
AOC_over_strain_range = integral(Fstrain, strain_range(1), strain_range(2)); %for our quadratic this in incorrect because it doesnt fit to zero

% This gives the area under the curve represented by ‘y’. Here ‘x’ is used to define the range or limits between which we want the area.
% A = trapz (x, y)
% A = cumtrapz (x, y) will compute the cumulative integration of Y wrt X.
% This is done using the trapezoidal integration and can be used to calculate the area under the curve for a portion.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PLOTTING %%%%
% Plot orginal, interpolated, curve fits + CI, local peaks, etc.
figure;
plot(Ringdat_SSsigep_UTSS(:,1), Ringdat_SSsigep_UTSS(:,2), 'ko','LineWidth',1.5 ); %original
hold on; plot( xq(strain_range_index(1):strain_range_index(2)), vq1(strain_range_index(1):strain_range_index(2)) ); %interpolated
hold on; plot( intpdat((strain_range_index(1):strain_range_index(2)), 1), intpdat((strain_range_index(1):strain_range_index(2)), 2), 'r--','LineWidth',4 ); %zero_interpolated
hold on; plot( intpdat(strain_range_index(1):strain_range_index(2) ,1), y_est, 'k','LineWidth',1.5 ); %fit
%hold on; plot( intpdat(strain_range_index(1):strain_range_index(2) ,1), y_est+2*delta, 'k--', intpdat(strain_range_index(1):strain_range_index(2) ,1), y_est-2*delta, 'k--' ); %95_CI
hold on; findpeaks(Ringdat_SSsigep_UTSS(:,2),Ringdat_SSsigep_UTSS(:,1));
hold on; plot(sec(:,1), sec(:,2), 'b-', 'LineWidth',1.5 ); %secant
hold on; plot(tang_1(:,1), tang_1(:,2), 'b-', 'LineWidth',1.5 ); %tangent_1
hold on; plot(tang_2(:,1), tang_2(:,2), 'b-', 'LineWidth',1.5 ); %tangent_2

hold on; % plot area of AOC (strain_range)
% area(x, y) % where y is Fstrain(x) func
% AOC_x is index of interpolated data at strain ranges
AOC_x = [ intpdat( strain_range_index(1):strain_range_index(2),  1) ]; 
FA1 = area(AOC_x, (p(1))*AOC_x.^2 + (p(2))*AOC_x + (p(3))) ;
FA1.FaceAlpha = 0.2;
FA1.EdgeAlpha = 0;

% Detail Plot %
xlabel('Strain ε'); ylabel('Stress σ [kPa]');
%xlim([ ]);
ylim([ 0, max(intpdat(:,2))+max(intpdat(:,2))*0.1 ])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store Data
data_cell{1,1} = 'original data'; data_cell{2,1} = Original_UTSS_Data;
data_cell{1,2} = 'interpolated data (zero-shifted)'; data_cell{2,2} = intpdat;
data_cell{1,3} = 'polynomial_order'; data_cell{2,3} = polynomial_order;
data_cell{1,4} = strcat('p_polyfits',mat2str(strain_range)); data_cell{2,4} = p;
data_cell{1,5} = 'y_est'; data_cell{2,5} = y_est;
data_cell{1,6} = 'R2'; data_cell{2,6} = R2;
data_cell{1,7} = 'delta_CI'; data_cell{2,7} = delta;
data_cell{1,8} = strcat('secant_E[',mat2str(strain_range(2)),']'); data_cell{2,8} = secant_slope;
data_cell{1,9} = strcat('tangent_1_E[',num2str(tangent_point_strain_1) ,']'); data_cell{2,9} = tangent_slope_1;
data_cell{1,10} = strcat('tangent_2_E[',num2str(tangent_point_strain_2) ,']'); data_cell{2,10} = tangent_slope_2;
data_cell{1,11} = strcat('strain energy',mat2str(strain_range)); data_cell{2,11} = AOC_over_strain_range;
data_cell{1,12} = 'yield stress'; data_cell{2,12} = Y_stress;
data_cell{1,13} = 'yield strain'; data_cell{2,13} = Y_strain;


%print: p1; p2; p3, r2; secant_E; tangent_E_1; tangent_E_2; AOC_over_strain_range; yield_stress; yield_strain
format shortG
disp([data_cell{2,4}(1); data_cell{2,4}(2); data_cell{2,4}(3);...
 data_cell{2,6}; data_cell{2,8}; data_cell{2,9}; data_cell{2,10}; data_cell{2,11}; data_cell{2,12}; data_cell{2,13}])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  MODIFICATION IDEAS  %%%
% Try Lo/Hi poly fit regions
% Try basing strain_range off Yield_stress
% Try higher order polyfits for UTSS data
% Integrate data to find AOC

end

