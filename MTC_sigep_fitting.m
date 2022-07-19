function data_cell = MTC_sigep_fitting(Ringdat_SSsigep_UTSS)
%        WT_T2_data_cell = MTC_sigep_fitting(WT_T2_UTSS)
%%% ============================ %%%
% add function description/summary
% include what each graph shows and how their output is calculated
% ref inputs and returns
%
%%%% Robert J. Wiener (c) Oct. 1, 2021 %%%%
%========================================

%%% INTERPOLATE stress-strain data from sigep_RR_avg output in 1% strain intervals
Ringdat_SSsigep_UTSS(end-1:end,:) = []; %interp error, end of broken tissue reads non-unique strain values
x = Ringdat_SSsigep_UTSS(:,1);
v = Ringdat_SSsigep_UTSS(:,2);
xq = [0:0.01:Ringdat_SSsigep_UTSS(end,1)];
vq1 = interp1(x,v,xq, 'linear'); %default='linear', also try 'spline' or 'cubic'
%use this to plot original and interpolated data
%figure; plot(x,v,'o',xq,vq1,':.');
xq = xq';
vq1 = vq1';
intpdat = horzcat(xq, vq1);
%shift interpolated data to stress starting at 0 for comparing power fits
intpdat(:,2) = intpdat(:,2) - intpdat(1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FITTING THE STRESS_STRAIN DATA
%Polynomial fit to calculated data ranges
%Fit curve to multiple data segments for stiffness interpretation
strain_range = [0.300, 0.650]; % set strain range (e.g. 150% strain_range=1.500)
intpdat(:,1) = round(intpdat(:,1),3); %round intpdat(strain) to 3 decimals for integer
strain_range_index = [ find( intpdat(:,1) == strain_range(1) ), find( intpdat(:,1) == strain_range(2) ) ];

%Linear fit "polyfit(~,~,1)" means p(1) is slope/stiffness(m)
[p, S] = polyfit(intpdat(strain_range_index(1):strain_range_index(2) ,1), intpdat(strain_range_index(1):strain_range_index(2) ,2), 2); %1=first order (y=mx+b)
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




%%%% PLOTTING %%%%
% Plot orginal, interpolated, curve fits + CI, and local peaks
figure; plot(Ringdat_SSsigep_UTSS(:,1), Ringdat_SSsigep_UTSS(:,2), 'ko','LineWidth',1.5)
hold on; plot( xq(strain_range_index(1):strain_range_index(2)), vq1(strain_range_index(1):strain_range_index(2)), '--','LineWidth',3);
hold on; plot( intpdat(strain_range_index(1):strain_range_index(2) ,1), y_est, 'k','LineWidth',1.5)
%hold on; plot( intpdat(strain_range_index(1):strain_range_index(2) ,1), y_est+2*delta, 'k--', intpdat(strain_range_index(1):strain_range_index(2) ,1), y_est-2*delta, 'k--' ) %plot 95% CI
hold on; findpeaks(Ringdat_SSsigep_UTSS(:,2),Ringdat_SSsigep_UTSS(:,1));
xlabel('Strain ε'); ylabel('Stress σ [kPa]');
%xlim([ ]); ylim([0, 10])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store Data
data_cell{1,1} = 'original data'; data_cell{2,1} = Ringdat_SSsigep_UTSS;
data_cell{1,2} = 'interpolated data'; data_cell{2,2} = intpdat;
data_cell{1,3} = 'p_polyfits'; data_cell{2,3} = p;
data_cell{1,4} = 'y_est'; data_cell{2,4} = y_est;
data_cell{1,5} = 'R2'; data_cell{2,5} = R2;
data_cell{1,6} = 'delta_CI'; data_cell{2,6} = delta;
data_cell{1,7} = 'yield stress'; data_cell{2,7} = Y_stress;
data_cell{1,8} = 'yield strain'; data_cell{2,8} = Y_strain;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  MODIFICATION IDEAS  %%%
% Try Lo/Hi poly fit regions
% Try basing strain_range off Yield_stress
% Try higher order polyfits for UTSS data

end

