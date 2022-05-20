
% store as [_SSsigep_UTSS_interp, _SSsigep_UTSS_c, _SSsigep_UTSS_yest]
function [intpdat, c, y_est] = MTC_sigep_fitting(Ringdat_SSsigep_UTSS)
%%% ============================ %%%
% add function description/summary
% include what each graph shows and how their output is calculated
% ref inputs and returns

%%%% Robert J. Wiener (c) Oct. 2021 %%%%
%========================================

%%% INTERPOLATE stress-strain data from sigep_RR_avg output in 5% strain intervals
Ringdat_SSsigep_UTSS(end-1:end,:) = []; %interp error, end of broken tissue reads non-unique strain values
x = Ringdat_SSsigep_UTSS(:,1);
v = Ringdat_SSsigep_UTSS(:,2);
xq = [0:0.05:Ringdat_SSsigep_UTSS(end,1)];
vq1 = interp1(x,v,xq, 'linear'); %default='linear', also try 'spline' or 'cubic'
%use this to plot original and interpolated data
figure; plot(x,v,'o',xq,vq1,':.');
%
xq = xq';
vq1 = vq1';
intpdat = horzcat(xq, vq1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% YOUNGS MODULUS E
%Linear regression for polynomial fit to calculated data ranges
%Fit line to data segments
strain_range = [0.100, 0.500]; % set strain range (e.g. to 150% strain_range=1.500)
intpdat(:,1) = round(intpdat(:,1),3); %round intpdat(strain) to 3 decimals for integer
strain_range_index = [ find( intpdat(:,1) == strain_range(1) ), find( intpdat(:,1) == strain_range(2) ) ];
[c, S] = polyfit(intpdat(strain_range_index(1):strain_range_index(2) ,1), intpdat(strain_range_index(1):strain_range_index(2) ,2), 1); %1=first order (y=mx+b)
y_est = polyval(c, intpdat(strain_range_index(1):strain_range_index(2), 1));%error


% Plot orginal, interpolated, and curve fits
figure; plot(Ringdat_SSsigep_UTSS(:,1), Ringdat_SSsigep_UTSS(:,2), 'ko','LineWidth',1.5)
hold on; plot( xq(strain_range_index(1):strain_range_index(2)), vq1(strain_range_index(1):strain_range_index(2)), '--','LineWidth',3);
hold on; plot( intpdat(strain_range_index(1):strain_range_index(2) ,1), y_est, 'k','LineWidth',1.5)
ylabel('Stress σ [kPa]');
xlabel('Strain ε');
%ylim([0, 10])
%xlim([ ])


%%% modify RR_avg fits to include various strain ranges
%%% include other fits for UTSS data

end


% Instead of this preprocessing only use final tissue stretch [UTSS]
function foo
%%% Preprocessing of RR runs into averaged tissue data
% SSsigep_RR_avg:
% col_1: strain_AVG
% col_2: stress_AVG
% col_3: strain_STDEV
% col_4: stress_STDEV
%
%%% Change to desired prefix (tissue ID)
for i = 1:length(KO_30_SSsigep_RR1)
    KO_30_SSsigep_RR_avg(i,1) = mean([KO_30_SSsigep_RR1(i,1), KO_30_SSsigep_RR2(i,1), KO_30_SSsigep_RR3(i,1)]);
    KO_30_SSsigep_RR_avg(i,2) = mean([KO_30_SSsigep_RR1(i,2), KO_30_SSsigep_RR2(i,2), KO_30_SSsigep_RR3(i,2)]);
    KO_30_SSsigep_RR_avg(i,3) = std([KO_30_SSsigep_RR1(i,1), KO_30_SSsigep_RR2(i,1), KO_30_SSsigep_RR3(i,1)]);
    KO_30_SSsigep_RR_avg(i,4) = std([KO_30_SSsigep_RR1(i,2), KO_30_SSsigep_RR2(i,2), KO_30_SSsigep_RR3(i,2)]);
end

end