function [c, y_est] = sigep_fitting(Ringdat_SS_sigep)
%%% ============================ %%%
% add function description/summary
% include what each graph shows and how their output is calculated
% ref inputs and returns

%%%% Robert J. Wiener (c) Oct. 2021 %%%%
%========================================

%%% INTERPOLATE stress-strain data from sigep output in 1% strain intervals
x = Ringdat_SS_sigep(:,1);
v = Ringdat_SS_sigep(:,2);
xq = [0:0.01:Ringdat_SS_sigep(end,1)];
vq1 = interp1(x,v,xq); %default='linear', also try 'spline' or 'cubic'
%use this to plot original and interpolated data
%plot(x,v,'o',xq,vq1,':.');

%%% SECTION stress-strain data into low-med-high ranges
%using a buffer of 5% strain i.e. -5
%segment remaining interplated data into 3 sections
% sec_interval = (length(xq)-5)/3;
% fl = floor(sec_interval);
% ce = ceil(sec_interval);
% 
% %create 2D arrays x=stain, y=stress(kPa)
% low = vertcat([xq(6:5+fl)], [vq1(6:5+fl)]);
% med = vertcat([xq(5+ce:5+fl+fl)], [vq1(5+ce:5+fl+fl)]);
% hi = vertcat([xq(5+fl+ce:5+fl+fl+fl)], [vq1(5+fl+ce:5+fl+fl+fl)]);
%%%%%%%%%%%%%%%%%%
intpdat = horzcat(xq', vq1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% YOUNGS MODULUS E
%Linear regression for polynomial fit to calculated data ranges
%Fit line to data segments
[c, S] = polyfit(intpdat(1:31,1), intpdat(1:31,2), 1); %1=first order (y=mx+b)
y_est = polyval(c, intpdat(1:31,1));%error


% Plot orginal, interpolated, and curve fits
plot(xq(1:31),vq1(1:31),'--','LineWidth',3);
hold on; plot(intpdat(1:31,1), y_est,'k','LineWidth',1.5)
ylabel('Stress σ [kPa]');
xlabel('Strain ε');
ylim([0,.30])
xlim([0,0.30])
end