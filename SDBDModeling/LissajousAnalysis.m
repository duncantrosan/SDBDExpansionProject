function [CapacitanceValues,Energy] = LissajousAnalysis(Time,Voltage,Charge);

[NFiles,NPoints] = size(Time);
t = Time(1,:);
VoltageMean = mean(Voltage);
ChargeMean = mean(Charge);
figure
plot(VoltageMean,ChargeMean);
[x,y] = getpts;
C_Cell1 = (y(2) - y(1))./ (x(2) - x(1));
C_d1 = (y(3) - y(2))./ (x(3) - x(2));
C_Cell2 = (y(4) - y(3))./ (x(4) - x(3));
C_d2 = (y(4) - y(1))./ (x(4) - x(1));
% for n = 1:NFiles
%     % Find Smooth Spline for Voltage and Charge 
%     % VoltageFit = fit(Time(n,:)', Voltage(n,:)', 'smoothingspline', 'SmoothingParam', 0.5);
%     % ChargeFit =  fit(Time(n,:)', Charge(n,:)', 'smoothingspline', 'SmoothingParam', 0.5);
% 
% end
CapacitanceValues = [C_Cell1;C_Cell2;C_d1;C_d2];
Energy = 1;