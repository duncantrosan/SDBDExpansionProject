function LissajousPlotting(Time,Voltage,Current,Charge,Resistance)

Ur = Voltage - Current*Resistance;
figure; 
plot(Voltage,Charge);hold on
plot(Ur,Charge)
end
