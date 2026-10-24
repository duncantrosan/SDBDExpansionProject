function [t,V,I,Info] = ReadScopeCSV(File)
% Read a two-channel Tektronix MSO64 CSV export: CH1 = HV probe output, CH2 = current probe.
% Returns rows of time (s), CH1 and CH2 exactly as stored (V), and header values:
% Info.SampleInterval (s), Info.VScale ([CH1 CH2] V/div), Info.File.
% Much faster than OscopeReadin for 1.25 M-point records, and reads one file at a time
% so a folder of captures never has to sit in memory at once.
M = readmatrix(File,'NumHeaderLines',21);
t = M(:,1)'; V = M(:,2)'; I = M(:,3)';
Hdr = fileread(File); Hdr = Hdr(1:min(2000,end));
Info.File = File;
Info.SampleInterval = str2double(regexp(Hdr,'Sample Interval,([^,\r\n]+)','tokens','once'));
Info.VScale = str2double(regexp(Hdr,'Vertical Scale,([^,]+),([^\r\n]+)','tokens','once'));
end
