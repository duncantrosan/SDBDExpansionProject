clear
close all
imagefiles = dir('/Users/duncantrosan/Documents/MATLAB/Data/2021/ICCDImages/09023c');      
nfiles = length(imagefiles);    % Number of files found

for n = 4:nfiles
currentfilename = fullfile(imagefiles(n).folder,imagefiles(n).name);
files{1,n-3} = currentfilename;
names{1,n-3} = imagefiles(n).name;
names{1,n-3} = names{1,n-3}(1:end-4);
end

names = natsort(names);
files = natsort(files);

nfiles = nfiles-3;
% writerObj = VideoWriter('ICCD.gif');
% writerObj.FrameRate = 1;
   % open the video writer
   %open(writerObj);
   
   
dis = 'discharge.mat';   
st = load(dis);
discharge = st.sdischarge;
dT = discharge(1,2) - discharge(1,1);
inital = round(-discharge(1,1)/dT);

 
filename = 'Testdis.gif';
for ii=1:nfiles
   h = figure();
   nexttile
   currentfilename = files{1,ii};
   currentimage = imread(currentfilename);
   images{ii} = currentimage;
   
   image(currentimage)
   colorMap = [linspace(0,1,5000)', linspace(0,1,5000)',linspace(0,1,5000)'];
   colormap(colorMap);
   colorbar
   set(gca,'visible','off')
   set(gca,'xtick',[])
   text(0,100,[names{ii}, 'μs'],'FontSize',40,'FontName','Arial','Color','white')
   
   % Discharge Plot
   nexttile
   plot(10^6*discharge(1,1:round(inital+(ii-1)*10^-6*10/dT)),10^3*discharge(2,1:round(inital+(ii-1)*10*10^-6/dT)),'Color','r');
   axis([-10,180,-1,20])
   xlabel('Time / μs')
   ylabel('Current / mA')
   
   frame = getframe(h);
   im = frame2im(frame);
   
  [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if ii == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','Delaytime',0.5); 
      end 
   
    clf
   
end
