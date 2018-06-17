%updated from plot_raster_2.m to read in from read_toe_2.m (14 May 2010, AK) 
%updated from plot_raster_2.m to read in from read_toe_2.m (14 Nov 2017, EL)

function [ax, psth, spec]=plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap, ~, ~)

stimpath = 'freq_shifted_song/stimulus';

%% read in the spike2 format toelist
[stimfile, subjectID, ~, site, sort1, ~, nreps, ...
    nspikes, toes, alltoes,~] = readtoe_2(toelist);
stimfile = strrep(stimfile, '.wav', '');
figname = strrep(sprintf('subj:%s  site:%s  sort1:%s  stim:%s', ...
    subjectID,site,sort1,stimfile), '_', '\_');

%% generate the raster
ax(1) = subplot(3,1,1);
for i=1:nreps                                                  
  plot(cell2mat(toes{:,i}), ones(nspikes(i),1) ...
      .* double(i), 'k.','MarkerSize',8 );
  if i==1                                                         
    hold;                                                           
  end                                                            
end
hold;
xmin=start;xmax=stop;
ymin=0;ymax=nreps+1;
axis([xmin xmax ymin ymax]);
set(ax(1), 'FontSize', 20);
ylabel('Rep');

% ax_pos = get(ax(1), 'Position');
set(ax(1),'XTickLabel', '');
% title(figname);

%% generate histogram
ax(2) = subplot(3,1,2);
set(ax(2), 'FontSize', 20);
ax2_pos = get(ax(1), 'Position');
ax2_pos(2)= ax2_pos(2)-ax2_pos(4);
set(ax(2), 'Position',ax2_pos); 

newtoes = alltoes(alltoes>=start);
newtoes = newtoes(newtoes<=stop);
%xbins = linspace(start, stop, (stop-start)*bps);
xbins = start:step:stop; %145128 for slow, 14513 for normal, 145135 for fast
% xbins = start:0.001:stop; %this is for 1 ms bins
if(isempty(newtoes))
    %do nothing
else
    n = histc(newtoes, xbins); % generate the spike count in each  bin
    ymax = max(n)+1;
end
 
nreps = double(nreps);
 
m = rdivide(n,nreps); %this thing did not work because integers can be combined only with integers not doubles, so now it works fine (AK) 
 
% if smooth == 1 %do smoothing
w=3;
ng=[0.2261 0.5478 0.2261];
tamw=(w-1)/2;
if (w>3)
    for i=1:tamw-1
     ng=conv(ng,ng);
    end
end
lng=length(ng);
limite=(lng+1)/2;
lc=length(n);
for i=limite:lc-limite
    n_int=0;
    for k=1:length(ng)-1
        n_int=(ng(k)*n(i-limite+k))+n_int;
    end
    n(i)=n_int;
end
psth=m;
plot(xbins, psth,'LineWidth',2);
% else
%     psth=m;
%     bar(xbins,psth,'LineWidth',2);
% end

axis([xmin xmax ymin ymax]);
ylabel('Spikes probability');
set(ax(2),'XTickLabel', '');
ylim([0 1]);

psth = psth(1:end-1, 1);

%% read in audio files
ax(3) = subplot(3,1,3);
ax3_pos = get(ax(3), 'Position');
ax3_pos(4)= ax2_pos(2)-ax3_pos(2);
set(ax(3), 'FontSize', 20);
set(ax(3), 'Position',ax3_pos);

fullstim = [stimpath '/' stimfile '.wav' ];
[Y,FS]=audioread(fullstim);
nsamples = length(Y);
%now pad the sound so that it is aligned with the raster and psth
if (start<0)
    Y = [zeros(abs(start)* FS,1); Y];
%     fprintf('adding %g samples to the start of the stim\n', abs(start)*FS);
elseif (xmin>0)
    Y = Y((start*FS)+1:nsamples); 
%     fprintf('subtracting %g samples from the front of the stimuli\n', start*FS);
end

if (stop*FS>nsamples)
    Y = [Y; zeros((stop*FS)-nsamples,1) ];
%     fprintf('adding %g samples to the end of the stimulus\n', stop*FS-nsamples);
else
    Y = Y(1:((stop-start)*FS));   
%     fprintf('subtracting %g samples from the end of the stimulus\n',...
%     nsamples-(stop*FS));
end

%% do the spectrogram
nlap = round(window*overlap);
[~,~,T,P] = spectrogram(Y,window,nlap,nfft,FS, 'yaxis');

newT = xmin:((xmax-xmin)/length(T)):xmax;
freqs=0:22050/(nfft/2):22050;
if(length(newT)~=length(newT))
    newT=newT(1:length(T));
end
clim = [-200  -65];
spec = 20*log10(P);
spec = spec(2:end, 1:end);
imagesc(newT,freqs,spec,clim);  
axis xy;
colormap(colormap(jet(256)));

axis([xmin xmax 0 15000]); %11025
xlabel('Time(sec)');
ylabel('Freq (Hz)');

set(gcf, 'Color', [1 1 1]);
%set(gcf, 'Position', [360 635 883 286]);
set(gcf,'Name',figname,'NumberTitle','off' );
set(gcf, 'PaperPositionMode', 'auto', 'Toolbar', 'none');

% if (doprintout)
%   set(gcf,'PaperOrientation', 'landscape');
%   print -noui 
% end

% Now do the SMI calculation
% 
% 
% newtoes_SMI = alltoes(alltoes>=tmin);
% newtoes_SMI = newtoes_SMI(newtoes_SMI<=tmax);
% 
% 
% number_of_spikes = length(newtoes_SMI);
% 
% disp(number_of_spikes)

end  %function end
