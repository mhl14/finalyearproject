function [convcoef] = self_cc(infile, stop, step)
%% input toe
stimpath = 'freq_shifted_song/stimulus';
% infile = ...
%     'response/concat_chan_23_24_electrode_7_4/10 reps/ss001m_165_s31_toe.txt';
start =0; 
% stop = stop; 

fid = fopen(infile);
if fid == -1
    error(['Cannot open file ' infile]);
end
toeinfo = textscan(fid, '%s %s', 1, 'delimiter', '\n');
 stimfile = cell2mat(toeinfo{1});
 if (strncmp(cell2mat(toeinfo{2}), 'StimSampleRate', 14))
     stimsampratetemp = cell2mat(toeinfo{2});
     stimsamprate = str2double(stimsampratetemp(17:21));
     toeinfoContinue = textscan(fid, '%s %s %s %s %s', 1, 'delimiter', '\n');
     subjectID = cell2mat(toeinfoContinue{1});
     pen = cell2mat(toeinfoContinue{2});
     site = cell2mat(toeinfoContinue{3});
     sort = cell2mat(toeinfoContinue{4});
     channels = cell2mat(toeinfoContinue{5});
     
 else
     Stimulus sample rate not available.
     stimsamprate = [];
     subjectID = cell2mat(toeinfo{2});
     toeinfoContinue = textscan(fid, '%s %s %s %s', 1, 'delimiter', '\n');
     pen = cell2mat(toeinfoContinue{1});
     site = cell2mat(toeinfoContinue{2});
     sort = cell2mat(toeinfoContinue{3});
     channels = cell2mat(toeinfoContinue{4});
 
 end
 nreps = cell2mat(textscan(fid, '%d', 1,'HeaderLines',5, 'delimiter', '\n') );
 %nreps = repetitions;
 nspikes = cell2mat(textscan(fid, '%d', nreps, 'delimiter', '\n')); % gets 32-bit ingegers nreps times with end-of-the-line step down each time

 toestart=ftell(fid); %indicates the current position within the toe file, i.e. the time of first AP

toes{:,1} = []; %in case there are NO toes in this toe file.
for i=1:nreps
 toes{:,i} = textscan(fid, '%.6f64 %*n', nspikes(i) ); % get into each toe the '%.6f64 %*n' stuff (APs times) taking them vertically as indicated by the current value of nspikes
end

fseek(fid, toestart, 'bof');

fclose(fid);

%% extract half the data
nreps=nreps/2;

toes1= toes {1,1}{1,1};toes2= toes {1,2}{1,1};toes3= toes {1,3}{1,1};
toes4= toes {1,4}{1,1};toes5= toes {1,5}{1,1};toes6= toes {1,6}{1,1};
toes7= toes {1,7}{1,1};toes8= toes {1,8}{1,1};toes9= toes {1,9}{1,1};
toes10= toes {1,10}{1,1};
toes11= toes {1,11}{1,1};toes12= toes {1,12}{1,1};
toes13= toes {1,13}{1,1};toes14= toes {1,14}{1,1};toes15= toes {1,15}{1,1};
toes16= toes {1,16}{1,1};toes17= toes {1,17}{1,1};toes18= toes {1,18}{1,1};
toes19= toes {1,19}{1,1};toes20= toes {1,20}{1,1};toes21= toes {1,21}{1,1};
toes22= toes {1,22}{1,1};toes23= toes {1,23}{1,1};toes24= toes {1,24}{1,1};
toes25= toes {1,25}{1,1};toes26= toes {1,26}{1,1};toes27= toes {1,27}{1,1};
toes28= toes {1,28}{1,1};toes29= toes {1,29}{1,1};toes30= toes {1,30}{1,1};

alltoesodd = [toes1; toes3; toes5; toes7; toes9; ... 
    toes11; toes13; toes15; ...
    toes17; toes19; toes21; toes23; toes25; toes27; toes29];
alltoeseven = [toes2; toes4; toes6; toes8; toes10; ...
    toes12; toes14; toes16; ...
    toes18; toes20; toes22; toes24; toes26; toes28; toes30];

%% calculate spike probability 1
stimfile = strrep(stimfile, '.wav', '');
figname = strrep(sprintf('subj:%s  site:%s  sort1:%s  stim:%s', ...
    subjectID,site,sort,stimfile), '_', '\_');

newtoes = alltoesodd(alltoesodd>=start);
newtoes = newtoes(newtoes<=stop);
%xbins = linspace(start, stop, (stop-start)*bps);
xbins = start:step:stop;
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
psthodd=m;
% figure; 
% plot(xbins, psthodd,'LineWidth',2);

%% calculate spike probability 2
stimfile = strrep(stimfile, '.wav', '');
figname = strrep(sprintf('subj:%s  site:%s  sort1:%s  stim:%s', ...
    subjectID,site,sort,stimfile), '_', '\_');

newtoes = alltoeseven(alltoeseven>=start);
newtoes = newtoes(newtoes<=stop);
%xbins = linspace(start, stop, (stop-start)*bps);
xbins = start:step:stop;
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
pstheven=m;
% figure; 
% plot(xbins, pstheven,'LineWidth',2);

%% calculate correlation coefficient
convcoef = corrcoef(pstheven, psthodd);

end