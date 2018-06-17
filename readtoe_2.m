% readtoe_2.m, which is updated from readtoe.m to read in the 48xchannel
% sorts (AK)
%
% Reads in TOE (Time Of Event) files of the format *_toe.txt and converts
% the data contained therein into matlab vectors.
%
% Modifications:
% 7/13/07 JMJ: toes output modified to be double precision.
% 6/17/08 JMJ: added stimSampleFrequency output term.
% 14 May 2010 AK: converted to readtoe_2.m for the new sort format (48xchannels-just different header).
% 14 Nov 2017 EL: minor syntax update. 

function [stimfile, subjectID, pen, site, sort, channels, nreps, nspikes, toes, alltoes, stimsamprate] = readtoe_2(infile)

fid = fopen(infile);
if fid == -1
    error(['Cannot open file ' infile]);
end
toeinfo = textscan(fid, '%s %s', 1, 'delimiter', '\n');
 stimfile = cell2mat(toeinfo{1});
 if (strncmp(cell2mat(toeinfo{2}), 'StimSampleRate', 14))
     %Stimulus sample rate is available. Incorporate it into output
     stimsampratetemp = cell2mat(toeinfo{2});
     stimsamprate = str2double(stimsampratetemp(17:21));
%      stimsamprate = 41000;
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
alltoes = cell2mat(textscan(fid, '%.6f64 %*n'));

fclose(fid);
