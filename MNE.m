function [A_mean,J_mean,h_mean] = MNE(stimulus, response, MNE_params)

% Run the MNE analysis with given stimulus and response

% MNE_params is a struct with fields Ndim, Nlags,order, fittype, Njack

% if stimuli are images, then Ndim is the total amount of pixels as the image will have to be in a vector.
% For spectrograms, Ndim will be number of frequency bands

    Ndim = MNE_params.Ndim; % Number of frequency bands in STRF or pixels in image
    Nlags = MNE_params.Nlags; % Nlags will be the number of timepoints in the STRF for spectrograms
    % Set to 1 if using images (=spatial model)
    order   = MNE_params.order;   % order of MNE model to fit
    fittype = MNE_params.fittype;   % 0 for regular fitting, 1 for random fitting
    njack   = MNE_params.Njack;   % # jackknives to run (also determines the size of each jackknife)
    Nd = sqrt(Ndim); % no. of pixels per side of image
    [~,Nsamples] = size(stimulus);
    
    % redefine stimulus to include time lags
    if Nlags>1
        Nsamples = Nsamples - (Nlags-1); 
        Ndimtotal = Ndim*Nlags;
        stim_ = zeros(Ndimtotal, Nsamples);
        for i=1:Nlags
            stim_(Ndim*(i-1)+1:Ndim*i,:) = stimulus(:,i:Nsamples+i-1);
        end
    else
        stim_ = stimulus;
    end
        
    %Response
    resp = response;
    resp_ = resp(Nlags:length(resp));
    
    h=[];
    J=[];
    pfinal_all_jacks = [];
    
    % MNE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % break into training and test sets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for jack = 1:njack %loop over all njacks
        stim = stim_;
        resp = resp_;
        Ntest = floor(size(stim)/njack);
        teststim = stim(:, 1+(jack-1)*Ntest:jack*Ntest);
        testresp = resp(1+(jack-1)*Ntest:jack*Ntest);
        stim(:, 1+(jack-1)*Ntest:jack*Ntest) = [];
        resp(1+(jack-1)*Ntest:jack*Ntest) = [];

        [Ndim_total, Nsamples] = size(stim);
        psp = mean(mean(resp));   % spike probabilisty
        avg = (stim*resp')/Nsamples;  
        avg = mean(avg,2);  % Ndim x 1
        avgs = [psp;avg];
        if order>1
            avgsqrd = stim*(repmat(resp',[1,Ndim_total]).*stim')/Nsamples;  % Ndim x Ndim
            avgsqrd = reshape(avgsqrd,[Ndim_total^2,1]);
            avgs = [avgs;avgsqrd];
        end
        
        % Initialize parameters
        pstart = log(1/avgs(1)-1);
        pstart(2:Ndim_total+1) = .001*(2*rand([1,Ndim_total])-1);
        if order>1
            temp = .0005*(2*rand([Ndim_total,Ndim_total])-1);
            pstart(Ndim_total+2:length(pstart)+Ndim_total^2) = reshape((temp+temp'),[1,Ndim_total^2]);

            clear temp;
        end
        
        disp('Starting optimization');
        tic
        % Run conjugate gradient algorithm
        pfinal = frprmn(pstart, @logloss, @dlogloss, stim', resp', teststim', testresp', order, avgs, Nd, fittype);
        disp(['Optimization took ' num2str(toc/60) ' minutes']);
        %And now plot the results and save the figures
        close all

        a=pfinal(1);
        h=[h;pfinal(2:Nlags*Ndim+1)];
        pfinal_all_jacks = [pfinal_all_jacks;pfinal];
        J=[J;pfinal(Nlags*Ndim+2:end)];
    end
    
    A_mean = mean(a);
    J_mean = mean(J);
    h_mean = mean(h);
end
