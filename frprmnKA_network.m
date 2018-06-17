function [pbest, flist,ftestlist, its,tally, min_eigenvals_list,max_eigenvals_list] = frprmnKA(p, func, dfunc, stim, resp, teststim, testresp, order, avgs, Nd, fittype,Nlags,Nf,jack,njack)
min_eigenvals_list = []; max_eigenvals_list = [];
ITMAX = 1000;
fp = feval(func, p, stim, resp, order);
xi = feval(dfunc, p, stim, avgs, order);
exitCondition = 0;
g  = -xi;
h  = g;
xi = g;
besttest = 1000;
flist=[];
ftestlist=[];
tally = 0;

% Loop over iterations of minimization
for its=1:ITMAX,
    %disp(['Iteration ' num2str(its)]);
    disp(num2str(its)) ;
    [p, xi, fret] = dlinmin(p, xi, func, dfunc, stim, resp, order, avgs);
    flist(its)=fret;
    if fittype==0
        ftestlist(its)=feval(func, p, teststim, testresp, order);
    end
    
    
    figure(1)
    plot(flist)
    if fittype==0
        hold on
        plot(ftestlist,'r')
        hold off
    end
    drawnow

    
    if fittype==0
        
        if ftestlist(its)<besttest*1  || its<=2 % changed 0.9 to 1
            besttest = ftestlist(its);
            pbest = p;
            tally=0;
        else
            tally = tally+1;
            disp(num2str(tally)) ;
            
        end
        
        J=p(Nlags*Nf+2:end);
 [V,D] = eig(reshape(J,Nlags*Nf,Nlags*Nf));
 

        [eigenvals,inds] = sort((diag(D)));
        min_eigenval = min(eigenvals); max_eigenval = max(eigenvals); 
        disp(num2str([min_eigenval , max_eigenval]));
       
                min_eigenvals_list = [min_eigenvals_list;min_eigenval];
                                max_eigenvals_list = [max_eigenvals_list;max_eigenval];
                             
                                
        figure(2)
                plot(eigenvals,'o');
%         plot(eigenvals,'o','MarkerEdgeColor','r');
        axis([-10 330,-0.08,0.08]);
        hold on

%         text(50,-0.032,num2str(inds(1)));  text(50,-0.028,num2str(inds(2))); text(50,-0.024,num2str(inds(3))); text(50,-0.020,num2str(inds(4))); text(50,-0.016,num2str(inds(5)));  
%         text(250,0.032,num2str(inds(end))); text(250,0.028,num2str(inds(end-1))); text(250,0.024,num2str(inds(end-2))); text(250,0.020,num2str(inds(end-3))); text(250,0.016,num2str(inds(end-4)))
          
text(50,-0.068,num2str(inds(1)));  text(50,-0.060,num2str(inds(2))); text(50,-0.052,num2str(inds(3))); text(50,-0.044,num2str(inds(4))); text(50,-0.036,num2str(inds(5)));  
text(250,0.068,num2str(inds(end))); text(250,0.060,num2str(inds(end-1))); text(250,0.052,num2str(inds(end-2))); text(250,0.044,num2str(inds(end-3))); text(250,0.036,num2str(inds(end-4)))
        
                                
% refline(0,0.0121)
% refline(0,-0.0120)
refline(0,0.0147)
refline(0,-0.0148)

hline = refline(0,0.0305);
set(hline,'Color','r')
hline = refline(0,-0.0302);
set(hline,'Color','r')


hold off

figure_name = ['eigenvalues_noise_to_data_step_' , num2str(its), '_jack_' num2str(jack) '_of_' num2str(njack)];

hgsave(figure_name)
 drawnow
       
        
        if tally== 100 || its==400 %THIS IS THEY KEY EXIT CONDITION (tally default set was 10)
            disp('min of test set found');
            exitCondition = 1;
            break;
        end
        
    else
%         J = reshape(p(52:2551),50,50);
        [evecs,evals]=eig(J);
        [EV,inds] = sort((diag(evals)));
        disp(num2str([min(EV) , max(EV)]));
        if its==200
            pbest = p;
            disp('stopping algorithm');
            exitCondition = 1;
            break;
        end
    end
    
    xi = feval(dfunc, p, stim, avgs, order);
    gg = sum(g.^2);
    dgg = sum( (xi + g).*xi );   % This statement for Polak-Ribiere
    % dgg = sum( xi.^2);         % This statement for Fletcher-Reeves
    if gg == 0,            % Unlikely.  If gradient is exactly zero then
        exitCondition = 2;   % we are already done.
        disp('Gradient equal to zero, exiting frprmn.');
        break;
    end
    gam = dgg/gg;
    g = -xi;
    h = g + gam.*h;
    xi = h;
end
if exitCondition == 0,
    disp('Too many iterations in frprmn');
end