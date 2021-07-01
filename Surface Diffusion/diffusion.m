%% 1D
% this plots the probability distribution from distance traveled after insertion
% it uses 2* the gaussian pdf with a variance (sigma^2) for 1D diffusion =
% 2*D*T to account for both negative and positive from the origin and
% integrates from 0 - infinity.


binsz = 0.001; % size for estimation - smaller is better - could also use the trapezoid rule to be more accurate
d = 0.06; % diffusion coefficient in um^2/s - from DIV12-17DIV12-17 mobile receptor pool from Borgdorff and Choquet, as well as the Petrini Neuron paper i.e. average values for extrasynaptic receptors in Fig. 1D mobile receptor pool from Borgdorff and Choquet, as well as the Petrini Neuron paper
delT = 69.3; % time for 50% of receptors to be bound to antibody :  binding rate of the x-linking antibody (0.010 sec-1) from the tau value (1.66 min) of the single exponential fit in Fig. S2B. Thus, the mean squared displacement (MSD) of receptors over the time interval between insertion and antibody binding to 50% of the receptors (69.3 sec) 
delT = [30,60,90,120]; % or plot a range of delta T based on range of possible values
nsigm = 4; % number of sigma out to make plot

figure; hold on;
for ll = 1:numel(delT)
    sigm = sqrt(2*d*delT(ll));
    N = ceil(sigm*nsigm/binsz);
    x = [0:N]*binsz;
    
    px = (2/(sigm*sqrt(2*pi))) * exp(-x.^2/(2*sigm^2));
    displayname_curve = ['D = ' num2str(d) ' \mum^2/s, T=' num2str(delT(ll)) 'sec'];
    plot(x,px,'LineWidth',4,'DisplayName',displayname_curve);
    set(gca,'FontSize',16)
end
title('1D Diffusion');
xlabel('Distance from insertion point (\mum)');
ylabel('Probability');
legend
check_integral = sum(px)*binsz %should be 1
%% use this for single delT
expected_displacement = sum(px*binsz .* x)
displayname_expval = ['Expected Displacement = ' num2str(expected_displacement) '\mum'];
plot([expected_displacement,expected_displacement],[0, max(px)],'k','DisplayName',displayname_expval);
legend()