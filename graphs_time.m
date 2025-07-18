% graphs_time Script
 
% times_index (in initialvariables) takes in the INDEX of the time vector. 
% thus, times_index = 1 corresponds to the first value in the t_values 
% vector. i may change this so that you can put in any t value and it will 
% choose the t-value closest to generate. 
% anyway times_index takes in a vector and will generate mulitple graphs
% depending on how many values are in the vector.

drawnow
for i = 1:length(times_index)
    fig=figure(i);
    
    [OCu_index, OCp_index] = find(YNovercompressive_ind(:,:,i)==1);

    % then need to use those indicies to get the correct points that are
    % overcompressive
    pOvercompressive = zeros(length(OCp_index));
    uOvercompressive = zeros(length(OCu_index));
    
    % for loop is needed because if this isn't used, matlab just calculates the
    % combination of every index possible, giving us a square... lol
    for j = 1:length(OCu_index)
        pOvercompressive(j) = puR(OCu_index(j), OCp_index(j), 1);
        uOvercompressive(j) = puR(OCu_index(j), OCp_index(j), 2);
    end

    scatter(pOvercompressive, uOvercompressive, 5, 'filled', 'MarkerFaceColor', '#C05780')
    xlabel('\rho'), ylabel('u')
    title(['Case ', num2str(caseNum), ': Overcompressive points at t = ', num2str(t_values(times_index(i))), ', a = ', num2str(aexp)])


    % create figure and get curves for shocks
    plotX = linspace(0,20,200);
    plotX1 = linspace(0,pbar,200*pbar/20);
    plotX2 = linspace(pbar,20,200-(200*pbar/20));
    
    u1 = uL .* ones(size(plotX));
    u2 = (((pL^aexp) - (plotX.^aexp))./((plotX.^aexp) - (pbar^aexp))) .* (uL + (a_t.*t_values(times_index(i)))) + uL;

    u2first = (((pL^aexp) - (plotX1.^aexp))./((plotX1.^aexp) - (pbar^aexp))) .* (uL + (a_t.*t_values(times_index(i)))) + uL;
    u2second = (((pL^aexp) - (plotX2.^aexp))./((plotX2.^aexp) - (pbar^aexp))) .* (uL + (a_t.*t_values(times_index(i)))) + uL;

    % plot up curves and points
    hold on;
    
    % u2 becomes a vertical line when pbar = pL
    plot(plotX,u1,'-k', 'LineWidth',1);
    if pbar == pL
        xline(pbar, '-k', 'LineWidth',1);
    else
        if caseNum ==4
            plot(plotX1, u2first,'-k', 'LineWidth',1);
            plot(plotX2, u2second,'--k', 'LineWidth',1);
        elseif caseNum == 1
            plot(plotX1, u2first,'-k', 'LineWidth',1);
            plot(plotX2, u2second,'--k', 'LineWidth',1);
        elseif caseNum == 6
            plot(plotX1, u2first,'--k', 'LineWidth',1);
            plot(plotX2, u2second,'-k', 'LineWidth',1);
        elseif caseNum == 3
            plot(plotX1, u2first,'--k', 'LineWidth',1);
            plot(plotX2, u2second,'-k', 'LineWidth',1);
        else
            plot(plotX,u2,'-k', 'LineWidth',1);
        end
    end

    % asymptote
    xline(pbar, '--k', 'LineWidth',1);
    
    %y-axis and x-axis are solid black lines
     xline(0);
     yline(0);
    
    ylim([-10,10])


    axis square

    hold off

    % saving the image
    nameFile = ['Case' num2str(caseNum) 'OvercompressiveRegion_a_t' num2str(a_t) '_t' num2str(t_values(times_index(i)))];

    saveas(fig,[nameFile '.png']);
end