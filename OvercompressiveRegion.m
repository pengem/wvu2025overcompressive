% OvercompressiveRegion Script

% this script runs for a(t) = b where b is any constant. 
% it computes many arrays. we start with a fixed left
% state (uL and pL) and then generate an array that has different values of
% uR and pR (so nxmx2 where nxmx1 are pR and behind it is uR). then, for
% each of these nxn elements, we then do a fun euler thing. 

clear all;
clc; 
close all;
tic
% call initial varibales
initialvariables;

for i = 1:Nu

    for j = 1:Np 

        pR = puR(i, j, 1);
        uR = puR(i, j, 2);

        % solve for roots by eta, aka udelta(0). 
        a = pL - pR;
        b = ((pR*uR) - (pR*uR*((pR/pbar)^aexp))) - ((pL*uL) - (pL*uL*((pL/pbar)^aexp))) + (pR*uR) - (pL*uL);
        c = ((pL*(uL^2)) - (pL*(uL^2)*((pL/pbar)^aexp))) - ((pR*(uR^2)) - (pR*(uR^2)*((pR/pbar)^aexp)));
        % twoeta = roots([a, b, c]);

        % only take the plus sign thing in the quad formula
        % note that this should be changed into udelta first thing bc t=0
        % quad becomes a linear function
        if a == 0
            eta = -c/b;
        else
            eta = (-b - sqrt((b^2) - (4*a*c)))/(2*a);
        end

        % matlab will say that a complex value is in between a real bound
        % so we need to get rid of the complex numbers
        if ~isreal(eta)
            eta = NaN;
        end
        
        % store udelta(0) into the first slice of the udelta array
        udelta(i,j,1) = eta;
        
        % calculating lambda1 and lambda2 for right state
        lambdaR1 = -((uR)*(((1+aexp)*((pR/pbar)^aexp))-1));
        lambdaR2 = -((uR)*(((pR/pbar)^aexp)-1));

        % note that upperbound does not change, so it's already calculated
        % in initial variables. might move to main script tho, idk what the
        % convention is lol
        lower_bound = max(lambdaR1, lambdaR2);

        % calculating lambda1 and lambda2 for unchanging left state
        lambdaL1 = -(uL*( ( (1+aexp)*( (pL/pbar)^aexp ) ) - 1) );
        lambdaL2 = -(uL*( ( (pL/pbar)^aexp ) - 1) );

        upper_bound = min(lambdaL1, lambdaL2);

        if (udelta(i,j,1) > lower_bound) && (udelta(i,j,1) < upper_bound)
            YNovercompressive_ind(i,j,1) = 1;
        else
            YNovercompressive_ind(i,j,1) = 0;
        end

        % calculating times now!
        for t = 1:length(t_values)-1
        % this is going to get confusing because i will be using stuff like
        % udelta at time 0, but that is stored in udelta(1). indexing in
        % matlab is kinda stupid!
           
        % also note that i have embedded the euler approximation into this
        % script instead of creating a seperate function for both omega and
        % udelta. this is just because i got lazy. their interdependence on
        % each other kind of complicates it a bit but maybe I will change 
        % it later. 

        % euler approximation for omega
            w_prime(i,j,t) = ((pR-pL) * udelta(i,j,t)) + ((pL * (uL + a_t*t_values(t)))*(1 - (pL/pbar)^aexp)) - ((pR * (uR + a_t*t_values(t)))*(1 - (pR/pbar)^aexp)); 
            w(i,j,t+1) = w(i,j,t) + (t_step * w_prime(i,j,t));
            
        % euler approximation for udelta
            if w(i,j,t) == 0
                udelta_prime_temp = 0;
            %elseif  NEED TO ADD THE THIRD CONDITION HERE

            else
                udelta_prime_temp = (1/w(i,j,t)) * (((pR*uR - pL*uL)*udelta(i,j,t)) + (pL*uL*(uL + a_t*t_values(t))*(1-(pL/pbar)^aexp)) - (pR*uR * (uR + a_t*t_values(t))*(1-(pR/pbar)^aexp)) - (w_prime(i,j,t)*udelta(i,j,t)));
            end

            udelta(i,j,t+1) = udelta(i,j,t) + t_step * udelta_prime_temp;


            % calculating lambda1 and lambda2 for right state
            lambdaR1 = -((uR+(a_t*t_values(t)))*(((1+aexp)*((pR/pbar)^aexp))-1));
            lambdaR2 = -((uR+(a_t*t_values(t)))*(((pR/pbar)^aexp)-1));
            lower_bound = max(lambdaR1, lambdaR2);

            % calculating lambda1 and lambda2 for unchanging left state
            lambdaL1 = -((uL+(a_t*t_values(t)))*( ( (1+aexp)*( (pL/pbar)^aexp ) ) - 1) );
            lambdaL2 = -((uL+(a_t*t_values(t)))*( ( (pL/pbar)^aexp ) - 1) );
            upper_bound = min(lambdaL1, lambdaL2);

            if (udelta(i,j,t+1) > lower_bound) && (udelta(i,j,t+1) < upper_bound)
            YNovercompressive_ind(i,j,t+1) = 1;
            else
            YNovercompressive_ind(i,j,t+1) = 0;
            end

        end

    end

end


% call graphs
graphs_time;
toc