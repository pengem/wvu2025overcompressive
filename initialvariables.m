% this script provides the initial variables for the OvercompressiveRegion
% script. it also constructs a bunch of arrays. 

%% initial variables (scalar)
pbar = 3;
pL = 2;
uL = 3;
aexp = -1.5;


a_t = 0.1;

times_index = [2201];

t_0 = 0; % idk if this is actually necessary lol but whatever
t_end = 22; % idk lmfao
t_step = 0.01; % idk either lol

% idk if this should go here or if everything below this should go in the
% main script? but everything below here MUST go together.
t_values = 0:t_step:t_end;

%% initialize right state array

p_min = 0;
p_max = 20;
dp = 0.2;

u_max = 10;
u_min = -10;
du = 0.2;

p_vals_all = p_min:dp:p_max; 
u_vals_all = u_max:-du:u_min;

[P, U] = meshgrid(p_vals_all, u_vals_all);  % grid of p and u

puRAll = cat(3, P, U);  % % initializing the right states

% get points that aren't on the axes
u_vals_index = find(u_vals_all ~= 0);
p_vals_index = find(p_vals_all ~= 0);

% recombine for points not on axes. note that it the up switches because of
% matlab's row first operations. 
puR = puRAll(u_vals_index, p_vals_index, :);

% the actual u and p values minus the axes
p_vals = puRAll(:,:,1);
u_vals = puRAll(:,:,2);

Nu = length(u_vals_index);  % rows
Np = length(p_vals_index);  % columns

%% initialize array to store values of udelta
% basically the three nested for loop will produce a udelta value at a
% specific t at a specific right state. 

udelta = zeros(Nu, Np, length(t_values));

%% upperbound calculation
% calculating lambda1 and lambda2 for unchanging left state
lambdaL1 = -(uL*( ( (1+aexp)*( (pL/pbar)^aexp ) ) - 1) );
lambdaL2 = -(uL*( ( (pL/pbar)^aexp ) - 1) );

upper_bound = min(lambdaL1, lambdaL2);


%% bounds
lowerbounds = zeros(Nu, Np);

%% initializing array of 0 and 1s for overcompressive region
% same size as udelta. there's a 1 for points that satisfy the
% overcompressive region and 0 for points that don't

YNovercompressive_ind = zeros(Nu, Np, length(t_values));

%% grid for all omega values and omega prime
% already presets t = 0 to be 0 for all points
w = zeros(Nu, Np, length(t_values));
w_prime = zeros(Nu, Np, length(t_values)-1);


% this is just for the case number for naming files and such lol
caseNum = 0;

if aexp < -1
    if uL < 0
        if pL < pbar
            caseNum = 1; % case 1
        elseif pL == pbar
            caseNum = 2; % case 2
        else
            caseNum = 3; % case 3
        end % pL > pbar
    else % uL > 0
        if pL < pbar
            caseNum = 4; % case 4
        elseif pL == pbar
            caseNum = 5; % case 5
        else % pL > pbar
            caseNum = 6; % case 6
        end
    end
elseif aexp == -1
    if uL < 0
        if pL < pbar
            caseNum = 7; % case 7
        elseif pL == pbar
            caseNum = 8; % case 8
        else % pL > pbar
            caseNum = 9; % case 9
        end
    else % uL > 0
        if pL < pbar
            caseNum = 10; % case 10
        elseif pL == pbar
            caseNum = 11; % case 11
        else % pL > pbar
            caseNum = 12; % case 12
        end
    end
elseif aexp > -1 && aexp < 0
    if uL < 0
        if pL < pbar
            caseNum = 13; % case 13
        elseif pL == pbar
            caseNum = 14; % case 14
        else % pL > pbar
            caseNum = 15; % case 15
        end
    else % uL > 0
        if pL < pbar
            caseNum = 16; % case 16
        elseif pL == pbar
            caseNum = 17; % case 17
        else % pL > pbar
            caseNum = 18; % case 18
        end
    end
elseif aexp > 0
    if uL < 0
        if pL < pbar
            caseNum = 19; % case 19
        elseif pL == pbar
            caseNum = 20; % case 20
        else % pL > pbar
            caseNum = 21; % case 21
        end
    else % uL > 0
        if pL < pbar
            caseNum = 22; % case 22
        elseif pL == pbar
            caseNum = 23; % case 23
        else % pL > pbar
            caseNum = 24; % case 24
        end
    end
else
    error('error lol');
end