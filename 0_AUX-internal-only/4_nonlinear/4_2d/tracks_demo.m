close all
clear all
clc



% the code requires the following structure for a 2-sided Feynman diagram

number_of_mu = 4;

% make empty fields, empty = NaN
track.states = zeros(number_of_mu +1,2) + NaN;
track.interaction.column = zeros(number_of_mu +1,1) + NaN;
track.interaction.direction = zeros(number_of_mu +1,1) + NaN;
track.phasematch = zeros(number_of_mu +1,1) + NaN;

% index in first dimension is "running number of interaction + 1"
% index in second dimension is column: 1= left = ket, 2 = right = bra
% value is number of state = quanmtum number, but >=1 
% this number is reduced by 1 for display !

% we start from ground state 
track.states(1,:) = 1;  %


% 1. interaction: absorption into left
track.interaction.column(2)  = 1;  % = left
track.interaction.direction(2)  = 1;  % = absorption
track.phasematch(2)  = 1;  % left-to-right
track.states(2,1) = 2;  % incremented quantum number
track.states(2,2) = 1;

% 2. interaction: absorption into right
track.interaction.column(3)  = 2;  % = right
track.interaction.direction(3)  = 1;  % = absorption
track.phasematch(3)  = -1;  % right-to-left
track.states(3,1) = 2;  % incremented quantum number
track.states(3,2) = 2;

% 3. interaction: emission from right 
track.interaction.column(4)  = 2;  % = right
track.interaction.direction(4)  = -1;  % = emission
track.phasematch(4)  = 1;  % left-to-right
track.states(4,1) = 2;  
track.states(4,2) = 1;  % decremented quantum number

% 3. interaction = radiation polarization: emission from left 
track.interaction.column(5)  = 1;  % = left
track.interaction.direction(5)  = -1;  % = emission
track.phasematch(5)  = -1;  % right-to-left
track.states(5,1) = 1;  % decremented quantum number
track.states(5,2) = 1; 

% of course on better does this with a script. Your turn ...
% we get an array of possible diagrams
% here for demo only one

tracks = [track];


% display all tracks, grouped by phase matching
for id1=1:2
    for id2=1:2
        for id3=1:2
            pm = ([id1; id2; id3] - 1.5) .* 2;
            display_tracks_by_phasematch(tracks,pm)
        end
    end
end




% calculate spectra in frequency domain

% define frequency, damping and transition dipole moment of transitions
% between states
% mirroring of matrix sub-traingles is taken into account by formulas

omega = zeros(3) + NaN;
omega(1,2) = 2.370;  % 2 pi nu in 1/fs
omega(1,3) = 2.42;
omega(2,3) = -1 .* (omega(1,2) - omega(1,3));

gamma0 = 0.005;
gamma = zeros(3) + gamma0;


mu = zeros(3) + NaN;
mu(1,2) = 2.99;
mu(1,3) = 4.227 .* sqrt(1.15);

t2 = 100; % fs

w1f = linspace(2.32,2.46,500);
w3f = linspace(2.32,2.46,500);

[omega3, omega1] = meshgrid(w3f, w1f);


%---------------------------------------
% non-rephasing

pm_nr = [+1; -1 ; +1];


trs = tracks(select_by_phasematch(tracks, pm_nr));
for id = 1:length(trs)
    formula = make_frequency_domain_formula(trs(id));
    frequency_nr_h(id,:,:) = eval(formula);
end

if (length(trs) > 1) % we have more than one diagram: sum over it
    spec_f_nr =   squeeze(sum(frequency_nr_h(:,:,:)));
else
    spec_f_nr =   squeeze((frequency_nr_h(:,:,:)));
end

imagesc(w3f, w1f, real(spec_f_nr))
set(gca,'YDir','normal')
title('real non-rephasing ')

