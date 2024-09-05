%%This script allows you to run the simulink model of the fish
%%
NumOfSegs=5;                    %Number of segments you want to split the fish into
bsi=NumOfSegs-1;                %The head segment is unchanged through this whole process
head_ogn=57;                    %Original head length of the fish model
fishlen= 363.5;                  %Body length. Total length of the fish
fish_model_from_head=59.1375;   %Calculated by calculating the length from the tip of the fish's head to the second equidistant point and scaling it to the model
head=fishlen*headsize_factor;    %How much of the image qualifies as the head
bodlen=fishlen-head;             %Just the body length, no head
seglen=bodlen/bsi;              %Segment length calculation

%% Generic Params
ground_damping = 1e3;   
ground_mu_k = 0.3;
ground_mu_k_legs = 0.5*20;
ground_mu_s= 0.5*0.5;
ground_mu_s_legs = 0.7*20;
ground_stiffness= 1e3;
critical_velocity=1e-3; 
leg_critical_velocity=1e-3;

%% Mass Handling, all variables relate to the fish in the video
fish_length_mm=109;          %fish length in mm.
fish_mass_g=7.5;            %in grams
fish_width_half_trunk=(fish_length_mm*0.6/7.2)/2;  %width in mm
fish_height_half=(fish_length_mm/8.9)/2;        %height in mm
fish_cross_area_trunk= fish_width_half_trunk*fish_height_half*pi; %in mm^2

%% Handling of the body vars: This code creates a hexagonal cross-section for the robot in attempt to scale it evenly in all directions
lsf= (fishlen)/fish_length_mm;
ab_factor=1/3;
var_a= fish_width_half_trunk*2*lsf;
var_c= fish_height_half*2*lsf;
var_b= var_a- (fish_cross_area_trunk)*(lsf^2)/var_c;

%Co-ordinates for the hexagons, inside each fish block: 
%[0 var_c/2; var_b var_c; var_a-var_b var_c; a var_c/2; var_a-var_b 0; var_b 0]
area=var_c*(var_a-var_b);

%% Density per unit length calculation
vid_bod_to_mass= (fish_mass_g/fish_length_mm); %g per mm or kg/m
model_whole_mass= (fishlen)*vid_bod_to_mass; %g per mm

%----------
%Allocation
%-----------
%UDM: uniformly distributed mass
head_dense= vid_bod_to_mass*area/fish_cross_area_trunk;
trunk_dense=vid_bod_to_mass*area/fish_cross_area_trunk;
tail_dense=vid_bod_to_mass*area/fish_cross_area_trunk;

%In the event that a user would like to alter the density distribution of the fish, please find attached, below, the following function (you can define your
%own boundaries as to which segment qualifies as the start of the head, trunk, tail)

%[head_dense,trunk_dense,tail_dense] = DensityFinder(segnum, headportion,trunkportion,tailportion, vid_bod_to_mass, head, model_whole_mass, seglen);

%% Extras
%Fish input torque gains (the amplitude is=2)
Fish_Gain=-[1,1,1,1,1,1,1,1,1,1, 1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
%Segment:
%1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30    -For ease of view

%Set the initial bias of each segment if desired, else set it to 0
D=zeros(NumOfSegs);

%Calculate the simulation time to match your video
t_diff=1/frameRate;
sim_time=t_diff*numel(midlineFiles);
disp('run sim for');
disp(sim_time);





