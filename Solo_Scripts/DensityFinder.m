function [head_dense,trunk_dense,tail_dense] = DensityFinder(segnum, headportion,trunkportion,tailportion, vid_bod_to_mass, head, model_whole_mass, seglen)
%   Returns the density variables for a designated three sub-sections of the fish
%'__' portion means the amount in the ratio that would be allocated to that section's density, e.g. if the head is three times as heavy as the other
%two sections, the ratio is 3:1:1 which means headportion=3, trunkportion=1, tailportion=1

switch segnum
    case 30
        l_trunk= 19*seglen; %length of the 'trunk' of the fish, assumed to be 19 segments long, change as needed
        l_tail=10*seglen;   %length of the 'tail' of the fish, assumed to be 10 segments long, change as needed
    
    case 15
        l_trunk= []*seglen;
        l_tail=[]*seglen;
    
    case 10
        l_trunk= []*seglen;
        l_tail=[]*seglen;
    
    case 5
        l_trunk= []*seglen;
        l_tail=[]*seglen;
    
    case 3
        l_trunk= []*seglen;
        l_tail=[]*seglen;
    
end

denom=headportion+trunkportion+tailportion;
head_dense=vid_bod_to_mass*(headportion/denom);
trunk_dense=vid_bod_to_mass*(trunkportion/denom);
tail_dense=vid_bod_to_mass*(tailportion/denom);
mass_new= head*head_dense+l_trunk*trunk_dense+l_tail*tail_dense;
re_ratio= model_whole_mass/mass_new;
head_dense=head_dense*re_ratio;
trunk_dense=trunk_dense*re_ratio;
tail_dense=tail_dense*re_ratio;

end