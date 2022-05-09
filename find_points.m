function [main_sink_loc, main_source_loc] = find_points(P, Tracings)

% Used to find initial reference points for each video

init_refpoints = [Tracings(1,2), Tracings(1,3);
    Tracings(end/2,2), Tracings(end/2,3)]; 
init_refpoints  = table2array(init_refpoints); 
init_refpoints_distx = abs(init_refpoints(1,1)-init_refpoints(2,1)); 
init_refpoints_disty = abs(init_refpoints(1,2)-init_refpoints(2,2));

main_sink_loc = [init_refpoints(1,1)+P(11)*init_refpoints_distx, 
    init_refpoints(1,2)+P(11)*init_refpoints_disty]; 
main_source_loc = [init_refpoints(1,1)+P(12)*init_refpoints_distx, 
    init_refpoints(1,2)+P(12)*init_refpoints_disty];