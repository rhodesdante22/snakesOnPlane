%% Source Sink Method
% Preston Zimprich

clc; clearvars; 
%%

Tables = load("fileListTables.mat");
Tracings = load("tracingsListTable.mat");
Tracings = Tracings.tracingsListTables;

%P = [main intensity threshold, keep threshold, source
%strength, main sink strength, bottom sup sinks, mid sup sinks, top sup
%sinks, v_inf coefficient, top gap, bottom gap, spacing for top, spacing for bottom]

%P = [12, 0.75, 30, 80, 550, 50, 40, 40, 0.65, 6/7, 3/7, 1/4, 3/4];
P = [0.039, 0.3, 7, 15, 3, 3, 3, 0.009, 4/7, 3/7, 2/5, 3/4];
%P = [-1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000];

for i = 1:10   
    [main_sink_loc{i}, main_source_loc{i}] = find_points(P, Tracings{i}); 
end

vid = VideoReader("1.avi");  
main_sink_loc = main_sink_loc{1}; 
main_source_loc = main_source_loc{1};

%%
FN = 1:1:50; 
M = zeros(vid.NumFrames*50, 4); 



for m = 1:vid.NumFrames
m

%vid = VideoReader("00001.avi");
frame = read(vid,m); 

%main_sink_loc, main_source_loc

[contour] = source_sink(P, m, frame, main_sink_loc, main_source_loc);
outArr = xInterp(contour, 50);
ImN = zeros(50,1) + m; 

M = [M; ImN, FN', outArr]; 

area = find_area(outArr);
longAxis = findLongAxis(outArr);
length = sqrt((longAxis(1,1)-longAxis(2,1))^2+(longAxis(1,2)-longAxis(2,2))^2);

vol(m) = dodgeVolume(area, length);

end


%%

figure
plot(vol)
ylim([0 80000]); 
xlabel('Frame'); ylabel('Area');

%%
% trace1 = Tracings{2}; 
% trace1 = [trace1(:,2), trace1(:,3)]; 
% trace1 = table2array(trace1); 
% 
% figure 
% imshow(frame)
% hold on
% plot(trace1(:,1), trace1(:,2)); 
% plot(traces)
% hold off