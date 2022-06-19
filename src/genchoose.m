function [x,y, cs] = genchoose(NN)
%-----------------------------------------------------------------------------------------
% funtion [x,y] = genChoose(NN) allows you to choose NN number of boxes on a field of view
% and plots squares on whatever imagesc has been displayed before
% NN can go up to 14 and colors rotate after 7 colors. 
%-----------------------------------------------------------------------------------------

cmap = jet;
temp = linspace(1,size(cmap,1),NN);
cs = cmap(floor(temp),:);
% cs = [98 103 114 99 109 121  107 98 103 114 99 109 121  107 98 103 114 99 109 121  107 98 114 103 121 98 107 99 114 98 114 103 121 109 107 99 98 114 103 121 109 107 99 98 114 103 121 109 107 99 98 114 103 121 109 107 99 98 114 103 121 109 107 99 98 114];
% cs = [98 114 103 121 109 99 107 98 114 103 121 98 107 99 114 98 114 103 121 109 107 99 98 114 ];
%figure;
for i = 1:NN
    [x(2*(i-1)+1:i*2), y(2*(i-1)+1:i*2)] = ginput(2);
    hold on;
    plot([x(2*(i-1)+1) x(2*(i-1)+1)],[y(2*(i-1)+1) y(i*2)],'LineWidth',2,'color',cs(i,:))
    plot([x(i*2) x(i*2)],[y(2*(i-1)+1) y(i*2)],'LineWidth',2,'color',cs(i,:))
    plot([x(2*(i-1)+1) x(i*2)],[y(2*(i-1)+1) y(2*(i-1)+1)],'LineWidth',2,'color',cs(i,:))
    plot([x(2*(i-1)+1) x(i*2)],[y(i*2) y(i*2)],'LineWidth',2,'color',cs(i,:))   
end

