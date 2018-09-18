function fig = LG_6DoFAnimation_Tra( varargin )
%LG_6DOFANIMATION Summary of this function goes here
%   Detailed explanation goes here
%% Inputs:
pQ_total = varargin{1};
R_total  = varargin{2};
pP_total = varargin{3};
SamplePlotFreq = varargin{4};
Trail = varargin{5};
CreateAVI = varargin{6};
isFixView = varargin{7};
%isTensionTotal = varargin{8};

[numSamples dummy] = size(pQ_total);% total number of samples
%% Default configuration parameters:
FullScreen = false;
% Trail = 'DotsOnly';
AVIfileName = 'Demo';
AVIfileNameEnum = true; 

ShowArrowHead = true;
ShowLegend = false;
Title = 'Quadcopter --';
Position = [50 50 866 600];
Spin = 120;
View = [(120:(Spin/(length(pQ_total)-1)):(120+Spin))', 25*ones(length(pQ_total), 1)];
AxisLength = 0.6;
LimitRatio = 1;
Xlabel = 'X(m)';
Ylabel = 'Y(m)';
Zlabel = 'Z(m)';

%isTension = isTensionTotal(1:SamplePlotFreq:numSamples, :);
pQ = pQ_total(1:SamplePlotFreq:numSamples, :);
R = R_total(:, :, 1:SamplePlotFreq:numSamples) * AxisLength;
pP = pP_total(1:SamplePlotFreq:numSamples, :);
if(numel(View) > 2)
    View = View(1:SamplePlotFreq:numSamples, :);
end
[numPlotSamples dummy] = size(pQ);
    
%% Setup AVI file

    aviobj = [];                                                            	% create null object
    if(CreateAVI)
        fileName = strcat(AVIfileName, '.avi');
        if(exist(fileName, 'file'))
            if(AVIfileNameEnum)                                              	% if file name exists and enum enabled
                i = 0;
                while(exist(fileName, 'file'))                                  % find un-used file name by appending enum
                    fileName = strcat(AVIfileName, sprintf('%i', i), '.avi');
                    i = i + 1;
                end
            else                                                                % else file name exists and enum disabled
                fileName = [];                                                  % file will not be created
            end
        end
        if(isempty(fileName))
            sprintf('AVI file not created as file already exists.')
        else
%             aviobj = avifile(fileName, 'fps', AVIfps, 'compression', 'Cinepak', 'quality', 100);
              aviobj = VideoWriter(fileName);
              open(aviobj);
        end
    end    
%%
fig = figure('Name', '6DOF Animation(Email:gmh_njust@163.com)');
    if(FullScreen)
        screenSize = get(0, 'ScreenSize');
        set(fig, 'Position', [0 0 screenSize(3) screenSize(4)]);
    elseif(~isempty(Position))
        set(fig, 'Position', Position);
    end
    %set(gca, 'drawmode', 'fast');
    lighting phong;
    set(gcf, 'Renderer', 'zbuffer');
    hold on;
    axis equal;
    grid on;
    view([View(1,1),View(1,2)]);
%     title(i);
    xlabel(Xlabel);
    ylabel(Ylabel);
    zlabel(Zlabel);
    
   % Create plot data arrays
    if(strcmp(Trail, 'DotsOnly') || strcmp(Trail, 'All'))
        xQ = zeros(numPlotSamples, 1);
        yQ = zeros(numPlotSamples, 1);
        zQ = zeros(numPlotSamples, 1);
        xP = zeros(numPlotSamples, 1);
        yP = zeros(numPlotSamples, 1);
        zP = zeros(numPlotSamples, 1); 
    end
    if(strcmp(Trail, 'All'))
        xQ = zeros(numPlotSamples, 1);
        yQ = zeros(numPlotSamples, 1);
        zQ = zeros(numPlotSamples, 1);
        xP = zeros(numPlotSamples, 1);
        yP = zeros(numPlotSamples, 1);
        zP = zeros(numPlotSamples, 1);
        ux = zeros(numPlotSamples, 1);
        vx = zeros(numPlotSamples, 1);
        wx = zeros(numPlotSamples, 1);
        uy = zeros(numPlotSamples, 1);
        vy = zeros(numPlotSamples, 1);
        wy = zeros(numPlotSamples, 1);
        uz = zeros(numPlotSamples, 1);
        vz = zeros(numPlotSamples, 1);
        wz = zeros(numPlotSamples, 1);
    end
    xQ(1) = pQ(1,1);
    yQ(1) = pQ(1,2);
    zQ(1) = pQ(1,3);  
    xP(1) = pP(1,1);
    yP(1) = pP(1,2);
    zP(1) = pP(1,3);  
    
    oxQ(1) = pQ(1,1);% Quadcopter initial position
    oyQ(1) = pQ(1,2);
    ozQ(1) = pQ(1,3);
    oxP(1) = pP(1,1);% payload initial position
    oyP(1) = pP(1,2);
    ozP(1) = pP(1,3);
    
  
%     ratio_quiver = 1;
%     R = ratio_quiver.*R;
    ux(1) = R(1,1,1:1);% each column vector of R is a direction vector
    vx(1) = R(2,1,1:1);
    wx(1) = R(3,1,1:1);
    uy(1) = R(1,2,1:1);
    vy(1) = R(2,2,1:1);
    wy(1) = R(3,2,1:1);
    uz(1) = R(1,3,1:1);
    vz(1) = R(2,3,1:1);
    wz(1) = R(3,3,1:1); 

    % Create graphics handles
%     subplot(2,1,1);
    orgHandle = plot3(oxQ, oyQ, ozQ, 'k.'); % mass center of the quadcopter
    PorgHandle = plot3(oxP,oyP,ozQ,'g-');
    if(ShowArrowHead)
        ShowArrowHeadStr = 'on';
    else
        ShowArrowHeadStr = 'off';
    end
    
    
    quivXhandle = quiver3(oxQ, oyQ, ozQ, ux, vx, wx,  'r', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');
    quivYhandle = quiver3(oxQ, oyQ, ozQ, uy, vy, wy,  'g', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');
    quivZhandle = quiver3(oxQ, oyQ, ozQ, uz, vz, wz,  'b', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');
    
    len = 0.34; % distance from the center of rotor 1 to rotor 3; 
    x1 = [oxQ,oyQ,ozQ] + 0.5*len*[ux,vx,wx];
    x3 = [oxQ,oyQ,ozQ] - 0.5*len*[ux,vx,wx];
    x2 = [oxQ,oyQ,ozQ] + 0.5*len*[uy,vy,wy];
    x4 = [oxQ,oyQ,ozQ] - 0.5*len*[uy,vy,wy];
    QX1 = [x1(1),x3(1)];
    QX2 = [x1(2),x3(2)];
    QX3 = [x1(3),x3(3)];
    QY1 = [x2(1),x4(1)];
    QY2 = [x2(2),x4(2)];
    QY3 = [x2(3),x4(3)];
    
    
    eb1 = [ux,vx,wx];
    eb2 = [uy,vy,wy];
    eb3 = [uz,vz,wz];
    temp = 0:(2*pi)/15 :2*pi;
    r0 = 0.04;
    xr1 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x1;
    xr2 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x2;
    xr3 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x3;
    xr4 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x4;
    
    % Create legend
    if(ShowLegend)
        legend('Origin', 'X', 'Y', 'Z');
    end
    lineQXhandle = line('XData',QX1,'YData',QX2,'ZData',QX3,'LineWidth',2,'Color',[1 0 0]);
    lineQYhandle = line('XData',QY1,'YData',QY2,'ZData',QY3,'LineWidth',2,'Color',[1 0 0]);
    circle1handle = line('XData',xr1(:,1),...
        'YData',xr1(:,2),...
        'ZData',xr1(:,3),...
        'LineWidth',2,'Color',[1 0 0]);
    
        circle2handle = line('XData',xr2(:,1),...
        'YData',xr2(:,2),...
        'ZData',xr2(:,3),...
        'LineWidth',2,'Color',[1 0 0]);
    
        circle3handle = line('XData',xr3(:,1),...
        'YData',xr3(:,2),...
        'ZData',xr3(:,3),...
        'LineWidth',2,'Color',[1 0 0]);
    
        circle4handle = line('XData',xr4(:,1),...
        'YData',xr4(:,2),...
        'ZData',xr4(:,3),...
        'LineWidth',2,'Color',[1 0 0]);
    [sx,sy,sz] = sphere;
    SPx = 0.02*sx + oxP;
    SPy = 0.02*sy + oyP;
    SPz = 0.02*sz + ozP;
    payloadhandle = surf(SPx,SPy,SPz,'FaceColor',[1 0 0],'EdgeColor','none');
    cablehandle = line('Xdata',[oxQ,oxP],...
        'Ydata',[oyQ,oyP],...
        'Zdata',[ozQ,ozP],...
        'LineWidth',2,'Color','k');
    
    % Set initial limits
    Xlim = [xQ(1)-AxisLength xQ(1)+AxisLength] * LimitRatio;
    Ylim = [yQ(1)-AxisLength yQ(1)+AxisLength] * LimitRatio;
    Zlim = [zQ(1)-AxisLength zQ(1)+AxisLength] * LimitRatio;
    set(gca, 'Xlim', Xlim, 'Ylim', Ylim, 'Zlim', Zlim);
    
        % Set initial view
    view(View(1, :));
%%
    for i = 1:numPlotSamples
        % Update graph title
        if(strcmp(Title, ''))
            titleText = sprintf('Sample %i of %i', 1+((i-1)*SamplePlotFreq), numSamples);
        else
            titleText = strcat(Title, ' (', sprintf('Sample %i of %i', 1+((i-1)*SamplePlotFreq), numSamples), ')');
        end
        title(titleText);
        
        % Plot body x y z axes
        if(strcmp(Trail, 'DotsOnly') || strcmp(Trail, 'All'))
            xQ(1:i) = pQ(1:i,1);
            yQ(1:i) = pQ(1:i,2);
            zQ(1:i) = pQ(1:i,3);
            oxQ = pQ(i,1);
            oyQ = pQ(i,2);
            ozQ = pQ(i,3);
            
            xP(1:i) = pP(1:i,1);
            yP(1:i) = pP(1:i,2);
            zP(1:i) = pP(1:i,3);
            oxP = pP(i,1);
            oyP = pP(i,2);
            ozP = pP(i,3);
            
            SPx = 0.02*sx + oxP;
            SPy = 0.02*sy + oyP;
            SPz = 0.02*sz + ozP;
        else
            xQ = pQ(i,1);
            yQ = pQ(i,2);
            zQ = pQ(i,3);
            oxQ = pQ(i,1);
            oyQ = pQ(i,2);
            ozQ = pQ(i,3);
            
            xP = pP(i,1);
            yP = pP(i,2);
            zP = pP(i,3);
            oxP = pP(i,1);
            oyP = pP(i,2);
            ozP = pP(i,3);
            
            SPx = 0.06*sx + oxP;
            SPy = 0.06*sy + oyP;
            SPz = 0.06*sz + ozP;
            
        end
        if(strcmp(Trail, 'All'))
            ox(1:i) = pQ(1:i,1);
            oy(1:i) = pQ(1:i,2);
            oz(1:i) = pQ(1:i,3);
            ux(1:i) = R(1,1,1:i);
            vx(1:i) = R(2,1,1:i);
            wx(1:i) = R(3,1,1:i);
            uy(1:i) = R(1,2,1:i);
            vy(1:i) = R(2,2,1:i);
            wy(1:i) = R(3,2,1:i);
            uz(1:i) = R(1,3,1:i);
            vz(1:i) = R(2,3,1:i);
            wz(1:i) = R(3,3,1:i);
        else
            ox = pQ(i,1);
            oy = pQ(i,2);
            oz = pQ(i,3);
            ux = R(1,1,i);
            vx = R(2,1,i);
            wx = R(3,1,i);
            uy = R(1,2,i);
            vy = R(2,2,i);
            wy = R(3,2,i);
            uz = R(1,3,i);
            vz = R(2,3,i);
            wz = R(3,3,i);
        end
        x1 = [oxQ,oyQ,ozQ] + 0.5*len*[ux,vx,wx];
        x3 = [oxQ,oyQ,ozQ] - 0.5*len*[ux,vx,wx];
        x2 = [oxQ,oyQ,ozQ] + 0.5*len*[uy,vy,wy];
        x4 = [oxQ,oyQ,ozQ] - 0.5*len*[uy,vy,wy];

        QX1 = [x1(1),x3(1)];
        QX2 = [x1(2),x3(2)];
        QX3 = [x1(3),x3(3)];
        QY1 = [x2(1),x4(1)];
        QY2 = [x2(2),x4(2)];
        QY3 = [x2(3),x4(3)];
        
        
        eb1 = [ux,vx,wx];
        eb2 = [uy,vy,wy];
        eb3 = [uz,vz,wz];
        temp = 0:(2*pi)/10 :2*pi;
        r0 = 0.04;
        xr1 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x1;
        xr2 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x2;
        xr3 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x3;
        xr4 = r0*cos(temp')*eb1 + r0*sin(temp')*eb2 + ones(size(temp'))*x4;
        ratio = 0.4;
        set(orgHandle, 'xdata', xQ, 'ydata', yQ, 'zdata', zQ,'Color','b','LineWidth',0.1);% 原点的句柄设置
        set(quivXhandle, 'xdata', oxQ, 'ydata', oyQ, 'zdata', ozQ,'udata', ratio*ux, 'vdata', ratio*vx, 'wdata', ratio*wx);
        set(quivYhandle, 'xdata', oxQ, 'ydata', oyQ, 'zdata', ozQ,'udata', ratio*uy, 'vdata', ratio*vy, 'wdata', ratio*wy);
        set(quivZhandle, 'xdata', oxQ, 'ydata', oyQ, 'zdata', ozQ,'udata', ratio*uz, 'vdata', ratio*vz, 'wdata', ratio*wz);
        set(lineQXhandle,'XData',QX1,'YData',QX2,'ZData',QX3);
        set(lineQYhandle,'XData',QY1,'YData',QY2,'ZData',QY3);
        set(circle1handle,'XData',xr1(:,1),'YData',xr1(:,2),'ZData',xr1(:,3));
        set(circle2handle,'XData',xr2(:,1),'YData',xr2(:,2),'ZData',xr2(:,3));
        set(circle3handle,'XData',xr3(:,1),'YData',xr3(:,2),'ZData',xr3(:,3));
        set(circle4handle,'XData',xr4(:,1),'YData',xr4(:,2),'ZData',xr4(:,3));
        %set(PorgHandle,'xdata', xP, 'ydata', yP, 'zdata', zP,'Color','b','LineWidth',0.1);
        set(payloadhandle,'XData',SPx,'YData',SPy,'ZData',SPz,'FaceColor',[1 0 0],'EdgeColor','none');
%     cablehandle = line('Xdata',[oxQ,oxP],...
%         'Ydata',[oyQ,oyP],...
%         'Zdata',[ozQ,ozP],...
%         'LineWidth',2,'Color','k');
%     if(isTension(i)==0)
%         set(cablehandle,'Color','g');
%     else
%         set(cablehandle,'Color','k');
%     end
    set(cablehandle, 'XData',[oxQ,oxP],'YData',[oyQ,oyP],'ZData',[ozQ,ozP]);
        if(isFixView)
             Xlim(1) = pQ(i,1) - LimitRatio*AxisLength;
             Ylim(1) = pQ(i,2) - LimitRatio*AxisLength;
             Zlim(1) = pQ(i,3) - LimitRatio*AxisLength;
             Xlim(2) = pQ(i,1) + LimitRatio*AxisLength;
             Ylim(2) = pQ(i,2) + LimitRatio*AxisLength;
             Zlim(2) = pQ(i,3) + LimitRatio*AxisLength;
             set(gca, 'Xlim', Xlim, 'Ylim', Ylim, 'Zlim', Zlim);
        else
    %       % Adjust axes for snug fit and draw
            axisLimChanged = false;
            if((pQ(i,1) - AxisLength) < Xlim(1)), Xlim(1) = pQ(i,1) - LimitRatio*AxisLength; axisLimChanged = true; end
            if((pQ(i,2) - AxisLength) < Ylim(1)), Ylim(1) = pQ(i,2) - LimitRatio*AxisLength; axisLimChanged = true; end
            if((pQ(i,3) - AxisLength) < Zlim(1)), Zlim(1) = pQ(i,3) - LimitRatio*AxisLength; axisLimChanged = true; end
            if((pQ(i,1) + AxisLength) > Xlim(2)), Xlim(2) = pQ(i,1) + LimitRatio*AxisLength; axisLimChanged = true; end
            if((pQ(i,2) + AxisLength) > Ylim(2)), Ylim(2) = pQ(i,2) + LimitRatio*AxisLength; axisLimChanged = true; end
            if((pQ(i,3) + AxisLength) > Zlim(2)), Zlim(2) = pQ(i,3) + LimitRatio*AxisLength; axisLimChanged = true; end
            if(axisLimChanged), set(gca, 'Xlim', Xlim, 'Ylim', Ylim, 'Zlim', Zlim); end
        end
        drawnow;
% 
        %%% Adjust view
        if(CreateAVI)
%             if(numel(View) > 2)
%                 view(View(i, :));
%             end
        end
% 
        % Add frame to AVI object
        if(~isempty(aviobj))
            frame = getframe(fig);
            writeVideo(aviobj,frame);
%             aviobj = addframe(aviobj, frame);
        end

    end

end

