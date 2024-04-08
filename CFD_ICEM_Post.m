function CFD_ICEM_Post(Dir, FileName_Catia, Act_Save)


Dir_ICEM = Dir.ICEM;
Dir_Sub_ICEM = [Dir_ICEM,'\',FileName_Catia];

FontSize = 10;
LineWidth_Grid = 0.2;


%% 读取数据
File_Msh = fopen([Dir_Sub_ICEM,'\',FileName_Catia,'.msh'],'r');

Data_Msh = textscan(File_Msh, '%s', 'Delimiter','\n');
Data_Msh = Data_Msh{1};

fclose(File_Msh);


%% 提取数据 
%NodeCoord
Index = strcmp(Data_Msh, '(');
nRow_Init = find(Index == 1, 1) + 1;

Index = strcmp(Data_Msh, '))');
nRow_End = find(Index == 1, 1) - 1;
    
Data_NodeCoord = Data_Msh(nRow_Init:nRow_End);

Data_Cache = textscan(sprintf('%s\n', Data_NodeCoord{:}), '%f %f %f');
Data_Cache = [Data_Cache{1}, Data_Cache{2}, Data_Cache{3}];
Data_Mesh.NodeCoord = Data_Cache;

%ElementConnect_Fluid
ZoneName_Total = {'Plane','Symm'};

for nZone = 1:length(ZoneName_Total)
    
    Index = contains(Data_Msh, ['zone ',ZoneName_Total{nZone}], 'IgnoreCase', true);
    nRow_Init = find(Index == 1, 1) + 2;
    
    Index = contains(Data_Msh(nRow_Init:end), ')', 'IgnoreCase', true);
    nRow_End = find(Index == 1, 1) - 1 + (nRow_Init - 1);

    Data_EC_Fluid = Data_Msh(nRow_Init:nRow_End);

    Data_Cache = textscan(sprintf('%s\n', Data_EC_Fluid{:}), '%s %s %s %s %s %s');
    Data_Cache = [hex2dec(Data_Cache{1}), hex2dec(Data_Cache{2}), hex2dec(Data_Cache{3}),...
                  hex2dec(Data_Cache{4}), hex2dec(Data_Cache{5}), hex2dec(Data_Cache{6})];
              
    eval(['Data_Mesh.ElementConnect.',ZoneName_Total{nZone},'=Data_Cache;']);
end

%粗估机翼几何数据
[m,n] = size(Data_Mesh.ElementConnect.Plane(:,1:4));
Index = unique(reshape(Data_Mesh.ElementConnect.Plane(:,1:4), m*n, 1));
Data_NodeCoord_Plane = Data_Mesh.NodeCoord(Index,:);

Span_Half = max(abs(Data_NodeCoord_Plane(:,2)));

Index = abs(Data_NodeCoord_Plane(:,2)) < 1e-2;
Data_NodeCoord_Plane_Root = Data_NodeCoord_Plane(Index,:);

Length_Root = max(Data_NodeCoord_Plane_Root(:,1)) - min(Data_NodeCoord_Plane_Root(:,1));

%清除变量
clear Data_Msh;


%% 作图准备
hFigure = figure('NumberTitle','off', 'Name',FileName_Catia);
FigureSize = [1400,940];
set(hFigure, 'Position',[(1920-FigureSize(1))/2,(1200-FigureSize(2))/2,FigureSize]);
hAxis = tight_subplot(1, 1, 0, [-1.5,-1.4], [-0.8,-0.8]);
set(hAxis, 'FontSize',FontSize);
axis equal;
hold on;
grid on;
view(3);
view(-20,40);

hTextBox = annotation('textbox',[0.01,0.888,0.1,0.1],'String',...
                      replace(FileName_Catia,'_','\_'),...
                      'EdgeColor','none','FontSize',FontSize+1);
hTextBox.HorizontalAlignment = 'Left'; 
hTextBox.VerticalAlignment = 'Top';


%% 作图
patch('Vertices', Data_Mesh.NodeCoord, 'Faces', Data_Mesh.ElementConnect.Plane(:,1:4),...
      'FaceColor','w', 'EdgeColor','k', 'LineWidth',LineWidth_Grid);
patch('Vertices', Data_Mesh.NodeCoord, 'Faces', Data_Mesh.ElementConnect.Symm(:,1:4),...
      'FaceColor','w', 'EdgeColor','k', 'LineWidth',LineWidth_Grid);

%设置Axis
RangeLimit_x_Total(1) = min(Length_Root/2 - 3 * Length_Root, min(Data_NodeCoord_Plane(:,1))) - 0.4 * Span_Half;
RangeLimit_x_Total(2) = max(Length_Root/2 + 3 * Length_Root, max(Data_NodeCoord_Plane(:,1))) + 0.4 * Span_Half;

RangeLimit_y_Total = [-1.2, 0] * Span_Half;

RangeLimit_z_Total = [-1.7, 1.9] * Span_Half;

xlim(RangeLimit_x_Total);
ylim(RangeLimit_y_Total);
zlim(RangeLimit_z_Total);


%% 保存
if Act_Save == 1
    %保存.fig
    % saveas(hFigure, [Dir_Sub_ICEM,'\',FileName_Catia,'.fig']);

    %保存.png
    print(hFigure, [Dir_Sub_ICEM,'\',FileName_Catia,'.png'],'-dpng','-r150');

    %删除.fig
    % Dir_Fig = [Dir_Sub_ICEM,'\',FileName_Catia,'.fig'];
    % if exist(Dir_Fig,'file')
    %     dos(['del ',Dir_Fig]);
    % end
    
    %关闭figure
    close(hFigure);
end


end























