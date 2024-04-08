function CFD_ICEM_Gen(Dir, Para_CFD, FlightStateRef, Data_Geometry, FileName_Catia)


Dir_Catia = Dir.Catia;
Dir_ICEM = Dir.ICEM;

MeshNum_x = Para_CFD.ICEM.MeshNum_x;
MeshNum_y = Para_CFD.ICEM.MeshNum_y;
MeshNum_z = Para_CFD.ICEM.MeshNum_z;
Layer_FirstYPlus = Para_CFD.ICEM.Layer_FirstYPlus;

Rho =  FlightStateRef.Rho;
DynamicVisc = FlightStateRef.DynamicVisc;
Velocity = FlightStateRef.Velocity;

MAC_Length = Data_Geometry.MAC_Length;

ScaleRate = Data_Geometry.ScaleRate;

Dir_Sub_ICEM = [Dir_ICEM,'\',FileName_Catia];
if exist(Dir_Sub_ICEM,'dir')
    for i = 1:3
        dos(['Recycle.exe ',Dir_Sub_ICEM]);
    end
    mkdir(Dir_Sub_ICEM);
else
    mkdir(Dir_Sub_ICEM);
end


%% 求解_附面层
Re = Rho * Velocity * MAC_Length / DynamicVisc;

cf = (2 * log10(Re) - 0.65) ^ -2.3;
uf = sqrt(cf * 0.5 * Rho * Velocity^2 / Rho);

yp = Layer_FirstYPlus * DynamicVisc / (Rho * uf);
yH = 2 * yp;
% Layer_FirstHeight = yH * 1000;
Layer_FirstHeight = yH/2 * 1000;

Act_LayerAdjust = 0;
if Data_Geometry.WingPlan.Span(end) > 4000 && Layer_FirstYPlus == 1
    Act_LayerAdjust = 1;
end


%% 求解MeshSize
Data_Num.Span = [13, 6, 6];
Data_Num.Divide = [11, 11, 11];
Data_Num.Layer = [2, 4, 2];

Data_Num.Tip_x = [9, 9, 9];
Data_Num.Tip_y = [3, 3, 3];
Data_Num.Tip_z = Data_Num.Layer;

SublineName_Total = {'Subline_01','Subline_12','Subline_23'};

[Data_SN, Data_SectionPoint, Data_BlockNode, Data_EL] = CFD_ICEM_Gen_GetData(Dir_Catia, Dir_ICEM,...
                                                                             Data_Num, FileName_Catia, SublineName_Total);

Data_MeshSize = CFD_ICEM_Gen_SolveMesh(Data_SectionPoint, Data_BlockNode, Data_EL, Data_Num,...
                                       MeshNum_x, MeshNum_y, MeshNum_z, Layer_FirstHeight, Act_LayerAdjust);
                               

%% 定义Script_加载
%抬头
Str_Rpl = {'ic_chdir Input_Dir_Sub_ICEM'};

%加载.model
GeometryName_Total = {'PLANE', 'BOI_FLUID', 'IN', 'SYMM', 'SUBLINE_01', 'SUBLINE_12', 'SUBLINE_23'};

for nGeometry = 1:length(GeometryName_Total)
    
    Str_Rpl = [Str_Rpl; 'ic_set_meshing_params global 0 gttol 1e-7'];
    
    Str_Rpl = [Str_Rpl; ['ic_run_application_exec . bin c4tetin ',...
                         '{-tcl --overwrite --triangulation-tolerance=1e-4 ',...'
                         '--one-tetin-file="{Input_Dir_Sub_ICEM/Input_FileName_Catia_',GeometryName_Total{nGeometry},'.tin}" ',...
                         '"Input_Dir_Catia/Input_FileName_Catia_',GeometryName_Total{nGeometry},'.model"}']];
                      
    Str_Rpl = [Str_Rpl; ['ic_load_tetin Input_Dir_Sub_ICEM/Input_FileName_Catia_',GeometryName_Total{nGeometry},'.tin 0 1']];
    
    if contains(GeometryName_Total{nGeometry}, 'SUBLINE')
        Str_Rpl = [Str_Rpl; ['ic_geo_rename_family W3L0C0 ',GeometryName_Total{nGeometry},' 1']];
    else
        Str_Rpl = [Str_Rpl; ['ic_geo_rename_family W3L0C4 ',GeometryName_Total{nGeometry},' 1']];
    end

end

%对称几何
Str_Rpl = [Str_Rpl; 'ic_move_geometry point names [ic_geo_get_objects point] mirror {0 1 0} cent {0 0 0}'];
Str_Rpl = [Str_Rpl; 'ic_move_geometry curve names [ic_geo_get_objects curve] mirror {0 1 0} cent {0 0 0}'];
Str_Rpl = [Str_Rpl; 'ic_move_geometry surface names [ic_geo_get_objects surface] mirror {0 1 0} cent {0 0 0}'];

%缩放几何
for nGeometry = 1:length(GeometryName_Total)
    
    if strcmp(GeometryName_Total{nGeometry},'BOI_FLUID') ||...
       strcmp(GeometryName_Total{nGeometry},'IN') ||...
       strcmp(GeometryName_Total{nGeometry},'SYMM')
        
        Str_Rpl = [Str_Rpl; 'ic_move_geometry surface names [ic_geo_objects_in_family surface ',GeometryName_Total{nGeometry},']',...
                            ' scale {',num2str(1/ScaleRate),' ',num2str(1/ScaleRate),' ',num2str(1/ScaleRate),'} cent {0 0 0}'];
    end
    
    if strcmp(GeometryName_Total{nGeometry},'SUBLINE_12') ||...
       strcmp(GeometryName_Total{nGeometry},'SUBLINE_23')
        
        Str_Rpl = [Str_Rpl; 'ic_move_geometry point names [ic_geo_objects_in_family point ',GeometryName_Total{nGeometry},']',...
                            ' scale {',num2str(1/ScaleRate),' ',num2str(1/ScaleRate),' ',num2str(1/ScaleRate),'} cent {0 0 0}'];
        
        Str_Rpl = [Str_Rpl; 'ic_move_geometry curve names [ic_geo_objects_in_family curve ',GeometryName_Total{nGeometry},']',...
                            ' scale {',num2str(1/ScaleRate),' ',num2str(1/ScaleRate),' ',num2str(1/ScaleRate),'} cent {0 0 0}'];
   
    end
    
end

%加载.blk
Str_Rpl = [Str_Rpl; 'ic_hex_restore_blocking Input_Dir_ICEM/Template/ICEM_Block_Template.blk'];


%% 定义Script_映射
%点
for nSubline = 1:length(SublineName_Total)
    
    %Point_Main (nSpan/nDivide/nDU/nLayer)
    %Node_Main  (nSpan/nDivide/nDU/nLayer)
    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline)
                    
                    Str_Rpl = [Str_Rpl; 'ic_hex_move_node ',...
                                        Data_SN(nSubline).Main.Node{nSpan, nDU, nDivide, nLayer},' ',...
                                        Data_SN(nSubline).Main.Point{nSpan, nDU, nDivide, nLayer}];
                                    
                end
            end
        end
    end
    
    %Point_Tip (nTip_x/nTip_y/nTip_z)
    %Node_Tip  (nTip_x/nTip_y/nTip_z)
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline)
            
                Str_Rpl = [Str_Rpl; 'ic_hex_move_node ',...
                                    Data_SN(nSubline).Tip.Node{nTip_x, nTip_y, nTip_z},' ',...
                                    Data_SN(nSubline).Tip.Point{nTip_x, nTip_y, nTip_z}];
                                
            end                 
        end
    end
    
end

%线
for nSubline = 1:length(SublineName_Total)
    
    %Curve_x_Main (nSpan/nDU/           /nLayer)
    %Edge_x_Main  (nSpan/nDU/(nDivide-1)/nLayer)
    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline) - 1
                for nLayer = 1:Data_Num.Layer(nSubline)
                    
                    if ~strcmp(Data_SN(nSubline).Main.Edge_x{nSpan, nDU, nDivide, nLayer}, 'NaN') &&...
                       ~strcmp(Data_SN(nSubline).Main.Curve_x{nSpan, nDU, nLayer}, 'NaN')
                   
                        Str_Rpl = [Str_Rpl; 'ic_hex_set_edge_projection ',...
                                            Data_SN(nSubline).Main.Edge_x{nSpan, nDU, nDivide, nLayer},' 1 ',...
                                            Data_SN(nSubline).Main.Curve_x{nSpan, nDU, nLayer}];
                                        
                    end

                end
            end
        end
    end
    
    %Curve_y_Main (         /nDU/nDivide/nLayer)
    %Edge_y_Main  ((nSpan-1)/nDU/nDivide/nLayer)
    for nSpan = 1:Data_Num.Span(nSubline) - 1
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline)
                    
                    if ~strcmp(Data_SN(nSubline).Main.Edge_y{nSpan, nDU, nDivide, nLayer}, 'NaN') &&...
                       ~strcmp(Data_SN(nSubline).Main.Curve_y{nDU, nDivide, nLayer}, 'NaN')
                        Str_Rpl = [Str_Rpl; 'ic_hex_set_edge_projection ',...
                                            Data_SN(nSubline).Main.Edge_y{nSpan, nDU, nDivide, nLayer},' 1 ',...
                                            Data_SN(nSubline).Main.Curve_y{nDU, nDivide, nLayer}];
                    end

                end
            end
        end
    end
    
    %Curve_z_Main (nSpan/nDU/nDivide/          )
    %Edge_z_Main  (nSpan/nDU/nDivide/(nLayer-1))
    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline) - 1
                    
                    if ~strcmp(Data_SN(nSubline).Main.Edge_z{nSpan, nDU, nDivide, nLayer}, 'NaN') &&...
                       ~strcmp(Data_SN(nSubline).Main.Curve_z{nSpan, nDU, nDivide}, 'NaN')
                   
                        Str_Rpl = [Str_Rpl; 'ic_hex_set_edge_projection ',...
                                            Data_SN(nSubline).Main.Edge_z{nSpan, nDU, nDivide, nLayer},' 1 ',...
                                            Data_SN(nSubline).Main.Curve_z{nSpan, nDU, nDivide}];
                                        
                    end

                end
            end
        end
    end
    
    %Curve_x_Tip (          /nTip_y/nTip_z)
    %Edge_x_Tip  ((nTip_x-1)/nTip_y/nTip_z）
    for nTip_x = 1:Data_Num.Tip_x(nSubline) - 1
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline)
                
                if ~strcmp(Data_SN(nSubline).Tip.Edge_x{nTip_x, nTip_y, nTip_z}, 'NaN') &&...
                   ~strcmp(Data_SN(nSubline).Tip.Curve_x{nTip_y, nTip_z}, 'NaN')     
               
                    Str_Rpl = [Str_Rpl; 'ic_hex_set_edge_projection ',...
                                        Data_SN(nSubline).Tip.Edge_x{nTip_x, nTip_y, nTip_z},' 1 ',...
                                        Data_SN(nSubline).Tip.Curve_x{nTip_y, nTip_z}];
                                    
                end

            end
        end
    end
    
    %Curve_y_Tip (nTip_x/      /nTip_z)
    %Edge_y_Tip  (nTip_x/(nTip_y-1)/nTip_z）
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline) - 1
            for nTip_z = 1:Data_Num.Tip_z(nSubline)
                
                if ~strcmp(Data_SN(nSubline).Tip.Edge_y{nTip_x, nTip_y, nTip_z}, 'NaN') &&...
                   ~strcmp(Data_SN(nSubline).Tip.Curve_y{nTip_x, nTip_z}, 'NaN')   
               
                    Str_Rpl = [Str_Rpl; 'ic_hex_set_edge_projection ',...
                                        Data_SN(nSubline).Tip.Edge_y{nTip_x, nTip_y, nTip_z},' 1 ',...
                                        Data_SN(nSubline).Tip.Curve_y{nTip_x, nTip_z}];
                                    
                end
              
            end
        end
    end
    
    %Curve_z_Tip (nTip_x/nTip_y/          )
    %Edge_z_Tip  (nTip_x/nTip_y/(nTip_z-1)）
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline) - 1

                if ~strcmp(Data_SN(nSubline).Tip.Edge_z{nTip_x, nTip_y, nTip_z}, 'NaN') &&...
                   ~strcmp(Data_SN(nSubline).Tip.Curve_z{nTip_x, nTip_y}, 'NaN') 
               
                    Str_Rpl = [Str_Rpl; 'ic_hex_set_edge_projection ',...
                                        Data_SN(nSubline).Tip.Edge_z{nTip_x, nTip_y, nTip_z},' 1 ',...
                                        Data_SN(nSubline).Tip.Curve_z{nTip_x, nTip_y}];
                                    
                end
                                    
            end
        end
    end
    
end

%面
for nSubline = 1:length(SublineName_Total)
    
    if strcmp(SublineName_Total{nSubline}, 'Subline_01')
        
        Str_Cache = '';
        for i = 1:length(Data_SN(nSubline).Surface.PLANE)
            Str_Cache = [Str_Cache, ' { ',Data_SN(1).Surface.PLANE{i},' } '];
        end
        
        Str_Rpl = [Str_Rpl; 'ic_hex_project_face node_numbers', Str_Cache, 'PLANE'];
        
    end
    
    if strcmp(SublineName_Total{nSubline}, 'Subline_12')
        
        Str_Cache = '';
        for i = 1:length(Data_SN(nSubline).Surface.BOI_FLUID)
            Str_Cache = [Str_Cache, ' { ',Data_SN(2).Surface.BOI_FLUID{i},' } '];
        end
        
        Str_Rpl = [Str_Rpl; 'ic_hex_project_face node_numbers', Str_Cache, 'BOI_FLUID'];
        
    end
    
    if strcmp(SublineName_Total{nSubline}, 'Subline_23')
        
        Str_Cache = '';
        for i = 1:length(Data_SN(nSubline).Surface.IN)
            Str_Cache = [Str_Cache, ' { ',Data_SN(3).Surface.IN{i},' } '];
        end
        
        Str_Rpl = [Str_Rpl; 'ic_hex_project_face node_numbers', Str_Cache, 'IN'];
        
    end

    Str_Cache = '';
    for i = 1:length(Data_SN(nSubline).Surface.SYMM)
        Str_Cache = [Str_Cache, ' { ',Data_SN(nSubline).Surface.SYMM{i},' } '];
    end
    
    Str_Rpl = [Str_Rpl; 'ic_hex_project_face node_numbers', Str_Cache, 'SYMM'];
    
end


%% 定义Script_节点规律
for nSubline = 1:length(SublineName_Total)
    %Main_x
    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline) - 1
                for nLayer = 1:Data_Num.Layer(nSubline)
                    
                    Law = Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Law;
                    
                    if ~isnan(Law)
                        Nodes = Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Nodes;
                        h1 = Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Size(1);
                        h2 = Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Size(2);
                        r1 = Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Rate(1);
                        r2 = Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Rate(2);

                        Str_Rpl = [Str_Rpl; 'ic_hex_set_mesh ',...
                                            Data_SN(nSubline).Main.Edge_x{nSpan, nDU, nDivide, nLayer}(1:end-2),...
                                            ' n ', num2str(Nodes),...
                                            ' h1 ', num2str(h1), ' h2 ', num2str(h2),...
                                            ' r1 ', num2str(r1), ' r2 ', num2str(r2),...
                                            ' lmax 0 ', Law, ' unlocked'];
                    end

                end
            end
        end
    end
    
    %Main_y
    for nSpan = 1:Data_Num.Span(nSubline) - 1
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline)
                    
                    Law = Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Law;
                    
                    if ~isnan(Law)
                        Nodes = Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Nodes;
                        if nSubline == 1
                            h1 = Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(1);
                            h2 = Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(2);
                        else
                            h1 = Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(2);
                            h2 = Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(1);
                        end
                        r1 = Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Rate(1);
                        r2 = Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Rate(2);

                        Str_Rpl = [Str_Rpl; 'ic_hex_set_mesh ',...
                                            Data_SN(nSubline).Main.Edge_y{nSpan, nDU, nDivide, nLayer}(1:end-2),...
                                            ' n ', num2str(Nodes),...
                                            ' h1 ', num2str(h1), ' h2 ', num2str(h2),...
                                            ' r1 ', num2str(r1), ' r2 ', num2str(r2),...
                                            ' lmax 0 ', Law, ' unlocked'];
                    end

                end
            end
        end
    end
    
    %Main_z
    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline) - 1
                    
                    Law = Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Law;
                    
                    if ~isnan(Law)
                        Nodes = Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Nodes;
                        h1 = Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Size(1);
                        h2 = Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Size(2);
                        r1 = Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Rate(1);
                        r2 = Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Rate(2);

                        Str_Rpl = [Str_Rpl; 'ic_hex_set_mesh ',...
                                            Data_SN(nSubline).Main.Edge_z{nSpan, nDU, nDivide, nLayer}(1:end-2),...
                                            ' n ', num2str(Nodes),...
                                            ' h1 ', num2str(h1), ' h2 ', num2str(h2),...
                                            ' r1 ', num2str(r1), ' r2 ', num2str(r2),...
                                            ' lmax 0 ', Law, ' unlocked'];
                    end

                end
            end
        end
    end
    
    %Tip_x
    for nTip_x = 1:Data_Num.Tip_x(nSubline) - 1
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline)

                Law = Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Law;
                
                if ~isnan(Law)
                    Nodes = Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Nodes;
                    h1 = Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Size(1);
                    h2 = Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Size(2);
                    r1 = Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Rate(1);
                    r2 = Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Rate(2);


                    Str_Rpl = [Str_Rpl; 'ic_hex_set_mesh ',...
                                        Data_SN(nSubline).Tip.Edge_x{nTip_x, nTip_y, nTip_z}(1:end-2),...
                                        ' n ', num2str(Nodes),...
                                        ' h1 ', num2str(h1), ' h2 ', num2str(h2),...
                                        ' r1 ', num2str(r1), ' r2 ', num2str(r2),...
                                        ' lmax 0 ', Law, ' unlocked'];
                end

            end
        end
    end
    
    %Tip_y
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline) - 1
            for nTip_z = 1:Data_Num.Tip_z(nSubline)
                
                Law = Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Law;
                
                if ~isnan(Law)
                    Nodes = Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Nodes;
                    h1 = Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Size(1);
                    h2 = Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Size(2);
                    r1 = Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Rate(1);
                    r2 = Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Rate(2);

                    Str_Rpl = [Str_Rpl; 'ic_hex_set_mesh ',...
                                        Data_SN(nSubline).Tip.Edge_y{nTip_x, nTip_y, nTip_z}(1:end-2),...
                                        ' n ', num2str(Nodes),...
                                        ' h1 ', num2str(h1), ' h2 ', num2str(h2),...
                                        ' r1 ', num2str(r1), ' r2 ', num2str(r2),...
                                        ' lmax 0 ', Law, ' unlocked'];
                end

            end
        end
    end
    
    %Tip_z
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline) - 1

                Law = Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Law;
                
                if ~isnan(Law)
                    Nodes = Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Nodes;
                    h1 = Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Size(1);
                    h2 = Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Size(2);
                    r1 = Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Rate(1);
                    r2 = Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Rate(2);

                    Str_Rpl = [Str_Rpl; 'ic_hex_set_mesh ',...
                                        Data_SN(nSubline).Tip.Edge_z{nTip_x, nTip_y, nTip_z}(1:end-2),...
                                        ' n ', num2str(Nodes),...
                                        ' h1 ', num2str(h1), ' h2 ', num2str(h2),...
                                        ' r1 ', num2str(r1), ' r2 ', num2str(r2),...
                                        ' lmax 0 ', Law, ' unlocked'];
                end

            end
        end
    end
  
end
    
                
%% 定义Script_保存
Str_Rpl = [Str_Rpl; 'ic_save_tetin Input_Dir_Sub_ICEM/Input_FileName_Catia.tin'];

Str_Rpl = [Str_Rpl; 'ic_hex_save_blocking Input_Dir_Sub_ICEM/Input_FileName_Catia.blk'];

Str_Rpl = [Str_Rpl; 'ic_save_project_file Input_Dir_Sub_ICEM/Input_FileName_Catia.prj ',...
                    '{array\ set\ file_name\ \{ {tetin Input_Dir_Sub_ICEM/Input_FileName_Catia.tin} ',...
                    '{settings Input_Dir_Sub_ICEM/Input_FileName_Catia.prj} ',...
                    '{blocking Input_Dir_Sub_ICEM/Input_FileName_Catia.blk} \} }'];


%% 定义Script_重加载
Str_Rpl = [Str_Rpl; 'ic_hex_unload_blocking ';... 
                    'ic_unload_tetin ';... 
                    'ic_empty_tetin '];

Str_Rpl = [Str_Rpl; 'ic_load_tetin Input_Dir_Sub_ICEM/Input_FileName_Catia.tin ';...
                    'ic_hex_restore_blocking Input_Dir_Sub_ICEM/Input_FileName_Catia.blk '];
                
                
%% 定义Script_输出.uns
Str_Rpl = [Str_Rpl; 'ic_hex_create_mesh ',...
                    'PLANE BOI_FLUID IN SYMM SUBLINE_01 SUBLINE_12 SUBLINE_23 FLUID ',...
                    'proj 2 dim_to_mesh 3'];

Str_Rpl = [Str_Rpl; 'ic_hex_write_file Input_Dir_Sub_ICEM/Input_FileName_Catia.uns ',...
                    'PLANE BOI_FLUID IN SYMM SUBLINE_01 SUBLINE_12 SUBLINE_23 FLUID ',...
                    'proj 2 dim_to_mesh 3 no_boco'];
 
                
%% 定义Script_输出网格质量
Str_Rpl = [Str_Rpl; 'ic_uns_load Input_Dir_Sub_ICEM/Input_FileName_Catia.uns 3 0 {} 1'];
    
Str_Rpl = [Str_Rpl; 'ic_uns_subset_create smooth_do_map 0';...
                    'ic_uns_update_family_type smooth_do_map {ORFN PLANE BOI_FLUID IN SYMM SUBLINE_01 SUBLINE_12 SUBLINE_23 FLUID} {HEXA_8} update 1';...
                    'ic_uns_metric smooth_do_map Determinant eval_at_node_method 0'];


%% 定义Script_输出.msh
Str_Rpl = [Str_Rpl; 'ic_boco_solver {Ansys Fluent} ';...
                    'ic_solver_mesh_info {Ansys Fluent} ';...
                    'ic_solution_set_solver {Ansys Fluent} 1'];

Str_Rpl = [Str_Rpl; 'ic_boco_save Input_Dir_Sub_ICEM/Input_FileName_Catia.fbc ';...
                    'ic_boco_save_atr Input_Dir_Sub_ICEM/Input_FileName_Catia.atr '];

Str_Rpl = [Str_Rpl; 'ic_exec {F:/ANSYS Inc/v212/icemcfd/win64_amd/icemcfd/output-interfaces/fluent6} ',...
                    '-dom Input_Dir_Sub_ICEM/Input_FileName_Catia.uns -b Input_Dir_Sub_ICEM/Input_FileName_Catia.fbc ',...
                    '-dim2d Input_Dir_Sub_ICEM/Input_FileName_Catia '];

                
%% 修改Script
for nSubline = 1:size(Str_Rpl,1)
    Str_Rpl{nSubline} = replace(Str_Rpl{nSubline}, 'Input_Dir_Catia', replace(Dir_Catia,'\','/'));
    Str_Rpl{nSubline} = replace(Str_Rpl{nSubline}, 'Input_FileName_Catia', FileName_Catia);
    Str_Rpl{nSubline} = replace(Str_Rpl{nSubline}, 'Input_Dir_ICEM', replace(Dir_ICEM,'\','/'));
    Str_Rpl{nSubline} = replace(Str_Rpl{nSubline}, 'Input_Dir_Sub_ICEM', replace(Dir_Sub_ICEM,'\','/'));
end


%% 运行ICEM
%生成.rpl
File_Rpl_Generate = fopen([Dir_Sub_ICEM,'\',FileName_Catia,'.rpl'],'w');
for nSubline = 1:size(Str_Rpl,1)
    fprintf(File_Rpl_Generate,'%s\n',Str_Rpl{nSubline,1});
end
fclose(File_Rpl_Generate);

%运行
[~,~] = dos(['cd /D ',Dir_Sub_ICEM,' && ',...
             'call icemcfd -batch', ' -script ',FileName_Catia,'.rpl',' > LogFile.log']); 

%删除多余.tin
for nGeometry = 1:length(GeometryName_Total)
    % delete([Dir_Sub_ICEM,'\',FileName_Catia,'_',GeometryName_Total{nGeometry},'.tin']);
    dos(['del ',Dir_Sub_ICEM,'\',FileName_Catia,'_',GeometryName_Total{nGeometry},'.tin']);
end   


end



       
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         