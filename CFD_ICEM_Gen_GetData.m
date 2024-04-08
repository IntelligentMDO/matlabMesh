function [Data_SN, Data_SectionPoint, Data_BlockNode, Data_EL]...
         = CFD_ICEM_Gen_GetData(Dir_Catia, Dir_ICEM, Data_Num, FileName_Catia, SublineName_Total)


%% 编号
for nSubline = 1:length(SublineName_Total)
    Data_Excel(nSubline).SN = readcell([Dir_ICEM,'\Template\ICEM_SerialNumber_Template.xls'], 'Sheet',SublineName_Total{nSubline});
end

%点
for nSubline = 1:length(SublineName_Total)

    %Point_Main (nSpan/nDivide/nDU/nLayer)
    %Node_Main  (nSpan/nDivide/nDU/nLayer)
    Data_SN(nSubline).Main.Point = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Point_Main', nSubline); 
    Data_SN(nSubline).Main.Node = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Node_Main', nSubline);
    
    
    %Point_Tip (nTip_x/nTip_y/nTip_z)
    %Node_Tip  (nTip_x/nTip_y/nTip_z)
    Data_SN(nSubline).Tip.Point = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Point_Tip', nSubline);
    Data_SN(nSubline).Tip.Node = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Node_Tip', nSubline);
    
end

%线
for nSubline = 1:length(SublineName_Total)
    
    %Curve_x_Main (nSpan/nDU/           /nLayer)
    %Edge_x_Main  (nSpan/nDU/(nDivide-1)/nLayer)
    Data_SN(nSubline).Main.Curve_x = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Curve_x_Main', nSubline);
    
    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline) - 1
                for nLayer = 1:Data_Num.Layer(nSubline)

                    Data_SN(nSubline).Main.Edge_x(nSpan, nDU, nDivide, nLayer)...
                    = {[Data_SN(nSubline).Main.Node{nSpan, nDU, nDivide, nLayer},' ',...
                        Data_SN(nSubline).Main.Node{nSpan, nDU, nDivide + 1, nLayer},' 0']};
                    
                end
            end
        end
    end
    
    %Curve_y_Main (         /nDU/nDivide/nLayer)
    %Edge_y_Main  ((nSpan-1)/nDU/nDivide/nLayer)
    Data_SN(nSubline).Main.Curve_y= GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Curve_y_Main', nSubline);
    
    for nSpan = 1:Data_Num.Span(nSubline) - 1
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline)

                    Data_SN(nSubline).Main.Edge_y(nSpan, nDU, nDivide, nLayer)...
                    = {[Data_SN(nSubline).Main.Node{nSpan, nDU, nDivide, nLayer},' ',...
                        Data_SN(nSubline).Main.Node{nSpan + 1, nDU, nDivide, nLayer},' 0']};

                end
            end
        end
    end
    
    %Curve_z_Main (nSpan/nDU/nDivide/          )
    %Edge_z_Main  (nSpan/nDU/nDivide/(nLayer-1))
    Data_SN(nSubline).Main.Curve_z= GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Curve_z_Main', nSubline);
    
    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline) - 1

                    Data_SN(nSubline).Main.Edge_z(nSpan, nDU, nDivide, nLayer)...
                    = {[Data_SN(nSubline).Main.Node{nSpan, nDU, nDivide, nLayer},' ',...
                        Data_SN(nSubline).Main.Node{nSpan, nDU, nDivide, nLayer + 1},' 0']};
                    
                end
            end
        end
    end

    %Curve_x_Tip (          /nTip_y/nTip_z)
    %Edge_x_Tip  ((nTip_x-1)/nTip_y/nTip_z）
    Data_SN(nSubline).Tip.Curve_x = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Curve_x_Tip', nSubline);
    
    for nTip_x = 1:Data_Num.Tip_x(nSubline) - 1
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline)

                Data_SN(nSubline).Tip.Edge_x(nTip_x, nTip_y, nTip_z)...
                = {[Data_SN(nSubline).Tip.Node{nTip_x, nTip_y, nTip_z},' ',...
                    Data_SN(nSubline).Tip.Node{nTip_x + 1, nTip_y, nTip_z},' 0']};

            end
        end
    end
    
    %Curve_y_Tip (nTip_x/      /nTip_z)
    %Edge_y_Tip  (nTip_x/(nTip_y-1)/nTip_z）
    Data_SN(nSubline).Tip.Curve_y = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Curve_y_Tip', nSubline);
    
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline) - 1
            for nTip_z = 1:Data_Num.Tip_z(nSubline)

                Data_SN(nSubline).Tip.Edge_y(nTip_x, nTip_y, nTip_z)...
                = {[Data_SN(nSubline).Tip.Node{nTip_x, nTip_y, nTip_z},' ',...
                    Data_SN(nSubline).Tip.Node{nTip_x, nTip_y + 1, nTip_z},' 0']};

            end
        end
    end
    
    %Curve_z_Tip (nTip_x/nTip_y/          )
    %Edge_z_Tip  (nTip_x/nTip_y/(nTip_z-1)）
    Data_SN(nSubline).Tip.Curve_z = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Curve_z_Tip', nSubline);
    
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline) - 1

                Data_SN(nSubline).Tip.Edge_z(nTip_x, nTip_y, nTip_z)...
                = {[Data_SN(nSubline).Tip.Node{nTip_x, nTip_y, nTip_z},' ',...
                    Data_SN(nSubline).Tip.Node{nTip_x, nTip_y, nTip_z + 1},' 0']};

            end
        end
    end
    
end

%面
for nSubline = 1:length(SublineName_Total)
    
    if strcmp(SublineName_Total{nSubline}, 'Subline_01')
        Data_SN(nSubline).Surface.PLANE = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Surface_PLANE', nSubline); 
    end
    
    if strcmp(SublineName_Total{nSubline}, 'Subline_12')
        Data_SN(nSubline).Surface.BOI_FLUID = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Surface_BOI_FLUID', nSubline); 
    end
    
    if strcmp(SublineName_Total{nSubline}, 'Subline_23')
        Data_SN(nSubline).Surface.IN = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Surface_IN', nSubline); 
    end

    Data_SN(nSubline).Surface.SYMM = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Surface_SYMM', nSubline); 
    
end

%修改Edge编号
for nSubline = 1:length(SublineName_Total)
      
    %Edge_y_Main
    if strcmp(SublineName_Total{nSubline}, 'Subline_12') || strcmp(SublineName_Total{nSubline}, 'Subline_23')
        
        Data_SN(nSubline).Main.Edge_y = GetSerialNumber(Data_Excel(nSubline).SN, Data_Num, SublineName_Total, 'Edge_y_Main', nSubline);

    end
    
    if strcmp(SublineName_Total{nSubline}, 'Subline_12')
        
        Data_SN(nSubline).Main.Edge_y(:, :, :, 1) = {'NaN'};

    end

end

clear Data_Excel;


%% 翼型坐标
Data_Excel = readcell([Dir_Catia,'\',FileName_Catia, '.xls']);

Num_Section = cell2mat(Data_Excel(3,1));
Num_PointPreSpline = cell2mat(Data_Excel(3,2));

Num_Point = Num_Section * 2 * Num_PointPreSpline;

Row_Init = 7 + Num_Section + 1;
Row_End = 7 + Num_Section + Num_Point;

Data_SectionPoint = cell2mat(Data_Excel(Row_Init:Row_End, 1:3));


%% 节点坐标
Data_Excel = readcell([Dir_Catia,'\',FileName_Catia, '.xls']);

for nSubline = 1:length(SublineName_Total)
    
    Row_Init = 7 + Num_Section + 1;
    Row_End = Row_Init + (Data_Num.Span(nSubline) * 2 * Data_Num.Divide(nSubline) * Data_Num.Layer(nSubline) - 1);

    nPoint = 0;
    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline)

                    nPoint = nPoint + 1;
                    nColumn_Init = 5 + (nSubline-1) * 4;
                    
                    Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 1) = cell2mat(Data_Excel(Row_Init - 1 + nPoint, nColumn_Init));
                    Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 2) = cell2mat(Data_Excel(Row_Init - 1 + nPoint, nColumn_Init + 1));
                    Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 3) = cell2mat(Data_Excel(Row_Init - 1 + nPoint, nColumn_Init + 2));

                end
            end
        end
    end


    Row_Init = Row_End + 2;
    Row_End = Row_Init + (Data_Num.Tip_x(nSubline) * Data_Num.Tip_y(nSubline) * Data_Num.Tip_z(nSubline) - 1);

    nPoint = 0;
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline)

                nPoint = nPoint + 1;
                nColumn_Init = 5 + (nSubline-1) * 4;

                Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 1) = cell2mat(Data_Excel(Row_Init - 1 + nPoint, nColumn_Init));
                Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 2) = cell2mat(Data_Excel(Row_Init - 1 + nPoint, nColumn_Init + 1));
                Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 3) = cell2mat(Data_Excel(Row_Init - 1 + nPoint, nColumn_Init + 2));

            end
        end
    end

end


%% 节点间距
%x
for nSubline = 1:length(SublineName_Total)

    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline) - 1
                for nLayer = 1:Data_Num.Layer(nSubline)
                    
                    XYZ_1 = [Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 1),...
                             Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 2),...
                             Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 3)];

                    XYZ_2 = [Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide + 1, nLayer, 1),...
                             Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide + 1, nLayer, 2),...
                             Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide + 1, nLayer, 3)];

                    Data_EL(nSubline).Main.x(nSpan, nDU, nDivide, nLayer) = pdist([XYZ_1; XYZ_2]);

                end

            end
        end
    end

    for nTip_x = 1:Data_Num.Tip_x(nSubline) - 1
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline)
                
                XYZ_1 = [Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 1),...
                         Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 2),...
                         Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 3)];

                XYZ_2 = [Data_BlockNode(nSubline).Tip(nTip_x + 1, nTip_y, nTip_z, 1),...
                         Data_BlockNode(nSubline).Tip(nTip_x + 1, nTip_y, nTip_z, 2),...
                         Data_BlockNode(nSubline).Tip(nTip_x + 1, nTip_y, nTip_z, 3)];

                Data_EL(nSubline).Tip.x(nTip_x, nTip_y, nTip_z) = pdist([XYZ_1; XYZ_2]);
                   
            end
        end
    end

end

%y
for nSubline = 1:length(SublineName_Total)

    for nSpan = 1:Data_Num.Span(nSubline) - 1
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline)
                    
                    XYZ_1 = [Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 1),...
                             Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 2),...
                             Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 3)];

                    XYZ_2 = [Data_BlockNode(nSubline).Main(nSpan + 1, nDU, nDivide, nLayer, 1),...
                             Data_BlockNode(nSubline).Main(nSpan + 1, nDU, nDivide, nLayer, 2),...
                             Data_BlockNode(nSubline).Main(nSpan + 1, nDU, nDivide, nLayer, 3)];

                    Data_EL(nSubline).Main.y(nSpan, nDU, nDivide, nLayer) = pdist([XYZ_1; XYZ_2]);

                end
            end
        end
    end
    
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline) - 1
            for nTip_z = 1:Data_Num.Tip_z(nSubline)
                
                XYZ_1 = [Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 1),...
                         Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 2),...
                         Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 3)];

                XYZ_2 = [Data_BlockNode(nSubline).Tip(nTip_x, nTip_y + 1, nTip_z, 1),...
                         Data_BlockNode(nSubline).Tip(nTip_x, nTip_y + 1, nTip_z, 2),...
                         Data_BlockNode(nSubline).Tip(nTip_x, nTip_y + 1, nTip_z, 3)];

                Data_EL(nSubline).Tip.y(nTip_x, nTip_y, nTip_z) = pdist([XYZ_1; XYZ_2]);
                   
            end
        end
    end

end

%z
for nSubline = 1:length(SublineName_Total)

    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline) - 1
                    
                    XYZ_1 = [Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 1),...
                             Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 2),...
                             Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer, 3)];

                    XYZ_2 = [Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer + 1, 1),...
                             Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer + 1, 2),...
                             Data_BlockNode(nSubline).Main(nSpan, nDU, nDivide, nLayer + 1, 3)];

                    Data_EL(nSubline).Main.z(nSpan, nDU, nDivide, nLayer) = pdist([XYZ_1; XYZ_2]);

                end

            end
        end
    end
    
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline) - 1
                
                XYZ_1 = [Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 1),...
                         Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 2),...
                         Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z, 3)];

                XYZ_2 = [Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z + 1, 1),...
                         Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z + 1, 2),...
                         Data_BlockNode(nSubline).Tip(nTip_x, nTip_y, nTip_z + 1, 3)];

                Data_EL(nSubline).Tip.z(nTip_x, nTip_y, nTip_z) = pdist([XYZ_1; XYZ_2]);
                   
            end
        end
    end

end


end


%%
function Data_SN = GetSerialNumber(Data_Excel, Data_Num, SublineName_Total, IndexName, nSubline)  

for nColumn = 1:size(Data_Excel,2)
    if strcmp(Data_Excel{1, nColumn}, [SublineName_Total{nSubline},'_',IndexName])
        nColumn_Local = nColumn;
        break;
    end
end

nGroup = 0;
for nRow = 5:size(Data_Excel,1)
    if ~ismissing(Data_Excel{nRow, nColumn_Local})
        if contains(Data_Excel{nRow, nColumn_Local}, '-')
            
            nGroup = nGroup + 1;
            InitRow_PreGroup_Total(nGroup, 1) = nRow;
            
            Index = strfind(Data_Excel{nRow, nColumn_Local}, '-');
            
            if strcmp(IndexName, 'Point_Main') || strcmp(IndexName, 'Node_Main') ||...
               strcmp(IndexName, 'Edge_y_Main')
                
                Index_ArrayPos_Total(nGroup, 1) = str2num(Data_Excel{nRow, nColumn_Local}(1:Index(1)-1));
                Index_ArrayPos_Total(nGroup, 2) = str2num(Data_Excel{nRow, nColumn_Local}(Index(1)+1:Index(2)-1));
                Index_ArrayPos_Total(nGroup, 3) = str2num(Data_Excel{nRow, nColumn_Local}(Index(2)+1:end));

            elseif strcmp(IndexName, 'Curve_x_Main') || strcmp(IndexName, 'Curve_y_Main') || strcmp(IndexName, 'Curve_z_Main') ||...
                   strcmp(IndexName, 'Point_Tip') || strcmp(IndexName, 'Node_Tip')
                
                Index_ArrayPos_Total(nGroup, 1) = str2num(Data_Excel{nRow, nColumn_Local}(1:Index(1)-1));
                Index_ArrayPos_Total(nGroup, 2) = str2num(Data_Excel{nRow, nColumn_Local}(Index(1)+1:end));
            
            elseif strcmp(IndexName, 'Curve_x_Tip') || strcmp(IndexName, 'Curve_y_Tip') || strcmp(IndexName, 'Curve_z_Tip') ||...
                   strcmp(IndexName, 'Surface_PLANE') || strcmp(IndexName, 'Surface_BOI_FLUID') || strcmp(IndexName, 'Surface_IN') ||...
                   strcmp(IndexName, 'Surface_SYMM') 
                
                Index_ArrayPos_Total(nGroup, 1) = str2num(Data_Excel{nRow, nColumn_Local}(1:Index(1)-1));
                
            end
            
        end
    end
end

if strcmp(IndexName, 'Point_Main') || strcmp(IndexName, 'Node_Main')
    
    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nDivide = 1:Data_Num.Divide(nSubline)

            nRow_Init = InitRow_PreGroup_Total(nGroup);

            nSpan = Index_ArrayPos_Total(nGroup, 1);
            nDU = Index_ArrayPos_Total(nGroup, 2);
            nLayer = Index_ArrayPos_Total(nGroup, 3);
            
            if strcmp(IndexName, 'Point_Main')
                K = 2;
            elseif strcmp(IndexName, 'Node_Main')
                K = 1;
            end
            
            Data_SN(nSpan, nDU, nDivide, nLayer) = Data_Excel(nRow_Init + nDivide * K, nColumn_Local);

        end
    end
             
elseif strcmp(IndexName, 'Curve_x_Main')

    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nSpan = 1:Data_Num.Span(nSubline)
            
            nRow_Init = InitRow_PreGroup_Total(nGroup);
            
            nDU = Index_ArrayPos_Total(nGroup, 1);
            nLayer = Index_ArrayPos_Total(nGroup, 2);
                       
            Data_SN(nSpan, nDU, nLayer) = Data_Excel(nRow_Init + nSpan * 2, nColumn_Local);

        end
    end
    
elseif strcmp(IndexName, 'Curve_y_Main')

    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nDivide = 1:Data_Num.Divide(nSubline)
            
            nRow_Init = InitRow_PreGroup_Total(nGroup);
            
            nDU = Index_ArrayPos_Total(nGroup, 1);
            nLayer = Index_ArrayPos_Total(nGroup, 2);
                       
            Data_SN(nDU, nDivide, nLayer) = Data_Excel(nRow_Init + nDivide * 2, nColumn_Local);

        end
    end

elseif strcmp(IndexName, 'Curve_z_Main')

    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nDivide = 1:Data_Num.Divide(nSubline)
            
            nRow_Init = InitRow_PreGroup_Total(nGroup);
            
            nSpan = Index_ArrayPos_Total(nGroup, 1);
            nDU = Index_ArrayPos_Total(nGroup, 2);
                       
            Data_SN(nSpan, nDU, nDivide) = Data_Excel(nRow_Init + nDivide * 2, nColumn_Local);

        end
    end
 
elseif strcmp(IndexName, 'Edge_y_Main')

    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nSpan = 1:Data_Num.Span(nSubline) - 1
            
            nRow_Init = InitRow_PreGroup_Total(nGroup);
            
            nDU = Index_ArrayPos_Total(nGroup, 1);
            nDivide = Index_ArrayPos_Total(nGroup, 2);
            nLayer = Index_ArrayPos_Total(nGroup, 3);
                       
            Data_SN(nSpan, nDU, nDivide, nLayer) = Data_Excel(nRow_Init + nSpan, nColumn_Local);

        end
    end

end

if strcmp(IndexName, 'Point_Tip') || strcmp(IndexName, 'Node_Tip')

    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nTip_x = 1:Data_Num.Tip_x(nSubline)
            
            nRow_Init = InitRow_PreGroup_Total(nGroup);
            
            nTip_y = Index_ArrayPos_Total(nGroup, 1);
            nTip_z = Index_ArrayPos_Total(nGroup, 2);
            
            if strcmp(IndexName, 'Point_Tip')
                K = 2;
            elseif strcmp(IndexName, 'Node_Tip')
                K = 1;
            end
            
            Data_SN(nTip_x, nTip_y, nTip_z) = Data_Excel(nRow_Init + nTip_x * K, nColumn_Local);
            
        end
    end
    
elseif strcmp(IndexName, 'Curve_x_Tip')

    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            
            nRow_Init = InitRow_PreGroup_Total(nGroup);

            nTip_z = Index_ArrayPos_Total(nGroup, 1);

            Data_SN(nTip_y, nTip_z) = Data_Excel(nRow_Init + nTip_y * 2, nColumn_Local);
        
        end
    end
 
elseif strcmp(IndexName, 'Curve_y_Tip')

    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nTip_x = 1:Data_Num.Tip_x(nSubline)
            
            nRow_Init = InitRow_PreGroup_Total(nGroup);
            
            nTip_z = Index_ArrayPos_Total(nGroup, 1);

            Data_SN(nTip_x, nTip_z) = Data_Excel(nRow_Init + nTip_x * 2, nColumn_Local);
            
        end
    end
    
elseif strcmp(IndexName, 'Curve_z_Tip')

    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nTip_x = 1:Data_Num.Tip_x(nSubline)
            
            nRow_Init = InitRow_PreGroup_Total(nGroup);
            
            nTip_y = Index_ArrayPos_Total(nGroup, 1);

            Data_SN(nTip_x, nTip_y) = Data_Excel(nRow_Init + nTip_x * 2, nColumn_Local);
            
        end
    end
    
end

if strcmp(IndexName, 'Surface_PLANE')
    
    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nFace = 1:((Data_Num.Span(nSubline)-1)*(Data_Num.Divide(nSubline)-1)*2 + (Data_Num.Divide(nSubline)-3)*2)
            
            nRow_Init = InitRow_PreGroup_Total(nGroup);
            
            Data_SN(nFace, 1) = Data_Excel(nRow_Init + nFace, nColumn_Local);
            
        end
    end

elseif strcmp(IndexName, 'Surface_SYMM')
    
    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nFace = 1:((Data_Num.Divide(nSubline)-1)*2*(Data_Num.Layer(nSubline)-1))
            
            nRow_Init = InitRow_PreGroup_Total(nGroup);
            
            Data_SN(nFace, 1) = Data_Excel(nRow_Init + nFace, nColumn_Local);
            
        end
    end

elseif strcmp(IndexName, 'Surface_BOI_FLUID') || strcmp(IndexName, 'Surface_IN')
    
    for nGroup = 1:length(InitRow_PreGroup_Total)
        for nFace = 1:((Data_Num.Span(nSubline)-1)*(Data_Num.Divide(nSubline)-1)*2+8*2)
            
            nRow_Init = InitRow_PreGroup_Total(nGroup);
            
            Data_SN(nFace, 1) = Data_Excel(nRow_Init + nFace, nColumn_Local);
            
        end
    end

end

end
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         