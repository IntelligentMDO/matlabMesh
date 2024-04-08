function Data_MeshSize = CFD_ICEM_Gen_SolveMesh(Data_SectionPoint, Data_BlockNode, Data_EL, Data_Num,...
                                                MeshNum_x, MeshNum_y, MeshNum_z, Layer_FirstHeight, Act_LayerAdjust)


Rate_FrontBack = 1;


%% 求解_MeshNum_x_PerEdge
if Rate_FrontBack == 1
    MeshNum_x_PerEdge = MeshNum_x / 20;
elseif Rate_FrontBack == 2
    MeshNum_x_PerEdge = MeshNum_x / 24;
end


%% 求解_MeshNum_y_PerEdge
%Dihedral
Num_Section = 5;
Num_Spline = Num_Section * 2;
Num_PointPerSpline = size(Data_SectionPoint, 1) / (2 * Num_Section);

for nSection = 1:5
    X1_Lead = Data_SectionPoint((nSection - 1) * Num_PointPerSpline + 1, 1);
    Y1_Lead = Data_SectionPoint((nSection - 1) * Num_PointPerSpline + 1, 2);
    Z1_Lead = Data_SectionPoint((nSection - 1) * Num_PointPerSpline + 1, 3);
    
    X1_Trail = (Data_SectionPoint(nSection * Num_PointPerSpline, 1) + Data_SectionPoint((Num_Spline - nSection + 1) * Num_PointPerSpline, 1))/2;
    Y1_Trail = (Data_SectionPoint(nSection * Num_PointPerSpline, 2) + Data_SectionPoint((Num_Spline - nSection + 1) * Num_PointPerSpline, 2))/2;
    Z1_Trail = (Data_SectionPoint(nSection * Num_PointPerSpline, 3) + Data_SectionPoint((Num_Spline - nSection + 1) * Num_PointPerSpline, 3))/2;
        
    Dx = X1_Trail - X1_Lead;
    Dy = Y1_Trail - Y1_Lead;
    Dz = Z1_Trail - Z1_Lead;
    
    Data_QuaterPoint(nSection, 1) = X1_Lead + Dx / 4;
    Data_QuaterPoint(nSection, 2) = Y1_Lead + Dy / 4;
    Data_QuaterPoint(nSection, 3) = Z1_Lead + Dz / 4;
end

Dihedral(1) = 0;
for nSection = 2:Num_Section
    Dz = Data_QuaterPoint(nSection, 3) - Data_QuaterPoint(nSection - 1, 3);
    Dy = Data_QuaterPoint(nSection, 2) - Data_QuaterPoint(nSection - 1, 2);

    Dihedral(nSection) = atan(Dz / Dy) * 180 / pi;
end

%SpanSLength_BlockSection
for nSpan = 1:Data_Num.Span(1)
    
    nDU = 1;
    nDivide = (Data_Num.Divide(1) + 1) / 2;
    nLayer = 1;

    Span_BlockSection(nSpan) = Data_BlockNode(1).Main(nSpan, nDU, nDivide, nLayer, 2);
    
end

for nSpan = 1:13
    if nSpan <= 4

        SpanSLength_BlockSection(nSpan) = Span_BlockSection(nSpan) / cos(deg2rad(Dihedral(2)));

    elseif nSpan > 4 && nSpan <= 7

        SpanSLength_BlockSection(nSpan)...
        = SpanSLength_BlockSection(4) + (Span_BlockSection(nSpan) - Span_BlockSection(4)) / cos(deg2rad(Dihedral(3)));

    elseif nSpan > 7 && nSpan <= 10

        SpanSLength_BlockSection(nSpan)...
        = SpanSLength_BlockSection(7) + (Span_BlockSection(nSpan) - Span_BlockSection(7)) / cos(deg2rad(Dihedral(4)));

    elseif nSpan > 10 && nSpan <= 13

        SpanSLength_BlockSection(nSpan)...
        = SpanSLength_BlockSection(10) + (Span_BlockSection(nSpan) - Span_BlockSection(10)) / cos(deg2rad(Dihedral(5)));

    end  
end

%K_SpanSLength
DSLength = (SpanSLength_BlockSection(13) - SpanSLength_BlockSection(11)) * 2; 
SpanSLength_Half_Adjust = SpanSLength_BlockSection(11) + DSLength;
    
for nSpan = 1:Data_Num.Span(1) - 1
    if nSpan <= 10 
        K_PerBlock(nSpan) = (SpanSLength_BlockSection(nSpan + 1) - SpanSLength_BlockSection(nSpan)) / SpanSLength_Half_Adjust;
    elseif nSpan >= 11
        if nSpan == 11 
            K_PerBlock(nSpan) = (SpanSLength_Half_Adjust - SpanSLength_BlockSection(11)) / 2 * 1.5 / SpanSLength_Half_Adjust;
        elseif nSpan == 12
            K_PerBlock(nSpan) = (SpanSLength_Half_Adjust - SpanSLength_BlockSection(11)) / 2 * 0.5 / SpanSLength_Half_Adjust;
        end
    end 
end

Index_nSpan = [3,4; 6,7; 9,10];
for nIndex = 1:size(Index_nSpan,1)
    K_PerBlock(Index_nSpan(nIndex, 1)) = (K_PerBlock(Index_nSpan(nIndex, 1)) + K_PerBlock(Index_nSpan(nIndex, 2))) / 2;
    K_PerBlock(Index_nSpan(nIndex, 2)) = K_PerBlock(Index_nSpan(nIndex, 1));
end

%MeshNum_y_PerEdge_1
nSpan = Data_Num.Span(3) - 1;
nDU = 1;
nDivide = 1;
nLayer = Data_Num.Layer(3);

nTip_x = 1;
nTip_y = 2;
nTip_z = 2;

if Rate_FrontBack == 1
    K = 0.9;
elseif Rate_FrontBack == 2
    K = 1.1;
end
MeshNum_y_PerEdge_1(Data_Num.Span(1) - 1) = round(Data_EL(3).Main.y(nSpan, nDU, nDivide, nLayer) /...
                                                 (Data_EL(3).Tip.x(nTip_x, nTip_y, nTip_z) / MeshNum_x_PerEdge) * K);

MeshNum_y_Local = MeshNum_y - MeshNum_y_PerEdge_1(Data_Num.Span(1) - 1);

for nSpan = 1:Data_Num.Span(1) - 2
    K_Local = K_PerBlock(nSpan) / sum(K_PerBlock(1:end-1));
    MeshNum_y_PerEdge_1(nSpan) = round(K_Local * round(MeshNum_y_Local * 0.88));
end

nDU = 2;
nDivide = 5;
nLayer = 1;

Index_1 = [4,3,5, 6,7,5,8, 9,10,8];
Index_2 = [5,5,4, 5,8,6,7, 8,8,9];
for i = 1:length(Index_1)
    for nSpan = 3:Data_Num.Span(1) - 2
        MeshSize_Cache_1(nSpan) = Data_EL(1).Main.y(nSpan, nDU, nDivide, nLayer) / MeshNum_y_PerEdge_1(nSpan);
    end
    if 0.85 * MeshSize_Cache_1(Index_1(i)) > MeshSize_Cache_1(Index_2(i))
        MeshNum_y_PerEdge_1(Index_1(i)) = MeshNum_y_PerEdge_1(Index_1(i)) + 1;
    end
end

Index_nSpan = [3,6,9];
for nIndex = 1:length(Index_nSpan)
    nSpan_1 = Index_nSpan(nIndex);
    nSpan_2 = Index_nSpan(nIndex) + 1;
    
    MeshNum_Max = max(MeshNum_y_PerEdge_1(nSpan_1), MeshNum_y_PerEdge_1(nSpan_2));
    
    MeshNum_y_PerEdge_1(nSpan_1) = MeshNum_Max;
    MeshNum_y_PerEdge_1(nSpan_2) = MeshNum_Max;
end

%MeshNum_y_PerEdge_2
for nSpan = 1:Data_Num.Span(1) - 1
    MeshNum_y_PerEdge_2(nSpan) = MeshNum_y_PerEdge_1(nSpan);
end

nDU = 2;
nDivide = 5;
nLayer = 1;

for nTry = 1:20
    if nTry == 1
        L = Data_EL(1).Main.y(1, nDU, nDivide, nLayer) + Data_EL(1).Main.y(2, nDU, nDivide, nLayer);
        N = MeshNum_y - sum(MeshNum_y_PerEdge_2(3:end)) + 1;
        SpEnd_Expect = Data_EL(1).Main.y(3, nDU, nDivide, nLayer) / MeshNum_y_PerEdge_2(3);
        Sp1 = SolveSp_Exponential(SpEnd_Expect, L, N, 'Sp1');

        S = SolveSp_Exponential(Sp1, L, N, 'Total');
    
        MeshNum_y_PerEdge_2(1) = find(S >= Data_EL(1).Main.y(1, nDU, nDivide, nLayer), 1) - 1;
        MeshNum_y_PerEdge_2(2) = (N - 1) - MeshNum_y_PerEdge_2(1);
    else
        MeshNum_y_PerEdge_2(1) = MeshNum_y_PerEdge_2(1) - 1;
    end
    
    L = Data_EL(1).Main.y(2, nDU, nDivide, nLayer);
    N = MeshNum_y_PerEdge_2(2) + 1;
    SpEnd_Expect = Data_EL(1).Main.y(3, nDU, nDivide, nLayer) / MeshNum_y_PerEdge_2(3);
    Sp1 = SolveSp_Exponential(SpEnd_Expect, L, N, 'Sp1');
    
    L = Data_EL(1).Main.y(1, nDU, nDivide, nLayer);
    N = MeshNum_y_PerEdge_2(1) + 1;
    SpEnd_Expect = Sp1;
    Sp1 = SolveSp_Exponential(SpEnd_Expect, L, N, 'Sp1');
    
    if Sp1 < SpEnd_Expect/5 && MeshNum_y_PerEdge_2(1) > 1
        MeshNum_y_PerEdge_2(11) = MeshNum_y_PerEdge_2(11) + 1;
    else
        break;
    end
end


%% 求解_MeshNum_z_PerEdge
%SLength_z
nSpan = 1;
nDU = 1;
nDivide = (Data_Num.Divide(1) + 1) / 2;

SLength_z = 0;
nS = 0;
for nLayer = 1:Data_Num.Layer(1) - 1
    nS = nS + 1;
    SLength_z(nS) = Data_EL(1).Main.z(nSpan, nDU, nDivide, nLayer);
end
for nLayer = 1:Data_Num.Layer(2) - 1
    nS = nS + 1;
    SLength_z(nS) = SLength_z(nS-1) + Data_EL(2).Main.z(nSpan, nDU, nDivide, nLayer);
end
for nLayer = 1:Data_Num.Layer(3) - 1
    nS = nS + 1;
    SLength_z(nS) = SLength_z(nS-1) + Data_EL(3).Main.z(nSpan, nDU, nDivide, nLayer);
end

%MeshNum_z_PerEdge
L = SLength_z(end);
N = round(MeshNum_z * 9/10) + 1;
Sp1 = Layer_FirstHeight;
S = SolveSp_Exponential(Sp1, L, N, 'Total');

for nS = 1:length(SLength_z)
    if nS == 1  
        MeshNum_z_PerEdge(nS) = find(S >= SLength_z(nS), 1) - 1;
    elseif nS > 1 && nS < length(SLength_z)
        MeshNum_z_PerEdge(nS) = (find(S >= SLength_z(nS), 1) - 1) - (find(S >= SLength_z(nS-1), 1) - 1);
    elseif nS == length(SLength_z)
        MeshNum_z_PerEdge(nS) = N - 1 - sum(MeshNum_z_PerEdge(1:nS-1)); 
    end
end

if Rate_FrontBack == 1
    MeshNum_z_PerEdge(2) = MeshNum_z_PerEdge(2) + round((MeshNum_z - (N - 1)) * 5/10);
    MeshNum_z_PerEdge(3) = MeshNum_z_PerEdge(3) + round((MeshNum_z - (N - 1)) * 2/10);
    MeshNum_z_PerEdge(4) = MeshNum_z_PerEdge(4) + round((MeshNum_z - (N - 1)) * 1/10);
    MeshNum_z_PerEdge(5) = MeshNum_z_PerEdge(5) +...
                           ((MeshNum_z - (N - 1)) - round((MeshNum_z - (N - 1)) * 5/10) -...
                                                    round((MeshNum_z - (N - 1)) * 2/10) -...
                                                    round((MeshNum_z - (N - 1)) * 1/10));   
elseif Rate_FrontBack == 2                
    MeshNum_z_PerEdge(2) = MeshNum_z_PerEdge(2) + round((MeshNum_z - (N - 1)) * 5/10);
    MeshNum_z_PerEdge(3) = MeshNum_z_PerEdge(3) + round((MeshNum_z - (N - 1)) * 3/10);
    MeshNum_z_PerEdge(5) = MeshNum_z_PerEdge(5) +...
                           ((MeshNum_z - (N - 1)) - round((MeshNum_z - (N - 1)) * 5/10) - round((MeshNum_z - (N - 1)) * 3/10)); 
end


%% Subline_01_x
%等效求解曲线长度
for nSpan = 1:Data_Num.Span(1)
    for nDU = 1:2
        for nLayer = 1:Data_Num.Layer(1)

            SLength_x(nSpan, nDU, nLayer) = 0;
            
            for nDivide = 2:Data_Num.Divide(1) - 1
                
                SLength_x(nSpan, nDU, nLayer) = SLength_x(nSpan, nDU, nLayer) + Data_EL(1).Main.x(nSpan, nDU, nDivide, nLayer);

            end
            
        end
    end
end

%赋值
% F = @(x) 1 - cos(x * pi);

K = 1.3;
F = @(x) -(1 - cos( (1-x) * pi)).^K  + 2^K;

x0_Total = F(linspace(0,1,(MeshNum_x_PerEdge * 8 + 1)));
x0_Total = x0_Total ./ (max(x0_Total));

x0_Didive = F(linspace(0,1,9));
x0_Didive = x0_Didive ./ (max(x0_Didive));
    
for nSpan = 1:Data_Num.Span(1)
    for nDU = 1:2
        for nLayer = 1:Data_Num.Layer(1)
            for nDivide = 1:Data_Num.Divide(1) - 1
                
                Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Law = 'biexponential';
                
                if nDivide == 1 || nDivide == Data_Num.Divide(1) - 1
                    
                    Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Nodes = MeshNum_x_PerEdge * Rate_FrontBack + 1;
                    
                    Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Size(1)...
                    = Data_EL(1).Main.x(nSpan, nDU, nDivide, nLayer) / (MeshNum_x_PerEdge * Rate_FrontBack);
                
                    Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Size(2)...
                    = Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Size(1);
                
                else
                    
                    Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Nodes = MeshNum_x_PerEdge + 1;
                    
                    if nLayer == 1
                        
                        Index_1 = find(x0_Total == x0_Didive(nDivide - 1));
                        Index_2 = find(x0_Total == x0_Didive(nDivide - 1 + 1));

                        Dx0_1 = x0_Total(Index_1 + 1) - x0_Total(Index_1);
                        Dx0_2 = x0_Total(Index_2)     - x0_Total(Index_2 - 1);
                        
                        if nDivide == 2
                            Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Size(1)...
                            = Data_MeshSize(1).Main.x(nSpan, nDU, 1, nLayer).Size(2);
                        else
                            Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Size(1) = Dx0_1 * SLength_x(nSpan, nDU, nLayer);
                        end
                        
                        Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Size(2) = Dx0_2 * SLength_x(nSpan, nDU, nLayer);

                    elseif nLayer == 2
                        
                        Rate = Data_EL(1).Main.x(nSpan, nDU, nDivide, nLayer) / Data_EL(1).Main.x(nSpan, nDU, nDivide, nLayer - 1);

                        Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Size(1)...
                        = Rate * Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer - 1).Size(1);
                    
                        Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Size(2)...
                        = Rate * Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer - 1).Size(2);

                    end
                    
                end
                
                Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Rate(1) = 0;
                Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Rate(2) = 0;
                
            end
        end
    end
end

%调整
for nSpan = 1:Data_Num.Span(1)
    for nDU = 1:2
        for nLayer = 1:Data_Num.Layer(1)

            if nLayer == 1
                
                Data_MeshSize(1).Main.x(nSpan, nDU, 2, nLayer).Size(1) = Data_MeshSize(1).Main.x(nSpan, nDU, 1, nLayer).Size(2);
                
            elseif nLayer == 2

                Data_MeshSize(1).Main.x(nSpan, nDU, Data_Num.Divide(1) - 2, nLayer).Size(1)...
                = Data_MeshSize(1).Main.x(nSpan, nDU, Data_Num.Divide(1) - 3, nLayer).Size(2);
            
            end
            
            Data_MeshSize(1).Main.x(nSpan, nDU, Data_Num.Divide(1) - 2, nLayer).Size(2)...
            = Data_MeshSize(1).Main.x(nSpan, nDU, Data_Num.Divide(1) - 1, nLayer).Size(1);

        end
    end
end

for nSpan = 1:Data_Num.Span(1)
    for nLayer = 1:Data_Num.Layer(1)
        
        MidSize = (Data_MeshSize(1).Main.x(nSpan, 1, 1, nLayer).Size(1) + ...
                   Data_MeshSize(1).Main.x(nSpan, 2, 1, nLayer).Size(1)) / 2;
        
        Data_MeshSize(1).Main.x(nSpan, 1, 1, nLayer).Size(1) = MidSize;
        Data_MeshSize(1).Main.x(nSpan, 2, 1, nLayer).Size(1) = MidSize;
                
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(1) - 2
                
                MidSize = (Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Size(2) + ...
                           Data_MeshSize(1).Main.x(nSpan, nDU, nDivide + 1, nLayer).Size(1)) / 2;
                
                Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Size(2) = MidSize;
                Data_MeshSize(1).Main.x(nSpan, nDU, nDivide + 1, nLayer).Size(1) = MidSize;
                
            end
        end
        
        MidSize = (Data_MeshSize(1).Main.x(nSpan, 1, Data_Num.Divide(1) - 1, nLayer).Size(2) + ...
                   Data_MeshSize(1).Main.x(nSpan, 2, Data_Num.Divide(1) - 1, nLayer).Size(2)) / 2;
        
        Data_MeshSize(1).Main.x(nSpan, 1, Data_Num.Divide(1) - 1, nLayer).Size(2) = MidSize;
        Data_MeshSize(1).Main.x(nSpan, 2, Data_Num.Divide(1) - 1, nLayer).Size(2) = MidSize;
        
    end
end

for nSpan = 1:Data_Num.Span(1)
    for nDU = 1:2
        for nLayer = 1:Data_Num.Layer(1)
            
            MinSize = min(Data_EL(1).Main.x(nSpan, nDU, Data_Num.Divide(1) - 2, nLayer) / MeshNum_x_PerEdge,...
                          Data_EL(1).Main.x(nSpan, nDU, Data_Num.Divide(1) - 1, nLayer) / (MeshNum_x_PerEdge * Rate_FrontBack));
                      
            Data_MeshSize(1).Main.x(nSpan, nDU, Data_Num.Divide(1) - 2, nLayer).Size(2) = MinSize;
            Data_MeshSize(1).Main.x(nSpan, nDU, Data_Num.Divide(1) - 1, nLayer).Size(1) = MinSize;
            
        end
    end
end


%% Subline_01_y
%赋值
for nDU = 1:2
    for nDivide = 1:Data_Num.Divide(1)
        for nLayer = 1:Data_Num.Layer(1)
            for nSpan = 1:Data_Num.Span(1) - 1
                
                if ismember(nSpan, [3,4,6,7,9,10])
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Law = 'biexponential';
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Nodes = MeshNum_y_PerEdge_2(nSpan) + 1;
                    
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(1)...
                    = Data_EL(1).Main.y(nSpan, nDU, nDivide, nLayer) / MeshNum_y_PerEdge_2(nSpan);

                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                    = Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(1);
                
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Rate(1) = 0;
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Rate(2) = 0;
                end
                
            end
        end
    end
end
                 
nSpan = 12;
nDU = 1;
nDivide = (Data_Num.Divide(1) + 1) / 2;
nLayer = 1;

nTip_x = (Data_Num.Tip_x(1) + 1) / 2;
nTip_y = nDU;
nTip_z = nLayer;

L = Data_EL(1).Main.y(nSpan, nDU, nDivide, nLayer);
N = MeshNum_y_PerEdge_2(nSpan) + 1;
Sp1 = Data_EL(1).Tip.y(nTip_x, nTip_y, nTip_z) / (MeshNum_x_PerEdge * Rate_FrontBack);
SpEnd = SolveSp_Exponential(Sp1, L, N, 'SpEnd');
            
for nDU = 1:2
    for nLayer = 1:Data_Num.Layer(1)
        for nDivide = 1:Data_Num.Divide(1)

            Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Law = 'biexponential';
            Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Nodes = MeshNum_y_PerEdge_2(nSpan) + 1;

            Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(1) = SpEnd;
            
            if nDivide == 1 || nDivide == 2
                
                Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                = (Data_MeshSize(1).Main.x(Data_Num.Span(1), 1, 1, nLayer).Size(1) +...
                   Data_MeshSize(1).Main.x(Data_Num.Span(1), 2, 1, nLayer).Size(2)) / 2;
            
            elseif nDivide == Data_Num.Divide(1) || nDivide == Data_Num.Divide(1) - 1
                
                Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                = (Data_MeshSize(1).Main.x(Data_Num.Span(1), 1, Data_Num.Divide(1) - 1, nLayer).Size(1) +...
                   Data_MeshSize(1).Main.x(Data_Num.Span(1), 2, Data_Num.Divide(1) - 1, nLayer).Size(2)) / 2;
            
            else
                
                nTip_x = nDivide - 1;
                nTip_y = nDU;
                nTip_z = nLayer;

                Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                = Data_EL(1).Tip.y(nTip_x, nTip_y, nTip_z) / (MeshNum_x_PerEdge * Rate_FrontBack);
            
            end
            
            Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Rate(1) = 0;
            Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Rate(2) = 0;       
 
        end
    end
end

for nDU = 1:2
    for nDivide = 1:Data_Num.Divide(1)
        for nLayer = 1:Data_Num.Layer(1)
            for nSpan = 1:Data_Num.Span(1) - 1
                
                if ismember(nSpan, [5,8,11])
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Law = 'biexponential';
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Nodes = MeshNum_y_PerEdge_2(nSpan) + 1;

                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(1)...
                    = Data_MeshSize(1).Main.y(nSpan - 1, nDU, nDivide, nLayer).Size(2);

                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                    = Data_MeshSize(1).Main.y(nSpan + 1, nDU, nDivide, nLayer).Size(1);

                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Rate(1) = 0;
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Rate(2) = 0;
                end
                
            end
        end
    end
end

for nDU = 1:2
    for nDivide = 1:Data_Num.Divide(1)
        for nLayer = 1:Data_Num.Layer(1)
            for nSpan = 1:Data_Num.Span(1) - 1
                
                if ismember(nSpan, 2)
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Law = 'biexponential';
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Nodes = MeshNum_y_PerEdge_2(nSpan) + 1;
    
                    L = Data_EL(1).Main.y(nSpan, nDU, nDivide, nLayer);
                    N = Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Nodes;
                    SpEnd_Expect = Data_MeshSize(1).Main.y(nSpan + 1, nDU, nDivide, nLayer).Size(1);
                    Sp1 = SolveSp_Exponential(SpEnd_Expect, L, N, 'Sp1');
    
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(1)...
                    = Sp1;

                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                    = SpEnd_Expect;

                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Rate(1) = 0;
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Rate(2) = 0;
                end
                
            end
        end
    end
end

for nDU = 1:2
    for nDivide = 1:Data_Num.Divide(1)
        for nLayer = 1:Data_Num.Layer(1)
            for nSpan = 1:Data_Num.Span(1) - 1
                
                if ismember(nSpan, 1)
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Law = 'biexponential';
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Nodes = MeshNum_y_PerEdge_2(nSpan) + 1;
                    
                    L = Data_EL(1).Main.y(nSpan, nDU, nDivide, nLayer);
                    N = Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Nodes;
                    SpEnd_Expect = Data_MeshSize(1).Main.y(nSpan + 1, nDU, nDivide, nLayer).Size(1);
                    Sp1 = SolveSp_Exponential(SpEnd_Expect, L, N, 'Sp1');
                    
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(1)...
                    = Sp1;

                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                    = SpEnd_Expect;

                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Rate(1) = 0;
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Rate(2) = 0;
                    
                    if nLayer == 2 && N == 3
                        Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(1)...
                        = max(Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, 1).Size(1),...
                              Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, 2).Size(1));
                    end
                end

            end
        end
    end
end

%调整
for nDU = 1:2
    for nDivide = 1:Data_Num.Divide(1)
        for nLayer = 1:Data_Num.Layer(1)
            for nSpan = 1:Data_Num.Span(1) - 1
                
                if ismember(nSpan, [3,6,9])
                    MidSize = (Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(2) + ...
                               Data_MeshSize(1).Main.y(nSpan + 1, nDU, nDivide, nLayer).Size(1)) / 2;

                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, nLayer).Size(2) = MidSize;
                    Data_MeshSize(1).Main.y(nSpan + 1, nDU, nDivide, nLayer).Size(1) = MidSize;
                end
                
            end
        end
    end
end

for nDU = 1:2
    for nSpan = 1:Data_Num.Span(1) - 1

        if ismember(nSpan, [2,5,8,11])
            Data_MeshSize(1).Main.y(nSpan, nDU, Data_Num.Divide(1)-2, 1).Size...
            = Data_MeshSize(1).Main.y(nSpan, nDU, Data_Num.Divide(1)-1, 1).Size;
        end

    end
end

for nDU = 1:2
    for nDivide = 1:Data_Num.Divide(1)
        for nLayer = 1:Data_Num.Layer(1)
            for nSpan = 1:Data_Num.Span(1) - 1
                
                if ismember(nSpan, [2,5,8,11])
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, 2).Size(1)...
                    = Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, 1).Size(1);
                    
                    Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, 2).Size(2)...
                    = Data_MeshSize(1).Main.y(nSpan, nDU, nDivide, 1).Size(2);
                end
                
            end
        end
    end
end

for nDU = 1:2
    for nDivide = 1:Data_Num.Divide(1)
            
        if nDivide == 1
            
            nTip_x = 1;
            nTip_y = 2;
            nTip_z = Data_Num.Tip_z(1);
            
            Data_MeshSize(1).Main.y(Data_Num.Span(1) - 1, nDU, nDivide, Data_Num.Layer(1)).Size(2)...
            = Data_EL(1).Tip.x(nTip_x, nTip_y, nTip_z) / MeshNum_x_PerEdge;

        elseif nDivide == Data_Num.Divide(1)
            
            nTip_x = Data_Num.Tip_x(1) - 1;
            nTip_y = 2;
            nTip_z = Data_Num.Tip_z(1);
            
            Data_MeshSize(1).Main.y(Data_Num.Span(1) - 1, nDU, nDivide, Data_Num.Layer(1)).Size(2)...
            = Data_EL(1).Tip.x(nTip_x, nTip_y, nTip_z) / MeshNum_x_PerEdge;

        end

    end
end


%% Subline_01_z
%赋值
for nSpan = 1:Data_Num.Span(1)
    for nDU = 1:2
        for nDivide = 1:Data_Num.Divide(1)
            for nLayer = 1:Data_Num.Layer(1) - 1
                
                Data_MeshSize(1).Main.z(nSpan, nDU, nDivide, nLayer).Law = 'biexponential';
                Data_MeshSize(1).Main.z(nSpan, nDU, nDivide, nLayer).Nodes = MeshNum_z_PerEdge(nLayer) + 1;
                
                Data_MeshSize(1).Main.z(nSpan, nDU, nDivide, nLayer).Size(1) = Layer_FirstHeight;
                
                L = Data_EL(1).Main.z(nSpan, nDU, nDivide, nLayer);
                N = Data_MeshSize(1).Main.z(nSpan, nDU, nDivide, nLayer).Nodes;
                Sp1 = Data_MeshSize(1).Main.z(nSpan, nDU, nDivide, nLayer).Size(1);
                SpEnd = SolveSp_Exponential(Sp1, L, N, 'SpEnd');
                
                Data_MeshSize(1).Main.z(nSpan, nDU, nDivide, nLayer).Size(2) = SpEnd;
                
                Data_MeshSize(1).Main.z(nSpan, nDU, nDivide, nLayer).Rate(1) = 0;
                Data_MeshSize(1).Main.z(nSpan, nDU, nDivide, nLayer).Rate(2) = 0;

            end
        end
    end
end

%调整
if Act_LayerAdjust == 1
    for nSpan = 1:Data_Num.Span(1)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(1)
                for nLayer = 1:Data_Num.Layer(1) - 1

                    if nDivide == Data_Num.Divide(1) - 1 && nLayer == 1
                        
                        Data_MeshSize(1).Main.z(nSpan, nDU, nDivide, nLayer).Size(1) = Layer_FirstHeight * 3;
                        
                    elseif nDivide == Data_Num.Divide(1) && nLayer == 1
                        
                        Data_MeshSize(1).Main.z(nSpan, nDU, nDivide, nLayer).Size(1) = Layer_FirstHeight * 6;
                        
                    end

                end
            end
        end
    end
end


%% Subline_01_Tip_x
%赋值
for nTip_x = 1:Data_Num.Tip_x(1) - 1
    for nTip_y = 1:Data_Num.Tip_y(1)
        for nTip_z = 1:Data_Num.Tip_z(1)
            
            if nTip_y == 2
                
                Data_MeshSize(1).Tip.x(nTip_x, nTip_y, nTip_z).Law = 'biexponential';

                nSpan = Data_Num.Span(1);
                nDU = 1;
                nDivide = nTip_x + 1;
                nLayer = nTip_z;

                Data_MeshSize(1).Tip.x(nTip_x, nTip_y, nTip_z).Nodes = Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Nodes;
                
                Data_MeshSize(1).Tip.x(nTip_x, nTip_y, nTip_z).Size(1)...
                = (Data_MeshSize(1).Main.x(nSpan, 1, nDivide, nLayer).Size(1) + ...
                   Data_MeshSize(1).Main.x(nSpan, 2, nDivide, nLayer).Size(1)) / 2;

                Data_MeshSize(1).Tip.x(nTip_x, nTip_y, nTip_z).Size(2)...
                = (Data_MeshSize(1).Main.x(nSpan, 1, nDivide, nLayer).Size(2) + ...
                   Data_MeshSize(1).Main.x(nSpan, 2, nDivide, nLayer).Size(2)) / 2;

                Data_MeshSize(1).Tip.x(nTip_x, nTip_y, nTip_z).Rate(1) = 0;
                Data_MeshSize(1).Tip.x(nTip_x, nTip_y, nTip_z).Rate(2) = 0;
            
            else
                
                Data_MeshSize(1).Tip.x(nTip_x, nTip_y, nTip_z).Law = NaN;
                
            end
        
        end
    end
end

%调整
for nTip_x = 1:Data_Num.Tip_x(1) - 1
    for nTip_y = 1:Data_Num.Tip_y(1)
        for nTip_z = 1:Data_Num.Tip_z(1)
            
            if nTip_y == 2
                
                nSpan = Data_Num.Span(1) - 1;
                nLayer = nTip_z;
                
                if nTip_x == 1

                    Data_MeshSize(1).Tip.x(nTip_x, nTip_y, nTip_z).Size(1)...
                    = Data_MeshSize(1).Main.y(nSpan, 1, 1, nLayer).Size(2);
                    
                elseif nTip_x == Data_Num.Tip_x(1) - 1
                    
                    Data_MeshSize(1).Tip.x(nTip_x, nTip_y, nTip_z).Size(2)...
                    = Data_MeshSize(1).Main.y(nSpan, 1,  Data_Num.Divide(1), nLayer).Size(2);

                end
                
            end
            
            if nTip_y == 2 && nTip_z == Data_Num.Tip_z(1)
                
                if nTip_x == 2
                    
                    MidSize = (Data_EL(1).Tip.x(nTip_x - 1, nTip_y, nTip_z) / MeshNum_x_PerEdge + ...
                               Data_EL(1).Tip.x(nTip_x, nTip_y, nTip_z) / MeshNum_x_PerEdge) / 2;
                   
                    Data_MeshSize(1).Tip.x(nTip_x - 1, nTip_y, nTip_z).Size(2) = MidSize;
                    Data_MeshSize(1).Tip.x(nTip_x, nTip_y, nTip_z).Size(1) = MidSize;
                    
                elseif nTip_x == Data_Num.Tip_x(1) - 2
                    
                    MidSize = (Data_EL(1).Tip.x(nTip_x, nTip_y, nTip_z) / MeshNum_x_PerEdge + ...
                               Data_EL(1).Tip.x(nTip_x + 1, nTip_y, nTip_z) / MeshNum_x_PerEdge) / 2;
                    
                    Data_MeshSize(1).Tip.x(nTip_x, nTip_y, nTip_z).Size(2) = MidSize;
                    Data_MeshSize(1).Tip.x(nTip_x + 1, nTip_y, nTip_z).Size(1) = MidSize;
                    
                end
  
            end

        end
    end
end


%% Subline_01_Tip_y
%赋值
for nTip_x = 1:Data_Num.Tip_x(1)
    for nTip_y = 1:Data_Num.Tip_y(1) - 1
        for nTip_z = 1:Data_Num.Tip_z(1)
            
            if nTip_x ~= 1 && nTip_x ~= Data_Num.Tip_x(1)
                
                Data_MeshSize(1).Tip.y(nTip_x, nTip_y, nTip_z).Law = 'biexponential';

                nSpan = Data_Num.Span(1);
                nDU = 1;
                nDivide = 1;
                nLayer = nTip_z;

                Data_MeshSize(1).Tip.y(nTip_x, nTip_y, nTip_z).Nodes = Data_MeshSize(1).Main.x(nSpan, nDU, nDivide, nLayer).Nodes;

                Data_MeshSize(1).Tip.y(nTip_x, nTip_y, nTip_z).Size(1)...
                = Data_EL(1).Tip.y(nTip_x, nTip_y, nTip_z) / (Data_MeshSize(1).Tip.y(nTip_x, nTip_y, nTip_z).Nodes - 1);

                Data_MeshSize(1).Tip.y(nTip_x, nTip_y, nTip_z).Size(2)...
                = Data_MeshSize(1).Tip.y(nTip_x, nTip_y, nTip_z).Size(1);

                Data_MeshSize(1).Tip.y(nTip_x, nTip_y, nTip_z).Rate(1) = 0;
                Data_MeshSize(1).Tip.y(nTip_x, nTip_y, nTip_z).Rate(2) = 0;
            
            else
                
                Data_MeshSize(1).Tip.y(nTip_x, nTip_y, nTip_z).Law = NaN;
                
            end
        
        end
    end
end

%调整
for nTip_x = 1:Data_Num.Tip_x(1)
    for nTip_z = 1:Data_Num.Tip_z(1)
        
        if nTip_x ~= 1 && nTip_x ~= Data_Num.Tip_x(1)
            
            MidSize = (Data_MeshSize(1).Tip.y(nTip_x, 1, nTip_z).Size(2) + ...
                       Data_MeshSize(1).Tip.y(nTip_x, 2, nTip_z).Size(1)) / 2;

            Data_MeshSize(1).Tip.y(nTip_x, 1, nTip_z).Size(2) = MidSize;
            Data_MeshSize(1).Tip.y(nTip_x, 2, nTip_z).Size(1) = MidSize;
        
        end
        
    end
end


%% Subline_01_Tip_z
%赋值
for nTip_x = 1:Data_Num.Tip_x(1)
    for nTip_y = 1:Data_Num.Tip_y(1)
        for nTip_z = 1:Data_Num.Tip_z(1) - 1
            
            if nTip_x ~= 1 && nTip_x ~= Data_Num.Tip_x(1) && nTip_y == 2
                
                Data_MeshSize(1).Tip.z(nTip_x, nTip_y, nTip_z).Law = 'biexponential';

                nSpan = Data_Num.Span(1);
                nDU = 1;
                nDivide = 1;
                nLayer = 1;

                Data_MeshSize(1).Tip.z(nTip_x, nTip_y, nTip_z).Nodes = Data_MeshSize(1).Main.z(nSpan, nDU, nDivide, nLayer).Nodes;

                Data_MeshSize(1).Tip.z(nTip_x, nTip_y, nTip_z).Size(1) = Layer_FirstHeight;

                L = Data_EL(1).Tip.z(nTip_x, nTip_y, nTip_z);
                N = Data_MeshSize(1).Tip.z(nTip_x, nTip_y, nTip_z).Nodes;
                Sp1 = Data_MeshSize(1).Tip.z(nTip_x, nTip_y, nTip_z).Size(1);
                SpEnd = SolveSp_Exponential(Sp1, L, N, 'SpEnd');

                Data_MeshSize(1).Tip.z(nTip_x, nTip_y, nTip_z).Size(2) = SpEnd;

                Data_MeshSize(1).Tip.z(nTip_x, nTip_y, nTip_z).Rate(1) = 0;
                Data_MeshSize(1).Tip.z(nTip_x, nTip_y, nTip_z).Rate(2) = 0;
            
            else
                
                Data_MeshSize(1).Tip.z(nTip_x, nTip_y, nTip_z).Law = NaN;
                
            end
        
        end
    end
end


%% Subline_12&23_x
for nSubline = 2:3
    %赋值
    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline) - 1
                for nLayer = 1:Data_Num.Layer(nSubline)

                    if nLayer > 1

                        Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Law = 'biexponential';

                        if nDivide == 1 || nDivide == Data_Num.Divide(nSubline) - 1

                            Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Nodes = MeshNum_x_PerEdge * Rate_FrontBack + 1;

                            Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Size(1)...
                            = Data_EL(nSubline).Main.x(nSpan, nDU, nDivide, nLayer) / (MeshNum_x_PerEdge * Rate_FrontBack);

                        else

                            Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Nodes = MeshNum_x_PerEdge + 1;

                            Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Size(1)...
                            = Data_EL(nSubline).Main.x(nSpan, nDU, nDivide, nLayer) / MeshNum_x_PerEdge;

                        end

                        Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Size(2)...
                        = Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Size(1);

                        Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Rate(1) = 0;
                        Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Rate(2) = 0;

                    else

                        Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Law = NaN;

                    end

                end
            end
        end
    end

    %调整
    for nSpan = 1:Data_Num.Span(nSubline)
        for nLayer = 1:Data_Num.Layer(nSubline)

            if nLayer > 1

                MidSize = (Data_MeshSize(nSubline).Main.x(nSpan, 1, 1, nLayer).Size(1) + ...
                           Data_MeshSize(nSubline).Main.x(nSpan, 2, 1, nLayer).Size(1)) / 2;

                Data_MeshSize(nSubline).Main.x(nSpan, 1, 1, nLayer).Size(1) = MidSize;
                Data_MeshSize(nSubline).Main.x(nSpan, 2, 1, nLayer).Size(1) = MidSize;

                for nDU = 1:2
                    for nDivide = 1:Data_Num.Divide(nSubline) - 2

                        MidSize = (Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Size(2) + ...
                                   Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide + 1, nLayer).Size(1)) / 2;

                        Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Size(2) = MidSize;
                        Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide + 1, nLayer).Size(1) = MidSize;

                    end
                end

                MidSize = (Data_MeshSize(nSubline).Main.x(nSpan, 1, Data_Num.Divide(nSubline) - 1, nLayer).Size(2) + ...
                           Data_MeshSize(nSubline).Main.x(nSpan, 2, Data_Num.Divide(nSubline) - 1, nLayer).Size(2)) / 2;

                Data_MeshSize(nSubline).Main.x(nSpan, 1, Data_Num.Divide(nSubline) - 1, nLayer).Size(2) = MidSize;
                Data_MeshSize(nSubline).Main.x(nSpan, 2, Data_Num.Divide(nSubline) - 1, nLayer).Size(2) = MidSize;

            end

        end
    end
    
end


%% Subline_12&23_y
for nSubline = 2:3
    %赋值
    for nDU = 1:2
        for nDivide = 1:Data_Num.Divide(nSubline)
            for nLayer = 1:Data_Num.Layer(nSubline)
                for nSpan = 1:Data_Num.Span(nSubline) - 1

                    if ismember(nSpan, [1,2,3,4])
                        if nLayer > 1 

                            Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Law = 'biexponential';

                            if nSpan <= 3
                                
                                Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Nodes...
                                = (Data_MeshSize(1).Main.y((nSpan-1)*3+1, nDU, nDivide, 1).Nodes - 1) +...
                                  (Data_MeshSize(1).Main.y((nSpan-1)*3+2, nDU, nDivide, 1).Nodes - 1) +...
                                  (Data_MeshSize(1).Main.y((nSpan-1)*3+3, nDU, nDivide, 1).Nodes - 1) + 1;
                              
                            elseif nSpan == 4
                                
                                Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Nodes...
                                = (Data_MeshSize(1).Main.y(Data_Num.Span(1)-3, nDU, nDivide, 1).Nodes - 1) +...
                                  (Data_MeshSize(1).Main.y(Data_Num.Span(1)-2, nDU, nDivide, 1).Nodes - 1) + 1;
                              
                            end

                            Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(1)...
                            = Data_EL(nSubline).Main.y(nSpan, nDU, nDivide, nLayer) /...
                              (Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Nodes - 1);

                            Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                            = Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(1);

                            Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Rate(1) = 0;
                            Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Rate(2) = 0;

                        else

                            Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Law = NaN;

                        end
                    end

                end
            end
        end
    end

    nSpan = 5;
    for nDU = 1:2
        for nLayer = 1:Data_Num.Layer(nSubline)
            for nDivide = 1:Data_Num.Divide(nSubline)

                if nLayer > 1 

                    Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Law = 'biexponential';
                    Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Nodes = Data_MeshSize(1).Main.y(12, nDU, nDivide, 1).Nodes;

                    Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(1)...
                    = Data_EL(nSubline).Main.y(nSpan, nDU, nDivide, nLayer) /...
                      (Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Nodes - 1);

                    if nDivide == 1 || nDivide == 2

                        Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                        = (Data_MeshSize(nSubline).Main.x(Data_Num.Span(nSubline), 1, 1, nLayer).Size(1) +...
                           Data_MeshSize(nSubline).Main.x(Data_Num.Span(nSubline), 2, 1, nLayer).Size(2)) / 2;

                    elseif nDivide == Data_Num.Divide(nSubline) || nDivide == Data_Num.Divide(nSubline) - 1

                        Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                        = (Data_MeshSize(nSubline).Main.x(Data_Num.Span(nSubline), 1, Data_Num.Divide(nSubline) - 1, nLayer).Size(1) +...
                           Data_MeshSize(nSubline).Main.x(Data_Num.Span(nSubline), 2, Data_Num.Divide(nSubline) - 1, nLayer).Size(2)) / 2;

                    else

                        nTip_x = nDivide - 1;
                        nTip_y = nDU;
                        nTip_z = nLayer;

                        Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                        = Data_EL(nSubline).Tip.y(nTip_x, nTip_y, nTip_z) / (MeshNum_x_PerEdge * Rate_FrontBack);

                    end

                    Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Rate(1) = 0;
                    Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Rate(2) = 0;

                else

                    Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Law = NaN;

                end

            end
        end
    end

    %调整
    if nSubline == 2
        
        Index_nSpan = [1,3; 4,6; 7,8; 9,NaN];
        
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline)
                    for nSpan = 1:Data_Num.Span(nSubline) - 2
                        
                        if nLayer == 2
                            if nSpan == 1
                                
                                Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(1)...
                                = Data_MeshSize(1).Main.y(Index_nSpan(nSpan, 1), nDU, nDivide, nLayer).Size(1) * 2;

                                Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                                = Data_MeshSize(1).Main.y(Index_nSpan(nSpan, 2), nDU, nDivide, nLayer).Size(2);
                            
                            elseif nSpan > 1 && nSpan < Data_Num.Span(nSubline) - 2
                                
                                Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(1)...
                                = Data_MeshSize(1).Main.y(Index_nSpan(nSpan, 1), nDU, nDivide, nLayer).Size(1);

                                Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(2)...
                                = Data_MeshSize(1).Main.y(Index_nSpan(nSpan, 2), nDU, nDivide, nLayer).Size(2);
                            
                            else
                                
                                Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(1)...
                                = Data_MeshSize(1).Main.y(Index_nSpan(nSpan, 1), nDU, nDivide, nLayer).Size(1);
                            
                            end
                        end

                    end
                end
            end
        end
    
    end
    
    for nDU = 1:2
        for nDivide = 1:Data_Num.Divide(nSubline) 
            for nLayer = 1:Data_Num.Layer(nSubline)

                if nSubline == 2 && nLayer == 2 

                    Data_MeshSize(nSubline).Main.y(4, nDU, nDivide, nLayer).Size(2)...
                    = Data_MeshSize(nSubline).Main.y(5, nDU, nDivide, nLayer).Size(1);
                
                elseif (nSubline == 2 && nLayer > 2) || (nSubline == 3 && nLayer > 1)
                    
                    Data_MeshSize(nSubline).Main.y(5, nDU, nDivide, nLayer).Size(1)...
                    = Data_MeshSize(nSubline).Main.y(4, nDU, nDivide, nLayer).Size(2);

                end
                    
            end
        end
    end
    
    for nDU = 1:2
        for nDivide = 1:Data_Num.Divide(nSubline) 
            for nLayer = 1:Data_Num.Layer(nSubline)
                for nSpan = 1:Data_Num.Span(nSubline) - 2

                    if nLayer > 1 

                        MidSize = (Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(2) + ...
                                   Data_MeshSize(nSubline).Main.y(nSpan + 1, nDU, nDivide, nLayer).Size(1)) / 2;

                        Data_MeshSize(nSubline).Main.y(nSpan, nDU, nDivide, nLayer).Size(2) = MidSize;
                        Data_MeshSize(nSubline).Main.y(nSpan + 1, nDU, nDivide, nLayer).Size(1) = MidSize;

                    end

                end
            end
        end
    end

end


%% Subline_12&23_z
for nSubline = 2:3
    
    if nSubline == 2
        Index_nSpan = [1,4,7,10,12,13];
    else
        Index_nSpan = [1,2,3,4,5,6];
    end

    %赋值
    for nSpan = 1:Data_Num.Span(nSubline)
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)
                for nLayer = 1:Data_Num.Layer(nSubline) - 1

                    Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Law = 'biexponential';
                    
                    if nSubline == 2
                        Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Nodes...
                        = MeshNum_z_PerEdge(1 + nLayer) + 1;
                    else
                        Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Nodes...
                        = MeshNum_z_PerEdge(1 + (Data_Num.Layer(nSubline - 1) - 1) + nLayer) + 1;
                    end

                    if nLayer == 1
                        Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Size(1)...
                        = Data_MeshSize(nSubline - 1).Main.z(Index_nSpan(nSpan), nDU, nDivide, Data_Num.Layer(nSubline - 1) - 1).Size(2);
                    else
                        Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Size(1)...
                        = Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer - 1).Size(2);
                    end

                    L = Data_EL(nSubline).Main.z(nSpan, nDU, nDivide, nLayer);
                    N = Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Nodes;
                    Sp1 = Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Size(1);
                    SpEnd = SolveSp_Exponential(Sp1, L, N, 'SpEnd');

                    Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Size(2) = SpEnd;

                    Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Rate(1) = 0;
                    Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Rate(2) = 0;

                end
            end
        end
    end

    %调整
    if nSubline == 2
        for nDU = 1:2
            for nDivide = 1:Data_Num.Divide(nSubline)

                MidSize = (Data_MeshSize(nSubline).Main.z(1, nDU, nDivide, Data_Num.Layer(nSubline) - 1).Size(2) +...
                           Data_MeshSize(nSubline).Main.z(Data_Num.Span(nSubline), nDU, nDivide, Data_Num.Layer(nSubline) - 1).Size(2)) / 2;

                for nSpan = 1:Data_Num.Span(nSubline)
                    Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, Data_Num.Layer(nSubline) - 1).Size(2) = MidSize;
                end

            end
        end
    end

end


%% Subline_12&23_Tip_x
for nSubline = 2:3
    %赋值
    for nTip_x = 1:Data_Num.Tip_x(nSubline) - 1
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline)

                if nTip_y == 2 && nTip_z > 1

                    Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Law = 'biexponential';

                    nSpan = Data_Num.Span(nSubline);
                    nDU = 1;
                    nDivide = nTip_x + 1;
                    nLayer = nTip_z;

                    Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Nodes...
                    = Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Nodes;

                    Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Size(1)...
                    = Data_EL(nSubline).Tip.x(nTip_x, nTip_y, nTip_z) / (Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Nodes - 1);

                    Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Size(2)...
                    = Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Size(1);

                    Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Rate(1) = 0;
                    Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Rate(2) = 0;

                else

                    Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Law = NaN;

                end

            end
        end
    end

    %调整
    for nTip_x = 1:Data_Num.Tip_x(nSubline) - 1
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline)

                if nTip_y == 2 && nTip_z > 1

                    nSpan = Data_Num.Span(nSubline) - 1;
                    nLayer = nTip_z;

                    if nTip_x == 1

                        Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Size(1)...
                        = Data_MeshSize(nSubline).Main.y(nSpan, 1, 1, nLayer).Size(2);

                    elseif nTip_x == Data_Num.Tip_x(nSubline) - 1

                        Data_MeshSize(nSubline).Tip.x(nTip_x, nTip_y, nTip_z).Size(2)...
                        = Data_MeshSize(nSubline).Main.y(nSpan, 1,  Data_Num.Divide(1), nLayer).Size(2);

                    end

                end

            end
        end
    end

end


%% Subline_12&23_Tip_y
for nSubline = 2:3
    %赋值
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline) - 1
            for nTip_z = 1:Data_Num.Tip_z(nSubline)

                if nTip_x ~= 1 && nTip_x ~= Data_Num.Tip_x(nSubline) && nTip_z > 1

                    Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Law = 'biexponential';

                    nSpan = Data_Num.Span(nSubline);
                    nDU = 1;
                    nDivide = 1;
                    nLayer = nTip_z;

                    Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Nodes...
                    = Data_MeshSize(nSubline).Main.x(nSpan, nDU, nDivide, nLayer).Nodes;

                    Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Size(1)...
                    = Data_EL(nSubline).Tip.y(nTip_x, nTip_y, nTip_z) / (Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Nodes - 1);

                    Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Size(2)...
                    = Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Size(1);

                    Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Rate(1) = 0;
                    Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Rate(2) = 0;

                else

                    Data_MeshSize(nSubline).Tip.y(nTip_x, nTip_y, nTip_z).Law = NaN;

                end

            end
        end
    end

    %调整
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_z = 1:Data_Num.Tip_z(nSubline)

            if nTip_x ~= 1 && nTip_x ~= Data_Num.Tip_x(nSubline) && nTip_z > 1

                MidSize = (Data_MeshSize(nSubline).Tip.y(nTip_x, 1, nTip_z).Size(2) + ...
                           Data_MeshSize(nSubline).Tip.y(nTip_x, 2, nTip_z).Size(1)) / 2;

                Data_MeshSize(nSubline).Tip.y(nTip_x, 1, nTip_z).Size(2) = MidSize;
                Data_MeshSize(nSubline).Tip.y(nTip_x, 2, nTip_z).Size(1) = MidSize;

            end

        end
    end
end


%% Subline_12&23_Tip_z
for nSubline = 2:3
    %赋值
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline)
            for nTip_z = 1:Data_Num.Tip_z(nSubline) - 1

                if nTip_x ~= 1 && nTip_x ~= Data_Num.Tip_x(nSubline) && nTip_y == 2

                    Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Law = 'biexponential';

                    nSpan = Data_Num.Span(nSubline);
                    nDU = 1;
                    nDivide = 1;
                    nLayer = nTip_z;
        
                    Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Nodes...
                    = Data_MeshSize(nSubline).Main.z(nSpan, nDU, nDivide, nLayer).Nodes;

                    if nTip_z == 1
                        Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Size(1)...
                        = Data_MeshSize(nSubline - 1).Tip.z(nTip_x, nTip_y, Data_Num.Tip_z(nSubline - 1) - 1).Size(2);
                    else
                        Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Size(1)...
                        = Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z - 1).Size(2);
                    end

                    L = Data_EL(nSubline).Tip.z(nTip_x, nTip_y, nTip_z);
                    N = Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Nodes;
                    Sp1 = Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Size(1);
                    SpEnd = SolveSp_Exponential(Sp1, L, N, 'SpEnd');

                    Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Size(2) = SpEnd;

                    Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Rate(1) = 0;
                    Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Rate(2) = 0;

                else

                   Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Law = NaN;

                end

            end
        end
    end

    %调整
    for nTip_x = 1:Data_Num.Tip_x(nSubline)
        for nTip_y = 1:Data_Num.Tip_y(nSubline)

            if nTip_x ~= 1 && nTip_x ~= Data_Num.Tip_x(nSubline) && nTip_y == 2

                nSpan = Data_Num.Span(nSubline);
                nDivide = nTip_x + 1;
                nLayer = Data_Num.Layer(nSubline) - 1;

                MidSize = (Data_MeshSize(nSubline).Main.z(nSpan, 1, nDivide, nLayer).Size(2) +...
                           Data_MeshSize(nSubline).Main.z(nSpan, 2, nDivide, nLayer).Size(2)) / 2;

                Data_MeshSize(nSubline).Tip.z(nTip_x, nTip_y, nTip_z).Size(2) = MidSize;

            end

        end
    end

end


end


%%
function Sp_Output = SolveSp_Exponential(Sp_Input, L, N, Method)

if strcmp(Method,'Sp1')
    
    if N == 1
        error('*** Error：SolveSp_Exponential ***');
    elseif N == 2
        Sp_Output = Sp_Input;
    else
        SpEnd_Expect = Sp_Input;
        
        syms Sp1;
        R = -log((N-1)*Sp1/L)/(N-2);
        F_Sp1 = Sp1*(N-1)*exp(R*(N-2)) - Sp1*(N-2)*exp(R*(N-3)) - SpEnd_Expect;

        Sp_Output = double(solve(F_Sp1,Sp1));
    end

elseif strcmp(Method,'SpEnd')

    if N == 1
        error('*** Error：SolveSp_Exponential ***');
    elseif N == 2
        Sp_Output = Sp_Input;
    else
        Sp1 = Sp_Input;
        R = -log((N-1)*Sp1/L)/(N-2);
        
        Sp_Output = Sp1*(N-1).*exp(R*(N-2)) - Sp1*(N-2).*exp(R*(N-3));
    end

elseif strcmp(Method,'Total')
    
    if N == 1
        error('*** Error：SolveSp_Exponential ***');
    elseif N == 2
        Sp_Output = [0,Sp_Input];
    else
        Sp1 = Sp_Input;
        N = 2:N;

        R = -log((N(end)-1)*Sp1/L)/(N(end)-2);
        Sp_PerEdge = Sp1*(N-1).*exp(R*(N-2)) - Sp1*(N-2).*exp(R*(N-3));

        Sp_Output = [0,cumsum(Sp_PerEdge)];
    end
    
end
    
end

























