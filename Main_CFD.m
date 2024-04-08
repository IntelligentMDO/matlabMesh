function Timer = Main_CFD(Dir_Total, Mode, Num_Infill) %Manual_1,DOE_2,Infill_3


if Mode == 1
    FileName_AddOn = '_Manual';
elseif Mode == 2
    FileName_AddOn = '_DOE';
elseif Mode == 3
    FileName_AddOn = '_Infill';
end

if Mode == 1                                                                %计算节点数量
    Num_Node = 1;
elseif Mode == 2
    Num_Node = 4;
elseif Mode == 3
    Num_Node = Num_Infill;
end

% p = gcp('nocreate');
% if isempty(p)
%     parpool(8);
% elseif p.NumWorkers ~= 8
%     delete(p);
%     parpool(8);
% end


%% 文件目录
Dir.Total = Dir_Total; 
Dir.Catia = [Dir.Total,'\','Catia'];
Dir.ICEM = [Dir.Total,'\','ICEM'];


%% 参数_FlightState
Para_CFD.FlightState.Mode_PT = 2;                                           %赋值模式_压强&温度(手动_1,自动(依据高度)_2)
if Para_CFD.FlightState.Mode_PT == 1                                        %***手动***
    Para_CFD.FlightState.Ref_Pressure = 101325;                             %压强(Pa)
    Para_CFD.FlightState.Ref_Temperature = 288.15;                          %温度(K)
end

Para_CFD.FlightState.Mode_HVA = [1,1];                                    %赋值模式_高度&速度(手动_1,外部_2)(Mode_PT=1时高度无效)
if Mode == 1
    Para_CFD.FlightState.Mode_HVA = ones(1,3); 
end
if Para_CFD.FlightState.Mode_PT == 1
    Para_CFD.FlightState.Mode_HVA(1) = NaN; 
end

if Para_CFD.FlightState.Mode_PT == 2 &&...                                  %***手动_高度***
   Para_CFD.FlightState.Mode_HVA(1) == 1
    Para_CFD.FlightState.Ref_Height = 0;                                    %高度(m)
end

if Para_CFD.FlightState.Mode_HVA(2) == 1                                    %***手动_速度***
    Para_CFD.FlightState.Ref_Velocity = 0.6+1e10;                           %速度(m/s,定义Ma为A+1e10)
end


%% 参数_ICEM
Para_CFD.ICEM.MeshNum_x = 20 * 12;
Para_CFD.ICEM.MeshNum_y = Para_CFD.ICEM.MeshNum_x * (9/20);
Para_CFD.ICEM.MeshNum_z = Para_CFD.ICEM.MeshNum_x * (8/20);

Para_CFD.ICEM.Layer_FirstYPlus = 1;                                         %附面层初始Y+


%% 初始化
Para_CFD.Post.FileName_Load_DesignVar...                                    %Catia生成的DesignVar文件名称
= ['CatiaR_DesignVar',FileName_AddOn];

Para_CFD.Post.FileName_Load_Geometry...                                     %Catia生成的Geometry文件名称
= ['CatiaR_Geometry',FileName_AddOn];

Para_CFD.Post.FileName_LoSa_Aero = ['CFDR_Aero',FileName_AddOn];

if Mode == 1
    if exist([Para_CFD.Post.FileName_Load_DesignVar,'.mat'],'file')
        load([Para_CFD.Post.FileName_Load_DesignVar,'.mat'])
        eval(['Data_DesignVar=', Para_CFD.Post.FileName_Load_DesignVar,';']);
    end
    
    if exist([Para_CFD.Post.FileName_Load_Geometry,'.mat'],'file')
        load([Para_CFD.Post.FileName_Load_Geometry,'.mat'])
        eval(['Data_Geometry=', Para_CFD.Post.FileName_Load_Geometry,';']);
    end
    
    Data_DebugCase = '';
end

if Mode == 2
    load([Para_CFD.Post.FileName_Load_DesignVar,'.mat'])
    eval(['Data_DesignVar=', Para_CFD.Post.FileName_Load_DesignVar,';']);
    
    load([Para_CFD.Post.FileName_Load_Geometry,'.mat'])
    eval(['Data_Geometry=', Para_CFD.Post.FileName_Load_Geometry,';']);
end

if Mode == 3
    %加载DOE数据
    Para_CFD.Post.FileName_Load_DesignVar = replace(Para_CFD.Post.FileName_Load_DesignVar,'_Infill','_DOE');
    Para_CFD.Post.FileName_Load_Geometry = replace(Para_CFD.Post.FileName_Load_Geometry,'_Infill','_DOE');
    Para_CFD.DebugCase.FileName_LoSa_DebugCase = replace(Para_CFD.DebugCase.FileName_LoSa_DebugCase,'_Infill','_DOE');
    Para_CFD.Post.FileName_LoSa_Aero = replace(Para_CFD.Post.FileName_LoSa_Aero,'_Infill','_DOE');
    
    load([Para_CFD.Post.FileName_Load_DesignVar,'.mat'])
    eval(['Data_DesignVar_DOE=', Para_CFD.Post.FileName_Load_DesignVar,';']);
    
    load([Para_CFD.Post.FileName_Load_Geometry,'.mat'])
    eval(['Data_Geometry_DOE=', Para_CFD.Post.FileName_Load_Geometry,';']);
    
    %加载Infill数据
    Para_CFD.Post.FileName_Load_DesignVar = replace(Para_CFD.Post.FileName_Load_DesignVar,'_DOE','_Infill');
    Para_CFD.Post.FileName_Load_Geometry = replace(Para_CFD.Post.FileName_Load_Geometry,'_DOE','_Infill');
    Para_CFD.DebugCase.FileName_LoSa_DebugCase = replace(Para_CFD.DebugCase.FileName_LoSa_DebugCase,'_DOE','_Infill');
    Para_CFD.Post.FileName_LoSa_Aero = replace(Para_CFD.Post.FileName_LoSa_Aero,'_DOE','_Infill');
    
    if exist([Para_CFD.Post.FileName_Load_DesignVar,'.mat'],'file')
        load([Para_CFD.Post.FileName_Load_DesignVar,'.mat'])
        eval(['Data_DesignVar_Infill=', Para_CFD.Post.FileName_Load_DesignVar,';']);
    else
        Data_DesignVar_Infill = '';
    end
    
    if exist([Para_CFD.Post.FileName_Load_Geometry,'.mat'],'file')
        load([Para_CFD.Post.FileName_Load_Geometry,'.mat'])
        eval(['Data_Geometry_Infill=', Para_CFD.Post.FileName_Load_Geometry,';']);
    else
        Data_Geometry_Infill = '';
    end

    %合并
    Data_DesignVar = [Data_DesignVar_DOE, Data_DesignVar_Infill];
    Data_Geometry = [Data_Geometry_DOE, Data_Geometry_Infill];
end


%% 数据校验
if ~strcmp(pwd, Dir.Total) ||...
   ~exist(Dir_Total,'dir') || strcmp(Dir_Total(end),'\') || ~ischar(Dir_Total)
    error('*** Error：Dir_Total ***');
end

FileName_Catia_Error_Total = dir([Dir.Catia,'\','*.CATPart.Error']);
if ~isempty(FileName_Catia_Error_Total)
    error('*** Error：Errors in CATPart ***')
end

FileName_Catia_Total = dir([Dir.Catia,'\','*.CATPart']);
if mod(length(FileName_Catia_Total), Num_Node) ~= 0
    error('*** Error：Number of Catia is not a multiple of Num_Node ***');
end
if Mode == 2 || Mode == 3
    for nCatia = 1:length(FileName_Catia_Total)
        if ~contains(FileName_Catia_Total(nCatia).name, num2str(nCatia))
            error('*** Error：Order of Catia ***');
        end
    end
end
    
if Mode ~= 3 && nargin > 2
    error('*** Error：Mode & Input ***');
end


%% 主程序
FileName_Catia_Total = dir([Dir.Catia,'\','*.CATPart']);

if Mode == 1 || Mode == 2
    nStart_Catia = 1;
    nStart_Data = 1;
end
if Mode == 3
    nStart_Catia = length(FileName_Catia_Total) - Num_Infill + 1;
    nStart_Data = length(Data_Post_DOE) + 1;
end

for nRound = ((nStart_Catia-1)/Num_Node+1):length(FileName_Catia_Total)/Num_Node

    %提取Catia名称
    for nNode = 1:Num_Node
        nCatia = (nRound - 1) * Num_Node + nNode;
        FileName_Catia_PreRound{nNode} = FileName_Catia_Total(nCatia).name(1:end-8);
    end
    
    %飞行状态参数赋值
    clear FlightStateRef;
    for nNode = 1:Num_Node
        nCatia = (nRound - 1) * Num_Node + nNode;
        
        if Para_CFD.FlightState.Mode_HVA(1) == 2
            Para_CFD.FlightState.Ref_Height = Data_DesignVar(nCatia).FlightState.Height;
        end
        if Para_CFD.FlightState.Mode_HVA(2) == 2
            Para_CFD.FlightState.Ref_Velocity = Data_DesignVar(nCatia).FlightState.Velocity;
        end
        [FlightStateRef(nNode).Pressure,...
         FlightStateRef(nNode).Temperature,...
         FlightStateRef(nNode).Rho,...
         FlightStateRef(nNode).SoundSpeed,...
         FlightStateRef(nNode).DynamicVisc,...
         FlightStateRef(nNode).Velocity]...
         = F_FlightState(Para_CFD);
     
         FlightStateRef(nNode).Re = FlightStateRef(nNode).Rho * FlightStateRef(nNode).Velocity * Data_Geometry(nCatia).MAC_Length /...
                                    FlightStateRef(nNode).DynamicVisc;
    end
    
    if Mode == 2
        fprintf('剩余：%d\n',length(FileName_Catia_Total)-(nRound-1) * Num_Node);
    end
    fprintf('******************** ICEM *********************\n');
    for nNode = 1:Num_Node
        fprintf('Model_%d：%s\n', nNode, FileName_Catia_PreRound{nNode});
    end
    T0 = clock;

    fprintf('Run ICEM...');
    fprintf(' %dN..', Num_Node);
    
    % parfor nNode = 1:Num_Node
    for nNode = 1:Num_Node
        nCatia = (nRound - 1) * Num_Node + nNode;

        CFD_ICEM_Gen(Dir, Para_CFD, FlightStateRef(nNode), Data_Geometry(nCatia), FileName_Catia_PreRound{nNode});
    end
    
    for nNode = 1:Num_Node
        % parfeval(@CFD_ICEM_Post, 0, Dir, FileName_Catia_PreRound{nNode}, 1);
        % CFD_ICEM_Post(Dir, FileName_Catia_PreRound{nNode}, 1);
    end
    
    fprintf('\n');    

    fprintf('***********************************************\n');
    Timer = round(etime(clock,T0));
    fprintf('耗时：%d秒\n\n',Timer);
    
end


end


%%
function [Pressure, Temperature, Rho, SoundSpeed, DynamicVisc, Velocity] = F_FlightState(Para_CFD)

%压强(Pa),温度(K)
if Para_CFD.FlightState.Mode_PT == 1
    Pressure = Para_CFD.FlightState.Ref_Pressure;
    Temperature = Para_CFD.FlightState.Ref_Temperature;
elseif Para_CFD.FlightState.Mode_PT == 2 
    [Pressure, Temperature] = USAtmo(Para_CFD.FlightState.Ref_Height);
end

%密度(kg/m^3),声速(m/s),动力学粘度(kg/ms)
[Rho, SoundSpeed, DynamicVisc] = F_Atmosphere(Pressure, Temperature);

%速度(m/s)
if Para_CFD.FlightState.Ref_Velocity > 1e10
    Velocity = (Para_CFD.FlightState.Ref_Velocity - 1e10) * SoundSpeed;
else
    Velocity = Para_CFD.FlightState.Ref_Velocity;
end

end


%%
function [Rho, SoundSpeed, DynamicVisc] = F_Atmosphere(Pressure, Temperature)

Rho = Pressure / (287.0531 * Temperature);                                  %密度(kg/m3)

SoundSpeed = sqrt(1.4 * 287.0531 * Temperature);                            %声速(m/s)

DynamicVisc = 1.716e-05 * (Temperature/273.11)^(3/2) *...                   %动力学粘度,3参数方程(kg/ms)
              (273.11+110.56) / (Temperature+110.56);
 
end



















