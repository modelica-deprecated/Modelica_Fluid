package AST_BatchPlant 
  "Model of the experimental batch plant at Process Control Laboratory at University of Dortmund (Prof. Engell)" 
  model BatchPlant_StandardWater 
    parameter Real riseTime = 0.001;
    
    replaceable package BatchMedium = Modelica.Media.Water.StandardWater extends 
      Modelica.Media.Interfaces.PartialTwoPhaseMedium "Component media";
    
  /*
  replaceable package BatchMedium =Modelica.Media.Electrolytes.WaterNaCl extends 
    Modelica.Media.Interfaces.PartialTwoPhaseMedium "Component media";
*/
    
    BaseClasses.TankWith3InletOutletArraysWithEvaporatorCondensor B5(
      redeclare package Medium = BatchMedium,
      height=0.5,
      n_SidePorts=0,
      V0=0.001,
      n_BottomPorts=1,
      bottom_pipeArea={0.0001},
      top_pipeArea={0.0001},
      n_TopPorts=1,
      min_level_for_heating=0.0001,
      level_start=0.0009,
      area=0.05,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      side_pipeArea=fill(0, 0)) 
      annotation (extent=[-100,-60; -20,-20]);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V12(
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      m_flow_small=0) 
      annotation (extent=[-90,-8; -70,12],    rotation=90);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V15(
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      m_flow_small=0) 
      annotation (extent=[-70,-72; -90,-92],    rotation=90);
    Modelica.Thermal.HeatTransfer.PrescribedHeatFlow HeatB5 
      annotation (extent=[-134,-50; -114,-30], rotation=0);
    Modelica.Thermal.HeatTransfer.PrescribedHeatFlow CoolingB7 
      annotation (extent=[-140,-130; -120,-110],
                                               rotation=0);
    Modelica.Thermal.HeatTransfer.PrescribedHeatFlow CoolingB6 
      annotation (extent=[92,-70; 112,-50],    rotation=180);
    
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.Controller controller(
        Transition3(enableTimer=true, waitTime=30), Transition7(
        condition=true,
        enableTimer=true,
        waitTime=300))       annotation (extent=[80,60; 120,100]);
    
    annotation (Diagram, Coordsys(extent=[-300,-300; 300,300]),
      experiment(StopTime=3100),
      experimentSetupOutput,
      Commands(file=
            "../Scripts/Examples/AST_BatchPlant_StandardWater/plot level.mos" 
          "plot level"));
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V11(
      m_flow_small=0,
      Kv=0.01,
      redeclare package Medium = BatchMedium) 
      annotation (extent=[-60,78; -40,98],  rotation=0);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V8(
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      m_flow_small=0) 
      annotation (extent=[-90,152; -70,172],   rotation=90);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V9(
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      m_flow_small=0) 
      annotation (extent=[90,152; 70,172], rotation=90);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V2(
      m_flow_small=0,
      Kv=0.01,
      redeclare package Medium = BatchMedium) 
      annotation (extent=[-30,230; -50,250]);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V4(
      m_flow_small=0,
      Kv=0.01,
      redeclare package Medium = BatchMedium) 
      annotation (extent=[50,230; 30,250]);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V3(
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      m_flow_small=0,
      riseTime=riseTime,
      finiteRiseTime=false) 
      annotation (extent=[-114,210; -134,230]);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume2(
      redeclare package Medium = BatchMedium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V=0.001) annotation (extent=[-190,210; -170,230]);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V6(
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      m_flow_small=0) 
      annotation (extent=[112,210; 132,230]);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume8(
      redeclare package Medium = BatchMedium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V=0.001) annotation (extent=[150,210; 170,230]);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V23(
      m_flow_small=0,
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      riseTime=riseTime,
      finiteRiseTime=false) 
      annotation (extent=[-116,-240; -96,-260], rotation=180);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume3(
      redeclare package Medium = BatchMedium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V=0.001) annotation (extent=[-190,-260; -170,-240]);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V1(
      m_flow_small=0,
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      riseTime=riseTime,
      finiteRiseTime=false) 
      annotation (extent=[-170,100; -190,120],  rotation=90);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V22(
      m_flow_small=0,
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      riseTime=riseTime,
      finiteRiseTime=false) 
      annotation (extent=[-170,-66; -190,-46],  rotation=90);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume1(
      redeclare package Medium = BatchMedium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V=0.001) annotation (extent=[-190,60; -170,80]);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V5(
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      m_flow_small=0) 
      annotation (extent=[170,120; 150,100],    rotation=270);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V24(
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      m_flow_small=0,
      riseTime=riseTime,
      finiteRiseTime=false) 
      annotation (extent=[104,-240; 84,-260],   rotation=180);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V25(
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      m_flow_small=0,
      riseTime=riseTime,
      finiteRiseTime=false) 
      annotation (extent=[170,-10; 150,-30],    rotation=270);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume6(
      redeclare package Medium = BatchMedium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V=0.001) annotation (extent=[150,-258; 170,-238]);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume7(
      redeclare package Medium = BatchMedium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V=0.001) annotation (extent=[150,50; 170,70]);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V20(
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      m_flow_small=0,
      riseTime=riseTime,
      finiteRiseTime=false) 
      annotation (extent=[70,-210; 50,-190],    rotation=90);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V19(
      m_flow_small=0,
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      riseTime=riseTime,
      finiteRiseTime=false) 
      annotation (extent=[2,-210; -18,-190],    rotation=90);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume4(
      redeclare package Medium = BatchMedium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V=0.001) annotation (extent=[-38,-260; -18,-240]);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V10(
      m_flow_small=0,
      Kv=0.01,
      redeclare package Medium = BatchMedium) 
      annotation (extent=[30,60; 10,80],    rotation=90);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V21(
      m_flow_small=0,
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      riseTime=riseTime,
      finiteRiseTime=false) 
      annotation (extent=[42,-240; 22,-260],    rotation=180);
    Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume5(
      redeclare package Medium = BatchMedium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V=0.001) annotation (extent=[50,-260; 70,-240]);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ValveDiscrete V18(
      Kv=0.01,
      redeclare package Medium = BatchMedium,
      m_flow_small=0,
      riseTime=riseTime,
      finiteRiseTime=false) 
      annotation (extent=[-70,-242; -90,-222],  rotation=90);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.PumpWithAssertOfCavitation
      P1(
      redeclare package Medium = BatchMedium,
      M=0.01,
      pin_start=1e5,
      pout_start=1e5,
      m_flow_start=0.1,
      redeclare function flowCharacteristic = 
          Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticFlow (q_nom={0,
              0.001,0.0015}, head_nom={100,50,0}),
      redeclare package SatMedium = BatchMedium,
      N_nom=200,
      checkValve=false) 
      annotation (extent=[-128,-260; -148,-240]);
    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.PumpWithAssertOfCavitation
      P2(
      redeclare package Medium = BatchMedium,
      M=0.01,
      checkValve=false,
      pin_start=1e5,
      pout_start=1e5,
      m_flow_start=0.1,
      redeclare function flowCharacteristic = 
          Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticFlow(q_nom={0,
              0.001,0.0015}, head_nom={100,50,0}),
      N_nom=200,
      redeclare package SatMedium = BatchMedium) 
      annotation (extent=[112,-256; 132,-236]);
    Volumes.Tank B1(
      level_start=0.2,
      redeclare package Medium = BatchMedium,
      levelMax=0.5,
      area=0.05,
      V0=0.0001,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0)}) annotation (extent=[-100,180; -60,220]);
    inner Ambient ambient annotation (extent=[-172,250; -152,270]);
    Modelica.Blocks.Logical.TriggeredTrapezoid P1_on(amplitude=100, rising=0) 
      annotation (extent=[-158,-234; -138,-214]);
    Modelica.Blocks.Logical.TriggeredTrapezoid P2_on(amplitude=50, rising=0) 
      annotation (extent=[90,-228; 110,-208]);
    Volumes.Tank B2(
      level_start=0.2,
      redeclare package Medium = BatchMedium,
      levelMax=0.5,
      area=0.05,
      V0=0.0001,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0)}) annotation (extent=[60,180; 100,220]);
    Volumes.Tank B3(
      redeclare package Medium = BatchMedium,
      levelMax=0.5,
      area=0.05,
      V0=0.0001,
      nTopPorts=2,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0),Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=
           0.011, portLevel=0)},
      level_start=0.02)  annotation (extent=[-20,100; 20,140]);
    BaseClasses.CoolingTank B4(
      redeclare package Medium = BatchMedium,
      levelMax=0.5,
      area=0.05,
      V0=0.0001,
      level_start=0.015,
      nTopPorts=1,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0)}) annotation (extent=[-100,30; -60,70]);
    BaseClasses.CoolingTank B7(
      redeclare package Medium = BatchMedium,
      V0=0.0001,
      nTopPorts=1,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0)},
      level_start=0.009,
      T_start=298,
      alpha0=4.9,
      levelMax=0.5,
      area=0.05,
      stiffCharacteristicForEmptyPort=false) 
                         annotation (extent=[-100,-140; -60,-100]);
    PressureLosses.WallFrictionAndGravity PipeB1B2(
      redeclare package Medium = BatchMedium,
      length=1,
      diameter=0.1,
      height_ab=0) annotation (extent=[10,230; -10,250]);
    PressureLosses.WallFrictionAndGravity PipeB1B3(
      redeclare package Medium = BatchMedium,
      length=1,
      diameter=0.1,
      height_ab=0.1,
      from_dp=false,
      m_flow_small=1,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction) 
      annotation (extent=[-42,134; -62,154]);
    PressureLosses.WallFrictionAndGravity PipeB2B3(
      redeclare package Medium = BatchMedium,
      length=1,
      diameter=0.1,
      height_ab=0.1,
      from_dp=false,
      m_flow_small=1,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction) 
      annotation (extent=[36,134; 56,154]);
    PressureLosses.WallFrictionAndGravity PipeB1B1(
      redeclare package Medium = BatchMedium,
      length=1,
      diameter=0.1,
      from_dp=false,
      m_flow_small=1,
      height_ab=0.5,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction) 
      annotation (extent=[30,20; 10,40],   rotation=90);
    PressureLosses.WallFrictionAndGravity PipeB6Pump(
      redeclare package Medium = BatchMedium,
      m_flow_small=1,
      roughness=0,
      from_dp=true,
      length=0.5,
      diameter=0.1,
      height_ab=0.5,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed) 
                   annotation (extent=[70,-116; 50,-96],  rotation=90);
    PressureLosses.WallFrictionAndGravity PipeB7Pump(
      redeclare package Medium = BatchMedium,
      length=1,
      diameter=0.1,
      from_dp=false,
      m_flow_small=1,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction,
      height_ab=0.1,
      roughness=0) annotation (extent=[-70,-210; -90,-190],   rotation=90);
    PressureLosses.WallFrictionAndGravity PipePump1B1(
      redeclare package Medium = BatchMedium,
      length=1,
      diameter=0.1,
      from_dp=false,
      m_flow_small=1,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction,
      roughness=0,
      height_ab=3) annotation (extent=[-170,-14; -190,6],     rotation=90);
    PressureLosses.WallFrictionAndGravity PipePump2B2(
      redeclare package Medium = BatchMedium,
      length=1,
      diameter=0.1,
      from_dp=false,
      m_flow_small=1,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction,
      roughness=0,
      height_ab=3) annotation (extent=[170,0; 150,20],      rotation=90);
    BaseClasses.CoolingTank B6(
      redeclare package Medium = BatchMedium,
      V0=0.0001,
      nTopPorts=1,
      T_start=298,
      alpha0=4.9,
      levelMax=0.5,
      area=0.05,
      level_start=0.02,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0)},
      stiffCharacteristicForEmptyPort=false) 
                         annotation (extent=[80,-80; 40,-40]);
  equation 
    controller.sensors.LIS_301 = B3.level;
    controller.sensors.QI_302 = 0;//B3.medium.X[2];
    controller.sensors.LIS_501 = B5.level;
    controller.sensors.QIS_502 = 0;//B5.medium.X[2];
    controller.sensors.TI_503 = B5.medium.T;
    controller.sensors.LIS_601 = B6.level;
    controller.sensors.TIS_602 = B6.medium.T;
    controller.sensors.LIS_701 = B7.level;
    controller.sensors.TIS_702 = B7.medium.T;
    
    P1_on.u = controller.actuators.P1;
    P2_on.u = controller.actuators.P2;
    V1.open = controller.actuators.V1;
    V2.open = controller.actuators.V2;
    V3.open = controller.actuators.V3;
    V4.open = controller.actuators.V4;
    V5.open = controller.actuators.V5;
    V6.open = controller.actuators.V6;
    V8.open = controller.actuators.V8;
    V9.open = controller.actuators.V9;
    V10.open = controller.actuators.V10;
    V11.open = controller.actuators.V11;
    V12.open = controller.actuators.V12;
    V15.open = controller.actuators.V15;
    V18.open = controller.actuators.V18;
    V19.open = controller.actuators.V19;
    V20.open = controller.actuators.V20;
    V21.open = controller.actuators.V21;
    V22.open = controller.actuators.V22;
    V23.open = controller.actuators.V23;
    V24.open = controller.actuators.V24;
    V25.open = controller.actuators.V25;
    HeatB5.Q_flow = if controller.actuators.T5_Heater then 20000 else 0;
    CoolingB6.Q_flow = if controller.actuators.T6_Cooling then -2000 else 0;
    CoolingB7.Q_flow = if controller.actuators.T7_Cooling then -2000 else 0;
    
    connect(PortVolume2.port, V3.port_b) annotation (points=[-180,220; -134,220],
                style(color=69, rgbcolor={0,127,255}));
    connect(V2.port_b, PortVolume2.port) annotation (points=[-50,240; -180,240;
          -180,220],      style(color=69, rgbcolor={0,127,255}));
    connect(V6.port_b, PortVolume8.port) annotation (points=[132,220; 160,220],
        style(color=69, rgbcolor={0,127,255}));
    connect(V4.port_a, PortVolume8.port) annotation (points=[50,240; 160,240;
          160,220],                 style(color=69, rgbcolor={0,127,255}));
    connect(PortVolume1.port, V1.port_a) annotation (points=[-180,70; -180,100],
        style(color=69, rgbcolor={0,127,255}));
    connect(PortVolume2.port, V1.port_b) annotation (points=[-180,220; -180,120],
        style(color=69, rgbcolor={0,127,255}));
    connect(V22.port_a, PortVolume3.port) annotation (points=[-180,-66; -180,
          -250], style(color=69, rgbcolor={0,127,255}));
    connect(PortVolume5.port, V24.port_a) annotation (points=[60,-250; 84,-250],
        style(color=69, rgbcolor={0,127,255}));
    connect(V25.port_a, PortVolume6.port) annotation (points=[160,-30; 160,-248],
        style(color=69, rgbcolor={0,127,255}));
    connect(V5.port_a, PortVolume7.port) annotation (points=[160,100; 160,60],
        style(color=69, rgbcolor={0,127,255}));
    connect(V5.port_b, PortVolume8.port) annotation (points=[160,120; 160,220],
        style(color=69, rgbcolor={0,127,255}));
    connect(V19.port_a, PortVolume4.port) annotation (points=[-8,-210; -8,-250;
          -28,-250], style(color=69, rgbcolor={0,127,255}));
    connect(V23.port_a, PortVolume4.port) annotation (points=[-96,-250; -28,
          -250],     style(color=69, rgbcolor={0,127,255}));
    connect(V21.port_b, PortVolume5.port) annotation (points=[42,-250; 60,-250],
        style(color=69, rgbcolor={0,127,255}));
    connect(PortVolume4.port, V21.port_a) annotation (points=[-28,-250; 22,-250],
                 style(color=69, rgbcolor={0,127,255}));
    connect(V20.port_a, PortVolume5.port) annotation (points=[60,-210; 60,-250],
                 style(color=69, rgbcolor={0,127,255}));
    connect(V18.port_a, V23.port_a) annotation (points=[-80,-242; -80,-250; -96,
          -250],            style(color=69, rgbcolor={0,127,255}));
    connect(P1.outlet, PortVolume3.port) annotation (points=[-144,-246.8; -144,
          -250; -180,-250],      style(color=69, rgbcolor={0,127,255}));
    connect(P1.inlet, V23.port_b) annotation (points=[-130,-252; -114,-252;
          -114,-250; -116,-250],
        style(color=69, rgbcolor={0,127,255}));
    connect(V24.port_b, P2.inlet) annotation (points=[104,-250; 106,-250; 106,
          -248; 114,-248],
        style(color=69, rgbcolor={0,127,255}));
    connect(P2.outlet, PortVolume6.port) annotation (points=[128,-242.8; 128,
          -248; 160,-248],     style(color=69, rgbcolor={0,127,255}));
    connect(V15.port_a, B5.BottomFluidPort[1]) annotation (points=[-80,-72; -80,
          -60.4],          style(color=69, rgbcolor={0,127,255}));
    connect(V3.port_a, B1.topPorts[1]) annotation (points=[-114,220; -106,220;
          -106,230; -80,230; -80,220],   style(color=69, rgbcolor={0,127,255}));
    connect(B1.ports[1], V8.port_b) annotation (points=[-80,180; -80,172],
        style(color=69, rgbcolor={0,127,255}));
    connect(P1_on.y, P1.N_in) annotation (points=[-137,-224; -135.4,-224;
          -135.4,-245.6],
                   style(color=74, rgbcolor={0,0,127}));
    connect(P2_on.y, P2.N_in) annotation (points=[111,-218; 119.4,-218; 119.4,
          -241.6], style(color=74, rgbcolor={0,0,127}));
    connect(B2.topPorts[1], V6.port_a) annotation (points=[80,220; 80,228; 106,
          228; 106,220; 112,220],
                           style(color=69, rgbcolor={0,127,255}));
    connect(B2.ports[1], V9.port_b) annotation (points=[80,180; 80,172], style(
          color=69, rgbcolor={0,127,255}));
    connect(V11.port_b, B3.ports[1]) annotation (points=[-40,88; 0,88; 0,100],
                style(color=69, rgbcolor={0,127,255}));
    connect(V10.port_b, B3.ports[2]) annotation (points=[20,80; 20,88; 0,88; 0,
          100],          style(color=69, rgbcolor={0,127,255}));
    connect(B4.ports[1], V12.port_b) annotation (points=[-80,30; -80,12],   style(
          color=69, rgbcolor={0,127,255}));
    connect(CoolingB7.port, B7.heatPort) annotation (points=[-120,-120; -100,
          -120],
        style(color=42, rgbcolor={191,0,0}));
    connect(V2.port_a, PipeB1B2.port_b) annotation (points=[-30,240; -10,240],
        style(color=69, rgbcolor={0,127,255}));
    connect(PipeB1B2.port_a, V4.port_b) annotation (points=[10,240; 30,240],
        style(color=69, rgbcolor={0,127,255}));
    connect(PipeB1B3.port_b, V8.port_a) annotation (points=[-62,144; -80,144;
          -80,152],  style(color=69, rgbcolor={0,127,255}));
    connect(PipeB1B3.port_a, B3.topPorts[1]) annotation (points=[-42,144; 0,144;
          0,139],   style(color=69, rgbcolor={0,127,255}));
    connect(PipeB2B3.port_a, B3.topPorts[2]) annotation (points=[36,144; 0,144;
          0,141],   style(color=69, rgbcolor={0,127,255}));
    connect(PipeB2B3.port_b, V9.port_a) annotation (points=[56,144; 80,144; 80,
          152], style(color=69, rgbcolor={0,127,255}));
    connect(V11.port_a, B4.topPorts[1]) annotation (points=[-60,88; -80,88; -80,
          70],       style(color=69, rgbcolor={0,127,255}));
    connect(V10.port_a, PipeB1B1.port_b) annotation (points=[20,60; 20,40],
        style(color=69, rgbcolor={0,127,255}));
    connect(PipeB1B1.port_a, V21.port_a) annotation (points=[20,20; 20,-250; 22,
          -250],     style(color=69, rgbcolor={0,127,255}));
    connect(B5.TopFluidPort[1], V12.port_a) annotation (points=[-80,-19.6; -80,
          -8],style(color=3, rgbcolor={0,0,255}));
    connect(V15.port_b, B7.topPorts[1]) annotation (points=[-80,-92; -80,-100],
        style(color=69, rgbcolor={0,127,255}));
    connect(B7.ports[1], PipeB7Pump.port_b) annotation (points=[-80,-140; -80,
          -190], style(color=69, rgbcolor={0,127,255}));
    connect(PipeB7Pump.port_a, V18.port_b) annotation (points=[-80,-210; -80,
          -222], style(color=69, rgbcolor={0,127,255}));
    connect(PipePump1B1.port_a, V22.port_b) annotation (points=[-180,-14; -180,
          -46],  style(color=69, rgbcolor={0,127,255}));
    connect(PipePump1B1.port_b, PortVolume1.port) annotation (points=[-180,6;
          -180,70],  style(color=69, rgbcolor={0,127,255}));
    connect(PipePump2B2.port_b, PortVolume7.port) annotation (points=[160,20;
          160,60],  style(color=69, rgbcolor={0,127,255}));
    connect(V25.port_b, PipePump2B2.port_a) annotation (points=[160,-10; 160,0],
                           style(color=69, rgbcolor={0,127,255}));
    connect(B6.ports[1], PipeB6Pump.port_b) annotation (points=[60,-80; 60,-96],
        style(color=69, rgbcolor={0,127,255}));
    connect(B6.topPorts[1], B5.Condensed) annotation (points=[60,-40; 60,-28;
          -19.6,-28],            style(color=69, rgbcolor={0,127,255}));
    connect(CoolingB6.port, B6.heatPort) annotation (points=[92,-60; 80,-60],
                        style(color=42, rgbcolor={191,0,0}));
    connect(V19.port_b, PipeB6Pump.port_a) annotation (points=[-8,-190; -8,-140;
          60,-140; 60,-116],       style(color=69, rgbcolor={0,127,255}));
    connect(V20.port_b, PipeB6Pump.port_a) annotation (points=[60,-190; 60,-116],
        style(color=69, rgbcolor={0,127,255}));
    connect(HeatB5.port, B5.HeatPort) annotation (points=[-114,-40; -102,-40],
        style(
        color=42,
        rgbcolor={191,0,0},
        fillColor=30,
        rgbfillColor={215,215,215},
        fillPattern=10));
  end BatchPlant_StandardWater;
  
  package BaseClasses 
    extends Modelica_Fluid.Icons.BaseClassLibrary;
    block TriggeredTrapezoid "Triggered trapezoid generator" 
      extends Modelica.Blocks.Interfaces.partialBooleanBlockIcon;
      
      parameter Real amplitude=1 "Amplitude of trapezoid";
      parameter Modelica.SIunits.Time rising(final min=0)=0 
        "Rising duration of trapezoid";
      parameter Modelica.SIunits.Time falling(final min=0)=rising 
        "Falling duration of trapezoid";
      parameter Real offset=0 "Offset of output signal";
      
      Modelica.Blocks.Interfaces.BooleanInput u 
        "Connector of Boolean input signal" 
                                       annotation(extent=[-140,-20; -100,20]);
      Modelica.Blocks.Interfaces.RealOutput y "Connector of Real output signal"
        annotation (extent=[100,-10; 120,10]);
      
      annotation (
        Icon(
          Line(points=[-60, -70; -60, -70; -30, 40; 8, 40; 40, -70; 40, -70]),
          Line(points=[-90, -70; 82, -70], style(color=8)),
          Line(points=[-80, 68; -80, -80], style(color=8)),
          Polygon(points=[90, -70; 68, -62; 68, -78; 90, -70], style(
              color=8,
              fillColor=8,
              fillPattern=1)),
          Polygon(points=[-80, 90; -88, 68; -72, 68; -80, 90], style(color=8,
                fillColor=8)),
          Line(points=[-80, -70; -60, -70; -60, 24; 8, 24; 8, -70; 60, -70],
              style(color=5))),
        Diagram(
          Line(points=[-80, -20; -60, -20; -30, 40; 8, 40; 40, -20; 60, -20]),
          Line(points=[-90, -70; 82, -70], style(color=0)),
          Line(points=[-80, 68; -80, -80], style(color=0)),
          Polygon(points=[90, -70; 68, -62; 68, -78; 90, -70], style(color=0,
                fillColor=7)),
          Polygon(points=[-80, 90; -88, 68; -72, 68; -80, 90], style(color=0,
                fillColor=7)),
          Line(points=[-80, -68; -60, -68; -60, -42; 8, -42; 8, -68; 60, -68],
              style(color=5)),
          Line(points=[-60, 40; -60, -42], style(
              color=0,
              pattern=3,
              fillColor=7,
              fillPattern=1)),
          Line(points=[8, -42; 8, 40], style(
              color=0,
              pattern=3,
              fillColor=7,
              fillPattern=1)),
          Line(points=[-58,40; -28,40],   style(
              color=0,
              fillColor=0,
              fillPattern=1)),
          Line(points=[8, -20; 40, -20], style(
              color=0,
              fillColor=0,
              fillPattern=1)),
          Line(points=[-20, 40; -20, -20], style(
              color=0,
              fillColor=0,
              fillPattern=1)),
          Line(points=[-20, -20; -20, -70], style(
              color=0,
              fillColor=0,
              fillPattern=1)),
          Text(
            extent=[-42,48; -42,38],
            style(
              color=0,
              arrow=3,
              fillColor=0,
              fillPattern=1),
            string="rising"),
          Text(
            extent=[24, -10; 24, -20],
            style(
              color=0,
              arrow=3,
              fillColor=0,
              fillPattern=1),
            string="falling"),
          Polygon(points=[-58,40; -54,42; -54,38; -58,40],     style(color=0,
                fillColor=7)),
          Polygon(points=[-30, 40; -34, 42; -34, 38; -30, 40], style(color=0,
                fillColor=7)),
          Polygon(points=[8, -20; 12, -18; 12, -22; 8, -20], style(color=0,
                fillColor=7)),
          Polygon(points=[40, -20; 36, -18; 36, -22; 40, -20], style(color=0,
                fillColor=7)),
          Polygon(points=[-22, -24; -20, -20; -18, -24; -22, -24], style(color=0,
                  fillColor=7)),
          Polygon(points=[-18, -66; -22, -66; -20, -70; -18, -66], style(color=0,
                  fillColor=7)),
          Polygon(points=[-22, 36; -20, 40; -18, 36; -22, 36], style(color=0,
                fillColor=7)),
          Polygon(points=[-18, -16; -22, -16; -20, -20; -18, -16], style(color=0,
                  fillColor=7)),
          Rectangle(extent=[-40, 6; 0, -4], style(
              color=7,
              fillColor=7,
              fillPattern=1)),
          Text(
            extent=[-20, 6; -20, -4],
            style(
              color=0,
              arrow=3,
              fillColor=0,
              fillPattern=1),
            string="amplitude"),
          Rectangle(extent=[-40, -48; 0, -58], style(
              color=7,
              fillColor=7,
              fillPattern=1)),
          Text(
            extent=[-20, -48; -20, -58],
            style(
              color=0,
              arrow=3,
              fillColor=0,
              fillPattern=1),
            string="offset"),
          Text(
            extent=[60,-82; 94,-92],
            string="time",
            style(color=0, rgbcolor={0,0,0})),
          Text(
            extent=[-88,-4; -54,-14],
            style(color=0, rgbcolor={0,0,0}),
            string="y"),
          Text(
            extent=[-88,-46; -54,-56],
            style(color=0, rgbcolor={0,0,0}),
            string="u"),
          Polygon(points=[40,60; 36,62; 36,58; 40,60],         style(color=0,
                fillColor=7)),
          Line(points=[-56,60; 36,60],    style(
              color=0,
              fillColor=0,
              fillPattern=1)),
          Polygon(points=[-60,60; -56,62; -56,58; -60,60],     style(color=0,
                fillColor=7)),
          Text(
            extent=[-16,72; -16,62],
            style(
              color=0,
              arrow=3,
              fillColor=0,
              fillPattern=1),
            string="y_high")),
        Documentation(info="<HTML>
<p>The block TriggeredTrapezoid has a boolean input and a real
output signal and requires the parameters <i>amplitude</i>,
<i>rising</i>, <i>falling</i> and <i>offset</i>. The
output signal <b>y</b> represents a trapezoidal signal dependent on the
input signal <b>u</b>.
</p>
<p>The behaviour is as follows: Assume the initial input to be false. In this
case, the output will be <i>offset</i>. After a rising edge (i.e. the input
changes from false to true), the output is rising during <i>rising</i> to the
sum of <i>offset</i> and <i>amplitude</i>. In contrast, after a falling
edge (i.e. the input changes from true to false), the output is falling
during <i>falling</i> to a value of <i>offset</i>.
</p>
<p>Note, that the case of edges before expiration of rising or falling is
handled properly.</p>
</HTML>
"));
    protected 
      discrete Real endValue "Value of y at time of recent edge";
      discrete Real rate "Current rising/falling rate";
      discrete Modelica.SIunits.Time T 
        "Predicted time of output reaching endValue";
    public 
      Modelica.Blocks.Interfaces.BooleanOutput y_high 
        annotation (extent=[100,-90; 120,-70]);
    initial equation 
      /* A start value of y is set, because pre(y) is present
     to avoid a warning message from the compiler. However,
     this setting does not have an effect, because y is initialized
     correctly, before pre(y) is used
  */
      pre(y) = 0;
    equation 
        y_high = time < T;
        y = if y_high then endValue - (T - time)*rate else  endValue;
      
        when {initial(),u,not u} then
          endValue = if u then offset + amplitude else offset;
          rate = if u and (rising > 0) then amplitude/rising else 
            if not u and (falling > 0) then -amplitude/falling else 0;
          T = if u and not (rising > 0) or not u and not (falling
             > 0) or not abs(amplitude) > 0 or initial() then time else time
             + (endValue - pre(y))/rate;
        end when;
    end TriggeredTrapezoid;
    
    block setReal "Set output signal to a time varying Real expression" 
      
      Modelica.Blocks.Interfaces.RealInput u "Set value of Real input" 
        annotation (extent=[-140,-20; -100,20], Dialog(group=
              "Time varying input signal"));
      
      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Window(
          x=0.29,
          y=0.23,
          width=0.6,
          height=0.6),
        Icon(
          Rectangle(extent=[-100,40; 100,-40],   style(
              color=0,
              fillColor=30,
              fillPattern=11)),
          Text(
            extent=[-96,15; 96,-15],
            string="%u",
            style(
              color=0,
              fillColor=2,
              fillPattern=1)), Text(extent=[-150,90; 140,50],     string="%name")),
        Diagram,
        Documentation(info="<html>
 
</html>"));
      
    end setReal;
    
    model ValveDiscrete "Valve for water/steam flows with linear pressure drop" 
      extends Modelica_Fluid.Interfaces.PartialTwoPortTransport;
      parameter Modelica_Fluid.Types.HydraulicConductance Kv 
        "Hydraulic conductance at full opening";
      Modelica.Blocks.Interfaces.BooleanInput open annotation (extent=[10,94;
            -10,74], rotation=90);
      parameter Real m_flow_small = 1e-7 "massflow when valve is closed";
      parameter Boolean finiteRiseTime = false 
        "= true, if valve is opened/closed linearly in riseTime time";
      parameter Real riseTime = 0 "Time to open or close the valve" 
                                           annotation(Dialog(enable=finiteRiseTime));
    protected 
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.TriggeredTrapezoid 
        trapezoid(rising=riseTime) if finiteRiseTime 
        annotation (extent=[-10,40; 10,60], rotation=-90);
      Modelica.Blocks.Sources.RealExpression m_flow_trapezoid(y=if trapezoid.y_high then 
                  trapezoid.y*Kv*dp else Kv*dp*m_flow_small) if finiteRiseTime 
        annotation (extent=[-90,-40; 40,-20]);
      Modelica.Blocks.Sources.RealExpression m_flow_pulse(y=if open then Kv*dp else 
                  Kv*dp*m_flow_small) if not finiteRiseTime 
        annotation (extent=[-90,-60; 40,-40]);
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.setReal set(u=m_flow) 
        annotation (extent=[70,-50; 90,-30]);
    annotation (
      Icon(
          Line(points=[0,50; 0,0],   style(
              color=0,
              rgbcolor={0,0,0},
              fillPattern=1)),
          Rectangle(extent=[-20,60; 20,50],   style(
              color=0,
              fillColor=0,
              fillPattern=1)),
             Text(extent=[-145,-58; 146,-98],   string="%name"),
          Polygon(points=[-100,50; 100,-50; 100,50; 0,0; -100,-50; -100,50], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=DynamicSelect(7, if open > 0.5 then 2 else 7)))),
      Diagram,
      Documentation(info="<HTML>
<p>This very simple model provides a pressure drop which is proportional to the flowrate and to the <tt>opening</tt> signal, without computing any fluid property.
<p>A medium model must be nevertheless be specified, so that the fluid ports can be connected to other components using the same medium model.
</HTML>",
        revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"),
        Coordsys(grid=[1,1], scale=0));
    equation 
      connect(trapezoid.u, open) annotation (points=[-7.34764e-016,62; 0,62; 0,
            84],
          style(color=5, rgbcolor={255,0,255}));
      connect(m_flow_trapezoid.y, set.u) annotation (points=[46.5,-30; 60,-30; 60,
            -40; 68,-40], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_pulse.y, set.u) annotation (points=[46.5,-50; 60,-50; 60,-40;
            68,-40], style(color=74, rgbcolor={0,0,127}));
    end ValveDiscrete;
    
    model TankWith3InletOutletArraysWithEvaporatorCondensor 
      "Tank with Heating and Evaporation" 
      import Modelica.SIunits.Conversions.*;
      import Modelica_Fluid.Types.Init.*;
      replaceable package Medium = 
          Modelica.Media.Interfaces.PartialTwoPhaseMedium 
        extends Modelica.Media.Interfaces.PartialTwoPhaseMedium 
        "Medium in the component" 
        annotation (choicesAllMatching=true);
    // parameter for Tank
      parameter Modelica.SIunits.Area area "Tank area";
      parameter SI.Area top_pipeArea[n_TopPorts] "Area of outlet pipe";
      parameter SI.Area side_pipeArea[n_SidePorts] "Area of outlet pipe";
      parameter SI.Area bottom_pipeArea[n_BottomPorts] "Area of outlet pipe";
      parameter Modelica.SIunits.Height height(min=0) = 10 "Height of Tank";
      parameter SI.Volume V0=0 "Volume of the liquid when the level is zero";
      constant Modelica.SIunits.Acceleration g=Modelica.Constants.g_n;
      parameter Real side_heights[n_SidePorts]=zeros(n_SidePorts);
      parameter Real bottom_heights[n_BottomPorts]=zeros(n_BottomPorts);
      parameter Real top_heights[n_TopPorts]=fill(height, n_TopPorts);
      parameter SI.Height level_start(min=0) "Initial tank level" 
        annotation(Dialog(tab="Initialization"));
      parameter Modelica_Fluid.Types.InitTypes.Temp initType=NoInit 
        "Initialization option" 
        annotation(Dialog(tab = "Initialization"));
      parameter Boolean use_T_start=true 
        "Use T_start if true, otherwise h_start"                                    annotation(Dialog(tab = "Initialization"), Evaluate = true);
      parameter Medium.Temperature T_start=if use_T_start then 293.15 else 
          Medium.temperature_phX(p_ambient, h_start, X_start) 
        "Start value of temperature" 
        annotation(Dialog(tab = "Initialization", enable = use_T_start));
      parameter Medium.SpecificEnthalpy h_start=if use_T_start then Medium.specificEnthalpy_pTX(
          p_ambient, T_start, X_start[1:Medium.nXi]) else 1e4 
        "Start value of specific enthalpy" 
        annotation(Dialog(tab = "Initialization", enable = not use_T_start));
      parameter Medium.MassFraction X_start[Medium.nX]=Medium.reference_X 
        "Start value of mass fractions m_i/m" 
        annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
      parameter Medium.AbsolutePressure p_ambient=101325 
        "Tank surface pressure";
      parameter Medium.Temperature T_ambient=293.15 "Tank surface Temperature";
      parameter Integer n_TopPorts=1 "number of Top connectors";
      parameter Integer n_SidePorts=1 "number of side connectors";
      parameter Integer n_BottomPorts=1 "number of bootom connectors";
      Medium.BaseProperties medium(
        preferredMediumStates=true,
        p(start=p_ambient),
        T(start=T_start),
        Xi(start=X_start[1:Medium.nXi]));
      Modelica.SIunits.Height level(
        stateSelect=StateSelect.prefer,
        min=0,
        max=height) "Level height of tank";
      SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
      SI.Energy U "Internal energy of tank volume";
      Real m(quantity=Medium.mediumName, unit="kg") "Mass of tank volume";
      Real mXi[Medium.nXi](quantity=Medium.substanceNames, each unit="kg") 
        "Component masses of the independent substances";
    // Hilfsvariablen  
      Real H_flow_BottomPorts[n_BottomPorts];
      Real H_flow_SidePorts[n_SidePorts];
      Real H_flow_TopPorts[n_TopPorts];
      Real m_flow_BottomPorts[n_BottomPorts];
      Real m_flow_SidePorts[n_SidePorts];
      Real m_flow_TopPorts[n_TopPorts];
      
      Real m_flow_BottomPorts_pos[n_BottomPorts];
      Real m_flow_SidePorts_pos[n_SidePorts];
      Real m_flow_TopPorts_pos[n_TopPorts];
      Real m_flow_pos;
      Medium.MassFlowRate mXi_flow_topPorts[n_TopPorts,Medium.nXi];
      Medium.MassFlowRate mXi_flow_bottomPorts[n_BottomPorts,Medium.nXi];
      Medium.MassFlowRate mXi_flow_sidePorts[n_SidePorts,Medium.nXi];
      
    // Connectors and InnerTanks
      Modelica_Fluid.Interfaces.FluidPort_b BottomFluidPort[n_BottomPorts](
        redeclare package Medium = Medium,
        m_flow(each start=0),
        mXi_flow(each start=0)) 
        annotation (extent=[-110,-112; -90,-92],  rotation=90);
      Modelica_Fluid.Interfaces.FluidPort_a TopFluidPort[n_TopPorts](
        redeclare package Medium = Medium,
        m_flow(each start=0),
        mXi_flow(each start=0)) 
        annotation (extent=[-110,92; -90,112]);
      Modelica_Fluid.Interfaces.FluidPort_b SideFluidPort[n_SidePorts](
        redeclare package Medium = Medium,
        m_flow(each start=0),
        mXi_flow(each start=0)) 
        annotation (extent=[0,-10; 20,10]);
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.InnerTank InnerTankTop[n_TopPorts](
        each h=medium.h,
        each p_ambient=p_ambient,
        each d=medium.d,
        each Xi = medium.Xi,
        H_flow=H_flow_TopPorts,
        m_flow=m_flow_TopPorts,
        mXi_flow=mXi_flow_topPorts,
        aboveLevel={level - top_heights[i] for i in 1:n_TopPorts},
        pipeArea={top_pipeArea[i] for i in 1:n_TopPorts},
        redeclare package Medium = Medium) 
          annotation (extent=[-140,60; -120,80]);
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.InnerTank 
        InnerTankSide[                                                           n_SidePorts](
        each h=medium.h,
        each p_ambient=p_ambient,
        each d=medium.d,
        each Xi = medium.Xi,
        H_flow=H_flow_SidePorts,
        m_flow=m_flow_SidePorts,
        mXi_flow=mXi_flow_sidePorts,
        aboveLevel={level - side_heights[i] for i in 1:n_SidePorts},
        pipeArea={side_pipeArea[i] for i in 1:n_SidePorts},
        redeclare package Medium = Medium) 
          annotation (extent=[-20,0; 0,20]);
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.InnerTank 
        InnerTankBottom[                                                           n_BottomPorts](
        each h=medium.h,
        each p_ambient=p_ambient,
        each d=medium.d,
        each Xi = medium.Xi,
        H_flow=H_flow_BottomPorts,
        m_flow=m_flow_BottomPorts,
        mXi_flow=mXi_flow_bottomPorts,
        aboveLevel={level - bottom_heights[i] for i in 1:n_BottomPorts},
        pipeArea={bottom_pipeArea[i] for i in 1:n_BottomPorts},
        redeclare package Medium = Medium) 
          annotation (extent=[-80,-80; -60,-60]);
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a HeatPort 
        annotation (extent=[-220,-10; -200,10]);
      Modelica_Fluid.Interfaces.FluidPort_b Condensed(redeclare package Medium 
          =        Medium) 
        annotation (extent=[192,50; 212,70]);
    // parameter for Evaporator
      parameter Real min_level_for_heating;
      parameter Real k=4.9 "Wärmedurchgangskoeffizient";
      Real Q_lost "Wärmeverlust";
      
      Medium.SaturationProperties sat 
        "State vector to compute saturation properties";
      Medium.SpecificEnthalpy h_v=Medium.dewEnthalpy(sat) 
        "specific enthalpy of vapour";
      Medium.SpecificEnthalpy h_l=Medium.bubbleEnthalpy(sat) 
        "specific enthalpy of liquid";
      Medium.SpecificEnthalpy h "'is'specific enthalpy of liquid";
      Medium.Density rho_v=Medium.dewDensity(sat) "density in vapour phase";
      Medium.Density rho_l=Medium.bubbleDensity(sat) "density in liquid phase";
      Medium.Density rho "'is' density in liquid phase";
      
    equation 
      for i in 1:n_BottomPorts loop
        m_flow_BottomPorts_pos[i] = (if m_flow_BottomPorts[i] > 0 then 
          m_flow_BottomPorts[i] else 0);
      end for;
      for i in 1:n_SidePorts loop
        m_flow_SidePorts_pos[i] = if m_flow_SidePorts[i] > 0 then m_flow_SidePorts[
          i] else 0;
      end for;
      for i in 1:n_TopPorts loop
        m_flow_TopPorts_pos[i] = if m_flow_TopPorts[i] > 0 then m_flow_TopPorts[i] else 
                0;
      end for;
      for i in 1:n_BottomPorts loop
        connect(InnerTankBottom[i].port, BottomFluidPort[i]) annotation (points=[-70,-81;
              -70,-102; -100,-102],
                                style(color=3, rgbcolor={0,0,255}));
      end for;
      for i in 1:n_TopPorts loop
        connect(InnerTankTop[i].port, TopFluidPort[i])  annotation (points=[-130,59;
              -92,59; -92,102; -100,102],
                                        style(color=3, rgbcolor={0,0,255}));
      end for;
      for i in 1:n_SidePorts loop
        connect(InnerTankSide[i].port, SideFluidPort[i])  annotation (points=[-10,-1;
              24,-1; 24,0; 10,0],     style(color=3, rgbcolor={0,0,255}));
      end for;
      
      medium.p = p_ambient;
      medium.T = HeatPort.T;
    // Mass balance  
      der(m) = sum(m_flow_BottomPorts) + sum(m_flow_SidePorts) + sum(
        m_flow_TopPorts) + Condensed.m_flow;
    // Energy balance
      
      U = m*medium.h - p_ambient*V "Internal energy of fluid";
      Q_lost = -k*(2*area + 2*sqrt(Modelica.Constants.pi*area))*level*(medium.T -
        T_ambient);
      m = V*medium.d "Mass of fluid";
      V = area*level + V0 "Volume of fluid";
      mXi = m*medium.Xi "Mass of fluid components";
      sat.psat = medium.p;
      sat.Tsat = Medium.saturationTemperature(medium.p);
      
      if noEvent(medium.T < sat.Tsat) then
        if Medium.singleState then
          der(U) = sum(H_flow_BottomPorts) + sum(H_flow_SidePorts) + sum(
            H_flow_TopPorts) + Condensed.H_flow + HeatPort.Q_flow + Q_lost 
            "Mechanical work is neglected";
        else
          der(U) = sum(H_flow_BottomPorts) + sum(H_flow_SidePorts) + sum(
            H_flow_TopPorts) + Condensed.H_flow - p_ambient*der(V) + Q_lost +
            HeatPort.Q_flow;
        end if;
        Condensed.H_flow = 0;
        Condensed.m_flow = 0;
        rho = medium.d;
        h = medium.h;
      else
        if Medium.singleState then
                                 //Q_flow cooling = - (HeatPort.Q_flow-Q_lost)
          der(U) = sum(H_flow_BottomPorts) + sum(H_flow_SidePorts) + sum(H_flow_TopPorts)
             + Condensed.H_flow "Mechanical work is neglected";
        else
          der(U) = sum(H_flow_BottomPorts) + sum(H_flow_SidePorts) + sum(
            H_flow_TopPorts) + Condensed.H_flow - p_ambient*der(V);
        end if;
        Condensed.H_flow = Condensed.m_flow*medium.h;
        Condensed.m_flow = -(HeatPort.Q_flow - Q_lost)/(h_v - h_l);
        rho = rho_l;//Density = liquid Densety
        h = h_l;    //Enthalpy = liquid Enthalpy
        if noEvent(HeatPort.Q_flow > 0.0) then
          assert(noEvent(abs(m_flow_pos) <= 0.01), "Es wird beim verdampfen befüllt.");
        end if;
        
      end if;
      
      m_flow_pos = sum(m_flow_TopPorts_pos) + sum(m_flow_SidePorts_pos) + sum(
        m_flow_BottomPorts_pos);
      
      for i in 1:Medium.nXi loop
           der(mXi[i]) = sum(mXi_flow_bottomPorts[:,i]) +
                         sum(mXi_flow_sidePorts[:,i]) +
                         sum(mXi_flow_topPorts[:,i]);
      end for;
      
      assert(level < height, " 
    Tank ist überfüllt.
    ");
      
      assert(not (HeatPort.Q_flow > 0.0 and level <= min_level_for_heating), "
    Es wird leerer Tank bezeizt
  ");
      
    initial equation 
      if initType == NoInit then
        // no initial equations
      elseif initType == InitialValues then
        level = level_start;
        if use_T_start then
          medium.T = T_start;
        else
          medium.h = h_start;
        end if;
        medium.Xi = X_start[1:Medium.nXi];
      elseif initType == SteadyStateHydraulic then
        der(level) = 0;
        if use_T_start then
          medium.T = T_start;
        else
          medium.h = h_start;
        end if;
        medium.Xi = X_start[1:Medium.nXi];
      else
        assert(false, "Unsupported initialization option");
      end if;
      annotation (
        Icon(
          Rectangle(extent=[-200,100; 0,-90],
                                            style(color=7, fillColor=7)),
              Rectangle(extent=DynamicSelect([-200,-100; 0,0], [-200,-100; 0,(-100
                     + 200*level/levelMax)]), style(
                  color=69,
                  rgbcolor={0,127,255},
                  fillColor=71,
                  rgbfillColor={85,170,255},
                  fillPattern=1)),
          Line(points=[-200,100; -200,-100; 0,-100; 0,100], style(
              color=0,
              fillColor=69,
              fillPattern=1)),
          Text(
            extent=[-198,74; 0,38],
            string="%name",
            style(fillColor=69, fillPattern=1)),
          Text(
            extent=[-184,-64; -14,-86],
            style(color=0),
            string="%level_start"),
          Text(
            extent=[-192,-34; -12,-54],
            style(color=0),
            string="level_start ="),
            Line(points=[-200,100; 0,100],   style(
                color=0,
                rgbcolor={0,0,0},
                pattern=3)),
          Polygon(points=[0,100; 200,70; 200,50; 200,50; 0,80; 0,100], style(
              color=0,
              rgbcolor={0,0,0},
              gradient=2,
              fillColor=3,
              rgbfillColor={0,0,255})),
          Polygon(points=[20,98; 30,74; 52,84; 66,72; 86,78; 98,66; 118,74; 130,
                60; 144,70; 152,60; 168,66; 180,54; 196,74; 190,76; 180,64; 170,
                70; 156,66; 148,76; 132,68; 120,80; 100,74; 88,88; 70,78; 50,92;
                32,82; 28,100; 20,98; 20,98], style(
              color=0,
              rgbcolor={0,0,0},
              gradient=2,
              fillColor=67,
              rgbfillColor={170,255,255}))),
        Documentation(info="<HTML>
<p>
<p>This tank has the same geometric variables as TankWith3InletOutletArrays plus the feature of a HeatPort and the possibility of evaporation. 
(Assumption: The gas is condensed emidiatly afterwards so that a liquid boiling fluid is created.)
<p>The tank can be initialized with the following options:
<ul>
<li>NoInit: no explicit initial conditions
<li>InitialValues: initial values of temperature (or specific enthalpy), composition and level are specified
<li>SteadyStateHydraulic: initial values of temperature (or specific enthalpy) and composition are specified; the initial level is determined so that levels and pressure are at steady state.
</ul>
Full steady state initialization is not supported, because the corresponding intial equations for temperature/enthalpy are undetermined (the flow rate through the port at steady state is zero). 
</p>
</HTML>"),
        Diagram,
        Coordsys(extent=[-200,-100; 200,100]));
    equation 
      
    end TankWith3InletOutletArraysWithEvaporatorCondensor;
    
    model InnerTank 
      import Modelica_Fluid;
        replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching=true);
      
        Modelica_Fluid.Interfaces.FluidPort_a port(redeclare package Medium = 
            Medium) 
        annotation (extent=[-10, -120; 10, -100], rotation=90);
       // Real mXi_flow;
        Boolean m_flow_negative( start = true) "true= massflow out of tank";
        constant Modelica.SIunits.Acceleration g=Modelica.Constants.g_n;
       input Real aboveLevel;
        input Real d;
        input Real p_ambient;
        input Real h;
      input Medium.MassFraction Xi[Medium.nXi] 
        "Actual mass fractions of fluid in tank"                    annotation(Dialog);
        input Real pipeArea;
        output Real H_flow;
        output Real m_flow;
       output Medium.MassFlowRate mXi_flow[Medium.nXi] 
        "= port.mXi_flow (used to transform vector of connectors in vector of Real numbers)";
      
    equation 
    m_flow_negative = (pre(m_flow_negative) and not port.p>p_ambient) or (port.m_flow < -1e-6);
      
    if noEvent(aboveLevel > 0) then
      port.p = aboveLevel*g*d + p_ambient - smooth(2,noEvent(if noEvent(m_flow < 0) then m_flow^2/(2*d*pipeArea^2) else 0));
    else
     if pre(m_flow_negative) then
        port.m_flow = 0;
      else
        port.p = p_ambient;
      end if;
    end if;
      
      H_flow = port.H_flow;
      m_flow = port.m_flow;
      mXi_flow = port.mXi_flow;
      port.H_flow = semiLinear(port.m_flow, port.h, h);
        port.mXi_flow = semiLinear(port.m_flow, port.Xi, Xi);
      
    end InnerTank;
    
    model PumpWithAssertOfCavitation 
      "Centrifugal pump with ideally controlled speed" 
      extends Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.PartialPump(final 
          computeNPSHa =                                                                        true);
      import Modelica.SIunits.Conversions.NonSIunits.*;
      parameter AngularVelocity_rpm N_const = N_nom "Constant rotational speed";
      Modelica.Blocks.Interfaces.RealInput N_in "Prescribed rotational speed" 
        annotation (extent=[-36,34; -16,54],   rotation=-90);
    equation 
      
       assert(inlet.p >= pv,   " 
    wahrscheinlich ist ein Ventil zu oder ein Tank vor der Pumpe leer.
    ");
      
        N = N_in "Rotational speed";
      if cardinality(N_in)==0 then
        N_in = N_const "Rotational speed provided by parameter";
      end if;
      annotation (
        Icon(
          Text(extent=[-58,58; -30,38], string="n")),
        Diagram,
        Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with controlled speed, either fixed or provided by an external signal.
<p>The model extends <tt>PartialPump</tt>
<p>If the <tt>N_in</tt> input connector is wired, it provides rotational speed of the pumps (rpm); otherwise, a constant rotational speed equal to <tt>n_const</tt> (which can be different from <tt>N_nom</tt>) is assumed.</p>
</HTML>", revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
    end PumpWithAssertOfCavitation;
    
    partial model PartialPump "Base model for centrifugal pumps" 
      import Modelica.SIunits.Conversions.NonSIunits.*;
      import Modelica.Constants.*;
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model" annotation(choicesAllMatching=true);
      Medium.BaseProperties fluid(p(start=pin_start),h(start=h_start)) 
        "Fluid properties at the inlet";
      replaceable package SatMedium = 
          Modelica.Media.Interfaces.PartialTwoPhaseMedium 
        "Saturated medium model (required only for NPSH computation)" 
                                                                     annotation(choicesAllMatching=true);
      replaceable function flowCharacteristic = 
          Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.baseFlow 
        "Head vs. q_flow characteristic at nominal speed and density" 
        annotation(Dialog(group="Characteristics"), choicesAllMatching=true);
      parameter Boolean usePowerCharacteristic = false 
        "Use powerCharacteristic (vs. efficiencyCharacteristic)" 
         annotation(Dialog(group="Characteristics"));
        replaceable function powerCharacteristic = 
            Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticPower
          (q_nom={0,0,0},W_nom={0,0,0}) 
        "Power consumption vs. q_flow at nominal speed and density" 
          annotation(Dialog(group="Characteristics", enable = usePowerCharacteristic),
                     choicesAllMatching=true);
      replaceable function efficiencyCharacteristic = 
        Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.constantEfficiency
          ( eta_nom=0.8) extends 
        Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.baseEfficiency 
        "Efficiency vs. q_flow at nominal speed and density" 
        annotation(Dialog(group="Characteristics",enable = not usePowerCharacteristic),
                   choicesAllMatching=true);
      parameter AngularVelocity_rpm N_nom = 1500 "Nominal rotational speed" 
        annotation(Dialog(group="Characteristics"));
      parameter Medium.Density d_nom = 1000 "Nominal fluid density" 
        annotation(Dialog(group="Characteristics"));
      parameter Integer Np_nom(min=1) = 1 "Nominal number of pumps in parallel";
      parameter SI.Mass M = 0 "Fluid mass inside the pump";
      parameter Boolean checkValve=true "Reverse flow stopped";
      parameter Boolean allowFlowReversal = true 
        "Flow reversal at the ports is allowed by the equations";
      parameter Boolean computeNPSHa=false 
        "Compute NPSH Available at the inlet";
      parameter Medium.AbsolutePressure pin_start "Inlet Pressure Start Value" 
        annotation(Dialog(tab="Initialization"));
      parameter Medium.AbsolutePressure pout_start 
        "Outlet Pressure Start Value" 
        annotation(Dialog(tab="Initialization"));
      parameter Boolean use_T_start = true 
        "Use T_start if true, otherwise h_start" 
        annotation(Dialog(tab = "Initialization"), Evaluate = true);
      parameter Medium.Temperature T_start=
        if use_T_start then 293.15 else Medium.temperature_phX(pin_start,h_start,Medium.reference_X[1:Medium.nXi]) 
        "Start value of temperature" 
        annotation(Dialog(tab = "Initialization", enable = use_T_start));
      parameter Medium.SpecificEnthalpy h_start=
        if use_T_start then Medium.specificEnthalpy_pTX(pin_start, T_start, Medium.reference_X[1:Medium.nXi]) else 1e4 
        "Start value of specific enthalpy" 
        annotation(Dialog(tab = "Initialization", enable = not use_T_start));
      parameter SI.MassFlowRate m_flow_start = 0 
        "Start value of mass flow rate (total)" 
        annotation(Dialog(tab="Initialization"));
      constant SI.Acceleration g=Modelica.Constants.g_n;
    //  parameter Choices.Init.Options.Temp initOpt=Choices.Init.Options.noInit 
    //    "Initialisation option";
      Modelica_Fluid.Interfaces.FluidPort_a inlet(
        redeclare package Medium = Medium,
        p(start=pin_start),
        m_flow(start=m_flow_start, min=if allowFlowReversal and not checkValve then 
                    -inf else 0)) 
      annotation (extent=[-100,-40; -60,0]);
      Modelica_Fluid.Interfaces.FluidPort_b outlet(
        redeclare package Medium = Medium,
        p(start=pout_start),
        m_flow(start=-m_flow_start, max=if allowFlowReversal and not checkValve then 
                    +inf else 0)) 
      annotation (extent=[40,12; 80,52]);
      SI.Pressure dp = outlet.p - inlet.p "Pressure increase";
      SI.Height head = dp/(d*g) "Pump head";
      Medium.Density d "Liquid density at the inlet";
      Medium.SpecificEnthalpy h_out(start=h_start) 
        "Enthalpy of the liquid flowing out of the pump";
      Medium.Temperature Tin "Liquid inlet temperature";
      SI.MassFlowRate m_flow = inlet.m_flow "Mass flow rate (total)";
      SI.MassFlowRate m_flow_single = m_flow/Np "Mass flow rate (single pump)";
      SI.VolumeFlowRate q_flow = m_flow/d "Volume flow rate (total)";
      SI.VolumeFlowRate q_flow_single = q_flow/Np 
        "Volume flow rate (single pump)";
      AngularVelocity_rpm N "Shaft rotational speed";
      Integer Np(min=1) "Number of pumps in parallel";
      SI.Power W_single "Power Consumption (single pump)";
      SI.Power W_tot = W_single*Np "Power Consumption (total)";
      constant SI.Power W_eps=1e-8 
        "Small coefficient to avoid numerical singularities in efficiency computations";
      Real eta "Global Efficiency";
      SI.Length NPSHa "Net Positive Suction Head available";
      Medium.AbsolutePressure pv "Saturation pressure of inlet liquid";
      Real s(start = m_flow_start) 
        "Curvilinear abscissa for the flow curve in parametric form";
      Modelica.Blocks.Interfaces.IntegerInput in_Np 
        annotation (extent=[16,34; 36,54], rotation=-90);
    equation 
      // Number of pumps in parallel
      Np = in_Np;
      if cardinality(in_Np)==0 then
        in_Np = Np_nom "Number of pumps selected by parameter";
      end if;
      
      // Flow equations
      if noEvent(s > 0 or (not checkValve)) then
        // Flow characteristics when check valve is open
        q_flow_single = s;
        head = noEvent((((if abs(N) > 1e-6 then N else 1e-6))/N_nom)^2*flowCharacteristic(q_flow_single*N_nom/((if abs(N) > 1e-6 then N else 1e-6))));
      else
        // Flow characteristics when check valve is closed
        head = (N/N_nom)^2*flowCharacteristic(0) - s;
        q_flow_single = 0;
      end if;
      
      // Power consumption  
      if usePowerCharacteristic then
        W_single = (N/N_nom)^3*(d/d_nom)*powerCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-6))) 
          "Power consumption (single pump)";
        eta = (dp*q_flow_single)/(W_single + W_eps) "Hydraulic efficiency";
      else
        eta = efficiencyCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-10)));
        W_single = dp*q_flow/eta;
      end if;
      // Fluid properties
      fluid.p = inlet.p;
      fluid.h = inlet.h;
      fluid.Xi = inlet.Xi;
      d = fluid.d;
      Tin = fluid.T;
      
      // Mass and energy balances
      inlet.m_flow + outlet.m_flow = 0 "Mass balance";
      inlet.mXi_flow + outlet.mXi_flow = zeros(Medium.nXi) 
        "Substance mass balance";
      inlet.H_flow=semiLinear(inlet.m_flow,inlet.h,h_out) 
        "Enthalpy flow at the inlet";
      outlet.H_flow=semiLinear(outlet.m_flow,outlet.h,h_out) 
        "Enthalpy flow at the outlet";
      if M > 0 then
        M * der(h_out) = m_flow_single*(inlet.h - outlet.h) + W_single 
          "Dynamic energy balance (density variations neglected)";
      else
        inlet.H_flow + outlet.H_flow + W_single*Np = 0 "Static energy balance";
      end if;
      
      // NPSH computations
      if computeNPSHa then
          pv = SatMedium.saturationPressure(fluid.T);
        NPSHa = (inlet.p - pv)/(d*Modelica.Constants.g_n);
      else
        pv = 0;
        NPSHa = 0;
      end if;
    /*
initial equation 
  if initOpt == Choices.Init.Options.noInit then
    // do nothing
  elseif initOpt == Choices.Init.Options.steadyState then
    if ThermalCapacity then
      der(h)=0;
    end if;
  else
    assert(false, "Unsupported initialisation option");
  end if;
*/
      annotation (
        Icon(
          Polygon(points=[-40,-64; -60,-100; 60,-100; 40,-64; -40,-64],
              style(pattern=0, fillColor=74)),
          Ellipse(extent=[-60,40; 60,-80],   style(gradient=3)),
          Polygon(points=[-30,12; -30,-48; 48,-20; -30,12],   style(
              pattern=0,
              gradient=2,
              fillColor=7)),
          Text(extent=[-100,-110; 100,-136], string="%name"),
          Text(extent=[-10,60; 18,40],  string="Np")),
        Diagram,
        Documentation(info="<HTML>
<p>This is the base model for the <tt>Pump</tt> and <tt>
PumpMech</tt> pump models.
<p>The model describes a centrifugal pump, or a group of <tt>Np</tt> identical pumps in parallel. The pump model is based on the theory of kinematic similarity: the pump characteristics are given for nominal operating conditions (rotational speed and fluid density), and then adapted to actual operating condition, according to the similarity equations. 
<p><b>Modelling options</b></p>
<p> The nominal hydraulic characteristic (head vs. volume flow rate) is given by the the replaceable function <tt>flowCharacteristic</tt>. 
<p> The pump energy balance can be specified in two alternative ways:
<ul>
<li><tt>usePowerCharacteristic = false</tt> (default option): the replaceable function <tt>efficiencyCharacteristic</tt> (efficiency vs. volume flow rate in nominal conditions) is used to determine the efficiency, and then the power consumption. The default is a constant efficiency of 0.8.
<li><tt>usePowerCharacteristic = true</tt>: the replaceable function <tt>powerCharacteristic</tt> (power consumption vs. volume flow rate in nominal conditions) is used to determine the power consumption, and then the efficiency.
</ul>
<p>
Several functions are provided in the package <tt>PumpCharacteristics</tt> to specify the characteristics as a function of some operating points at nominal conditions.
<p>Depending on the value of the <tt>checkValve</tt> parameter, the model either supports reverse flow conditions, or includes a built-in check valve to avoid flow reversal.
<p>If the <tt>in_Np</tt> input connector is wired, it provides the number of pumps in parallel; otherwise,  <tt>Np_n</tt> parallel pumps are assumed.</p>
<p>It is possible to take into account the heat capacity of the fluid inside the pump by specifying its mass <tt>M</tt> at nominal conditions; this is necessary to avoid singularities in the computation of the outlet enthalpy in case of zero flow rate. If zero flow rate conditions are always avoided, this dynamic effect can be neglected by leaving the default value <tt>M = 0</tt>, thus avoiding a fast state variable in the model.
<p>If <tt>computeNPSHa = true</tt>, the available net positive suction head is also computed; this requires a two-phase medium model to provide the fluid saturation pressure.
</HTML>", revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
      
    end PartialPump;
    
    model Controller 
      
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ControllerUtilities.Port_Sensors
        sensors 
        annotation (extent=[-280,-40; -200,40]);
      annotation (Diagram, Coordsys(extent=[-200,-200; 200,200]),
        Icon(
          Rectangle(extent=[-200,200; 200,-200], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=7,
              rgbfillColor={255,255,255})),
          Text(
            extent=[-288,286; 262,208],
            style(color=3, rgbcolor={0,0,255}),
            string="%name"),
          Line(points=[-48,0; 0,0],         style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Rectangle(extent=[-170,60; -50,-60],  style(color=0, rgbcolor={0,0,0})),
          Line(points=[0,40; 0,-40],   style(color=0, rgbcolor={0,0,0})),
          Line(points=[0,0; 26,0],       style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Polygon(points=[26,10; 50,0; 26,-10; 26,10],         style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0})),
          Rectangle(extent=[50,60; 170,-60],  style(color=0, rgbcolor={0,0,0})),
          Polygon(points=[-24,10; 0,0; -24,-10; -24,10],       style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0}))));
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ControllerUtilities.Port_Actuators
        actuators 
        annotation (extent=[200,-20; 240,20]);
      
      parameter Real w_dilution=0.003;
      parameter Real w_concentrate=0.005;
      parameter Real startTime=1;
      parameter Real T5_batch_level=0.211;
      
      Modelica.StateGraph.InitialStep InitialStep1 
        annotation (extent=[-180,90; -160,110],  rotation=0);
      Modelica.StateGraph.Transition Transition1(enableTimer=true, waitTime=
            startTime) annotation (extent=[-150,90; -130,110],  rotation=0);
      Modelica.StateGraph.Step Step1 
        annotation (extent=[-120,90; -100,110],  rotation=0);
      Modelica.StateGraph.Transition Transition2(condition=LIS_301 >= 0.13) 
        annotation (extent=[-90,90; -70,110],  rotation=0);
      Modelica.StateGraph.Step Step2 annotation (extent=[-60,90; -40,110]);
      Modelica.StateGraph.Transition Transition3(
        condition=true,
        enableTimer=true,
        waitTime=500) 
        annotation (extent=[-30,90; -10,110]);
      Modelica.StateGraph.Step Step3 annotation (extent=[0,90; 20,110]);
      Modelica.StateGraph.Transition Transition4(condition=LIS_301 <= 0.01) 
        annotation (extent=[30,90; 50,110]);
      Modelica.StateGraph.Step Step4 annotation (extent=[60,90; 80,110]);
      Modelica.StateGraph.Transition Transition5(condition=T5_idle) 
        annotation (extent=[90,90; 110,110]);
      Modelica.StateGraph.Step Step5 annotation (extent=[120,90; 140,110]);
      Modelica.StateGraph.Transition Transition6(condition=LIS_501 >=
            T5_batch_level) annotation (extent=[150,90; 170,110]);
      Modelica.StateGraph.Step Step6 annotation (extent=[-120,30; -100,50]);
      Modelica.StateGraph.Transition Transition7(
        condition=true,
        enableTimer=true,
        waitTime=300) 
        annotation (extent=[-90,30; -70,50]);
      Modelica.StateGraph.Parallel Parallel1 annotation (extent=[-176,-100; 194,
            0]);
      Modelica.StateGraph.Step Step7 annotation (extent=[-122,-80; -102,-60]);
      Modelica.StateGraph.Step Step8 annotation (extent=[-62,-80; -42,-60]);
      Modelica.StateGraph.Step Step9 annotation (extent=[-2,-80; 18,-60]);
      Modelica.StateGraph.Step Step10 annotation (extent=[58,-80; 78,-60]);
      Modelica.StateGraph.Step Step11 annotation (extent=[118,-80; 138,-60]);
      Modelica.StateGraph.Step Step12 annotation (extent=[-62,-40; -42,-20]);
      Modelica.StateGraph.Step Step13 annotation (extent=[-2,-40; 18,-20]);
      Modelica.StateGraph.Step Step14 annotation (extent=[58,-40; 78,-20]);
      Modelica.StateGraph.Transition Transition8(condition=T7_idle) 
        annotation (extent=[-92,-80; -72,-60]);
      Modelica.StateGraph.Transition Transition9(condition=LIS_501 <= 0.01) 
        annotation (extent=[-32,-80; -12,-60]);
      Modelica.StateGraph.Transition Transition10(condition=TIS_702 <= 298) 
        annotation (extent=[28,-80; 48,-60]);
      Modelica.StateGraph.Transition Transition11(condition=LIS_701 <= 0.01) 
        annotation (extent=[88,-80; 108,-60]);
      Modelica.StateGraph.Transition Transition12(condition=TIS_602 <= 298) 
        annotation (extent=[-32,-40; -12,-20]);
      Modelica.StateGraph.Transition Transition13(condition=LIS_601 <= 0.01) 
        annotation (extent=[28,-40; 48,-20]);
      
      Real LIS_301;
      Real LIS_501;
      Real LIS_601;
      Real LIS_701;
      Real QI_302;
      Real QIS_502;
      Real TIS_602;
      Real TIS_702;
      Boolean T5_idle;
      Boolean T7_idle;
      Modelica.StateGraph.TransitionWithSignal TransitionWithSignal1 
        annotation (extent=[-12,-160; 8,-140],  rotation=180);
      Modelica.Blocks.Sources.BooleanExpression BooleanExpression1(y=time >
            2500) 
        annotation (extent=[-104,-148; -18,-116]);
    equation 
      LIS_301 = sensors.LIS_301;
      LIS_501 = sensors.LIS_501;
      LIS_601 = sensors.LIS_601;
      LIS_701 = sensors.LIS_701;
      QI_302 = sensors.QI_302;
      QIS_502 = sensors.QIS_502;
      TIS_602 = sensors.TIS_602;
      TIS_702 = sensors.TIS_702;
      T5_idle = not actuators.V12 and not actuators.V15 and not actuators.T5_Heater
         and sensors.LIS_501 < 0.01;
      T7_idle = not actuators.V15 and not actuators.V18 and not actuators.
        T7_Cooling and sensors.LIS_701 < 0.01;
      
      actuators.P1 = Step10.active;
      actuators.P2 = Step13.active;
      actuators.T5_Heater = Step6.active;
      actuators.T7_Cooling = Step9.active;
      actuators.T6_Cooling = Step12.active;
      actuators.V1 = Step10.active;
      actuators.V2 = false;
      actuators.V3 = Step10.active;
      actuators.V4 = false;
      actuators.V5 = Step13.active;
      actuators.V6 = Step13.active;
      actuators.V8 = Step1.active;
      actuators.V9 = Step2.active;
      actuators.V10 = false;
      actuators.V11 = Step3.active;
      actuators.V12 = Step5.active;
      actuators.V15 = Step8.active;
      actuators.V18 = Step10.active;
      actuators.V19 = false;
      actuators.V20 = Step13.active;
      actuators.V21 = false;
      actuators.V22 = Step10.active;
      actuators.V23 = Step10.active;
      actuators.V25 = Step13.active;
      actuators.V24 = Step13.active;
      
      connect(InitialStep1.outPort[1], Transition1.inPort) annotation (points=[-159.5,
            100; -144,100],        style(color=0, rgbcolor={0,0,0}));
      connect(Transition1.outPort, Step1.inPort[1]) annotation (points=[-138.5,
            100; -121,100],
                       style(color=0, rgbcolor={0,0,0}));
      connect(Step1.outPort[1], Transition2.inPort) 
        annotation (points=[-99.5,100; -84,100], style(color=0, rgbcolor={0,0,0}));
      connect(Transition2.outPort, Step2.inPort[1]) 
        annotation (points=[-78.5,100; -61,100], style(color=0, rgbcolor={0,0,0}));
      connect(Step2.outPort[1], Transition3.inPort) 
        annotation (points=[-39.5,100; -24,100], style(color=0, rgbcolor={0,0,0}));
      connect(Transition3.outPort, Step3.inPort[1]) 
        annotation (points=[-18.5,100; -1,100], style(color=0, rgbcolor={0,0,0}));
      connect(Step3.outPort[1], Transition4.inPort) 
        annotation (points=[20.5,100; 36,100], style(color=0, rgbcolor={0,0,0}));
      connect(Transition4.outPort, Step4.inPort[1]) 
        annotation (points=[41.5,100; 59,100], style(color=0, rgbcolor={0,0,0}));
      connect(Step4.outPort[1], Transition5.inPort) 
        annotation (points=[80.5,100; 96,100], style(color=0, rgbcolor={0,0,0}));
      connect(Transition5.outPort, Step5.inPort[1]) 
        annotation (points=[101.5,100; 119,100], style(color=0, rgbcolor={0,0,0}));
      connect(Step5.outPort[1], Transition6.inPort) 
        annotation (points=[140.5,100; 156,100], style(color=0, rgbcolor={0,0,0}));
      connect(Transition6.outPort, Step6.inPort[1]) annotation (points=[161.5,
            100; 184,100; 184,70; -160,70; -160,40; -121,40],style(color=0,
            rgbcolor={0,0,0}));
      connect(Step6.outPort[1], Transition7.inPort) 
        annotation (points=[-99.5,40; -84,40],   style(color=0, rgbcolor={0,0,0}));
      connect(Step12.inPort[1], Parallel1.split[1]) annotation (points=[-63,-30;
            -134.375,-30; -134.375,-25],
                                  style(color=0, rgbcolor={0,0,0}));
      connect(Step12.outPort[1], Transition12.inPort) 
        annotation (points=[-41.5,-30; -26,-30],
                                               style(color=0, rgbcolor={0,0,0}));
      connect(Transition12.outPort, Step13.inPort[1]) 
        annotation (points=[-20.5,-30; -3,-30],style(color=0, rgbcolor={0,0,0}));
      connect(Step13.outPort[1], Transition13.inPort) 
        annotation (points=[18.5,-30; 34,-30],
                                             style(color=0, rgbcolor={0,0,0}));
      connect(Transition13.outPort, Step14.inPort[1]) 
        annotation (points=[39.5,-30; 57,-30],
                                             style(color=0, rgbcolor={0,0,0}));
      connect(Step14.outPort[1], Parallel1.join[1]) annotation (points=[78.5,-30;
            152.375,-30; 152.375,-25],
                                style(color=0, rgbcolor={0,0,0}));
      connect(Step7.inPort[1], Parallel1.split[2]) annotation (points=[-123,-70;
            -138,-70; -138,-75; -134.375,-75],style(color=0, rgbcolor={0,0,0}));
      connect(Step7.outPort[1], Transition8.inPort) annotation (points=[-101.5,
            -70; -86,-70],
                      style(color=0, rgbcolor={0,0,0}));
      connect(Transition8.outPort, Step8.inPort[1]) 
        annotation (points=[-80.5,-70; -63,-70], style(color=0, rgbcolor={0,0,0}));
      connect(Step8.outPort[1], Transition9.inPort) 
        annotation (points=[-41.5,-70; -26,-70], style(color=0, rgbcolor={0,0,0}));
      connect(Transition9.outPort, Step9.inPort[1]) 
        annotation (points=[-20.5,-70; -3,-70],  style(color=0, rgbcolor={0,0,0}));
      connect(Step9.outPort[1], Transition10.inPort) 
        annotation (points=[18.5,-70; 34,-70], style(color=0, rgbcolor={0,0,0}));
      connect(Transition10.outPort, Step10.inPort[1]) 
        annotation (points=[39.5,-70; 57,-70], style(color=0, rgbcolor={0,0,0}));
      connect(Step10.outPort[1], Transition11.inPort) 
        annotation (points=[78.5,-70; 94,-70], style(color=0, rgbcolor={0,0,0}));
      connect(Transition11.outPort, Step11.inPort[1]) 
        annotation (points=[99.5,-70; 117,-70], style(color=0, rgbcolor={0,0,0}));
      connect(Step11.outPort[1], Parallel1.join[2]) annotation (points=[138.5,
            -70; 154,-70; 154,-75; 152.375,-75],
                                           style(color=0, rgbcolor={0,0,0}));
      connect(Transition7.outPort, Parallel1.inPort) annotation (points=[-78.5,40;
            -40,40; -40,20; -190,20; -190,-50; -181.55,-50],
                                                           style(color=0, rgbcolor=
              {0,0,0}));
      connect(TransitionWithSignal1.inPort, Parallel1.outPort) annotation (points=[2,-150;
            208,-150; 208,-50; 197.7,-50],       style(color=0, rgbcolor={0,0,0}));
      connect(TransitionWithSignal1.outPort, InitialStep1.inPort[1]) annotation (
          points=[-3.5,-150; -194,-150; -194,100; -181,100], style(color=0,
            rgbcolor={0,0,0}));
      connect(BooleanExpression1.y, TransitionWithSignal1.condition) annotation (
          points=[-13.7,-132; -2,-132; -2,-138],                   style(color=5,
            rgbcolor={255,0,255}));
    end Controller;
    
    package ControllerUtilities 
      extends Modelica.Icons.Library;
      class Adapter_Inference 
        Port_IdleTanks idleTanks;
      end Adapter_Inference;
      
      class Adapter_Superposition 
        Port_Actuators actuators;
      end Adapter_Superposition;
      
      class Block_Recipe_TBD 
        parameter Real startTime;
        parameter Real w_dilution=0.003;
        parameter Real w_concentrat=0.005;
        parameter Real T3_batch_level=0.1273;
        parameter Real T5_batch_level=0.211;
        Boolean trig;
        Boolean S0(start=true);
        Boolean S1;
        Boolean S2;
        Boolean S3;
        Boolean S4;
        Boolean S5;
        Boolean S6;
        Boolean S7;
        Boolean S8;
        Boolean S9;
        Boolean S10;
        Boolean S11;
        Boolean S12;
        Boolean S13;
        Boolean S14;
        Boolean tr0;
        Boolean tr1;
        Boolean tr2;
        Boolean tr3;
        Boolean tr4;
        Boolean tr5;
        Boolean tr6;
        Boolean tr7;
        Boolean tr8;
        Boolean tr9;
        Boolean tr10;
        Boolean tr11;
        Boolean tr12;
        Boolean tr13;
        Port_Actuators act annotation (extent=[-110, -10; -90, 10]);
      end Block_Recipe_TBD;
      
      class BlockMain 
        Boolean trig;
        
        Port_Actuators actuators annotation (extent=[90, -10; 110, 10]);
        Block_Recipe_TBD Recipe1 annotation (extent=[-50,10; -10,50]);
        Block_Recipe_TBD Recipe2 annotation (extent=[10,10; 50,50]);
        Adapter_Inference Inference annotation (extent=[-50,-50; -10,-10]);
        Adapter_Superposition Superposition annotation (extent=[10,-50; 50,-10]);
      end BlockMain;
      
      class Buffer_Recipe_TBD 
        Port_Actuators act;
        Boolean S0;
        Boolean S1;
        Boolean S2;
        Boolean S3;
        Boolean S4;
        Boolean S5;
        Boolean S6;
        Boolean S7;
        Boolean S8;
        Boolean S9;
        Boolean S10;
        Boolean S11;
        Boolean S12;
        Boolean S13;
        Boolean S14;
      end Buffer_Recipe_TBD;
      
      class BufferMain 
        Buffer_Recipe_TBD Recipe1;
        Buffer_Recipe_TBD Recipe2;
      end BufferMain;
      
      connector Port_Actuators 
        Boolean P1;
        Boolean P2;
        Boolean T5_Heater;
        Boolean T7_Cooling;
        Boolean T6_Cooling;
        Boolean V1;
        Boolean V2;
        Boolean V3;
        Boolean V4;
        Boolean V5;
        Boolean V6;
        Boolean V8;
        Boolean V9;
        Boolean V10;
        Boolean V11;
        Boolean V12;
        Boolean V15;
        Boolean V18;
        Boolean V19;
        Boolean V20;
        Boolean V21;
        Boolean V22;
        Boolean V23;
        Boolean V24;
        Boolean V25;
        
        annotation (Icon(
             Polygon(points=[-100,100; 100,0; -100,-100; -100,100], style(
                color=0,
                rgbcolor={0,0,0},
                thickness=2,
                fillColor=7,
                rgbfillColor={255,255,255}))), Diagram(
                Polygon(points=[0,50; 100,0; 0,-50; 0,50], style(
                color=0,
                rgbcolor={0,0,0},
                thickness=2,
                fillColor=7,
                rgbfillColor={255,255,255}))));
      end Port_Actuators;
      
      connector Port_IdleTanks 
        Boolean T5_idle;
        Boolean T7_idle;
      end Port_IdleTanks;
      
      connector Port_Sensors 
        Real LIS_301;
        Real QI_302;
        Real LIS_501;
        Real QIS_502;
        Real TI_503;
        Real LIS_601;
        Real TIS_602;
        Real LIS_701;
        Real TIS_702;
        annotation (Icon(Polygon(points=[-100,100; -100,-100; 100,0; -100,100],
                style(
                color=0,
                rgbcolor={0,0,0},
                thickness=2,
                fillColor=30,
                rgbfillColor={215,215,215}))), Diagram(Polygon(points=[0,50; 0,
                  -50; 100,0; 0,50], style(
                color=0,
                rgbcolor={0,0,0},
                thickness=2,
                fillColor=30,
                rgbfillColor={215,215,215}))));
      end Port_Sensors;
    end ControllerUtilities;
    
  model CoolingTank 
      "Open tank with top and bottom inlet/outlet ports at a defineable height, heat transfer through the walls and a cooling unit via heatPort connector" 
      import SI = Modelica.SIunits;
      import Modelica.Constants;
      import Modelica_Fluid.PressureLosses.BaseClasses.lossConstant_D_zeta;
      import Modelica_Fluid.Utilities.regRoot2;
      import Modelica_Fluid.Volumes.BaseClasses.TankPortData;
      
    replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
      
    SI.Height level(stateSelect=StateSelect.prefer, start=level_start) 
        "Fluid level in the tank";
      
  //Tank geometry  
      parameter SI.Height levelMax "Maximum level of tank before it overflows";
      parameter SI.Area area "Area of tank";
      parameter SI.Volume V0=0 "Volume of the liquid when level = 0";
      
  //Port definitions 
      parameter Integer nTopPorts(min=1) = 1 
        "Number of inlet ports above levelMax (>= 1)";
      
      Modelica_Fluid.Interfaces.FluidPort_a topPorts[nTopPorts](
      redeclare package Medium = Medium,
      m_flow(each start=0, each min=0),
      mXi_flow(each start=0, each min=0)) if   nTopPorts > 0 
        "Inlet ports over levelMax at top of tank (fluid flows only from the port in to the tank)"
      annotation (extent=[-10,90; 10,110]);
      
      parameter Modelica_Fluid.Volumes.BaseClasses.TankPortData portsData[:] = {TankPortData(diameter=0)} 
        "Data of inlet/outlet ports at side and bottom of tank";
      
      Modelica_Fluid.Interfaces.FluidPort_b ports[size(portsData,1)](
      redeclare package Medium = Medium,
      m_flow(each start=0),
      mXi_flow(each start=0)) 
        "inlet/outlet ports at bottom or side of tank (fluid flows in to or out of port; a port might be above the fluid level)"
      annotation (extent=[-10,-110; 10,-90]);
      
  // Heat transfer  
      parameter SI.CoefficientOfHeatTransfer alpha0=0 
        "Coefficient of heat transfer of tank wall";
     Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation (extent=[-110,-10; -90,10]);
      
  //Ambient  
     outer Modelica_Fluid.Ambient ambient "Ambient conditions";
     parameter Medium.AbsolutePressure p_ambient=ambient.default_p_ambient 
        "Tank surface pressure" 
      annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
     parameter Medium.Temperature T_ambient=ambient.default_T_ambient 
        "Tank surface Temperature" 
      annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
      
  //Initialization
      parameter Types.Init.Temp initType=Types.Init.InitialValues 
        "Initialization option" 
      annotation(Evaluate=true,Dialog(tab = "Ambient and Initialization", group = "Initialization"));
      parameter SI.Height level_start(min=0) "Start value of tank level" 
      annotation(Dialog(tab="Ambient and Initialization", group = "Initialization"));
      parameter Medium.Temperature T_start=T_ambient 
        "Start value of temperature" 
      annotation(Dialog(tab = "Ambient and Initialization", group = "Initialization"));
      parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
        "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Ambient and Initialization", group = "Initialization", enable=Medium.nXi > 0));
      
  // Advanced  
      parameter Real hysteresisFactor(min=0) = 0.1 
        "Hysteresis for empty pipe = diameter*hysteresisFactor" 
      annotation(Dialog(tab="Advanced", group="Numerical properties"));
      parameter SI.MassFlowRate m_flow_small(min=0) = 1e-5 
        "Regularization range at zero mass flow rate" 
      annotation(Dialog(tab="Advanced", group="Numerical properties"));
      parameter Boolean stiffCharacteristicForEmptyPort = true 
        "=true, if steep pressure loss characteristic for empty pipe port" 
      annotation(Dialog(tab="Advanced", group="Numerical properties"), Evaluate=true);
      parameter Real zetaLarge(min=0) = 1e5 
        "Large pressure loss factor if mass flows out of empty pipe port" 
      annotation(Dialog(tab="Advanced", group="Numerical properties", enable=stiffCharacteristicForEmptyPort));
      
  //Tank properties  
       final parameter Integer nPorts = size(ports,1) 
        "Number of inlet/outlet ports";
       final parameter Medium.SpecificEnthalpy h_start=Medium.specificEnthalpy_pTX(
          p_ambient,
          T_start,
          X_start) annotation(Hide=true);
      Medium.BaseProperties medium(
        preferredMediumStates=true,
        p(start=p_ambient),
        T(start=T_start),
        h(start=h_start),
        Xi(start=X_start[1:Medium.nXi]));
      SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
      SI.Energy U(stateSelect=StateSelect.never) 
        "Internal energy of tank volume";
      SI.Heat Q_lost "Heat lost through the walls";
      SI.Mass m(stateSelect=StateSelect.never) "Mass of fluid in tank";
      SI.Mass mXi[Medium.nXi](each stateSelect=StateSelect.never) 
        "Masses of independent components in the fluid";
      
    protected 
      parameter SI.Area bottomArea[nPorts]=Constants.pi*{(portsData[i].diameter/2)^2 for i in 1:nPorts};
      parameter SI.Diameter ports_emptyPipeHysteresis[nPorts] = portsData.diameter*hysteresisFactor;
      SI.Length levelAbovePort[nPorts] "Height of fluid over bottom ports";
      Boolean ports_m_flow_out[nPorts](each start = true, each fixed=true);
      Boolean aboveLevel[nPorts] "= true, if level >= ports[i].portLevel";
      Real zeta_out[nPorts];
      
  equation 
    assert(level <= levelMax, "Tank starts to overflow (level = levelMax = " + String(level) + ")");
    assert(m>=0, "Mass in tank is zero");
      
    // Total quantities
      medium.T = heatPort.T;
      medium.p = p_ambient;
      V = area*level + V0 "Volume of fluid";
      m = V*medium.d "Mass of fluid";
      mXi = m*medium.Xi "Mass of fluid components";
      U = m*medium.u "Internal energy of fluid";
      Q_lost = -alpha0*(area+2*sqrt(area*Modelica.Constants.pi)*level)*(medium.T-T_ambient) 
        "Q=-k*A*dT";
      
    // Mass balances
      der(m) = sum(topPorts.m_flow) + sum(ports.m_flow);
      for i in 1:Medium.nXi loop
        der(mXi[i]) = sum(topPorts.mXi_flow[i]) + sum(ports.mXi_flow[i]);
      end for;
      
    // Energy balance
      if Medium.singleState then
        der(U) = sum(topPorts.H_flow) + sum(ports.H_flow) + Q_lost + heatPort.Q_flow;
                                 //Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)
      else
        der(U) = sum(topPorts.H_flow) + sum(ports.H_flow) - p_ambient*der(V) + Q_lost + heatPort.Q_flow;
      end if;
      
    // Properties at top ports
      for i in 1:nTopPorts loop
         // It is assumed that fluid flows only into one of the top ports and never out of it 
         topPorts[i].H_flow   = semiLinear(topPorts[i].m_flow, topPorts[i].h, h_start);
         topPorts[i].mXi_flow = semiLinear(topPorts[i].m_flow, topPorts[i].Xi, X_start[1:Medium.nXi]);
         topPorts[i].p        = p_ambient;
  /*
       assert(topPorts[i].m_flow > -1, "Mass flows out of tank via topPorts[" + String(i) + "]\n" +
                                         "This indicates a wrong model");
*/
      end for;
      
    // Properties at bottom ports
      for i in 1:nPorts loop
         ports[i].H_flow = semiLinear(ports[i].m_flow, ports[i].h, medium.h);
         ports[i].mXi_flow = semiLinear(ports[i].m_flow, ports[i].Xi, medium.Xi);
         aboveLevel[i] = level >= (portsData[i].portLevel + ports_emptyPipeHysteresis[i])
                         or pre(aboveLevel[i]) and level >= (portsData[i].portLevel - ports_emptyPipeHysteresis[i]);
         levelAbovePort[i] = if aboveLevel[i] then level - portsData[i].portLevel else 0;
        
         if stiffCharacteristicForEmptyPort then
            // If port is above fluid level, use large zeta if fluid flows out of port (= small mass flow rate)
            zeta_out[i] = 1 + (if aboveLevel[i] then 0 else zetaLarge);
            ports[i].p = p_ambient + levelAbovePort[i]*ambient.g*medium.d
                                 + Modelica_Fluid.Utilities.regSquare2(ports[i].m_flow, m_flow_small,
                                    0, lossConstant_D_zeta(portsData[i].diameter, zeta_out[i])/medium.d);
            ports_m_flow_out[i] = false;
          
         else
            // Handling according to Remelhe/Poschlad
            ports_m_flow_out[i] = (pre(ports_m_flow_out[i]) and not ports[i].p>p_ambient)
                                       or ports[i].m_flow < -1e-6;
           if aboveLevel[i] then
               ports[i].p = p_ambient + levelAbovePort[i]*ambient.g*medium.d -
                                 smooth(2,noEvent(if ports[i].m_flow < 0 then ports[i].m_flow^2/
                                       (2*medium.d*bottomArea[i]^2) else 0));
           else
              if pre(ports_m_flow_out[i]) then
                 ports[i].m_flow = 0;
              else
                 ports[i].p = p_ambient;
              end if;
           end if;
            zeta_out[i] =0;
         end if;
       end for;
      
  initial equation 
      for i in 1:nPorts loop
         pre(aboveLevel[i]) = level_start >= portsData[i].portLevel;
      end for;
      
      if initType == Types.Init.NoInit then
      // no initial equations
      elseif initType == Types.Init.InitialValues then
        level = level_start;
        medium.T = T_start;
        medium.Xi = X_start[1:Medium.nXi];
      elseif initType == Types.Init.SteadyState then
        der(level) = 0;
        der(medium.T) = 0;
        der(medium.Xi) = zeros(Medium.nXi);
      elseif initType == Types.Init.SteadyStateHydraulic then
        der(level) = 0;
        medium.T = T_start;
        medium.Xi = X_start[1:Medium.nXi];
      else
        assert(false, "Unsupported initialization option");
      end if;
      
      annotation (
        Icon(
          Rectangle(extent=[-100,-100; 100,100], style(
              color=7,
              rgbcolor={255,255,255},
              fillColor=7,
              rgbfillColor={255,255,255})),
          Line(points=[-100,100; -100,-100; 100,-100; 100,100], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=69,
              rgbfillColor={0,127,255},
              fillPattern=1)),
            Rectangle(extent=DynamicSelect([-100,-100; 100,0], [-100,-100; 100,(-100
                   + 200*level/levelMax)]), style(
                color=69,
                rgbcolor={0,127,255},
                fillColor=71,
                rgbfillColor={85,170,255},
                fillPattern=1)),
          Text(
            extent=[-94,19; 96,-1],
          string=DynamicSelect(" ", realString(level, 1, 3)),
          style(color=0, rgbcolor={0,0,0})),
          Line(points=[-100,100; 100,100], style(
              color=0,
              rgbcolor={0,0,0},
              pattern=3)),
        Text(
            extent=[-94,90;95,60],
            style(color=3, rgbcolor={0,0,255}),
            string="%name"),
        Text(
          extent=[-95,-85; 95,-65],
          style(color=0),
            string="%level_start"),
        Text(
          extent=[-95,-55; 95,-35],
          style(color=0),
            string="level_start ="),
        Text(extent=[-95,50; 95,30], string="level =",
          style(color=0, rgbcolor={0,0,0}))),
        Documentation(info="<HTML>
<p> 
Model of a tank that is open to the environment at the fixed pressure
<tt>p_ambient</tt>. Heat transfer to the environment and to 
the tank walls is neglected.
The tank is filled with a single or multiple-substance liquid, 
assumed to have uniform temperature and mass fractions.
</p> 
 
<p>
At the top of the tank over the maximal fill level <b>levelMax</b> 
a vector of FluidPorts, called <b>topPorts</b>, is present.
The assumption is made that fluid flows always in to the tank via these
ports (and never back in to the connector).
If the tank has no top ports, set <b>nTopPorts</b> = 1, and do not
connect to this port (the default connection semantics of Modelica
leads to a behaviour as if the port would not be present; 
the reason is that some tools do currently no support zero sized
connectors).
</p>
 
<p>
The vector of connectors <b>ports</b> are fluid ports at the bottom
and side of the tank at a defineable height. Fluid can flow either out
of or in to this port. The fluid level of the tank may be below
one of these ports. This case is approximated by introducing a
large pressure flow coefficient so that the mass flow rate
through this port is very small in this case.
</p>
 
<p>
If the tank starts to over flow (i.e., level > levelMax), an
assertion is triggered.
</p>
 
<p>
When the diagram layer is open in the plot environment, the
level of the tank is dynamically visualized. Note, the speed
of the diagram animation in Dymola can be set via command
<b>animationSpeed</b>(), e.g., animationSpeed(speed = 10)
</p>
</HTML>",   revisions="<html>
<ul>
<li><i>Jul. 29, 2006</i> by Martin Otter (DLR):<br> 
   Improved handling of ports that are above the fluid level and
   simpler implementation.</li>
 
<li><i>Jan. 6, 2006</i> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund), 
   Martin Otter (DLR):<br> 
   Implementation based on former tank model but with several improvements
   (top, bottom, side ports; correctly treating kinetic energy for outlet
   and total dissipation for inlet; ports can be above the fluid level).</li>
</ul>
</html>"),
        Diagram,
        uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")),
        Coordsys(grid=[1,1], scale=0.2));
  end CoolingTank;
  end BaseClasses;
  annotation (Documentation(info="<html>
<p>
The process under consideration is an evaporation plant for a 
student lab at the Process Control Laboratory (AST) of the 
University of Dortmund that evaporates a water sodium chloride 
mixture so that a higher concentrated solution is produced. 
The task of the students is to learn how to program the process 
control system. A picture of the batch plant is shown in the figure
below.
</p>

<p>
The flow sheet diagram is shown in figure 2.
</p>

</html>"));
end AST_BatchPlant;
