within FluidSandbox.Examples.StandardConnectionSemantics;
model UncontrolledWaterSteamCycle 
  "Example inspired by PowerFluid: Simplified, without control" 
  extends Icons.Example;
  Sources.ControlledPump pump(   redeclare package Medium = 
        Modelica.Media.Water.StandardWater, redeclare package FluidInterface = 
        FluidInterface) 
    annotation (extent=[-80,10; -60,30]);
  HeatTransfer.EvaporatingVessel evaporator(
    cp_D=500,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_D=100e3,
    V_t=50,
    V_l_start=33,
    redeclare package FluidInterface = FluidInterface, 
    initType=Modelica_Fluid.Types.Init.InitialValues, 
    p_start=100000)     annotation (extent=[-40,10; -20,30]);
  Modelica.Blocks.Sources.TimeTable fuelTable(table=[0,0; 5400,1000; 7210,
        1000]) 
              annotation (extent=[-80,-60; -60,-40]);
public 
  Modelica.Blocks.Math.Gain MW2W(k=1e6) 
    annotation (extent=[-35,-37.5; -25,-26.5], rotation=90);
  Modelica.Thermal.HeatTransfer.PrescribedHeatFlow furnace 
    annotation (extent=[-40,-20; -20,0],     rotation=90);
  Turbomachinery.TurbineStage turbine(
                                  redeclare package Medium = 
        Modelica.Media.Water.StandardWater, K_t=0.001,
    redeclare package FluidInterface = FluidInterface) 
    annotation (extent=[9,32; 29,52]);
  Modelica.Mechanics.Rotational.ConstantSpeed load(w_fixed=-50) 
    annotation (extent=[87,38.5; 72,53.5]);
  HeatTransfer.EvaporatingVessel condenser(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    V_t=100,
    m_D=10e3,
    cp_D=500,
    redeclare package FluidInterface = FluidInterface, 
    p_start=10000) 
              annotation (extent=[20,-70; 40,-90]);
  Modelica.Thermal.HeatTransfer.PrescribedHeatFlow furnace1 
    annotation (extent=[20,-40; 40,-60],     rotation=90);
  Modelica.Blocks.Sources.Ramp const(
    startTime=0,
    height=180,
    duration=3600,
    offset=20) annotation (extent=[-100,60; -80,80]);
  Modelica.Blocks.Sources.Ramp const1(
    height=-4e8,
    duration=7200,
    offset=0,
    startTime=0) annotation (extent=[0,-20; 20,0]);
equation 
  connect(fuelTable.y,MW2W. u) 
                             annotation (points=[-59,-50; -30,-50; -30,-38.6],
                         style(color=74, rgbcolor={0,0,127}));
  connect(turbine.flange,load. flange) 
    annotation (points=[29,48; 50.5,48; 50.5,46; 72,46],
                                         style(color=0, rgbcolor={0,0,0}));
  annotation (Diagram, experiment(StopTime=1200));
  connect(MW2W.y, furnace.Q_flow) annotation (points=[-30,-25.95; -30,-20], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(evaporator.heatPort, furnace.port) annotation (points=[-30,10; -30,0],
      style(
      color=42,
      rgbcolor={191,0,0},
      smooth=0));
  connect(furnace1.port, condenser.heatPort) annotation (points=[30,-60; 30,-70],
      style(
      color=42,
      rgbcolor={191,0,0},
      smooth=0));
  connect(turbine.port_b, condenser.port_b) annotation (points=[29,42; 50,42;
        50,-80; 40,-80], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(condenser.port_a, pump.port_a) annotation (points=[20,-80; -90,-80;
        -90,20; -80,20], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(pump.port_b, evaporator.port_a) annotation (points=[-60,20; -40,20],
      style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(evaporator.port_b, turbine.port_a) annotation (points=[-20,20; 0,20;
        0,42; 9,42], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(const1.y, furnace1.Q_flow) annotation (points=[21,-10; 30,-10; 30,-40],
      style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(const.y, pump.massFlowRate) annotation (points=[-79,70; -70,70; -70,
        29], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
end UncontrolledWaterSteamCycle;
