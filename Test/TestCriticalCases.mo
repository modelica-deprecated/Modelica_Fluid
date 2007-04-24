package TestCriticalCases 
  "Collection of test cases which might be critical for the solvers" 
  
  model BranchingPipes1 
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      p_b_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-72,-4; -52,16]);
    
    ControlValves.ValveIncompressible valveIncompressible(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[10,36; 30,56]);
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[74,-20; 62,-8]);
    Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      p_b_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-40,36; -20,56]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      p_b_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-40,-50; -20,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (extent=[-20,70; 0,90]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 2,1;
          2.01,0; 100,0]) 
                       annotation (extent=[-20,-10; 0,10]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -72,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible1.port_b, sink.port) annotation (points=[28,-40;
          46,-40; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible.port_b, sink.port) annotation (points=[30,46;
          46,46; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (points=[-20,
          -40; 8,-40], style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valveIncompressible.port_a) annotation (points=[-20,
          46; 10,46], style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, pipe1.port_b) annotation (points=[-40,46; -46,46; -46,
          6; -52,6], style(color=69, rgbcolor={0,127,255}));
    connect(pipe1.port_b, pipe3.port_a) annotation (points=[-52,6; -46,6; -46,
          -40; -40,-40], style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (
        points=[1,80; 20,80; 20,55], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (
        points=[1,0; 18,0; 18,-31], style(color=74, rgbcolor={0,0,127}));
  end BranchingPipes1;
  
  model BranchingPipes2 
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-72,-4; -52,16]);
    
    ControlValves.ValveIncompressible valveIncompressible(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[10,36; 30,56]);
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[74,-20; 62,-8]);
    Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-40,36; -20,56]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-40,-50; -20,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,0; 1,0;
          1,1; 100,1]) annotation (extent=[-20,70; 0,90]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,0; 2,0;
          2,1; 100,1]) annotation (extent=[-20,-10; 0,10]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -72,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible1.port_b, sink.port) annotation (points=[28,-40;
          46,-40; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible.port_b, sink.port) annotation (points=[30,46;
          46,46; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (points=[-20,
          -40; 8,-40], style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valveIncompressible.port_a) annotation (points=[-20,
          46; 10,46], style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, pipe1.port_b) annotation (points=[-40,46; -46,46; -46,
          6; -52,6], style(color=69, rgbcolor={0,127,255}));
    connect(pipe1.port_b, pipe3.port_a) annotation (points=[-52,6; -46,6; -46,
          -40; -40,-40], style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (
        points=[1,80; 20,80; 20,55], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (
        points=[1,0; 18,0; 18,-31], style(color=74, rgbcolor={0,0,127}));
  end BranchingPipes2;
  
  model BranchingPipes3 
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-80,-4; -60,16]);
    
    ControlValves.ValveIncompressible valveIncompressible(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[10,36; 30,56]);
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[74,-20; 62,-8]);
    Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-40,36; -20,56]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-40,-50; -20,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (extent=[-20,70; 0,90]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 2,1;
          2,0; 100,0]) annotation (extent=[-20,-10; 0,10]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
    Junctions.Splitter splitter(redeclare package Medium = Medium) 
      annotation (extent=[-50,0; -36,12], rotation=90);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -80,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible1.port_b, sink.port) annotation (points=[28,-40;
          46,-40; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible.port_b, sink.port) annotation (points=[30,46;
          46,46; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (points=[-20,
          -40; 8,-40], style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valveIncompressible.port_a) annotation (points=[-20,
          46; 10,46], style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (
        points=[1,80; 20,80; 20,55], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (
        points=[1,0; 18,0; 18,-31], style(color=74, rgbcolor={0,0,127}));
    connect(pipe1.port_b, splitter.port_3) annotation (points=[-60,6; -50.7,6],
        style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, splitter.port_2) annotation (points=[-40,46; -43,46;
          -43,12.6], style(color=69, rgbcolor={0,127,255}));
    connect(splitter.port_1, pipe3.port_a) annotation (points=[-43,-0.6; -43,
          -40.3; -40,-40.3; -40,-40], style(color=69, rgbcolor={0,127,255}));
  end BranchingPipes3;
  
  model BranchingPipes4 
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-80,-4; -60,16]);
    
    ControlValves.ValveIncompressible valveIncompressible(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[10,36; 30,56]);
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Uses dynamic splitter. Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes. The simulation fails at t=0 due to lack of initialization of the splitter state variables.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[74,-20; 62,-8]);
    Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-40,36; -20,56]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (extent=[-40,-50; -20,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (extent=[-20,70; 0,90]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 2,1;
          2,0; 100,0]) annotation (extent=[-20,-10; 0,10]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
    Junctions.Junction_dynamic splitter(redeclare package Medium = Medium) 
      annotation (extent=[-50,0; -36,12], rotation=90);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -80,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible1.port_b, sink.port) annotation (points=[28,-40;
          46,-40; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible.port_b, sink.port) annotation (points=[30,46;
          46,46; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (points=[-20,
          -40; 8,-40], style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valveIncompressible.port_a) annotation (points=[-20,
          46; 10,46], style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (
        points=[1,80; 20,80; 20,55], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (
        points=[1,0; 18,0; 18,-31], style(color=74, rgbcolor={0,0,127}));
    connect(pipe1.port_b, splitter.port_3) annotation (points=[-60,6; -50.7,6],
        style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, splitter.port_2) annotation (points=[-40,46; -43,46;
          -43,12.6], style(color=69, rgbcolor={0,127,255}));
    connect(splitter.port_1, pipe3.port_a) annotation (points=[-43,-0.6; -43,
          -40.3; -40,-40.3; -40,-40], style(color=69, rgbcolor={0,127,255}));
  end BranchingPipes4;
  
  model SeriesPipes1 
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,-6; -88,6]);
    Pipes.LumpedPipe pipe1(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (extent=[-76,-10; -56,10]);
    
    ControlValves.ValveIncompressible valveIncompressible(
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=1000,
      redeclare package Medium = Medium) annotation (extent=[52,-10; 72,10]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with the valve open. At t=1, the valve is closed, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      T=300,
      p=1.0e5,
      redeclare package Medium = Medium) annotation (extent=[94,-6; 82,6]);
    Pipes.LumpedPipe pipe2(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (extent=[-14,-10; 6,10]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (extent=[-44,78; -24,98]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
    PressureLosses.SimpleGenericOrifice simpleGenericOrifice(
      zeta=0.4,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (extent=[-46,-10; -26,10]);
    Pipes.LumpedPipe pipe3(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (extent=[16,-10; 36,10]);
    
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,0; -76,0], style(
          color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible.port_b, sink.port) 
      annotation (points=[72,0; 82,0], style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (
        points=[-23,88; 62,88; 62,9], style(color=74, rgbcolor={0,0,127}));
    connect(pipe1.port_b, simpleGenericOrifice.port_a) annotation (points=[-56,
          0; -46,0], style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, simpleGenericOrifice.port_b) annotation (points=[-14,
          0; -26,0], style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, pipe3.port_a) 
      annotation (points=[6,0; 16,0], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valveIncompressible.port_a) 
      annotation (points=[36,0; 52,0], style(color=69, rgbcolor={0,127,255}));
  end SeriesPipes1;
  
  model SeriesPipes12 
    extends SeriesPipes1(
      pipe1(initType=Modelica_Fluid.Types.Init.SteadyState),
      pipe2(initType=Modelica_Fluid.Types.Init.SteadyState),
      pipe3(initType=Modelica_Fluid.Types.Init.SteadyState));
  equation 
    
    annotation (Documentation(info="<html>
Same as SeriesPipes1, but with steady-state initial conditions. Equal start attributes 
for pressures. Initialization fails.
</html>"));
  end SeriesPipes12;
  
  model SeriesPipes13 
    extends SeriesPipes1(
      pipe1(initType=Modelica_Fluid.Types.Init.SteadyState),
      pipe2(initType=Modelica_Fluid.Types.Init.SteadyState, p_start=4.95e5),
      pipe3(initType=Modelica_Fluid.Types.Init.SteadyState, p_start=4.9e5));
  equation 
    
    annotation (Documentation(info="<html>
Same as SeriesPipes1, but with steady-state initial conditions. Start attributes for 
pressure in order to get positive flow rates. Initialization succeeds, then the simulation
fails for zero flow rate.
</html>"));
  end SeriesPipes13;
  
  model IncompressibleFluidNetwork1 
    replaceable package Medium = 
        Modelica.Media.Incompressible.Examples.Essotherm650 
      extends Modelica.Media.Interfaces.PartialMedium;
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-98,4; -86,16]);
    Pipes.LumpedPipe pipe1(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (extent=[-78,0; -58,20]);
    
    Pipes.LumpedPipe pipe2(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      length=0.5)                        annotation (extent=[-60,26; -40,46],
        rotation=90);
    
    Pipes.LumpedPipe pipe3(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      length=0.5)                        annotation (extent=[-60,-26; -40,-6],
        rotation=-90);
    Pipes.LumpedPipe pipe4(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      length=2)                          annotation (extent=[-8,-46; 12,-26],
        rotation=0);
    Pipes.LumpedPipe pipe5(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      length=20)                         annotation (extent=[26,-60; 46,-40],
        rotation=0);
    Pipes.LumpedPipe pipe6(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      length=20)                         annotation (extent=[26,-32; 46,-12],
        rotation=0);
    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=3.0e4,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[-40,48; -26,64]);
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=3.0e4,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[-40,-28; -26,-44]);
    Pipes.LumpedPipe pipe7(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (extent=[-18,46; 2,66],
        rotation=0);
    ControlValves.ValveIncompressible valve3(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=3.0e4,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[62,2; 76,18]);
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) 
             annotation (extent=[98,4; 86,16]);
    inner Ambient ambient annotation (extent=[-98,80; -78,100]);
    Modelica.Blocks.Sources.Step valveOpening1(
      height=-0.2,
      offset=1,
      startTime=1) annotation (extent=[-66,80; -46,100]);
    Modelica.Blocks.Sources.Step valveOpening2(
      height=-0.2,
      offset=1,
      startTime=2) annotation (extent=[-82,-64; -62,-44]);
    Modelica.Blocks.Sources.Step valveOpening3(
      height=-0.2,
      offset=1,
      startTime=2) annotation (extent=[8,68; 28,88]);
    Pipes.LumpedPipe pipe8(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (extent=[-2,20; 18,40],
        rotation=-90);
    Pipes.LumpedPipe pipe9(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (extent=[12,46; 32,66],
        rotation=0);
    Pipes.LumpedPipe pipe10(
      p_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (extent=[18,-4; 38,16],
        rotation=0);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-86,10; -78,10], style(
          color=69, rgbcolor={0,127,255}));
    annotation (Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false));
    connect(pipe1.port_b, pipe3.port_a) annotation (points=[-58,10; -50,10; -50,
          -6], style(color=69, rgbcolor={0,127,255}));
    connect(pipe1.port_b, pipe2.port_a) annotation (points=[-58,10; -50,10; -50,
          26], style(color=69, rgbcolor={0,127,255}));
    connect(pipe4.port_b, pipe5.port_a) annotation (points=[12,-36; 18,-36; 18,
          -50; 26,-50], style(color=69, rgbcolor={0,127,255}));
    connect(pipe4.port_b, pipe6.port_a) annotation (points=[12,-36; 18,-36; 18,
          -22; 26,-22], style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valve1.port_a) annotation (points=[-50,46; -50,56;
          -40,56], style(color=69, rgbcolor={0,127,255}));
    connect(valve2.port_b, pipe4.port_a) annotation (points=[-26,-36; -8,-36],
        style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valve2.port_a) annotation (points=[-50,-26; -50,-36;
          -40,-36], style(color=69, rgbcolor={0,127,255}));
    connect(valve1.port_b, pipe7.port_a) annotation (points=[-26,56; -18,56],
        style(color=69, rgbcolor={0,127,255}));
    connect(pipe6.port_b, valve3.port_a) annotation (points=[46,-22; 54,-22; 54,
          10; 62,10], style(color=69, rgbcolor={0,127,255}));
    connect(pipe5.port_b, valve3.port_a) annotation (points=[46,-50; 54,-50; 54,
          10; 62,10], style(color=69, rgbcolor={0,127,255}));
    connect(valve3.port_b, sink.port) annotation (points=[76,10; 86,10], style(
          color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition) annotation (points=[-45,90;
          -33,90; -33,63.2], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition) annotation (points=[-61,-54;
          -33,-54; -33,-43.2], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening3.y, valve3.stemPosition) annotation (points=[29,78; 69,
          78; 69,17.2], style(color=74, rgbcolor={0,0,127}));
    connect(pipe7.port_b, pipe9.port_a) 
      annotation (points=[2,56; 12,56], style(color=69, rgbcolor={0,127,255}));
    connect(pipe7.port_b, pipe8.port_a) annotation (points=[2,56; 8,56; 8,40],
        style(color=69, rgbcolor={0,127,255}));
    connect(pipe9.port_b, valve3.port_a) annotation (points=[32,56; 62,56; 62,
          10], style(color=69, rgbcolor={0,127,255}));
    connect(pipe8.port_b, pipe10.port_a) annotation (points=[8,20; 8,6; 18,6],
        style(color=69, rgbcolor={0,127,255}));
    connect(pipe8.port_b, pipe4.port_b) annotation (points=[8,20; 8,-20; 12,-20;
          12,-36], style(color=69, rgbcolor={0,127,255}));
    connect(pipe10.port_b, valve3.port_a) annotation (points=[38,6; 50,6; 50,10;
          62,10], style(color=69, rgbcolor={0,127,255}));
  end IncompressibleFluidNetwork1;
  
  model BranchingPipes12 
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5) annotation (extent=[-78,-4; -58,16]);
    
    ControlValves.ValveIncompressible valveIncompressible(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[10,36; 30,56]);
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=1000) annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[74,-20; 62,-8]);
    Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (extent=[-40,36; -20,56]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (extent=[-40,-50; -20,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (extent=[-20,70; 0,90]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 2,1;
          2.01,1e-6; 100,0]) 
                       annotation (extent=[-20,-10; 0,10]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -78,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible1.port_b, sink.port) annotation (points=[28,-40;
          46,-40; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible.port_b, sink.port) annotation (points=[30,46;
          46,46; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (points=[-20,
          -40; 8,-40], style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valveIncompressible.port_a) annotation (points=[-20,
          46; 10,46], style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, pipe1.port_b) annotation (points=[-40,46; -46,46; -46,
          6; -58,6], style(color=69, rgbcolor={0,127,255}));
    connect(pipe1.port_b, pipe3.port_a) annotation (points=[-58,6; -46,6; -46,
          -40; -40,-40], style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (
        points=[1,80; 20,80; 20,55], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (
        points=[1,0; 18,0; 18,-31], style(color=74, rgbcolor={0,0,127}));
  end BranchingPipes12;
  
  model BranchingPipes13 
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5) annotation (extent=[-78,-4; -58,16]);
    
    ControlValves.ValveIncompressible valveIncompressible(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5)    annotation (extent=[10,36; 30,56]);
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5)    annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[74,-20; 62,-8]);
    Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (extent=[-34,36; -14,56]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (extent=[-30,-50; -10,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (extent=[-20,70; 0,90]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 2,1;
          2.01,1e-6; 100,0]) 
                       annotation (extent=[-20,-10; 0,10]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
    Junctions.JunctionIdeal junctionIdeal(redeclare package Medium = Medium,
        p_start=5.0e5) 
      annotation (extent=[-48,-4; -28,16], rotation=90);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -78,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible1.port_b, sink.port) annotation (points=[28,-40;
          46,-40; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible.port_b, sink.port) annotation (points=[30,46;
          46,46; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (points=[-10,-40;
          8,-40],      style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valveIncompressible.port_a) annotation (points=[-14,46;
          10,46],     style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (
        points=[1,80; 20,80; 20,55], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (
        points=[1,0; 18,0; 18,-31], style(color=74, rgbcolor={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (points=[-58,6; 
          -53.5,6; -53.5,6; -49,6],
                              style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (points=[-34,46; -38,
          46; -38,17], style(color=69, rgbcolor={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (points=[-38,-5; -38,
          -40; -30,-40], style(color=69, rgbcolor={0,127,255}));
  end BranchingPipes13;
  
  model BranchingPipes14 
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5) annotation (extent=[-78,-4; -58,16]);
    
    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5)    annotation (extent=[10,36; 30,56]);
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5)    annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[74,-20; 62,-8]);
    Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (extent=[-34,36; -14,56]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (extent=[-30,-50; -10,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          2,1e-2; 100,1e-2]) 
                       annotation (extent=[-20,72; 0,92]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 3,1;
          4,1e-2; 100,1e-2]) 
                       annotation (extent=[-20,-10; 0,10]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
    Junctions.JunctionVolume junctionIdeal(
                                          redeclare package Medium = Medium,
      V=1e-3,
      p_start=5.0e5,
      T_start=300,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
      annotation (extent=[-48,-4; -28,16], rotation=90);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -78,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(valve2.port_b, sink.port)               annotation (points=[28,-40;
          46,-40; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(valve1.port_b, sink.port)              annotation (points=[30,46;
          46,46; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valve2.port_a)               annotation (points=[-10,-40;
          8,-40],      style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valve1.port_a)              annotation (points=[-14,46;
          10,46],     style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition)              annotation (
        points=[1,82; 20,82; 20,55], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition)               annotation (
        points=[1,0; 18,0; 18,-31], style(color=74, rgbcolor={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (points=[-58,6;
          -53.5,6; -53.5,6; -49,6],
                              style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (points=[-34,46; -38,
          46; -38,17], style(color=69, rgbcolor={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (points=[-38,-5; -38,
          -40; -30,-40], style(color=69, rgbcolor={0,127,255}));
  end BranchingPipes14;
  
  model BranchingPipes15 
    replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    // replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5,
      dp_small=10)      annotation (extent=[-78,-4; -58,16]);
    
    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5)    annotation (extent=[10,36; 30,56]);
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5)    annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[74,-20; 62,-8]);
    Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (extent=[-34,36; -14,56]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (extent=[-30,-50; -10,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1; 2,
          1e-3; 100,1e-3]) 
                       annotation (extent=[-20,72; 0,92]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 3,1; 4,
          1e-3; 100,1e-3]) 
                       annotation (extent=[-20,-10; 0,10]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
    Junctions.JunctionVolume junctionIdeal(
                                          redeclare package Medium = Medium,
      V=1e-3,
      p_start=5.0e5,
      T_start=300,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
      annotation (extent=[-48,-4; -28,16], rotation=90);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -78,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(valve2.port_b, sink.port)               annotation (points=[28,-40;
          46,-40; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(valve1.port_b, sink.port)              annotation (points=[30,46;
          46,46; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valve2.port_a)               annotation (points=[-10,-40;
          8,-40],      style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valve1.port_a)              annotation (points=[-14,46;
          10,46],     style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition)              annotation (
        points=[1,82; 20,82; 20,55], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition)               annotation (
        points=[1,0; 18,0; 18,-31], style(color=74, rgbcolor={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (points=[-58,6;
          -53.5,6; -53.5,6; -49,6],
                              style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (points=[-34,46; -38,
          46; -38,17], style(color=69, rgbcolor={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (points=[-38,-5; -38,
          -40; -30,-40], style(color=69, rgbcolor={0,127,255}));
  end BranchingPipes15;
  
  model BranchingPipes16 
    replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    // replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5,
      dp_small=10)      annotation (extent=[-78,-4; -58,16]);
    
    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5)    annotation (extent=[10,36; 30,56]);
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5)    annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[74,-20; 62,-8]);
    Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (extent=[-34,36; -14,56]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (extent=[-30,-50; -10,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          2,0; 100,0]) annotation (extent=[-20,72; 0,92]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 3,1;
          4,0; 100,0]) annotation (extent=[-20,-12; 0,8]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
    Junctions.JunctionVolume junctionIdeal(
                                          redeclare package Medium = Medium,
      V=1e-3,
      p_start=5.0e5,
      T_start=300,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
      annotation (extent=[-48,-4; -28,16], rotation=90);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -78,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(valve2.port_b, sink.port)               annotation (points=[28,-40;
          46,-40; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(valve1.port_b, sink.port)              annotation (points=[30,46;
          46,46; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valve2.port_a)               annotation (points=[-10,-40;
          8,-40],      style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valve1.port_a)              annotation (points=[-14,46;
          10,46],     style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition)              annotation (
        points=[1,82; 20,82; 20,55], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition)               annotation (
        points=[1,-2; 18,-2; 18,-31],
                                    style(color=74, rgbcolor={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (points=[-58,6;
          -53.5,6; -53.5,6; -49,6],
                              style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (points=[-34,46; -38,
          46; -38,17], style(color=69, rgbcolor={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (points=[-38,-5; -38,
          -40; -30,-40], style(color=69, rgbcolor={0,127,255}));
  end BranchingPipes16;
  
  model BranchingPipes17 
    replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    // replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5,
      dp_small=10)      annotation (extent=[-78,-4; -58,16]);
    
    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5, 
      dp(start=10)) 
                  annotation (extent=[10,36; 30,56]);
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5)    annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5, Tolerance=1e-007),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[94,-18; 82,-6]);
    Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (extent=[-34,36; -14,56]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (extent=[-30,-50; -10,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          2,0; 100,0]) annotation (extent=[-20,72; 0,92]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 3,1;
          4,0; 100,0]) annotation (extent=[-18,-12; 2,8]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
    Junctions.JunctionVolume junctionIdeal(
                                          redeclare package Medium = Medium,
      V=1e-3,
      p_start=5.0e5,
      T_start=300,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
      annotation (extent=[-50,-4; -30,16], rotation=90);
    Junctions.JunctionVolume junctionVolume(redeclare package Medium = Medium,
        V=1e-3) annotation (extent=[66,-22; 46,-2], rotation=90);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -78,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valve2.port_a)               annotation (points=[-10,-40;
          8,-40],      style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valve1.port_a)              annotation (points=[-14,46;
          10,46],     style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition)              annotation (
        points=[1,82; 20,82; 20,55], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition)               annotation (
        points=[3,-2; 18,-2; 18,-31],
                                    style(color=74, rgbcolor={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (points=[-58,6; 
          -54.5,6; -54.5,6; -51,6],
                              style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (points=[-34,46; -40,
          46; -40,17], style(color=69, rgbcolor={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (points=[-40,-5; -40,
          -40; -30,-40], style(color=69, rgbcolor={0,127,255}));
    connect(junctionVolume.port_3, sink.port) annotation (points=[67,-12; 82,
          -12], style(color=69, rgbcolor={0,127,255}));
    connect(valve2.port_b, junctionVolume.port_1) annotation (points=[28,-40;
          56,-40; 56,-23], style(color=69, rgbcolor={0,127,255}));
    connect(valve1.port_b, junctionVolume.port_2) annotation (points=[30,46; 56,
          46; 56,-1], style(color=69, rgbcolor={0,127,255}));
  end BranchingPipes17;
  
  model BranchingPipes18 
    // replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5,
      dp_small=10)      annotation (extent=[-78,-4; -58,16]);
    
    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5)    annotation (extent=[10,36; 30,56]);
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5)    annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[94,-18; 82,-6]);
    Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (extent=[-34,36; -14,56]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (extent=[-30,-50; -10,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1; 2,
          0; 100,0])   annotation (extent=[-20,72; 0,92]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 3,1; 4,
          0; 100,0])   annotation (extent=[-18,-12; 2,8]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
    Junctions.JunctionVolume junctionIdeal(
                                          redeclare package Medium = Medium,
      V=1e-3,
      p_start=5.0e5,
      T_start=300,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
      annotation (extent=[-48,-4; -28,16], rotation=90);
    Junctions.JunctionVolume junctionVolume(redeclare package Medium = Medium, V=
          1e-3) annotation (extent=[66,-22; 46,-2], rotation=90);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -78,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valve2.port_a)               annotation (points=[-10,-40;
          8,-40],      style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_b, valve1.port_a)              annotation (points=[-14,46;
          10,46],     style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition)              annotation (
        points=[1,82; 20,82; 20,55], style(color=74, rgbcolor={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition)               annotation (
        points=[3,-2; 18,-2; 18,-31],
                                    style(color=74, rgbcolor={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (points=[-58,6; 
          -53.5,6; -53.5,6; -49,6],
                              style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (points=[-34,46; -38,
          46; -38,17], style(color=69, rgbcolor={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (points=[-38,-5; -38,
          -40; -30,-40], style(color=69, rgbcolor={0,127,255}));
    connect(junctionVolume.port_3, sink.port) annotation (points=[67,-12; 82,
          -12], style(color=69, rgbcolor={0,127,255}));
    connect(valve2.port_b, junctionVolume.port_1) annotation (points=[28,-40;
          56,-40; 56,-23], style(color=69, rgbcolor={0,127,255}));
    connect(valve1.port_b, junctionVolume.port_2) annotation (points=[30,46; 56,
          46; 56,-1], style(color=69, rgbcolor={0,127,255}));
  end BranchingPipes18;

  model BranchingPipes131 
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (extent=[-100,0; -88,12]);
    Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5, 
      frictionAndGravity2(dp(start=1000)), 
      frictionAndGravity1(dp(start=1000))) 
                        annotation (extent=[-78,-4; -58,16]);
    
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=5.0e5,
      dp_nom=4.0e5,
      m_flow_nom=1,
      d_nom=5, 
      dp(start=4.0e5)) 
                  annotation (extent=[8,-50; 28,-30]);
    annotation (
      Diagram,
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (extent=[74,-20; 62,-8]);
    
    Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5, 
      frictionAndGravity1(dp(start=1000)), 
      frictionAndGravity2(dp(start=1000))) 
                        annotation (extent=[-30,-50; -10,-30]);
    
    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (extent=[-20,70; 0,90]);
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 2,1; 
          2.01,1e-6; 100,0]) 
                       annotation (extent=[-20,-8; 0,12]);
    inner Ambient ambient annotation (extent=[-100,60; -80,80]);
  equation 
    connect(source.port, pipe1.port_a) annotation (points=[-88,6; -78,6], style(
          color=69, rgbcolor={0,127,255}));
    connect(valveIncompressible1.port_b, sink.port) annotation (points=[28,-40;
          46,-40; 46,-14; 62,-14], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (points=[-10,-40;
          8,-40],      style(color=69, rgbcolor={0,127,255}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (
        points=[1,2; 18,2; 18,-31], style(color=74, rgbcolor={0,0,127}));
    connect(pipe1.port_b, pipe3.port_a) annotation (points=[-58,6; -44,6; -44,
          -40; -30,-40], style(color=69, rgbcolor={0,127,255}));
  end BranchingPipes131;
end TestCriticalCases;
