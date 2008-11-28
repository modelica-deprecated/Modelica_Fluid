within Modelica_Fluid.Test;
package TestCriticalCases
  "Collection of test cases which might be critical for the solvers"
  model IdealMixing1 "Test properties of ideal mixing"
    // package Medium =  Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater;
    // Modelica.Media.IdealGases.MixtureGases.FlueGasSixComponents,package Medium = Modelica.Media.Air.DryAirNasa;
    package Medium = 
        Modelica.Media.IdealGases.MixtureGases.FlueGasSixComponents;

    PressureLosses.WallFrictionAndGravity pipeFriction1(
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
      length=1,
      diameter=0.2,
      redeclare package Medium = Medium) annotation (Placement(transformation(
            extent={{-32,-40},{-12,-20}}, rotation=0)));
    PressureLosses.WallFrictionAndGravity pipeFriction2(
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
      length=1,
      diameter=0.2,
      redeclare package Medium = Medium) annotation (Placement(transformation(
            extent={{12,-40},{32,-20}}, rotation=0)));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics),
      experiment(StopTime=10),
      experimentSetupOutput);
    PressureLosses.WallFrictionAndGravity pipeFriction3(
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
      length=1,
      diameter=0.2,
      redeclare package Medium = Medium) 
      annotation (Placement(transformation(
          origin={0,0},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Sources.PrescribedBoundary_pTX boundary1(
      usePressureInput=true,
      useTemperatureInput=true,
      redeclare package Medium = Medium) annotation (Placement(transformation(
            extent={{-68,-40},{-48,-20}}, rotation=0)));
    Sources.PrescribedBoundary_pTX boundary2(
      usePressureInput=false,
      useTemperatureInput=false,
      p=101000,
      T=320,
      redeclare package Medium = Medium) annotation (Placement(transformation(
            extent={{66,-40},{46,-20}}, rotation=0)));
    Sources.PrescribedBoundary_pTX boundary3(
      usePressureInput=true,
      useTemperatureInput=false,
      T=340,
      redeclare package Medium = Medium) 
      annotation (Placement(transformation(
          origin={0,30},
          extent={{-10,10},{10,-10}},
          rotation=270)));
    Modelica.Blocks.Sources.Sine sine1(
      amplitude=0.05e5,
      freqHz=2,
      offset=1e5,
      phase=0.013962634015955) annotation (Placement(transformation(extent={{
              -100,-20},{-80,0}}, rotation=0)));
    Modelica.Blocks.Sources.Sine sine2(
      amplitude=10,
      freqHz=1,
      phase=0.0017453292519943,
      offset=300) annotation (Placement(transformation(extent={{-100,-58},{-80,
              -38}}, rotation=0)));
    Modelica.Blocks.Sources.Sine sine3(
      amplitude=0.05e5,
      freqHz=2,
      offset=1e5) annotation (Placement(transformation(
          origin={0,70},
          extent={{10,-10},{-10,10}},
          rotation=90)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-88,60},
              {-68,80}}, rotation=0)));
    Sensors.TemperatureOnePort temperature(redeclare package Medium = Medium) 
      annotation (Placement(transformation(extent={{-10,-60},{10,-80}},
            rotation=0)));
  equation
    connect(pipeFriction1.port_b, pipeFriction2.port_a) annotation (Line(
        points={{-12,-30},{12,-30}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipeFriction3.port_a, pipeFriction1.port_b) annotation (Line(
        points={{-6.12323e-016,-10},{0,-10},{0,-30},{-12,-30}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(boundary1.ports[1], pipeFriction1.port_a) 
                                                  annotation (Line(
        points={{-48,-30},{-32,-30}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(boundary2.ports[1], pipeFriction2.port_b) 
                                                  annotation (Line(
        points={{46,-30},{32,-30}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(boundary3.ports[1], pipeFriction3.port_b) 
                                                  annotation (Line(
        points={{-1.83697e-015,20},{-1.83697e-015,10},{6.12323e-016,10}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(sine1.y, boundary1.p_in) annotation (Line(
        points={{-79,-10},{-76,-10},{-76,-24},{-70,-24}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sine2.y, boundary1.T_in) annotation (Line(
        points={{-79,-48},{-76,-48},{-76,-30},{-70,-30}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(sine3.y, boundary3.p_in) annotation (Line(
        points={{-6.73556e-016,59},{-6.73556e-016,50.5},{-6,50.5},{-6,42}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(temperature.port, pipeFriction3.port_a) annotation (Line(
        points={{0,-60},{0,-10},{-6.12323e-016,-10}},
        color={0,127,255},
        smooth=Smooth.None));
  end IdealMixing1;

  model BranchingPipes1
    //replaceable package Medium = Modelica.Media.Water.StandardWater;
    replaceable package Medium = 
        Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=200000,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      initType=Modelica_Fluid.Types.Init.SteadyStateHydraulic,
      p_a_start=200000,
      p_b_start=100000,
      T_start=300)      annotation (Placement(transformation(extent={{-72,-4},{
              -52,16}}, rotation=0)));

    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      m_flow_nominal=1,
      d_nominal=1000,
      dp_nominal=200000) 
                  annotation (Placement(transformation(extent={{10,36},{30,56}},
            rotation=0)));
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      m_flow_nominal=1,
      d_nominal=1000,
      dp_nominal=200000) 
                  annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      nPorts=2,
      p=100000,
      T=300)   annotation (Placement(transformation(extent={{74,-20},{62,-8}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      initType=Modelica_Fluid.Types.Init.SteadyStateHydraulic,
      p_a_start=200000,
      p_b_start=200000,
      T_start=300)      annotation (Placement(transformation(extent={{-40,36},{
              -20,56}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      initType=Modelica_Fluid.Types.Init.SteadyStateHydraulic,
      p_a_start=200000,
      p_b_start=200000,
      T_start=300)      annotation (Placement(transformation(extent={{-40,-50},
              {-20,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (Placement(transformation(extent={{-20,70},{0,
              90}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 2,1;
          2.01,0; 100,0]) 
                       annotation (Placement(transformation(extent={{-20,-10},{
              0,10}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-72,6}},
          color={0,127,255}));
    connect(valve2.port_b, sink.ports[2])               annotation (Line(points={{28,-40},
            {46,-40},{46,-15.2},{62,-15.2}},     color={0,127,255}));
    connect(valve1.port_b, sink.ports[1])              annotation (Line(points={{30,46},
            {46,46},{46,-12.8},{62,-12.8}}, color={0,127,255}));
    connect(pipe3.port_b, valve2.port_a)               annotation (Line(points=
            {{-20,-40},{8,-40}}, color={0,127,255}));
    connect(pipe2.port_b, valve1.port_a)              annotation (Line(points={
            {-20,46},{10,46}}, color={0,127,255}));
    connect(pipe2.port_a, pipe1.port_b) annotation (Line(points={{-40,46},{-46,
            46},{-46,6},{-52,6}}, color={0,127,255}));
    connect(pipe1.port_b, pipe3.port_a) annotation (Line(points={{-52,6},{-46,6},
            {-46,-40},{-40,-40}}, color={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition)              annotation (Line(
          points={{1,80},{20,80},{20,54}}, color={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition)               annotation (Line(
          points={{1,0},{18,0},{18,-32}}, color={0,0,127}));
  end BranchingPipes1;

  model BranchingPipes2
    replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=5.0e5,
      p_b_start=5.0e5)  annotation (Placement(transformation(extent={{-72,-4},{
              -52,16}}, rotation=0)));

    ControlValves.ValveIncompressible valveIncompressible(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=1000) annotation (Placement(transformation(extent={{10,36},{30,56}},
            rotation=0)));
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=1000) annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      nPorts=2,
      p=100000,
      T=300)   annotation (Placement(transformation(extent={{74,-20},{62,-8}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=5.0e5,
      p_b_start=5.0e5)  annotation (Placement(transformation(extent={{-40,36},{
              -20,56}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=5.0e5,
      p_b_start=5.0e5)  annotation (Placement(transformation(extent={{-40,-50},
              {-20,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,0; 1,0;
          1,1; 100,1]) annotation (Placement(transformation(extent={{-20,70},{0,
              90}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,0; 2,0;
          2,1; 100,1]) annotation (Placement(transformation(extent={{-20,-10},{
              0,10}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-72,6}},
          color={0,127,255}));
    connect(valveIncompressible1.port_b, sink.ports[2]) annotation (Line(points={{28,-40},
            {46,-40},{46,-15.2},{62,-15.2}},     color={0,127,255}));
    connect(valveIncompressible.port_b, sink.ports[1]) annotation (Line(points={{30,46},
            {46,46},{46,-12.8},{62,-12.8}}, color={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (Line(points=
            {{-20,-40},{8,-40}}, color={0,127,255}));
    connect(pipe2.port_b, valveIncompressible.port_a) annotation (Line(points={
            {-20,46},{10,46}}, color={0,127,255}));
    connect(pipe2.port_a, pipe1.port_b) annotation (Line(points={{-40,46},{-46,
            46},{-46,6},{-52,6}}, color={0,127,255}));
    connect(pipe1.port_b, pipe3.port_a) annotation (Line(points={{-52,6},{-46,6},
            {-46,-40},{-40,-40}}, color={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (Line(
          points={{1,80},{20,80},{20,54}}, color={0,0,127}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (Line(
          points={{1,0},{18,0},{18,-32}}, color={0,0,127}));
  end BranchingPipes2;

  model BranchingPipes3
    replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      p_b_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (Placement(transformation(extent={{-80,-4},{
              -60,16}}, rotation=0)));

    ControlValves.ValveIncompressible valveIncompressible(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=1000) annotation (Placement(transformation(extent={{10,36},{30,56}},
            rotation=0)));
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=1000) annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      nPorts=2,
      p=100000,
      T=300)   annotation (Placement(transformation(extent={{74,-20},{62,-8}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      p_b_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (Placement(transformation(extent={{-40,36},{
              -20,56}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      p_b_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2) annotation (Placement(transformation(extent={{-40,-50},
              {-20,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,0; 1,0;
          1,1; 100,1]) annotation (Placement(transformation(extent={{-20,70},{0,
              90}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,0; 2,0;
          2,1; 100,1]) annotation (Placement(transformation(extent={{-20,-10},{
              0,10}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Junctions.JunctionIdeal splitter(redeclare package Medium = Medium) 
      annotation (Placement(transformation(
          origin={-43,6},
          extent={{-6,-7},{6,7}},
          rotation=90)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-80,6}},
          color={0,127,255}));
    connect(valveIncompressible1.port_b, sink.ports[2]) annotation (Line(points={{28,-40},
            {46,-40},{46,-15.2},{62,-15.2}},     color={0,127,255}));
    connect(valveIncompressible.port_b, sink.ports[1]) annotation (Line(points={{30,46},
            {46,46},{46,-12.8},{62,-12.8}}, color={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (Line(points=
            {{-20,-40},{8,-40}}, color={0,127,255}));
    connect(pipe2.port_b, valveIncompressible.port_a) annotation (Line(points={
            {-20,46},{10,46}}, color={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (Line(
          points={{1,80},{20,80},{20,54}}, color={0,0,127}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (Line(
          points={{1,0},{18,0},{18,-32}}, color={0,0,127}));
    connect(pipe1.port_b, splitter.port_3) annotation (Line(points={{-60,6},{
            -55,6},{-55,6},{-50,6}},
                     color={0,127,255}));
    connect(pipe2.port_a, splitter.port_2) annotation (Line(points={{-40,46},{
            -43,46},{-43,12}}, color={0,127,255}));
    connect(splitter.port_1, pipe3.port_a) annotation (Line(points={{-43,0},{
            -43,-40.3},{-40,-40.3},{-40,-40}}, color={0,127,255}));
  end BranchingPipes3;

  model BranchingPipes4
    replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=5.0e5,
      p_b_start=5.0e5)  annotation (Placement(transformation(extent={{-80,-4},{
              -60,16}}, rotation=0)));

    ControlValves.ValveIncompressible valveIncompressible(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=1000) annotation (Placement(transformation(extent={{10,36},{30,56}},
            rotation=0)));
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=1000) annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Uses dynamic splitter. Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes. The simulation fails at t=0 due to lack of initialization of the splitter state variables.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      nPorts=2,
      p=100000,
      T=300)   annotation (Placement(transformation(extent={{74,-20},{62,-8}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=5.0e5,
      p_b_start=5.0e5)  annotation (Placement(transformation(extent={{-40,36},{
              -20,56}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=5.0e5,
      p_b_start=5.0e5)  annotation (Placement(transformation(extent={{-40,-50},
              {-20,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,0; 1,0;
          1,1; 100,1]) annotation (Placement(transformation(extent={{-20,70},{0,
              90}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,0; 2,0;
          2,1; 100,1]) annotation (Placement(transformation(extent={{-20,-10},{
              0,10}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Junctions.JunctionVolume splitter(redeclare package Medium = Medium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V=0.0002,
      p_start=500000) 
      annotation (Placement(transformation(
          origin={-43,6},
          extent={{-6,-7},{6,7}},
          rotation=90)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-80,6}},
          color={0,127,255}));
    connect(valveIncompressible1.port_b, sink.ports[2]) annotation (Line(points={{28,-40},
            {46,-40},{46,-15.2},{62,-15.2}},     color={0,127,255}));
    connect(valveIncompressible.port_b, sink.ports[1]) annotation (Line(points={{30,46},
            {46,46},{46,-12.8},{62,-12.8}}, color={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (Line(points=
            {{-20,-40},{8,-40}}, color={0,127,255}));
    connect(pipe2.port_b, valveIncompressible.port_a) annotation (Line(points={
            {-20,46},{10,46}}, color={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (Line(
          points={{1,80},{20,80},{20,54}}, color={0,0,127}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (Line(
          points={{1,0},{18,0},{18,-32}}, color={0,0,127}));
    connect(pipe1.port_b, splitter.port_3) annotation (Line(points={{-60,6},{
            -55,6},{-55,6},{-50,6}},
                     color={0,127,255}));
    connect(pipe2.port_a, splitter.port_2) annotation (Line(points={{-40,46},{
            -43,46},{-43,12}}, color={0,127,255}));
    connect(splitter.port_1, pipe3.port_a) annotation (Line(points={{-43,0},{
            -43,-40.3},{-40,-40.3},{-40,-40}}, color={0,127,255}));
  end BranchingPipes4;

  model SeriesPipes1
    replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=500000,
      T=300) annotation (Placement(transformation(extent={{-100,-6},{-88,6}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      p_a_start=500000,
      p_b_start=500000)                  annotation (Placement(transformation(
            extent={{-76,-10},{-56,10}}, rotation=0)));

    ControlValves.ValveIncompressible valveIncompressible(
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      m_flow_nominal=1,
      d_nominal=1000,
      redeclare package Medium = Medium,
      dp_nominal=400000,
      minStemPosition=0.01)              annotation (Placement(transformation(
            extent={{52,-10},{72,10}}, rotation=0)));
    annotation (
      Diagram(graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with the valve open. At t=1, the valve is closed, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      p=100000,
      T=300)                             annotation (Placement(transformation(
            extent={{94,-6},{82,6}}, rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      p_a_start=500000,
      p_b_start=500000)                  annotation (Placement(transformation(
            extent={{-14,-10},{6,10}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (Placement(transformation(extent={{-20,70},{0,
              90}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    PressureLosses.SimpleGenericOrifice simpleGenericOrifice(
      zeta=0.4,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (Placement(transformation(
            extent={{-46,-10},{-26,10}}, rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      p_a_start=500000,
      p_b_start=500000)                  annotation (Placement(transformation(
            extent={{16,-10},{36,10}}, rotation=0)));

  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,0},{-76,0}},
          color={0,127,255}));
    connect(valveIncompressible.port_b, sink.ports[1]) 
      annotation (Line(points={{72,0},{82,0}}, color={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (Line(
          points={{1,80},{62,80},{62,9}}, color={0,0,127}));
    connect(pipe1.port_b, simpleGenericOrifice.port_a) annotation (Line(points=
            {{-56,0},{-46,0}}, color={0,127,255}));
    connect(pipe2.port_a, simpleGenericOrifice.port_b) annotation (Line(points=
            {{-14,0},{-26,0}}, color={0,127,255}));
    connect(pipe2.port_b, pipe3.port_a) 
      annotation (Line(points={{6,0},{16,0}}, color={0,127,255}));
    connect(pipe3.port_b, valveIncompressible.port_a) 
      annotation (Line(points={{36,0},{52,0}}, color={0,127,255}));
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
      pipe2(initType=Modelica_Fluid.Types.Init.SteadyState, p_a_start=4.95e5, p_b_start=4.95e5),
      pipe3(initType=Modelica_Fluid.Types.Init.SteadyState, p_a_start=4.9e5, p_b_start=4.9e5));
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
      constrainedby Modelica.Media.Interfaces.PartialMedium;
    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-98,4},{-86,16}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (Placement(transformation(
            extent={{-78,0},{-58,20}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      length=0.5)                        annotation (Placement(transformation(
          origin={-50,36},
          extent={{-10,-10},{10,10}},
          rotation=90)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      length=0.5)                        annotation (Placement(transformation(
          origin={-50,-16},
          extent={{-10,-10},{10,10}},
          rotation=270)));
    Modelica_Fluid.Pipes.LumpedPipe pipe4(
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      length=2)                          annotation (Placement(transformation(
            extent={{-16,-46},{4,-26}}, rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe6(
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      length=20)                         annotation (Placement(transformation(
            extent={{26,-32},{46,-12}}, rotation=0)));
    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      m_flow_nominal=1,
      d_nominal=1000,
      dp_nominal=30000,
      minStemPosition=0.01) 
                  annotation (Placement(transformation(extent={{-40,48},{-26,64}},
            rotation=0)));
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      m_flow_nominal=1,
      d_nominal=1000,
      dp_nominal=30000,
      minStemPosition=0.01) 
                  annotation (Placement(transformation(extent={{-40,-28},{-26,
              -44}}, rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe7(
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (Placement(transformation(
            extent={{-18,46},{2,66}}, rotation=0)));
    ControlValves.ValveIncompressible valve3(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      m_flow_nominal=1,
      d_nominal=1000,
      dp_nominal=30000,
      minStemPosition=0.01) 
                  annotation (Placement(transformation(extent={{62,2},{76,18}},
            rotation=0)));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) 
             annotation (Placement(transformation(extent={{98,4},{86,16}},
            rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-98,80},
              {-78,100}}, rotation=0)));
    Modelica.Blocks.Sources.Step valveOpening1(
      height=-0.2,
      offset=1,
      startTime=1) annotation (Placement(transformation(extent={{-66,80},{-46,
              100}}, rotation=0)));
    Modelica.Blocks.Sources.Step valveOpening2(
      height=-0.2,
      offset=1,
      startTime=2) annotation (Placement(transformation(extent={{-82,-64},{-62,
              -44}}, rotation=0)));
    Modelica.Blocks.Sources.Step valveOpening3(
      height=-0.2,
      offset=1,
      startTime=2) annotation (Placement(transformation(extent={{8,68},{28,88}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe8(
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (Placement(transformation(
          origin={10,30},
          extent={{-10,-10},{10,10}},
          rotation=270)));
    Modelica_Fluid.Pipes.LumpedPipe pipe9(
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (Placement(transformation(
            extent={{16,46},{36,66}}, rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe10(
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.5e-2,
      redeclare package Medium = Medium) annotation (Placement(transformation(
            extent={{18,-4},{38,16}}, rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe5(
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      diameter=2.5e-2,
      redeclare package Medium = Medium,
      length=20)                         annotation (Placement(transformation(
            extent={{24,-60},{44,-40}}, rotation=0)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-86,10},{-78,
            10}}, color={0,127,255}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{
              -100,-100},{100,100}}),
                        graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false));
    connect(pipe1.port_b, pipe3.port_a) annotation (Line(points={{-58,10},{-50,
            10},{-50,-6}}, color={0,127,255}));
    connect(pipe1.port_b, pipe2.port_a) annotation (Line(points={{-58,10},{-50,
            10},{-50,26}}, color={0,127,255}));
    connect(pipe2.port_b, valve1.port_a) annotation (Line(points={{-50,46},{-50,
            56},{-40,56}}, color={0,127,255}));
    connect(valve2.port_b, pipe4.port_a) annotation (Line(points={{-26,-36},{
            -16,-36}},
                   color={0,127,255}));
    connect(pipe3.port_b, valve2.port_a) annotation (Line(points={{-50,-26},{
            -50,-36},{-40,-36}}, color={0,127,255}));
    connect(valve1.port_b, pipe7.port_a) annotation (Line(points={{-26,56},{-18,
            56}}, color={0,127,255}));
    connect(pipe6.port_b, valve3.port_a) annotation (Line(points={{46,-22},{54,
            -22},{54,10},{62,10}}, color={0,127,255}));
    connect(valve3.port_b, sink.ports[1]) annotation (Line(points={{76,10},{86,10}},
          color={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition) annotation (Line(points={{-45,
            90},{-33,90},{-33,63.2}}, color={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition) annotation (Line(points={{-61,
            -54},{-33,-54},{-33,-43.2}}, color={0,0,127}));
    connect(valveOpening3.y, valve3.stemPosition) annotation (Line(points={{29,
            78},{69,78},{69,17.2}}, color={0,0,127}));
    connect(pipe7.port_b, pipe9.port_a) 
      annotation (Line(points={{2,56},{16,56}}, color={0,127,255}));
    connect(pipe7.port_b, pipe8.port_a) annotation (Line(points={{2,56},{10,56},
            {10,40}},color={0,127,255}));
    connect(pipe9.port_b, valve3.port_a) annotation (Line(points={{36,56},{54,
            56},{54,10},{62,10}},
                          color={0,127,255}));
    connect(pipe8.port_b, pipe10.port_a) annotation (Line(points={{10,20},{10,6},
            {18,6}}, color={0,127,255}));
    connect(pipe10.port_b, valve3.port_a) annotation (Line(points={{38,6},{50,6},
            {50,10},{62,10}}, color={0,127,255}));
    connect(pipe4.port_b, pipe6.port_a) annotation (Line(
        points={{4,-36},{10,-36},{10,-22},{26,-22}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipe8.port_b, pipe4.port_b) annotation (Line(
        points={{10,20},{10,-36},{4,-36}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipe5.port_a, pipe4.port_b) annotation (Line(
        points={{24,-50},{10,-50},{10,-36},{4,-36}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipe5.port_b, valve3.port_a) annotation (Line(
        points={{44,-50},{54,-50},{54,10},{62,10}},
        color={0,127,255},
        smooth=Smooth.None));
  end IncompressibleFluidNetwork1;

  model BranchingPipes12
    replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5) annotation (Placement(transformation(extent={{-78,-4},{
              -58,16}}, rotation=0)));

    ControlValves.ValveIncompressible valveIncompressible(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=1000) annotation (Placement(transformation(extent={{10,36},{30,56}},
            rotation=0)));
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=1000) annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      nPorts=2,
      p=100000,
      T=300)   annotation (Placement(transformation(extent={{74,-20},{62,-8}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (Placement(transformation(extent={{-40,36},{
              -20,56}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (Placement(transformation(extent={{-40,-50},
              {-20,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (Placement(transformation(extent={{-20,70},{0,
              90}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 2,1;
          2.01,1e-6; 100,0]) 
                       annotation (Placement(transformation(extent={{-20,-10},{
              0,10}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-78,6}},
          color={0,127,255}));
    connect(valveIncompressible1.port_b, sink.ports[2]) annotation (Line(points={{28,-40},
            {46,-40},{46,-15.2},{62,-15.2}},     color={0,127,255}));
    connect(valveIncompressible.port_b, sink.ports[1]) annotation (Line(points={{30,46},
            {46,46},{46,-12.8},{62,-12.8}}, color={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (Line(points=
            {{-20,-40},{8,-40}}, color={0,127,255}));
    connect(pipe2.port_b, valveIncompressible.port_a) annotation (Line(points={
            {-20,46},{10,46}}, color={0,127,255}));
    connect(pipe2.port_a, pipe1.port_b) annotation (Line(points={{-40,46},{-46,
            46},{-46,6},{-58,6}}, color={0,127,255}));
    connect(pipe1.port_b, pipe3.port_a) annotation (Line(points={{-58,6},{-46,6},
            {-46,-40},{-40,-40}}, color={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (Line(
          points={{1,80},{20,80},{20,54}}, color={0,0,127}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (Line(
          points={{1,0},{18,0},{18,-32}}, color={0,0,127}));
  end BranchingPipes12;

  model BranchingPipes13
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5) annotation (Placement(transformation(extent={{-78,-4},{
              -58,16}}, rotation=0)));

    ControlValves.ValveIncompressible valveIncompressible(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      m_flow_nominal=1,
      d_nominal=5,
      dp_nominal=400000,
      minStemPosition=0.001) 
                  annotation (Placement(transformation(extent={{10,36},{30,56}},
            rotation=0)));
    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      m_flow_nominal=1,
      d_nominal=5,
      dp_nominal=400000,
      minStemPosition=0.001) 
                  annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      nPorts=2,
      p=100000,
      T=300)   annotation (Placement(transformation(extent={{74,-20},{62,-8}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (Placement(transformation(extent={{-34,36},{
              -14,56}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (Placement(transformation(extent={{-30,-50},
              {-10,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (Placement(transformation(extent={{-20,70},{0,
              90}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 2,1;
          2.01,1e-6; 100,0]) 
                       annotation (Placement(transformation(extent={{-20,-10},{
              0,10}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Junctions.JunctionIdeal junctionIdeal(redeclare package Medium = Medium) 
      annotation (Placement(transformation(
          origin={-38,6},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-78,6}},
          color={0,127,255}));
    connect(valveIncompressible1.port_b, sink.ports[2]) annotation (Line(points={{28,-40},
            {46,-40},{46,-15.2},{62,-15.2}},     color={0,127,255}));
    connect(valveIncompressible.port_b, sink.ports[1]) annotation (Line(points={{30,46},
            {46,46},{46,-12.8},{62,-12.8}}, color={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (Line(points=
            {{-10,-40},{8,-40}}, color={0,127,255}));
    connect(pipe2.port_b, valveIncompressible.port_a) annotation (Line(points={
            {-14,46},{10,46}}, color={0,127,255}));
    connect(valveOpening1.y, valveIncompressible.stemPosition) annotation (Line(
          points={{1,80},{20,80},{20,54}}, color={0,0,127}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (Line(
          points={{1,0},{18,0},{18,-32}}, color={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (Line(points={{-58,6},
            {-53.5,6},{-53.5,6},{-48,6}}, color={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (Line(points={{-34,
            46},{-38,46},{-38,16}}, color={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (Line(points={{-38,
            -4},{-38,-40},{-30,-40}}, color={0,127,255}));
  end BranchingPipes13;

  model BranchingPipes14
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5) annotation (Placement(transformation(extent={{-78,-4},{
              -58,16}}, rotation=0)));

    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
       dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=5)    annotation (Placement(transformation(extent={{10,36},{30,56}},
            rotation=0)));
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=5)    annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      nPorts=2,
      p=100000,
      T=300)   annotation (Placement(transformation(extent={{74,-20},{62,-8}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (Placement(transformation(extent={{-34,36},{
              -14,56}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (Placement(transformation(extent={{-30,-50},
              {-10,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          2,1e-2; 100,1e-2]) 
                       annotation (Placement(transformation(extent={{-20,72},{0,
              92}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 3,1;
          4,1e-2; 100,1e-2]) 
                       annotation (Placement(transformation(extent={{-20,-10},{
              0,10}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Junctions.JunctionVolume junctionIdeal(
                                          redeclare package Medium = Medium,
      V=1e-3,
      p_start=5.0e5,
      T_start=300,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
      annotation (Placement(transformation(
          origin={-38,6},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-78,6}},
          color={0,127,255}));
    connect(valve2.port_b, sink.ports[2])               annotation (Line(points={{28,-40},
            {46,-40},{46,-15.2},{62,-15.2}},     color={0,127,255}));
    connect(valve1.port_b, sink.ports[1])              annotation (Line(points={{30,46},
            {46,46},{46,-12.8},{62,-12.8}}, color={0,127,255}));
    connect(pipe3.port_b, valve2.port_a)               annotation (Line(points=
            {{-10,-40},{8,-40}}, color={0,127,255}));
    connect(pipe2.port_b, valve1.port_a)              annotation (Line(points={
            {-14,46},{10,46}}, color={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition)              annotation (Line(
          points={{1,82},{20,82},{20,54}}, color={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition)               annotation (Line(
          points={{1,0},{18,0},{18,-32}}, color={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (Line(points={{-58,6},
            {-53.5,6},{-53.5,6},{-48,6}}, color={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (Line(points={{-34,
            46},{-38,46},{-38,16}}, color={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (Line(points={{-38,
            -4},{-38,-40},{-30,-40}}, color={0,127,255}));
  end BranchingPipes14;

  model BranchingPipes15
    replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    // replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-78,-4},{
              -58,16}}, rotation=0)));

    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=5)    annotation (Placement(transformation(extent={{10,36},{30,56}},
            rotation=0)));
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=5)    annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      nPorts=2,
      p=100000,
      T=300)   annotation (Placement(transformation(extent={{74,-20},{62,-8}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-34,36},{
              -14,56}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-30,-50},
              {-10,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1; 2,
          1e-3; 100,1e-3]) 
                       annotation (Placement(transformation(extent={{-20,72},{0,
              92}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 3,1; 4,
          1e-3; 100,1e-3]) 
                       annotation (Placement(transformation(extent={{-20,-10},{
              0,10}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Junctions.JunctionVolume junctionIdeal(
                                          redeclare package Medium = Medium,
      V=1e-3,
      p_start=5.0e5,
      T_start=300,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
      annotation (Placement(transformation(
          origin={-38,6},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-78,6}},
          color={0,127,255}));
    connect(valve2.port_b, sink.ports[2])               annotation (Line(points={{28,-40},
            {46,-40},{46,-15.2},{62,-15.2}},     color={0,127,255}));
    connect(valve1.port_b, sink.ports[1])              annotation (Line(points={{30,46},
            {46,46},{46,-12.8},{62,-12.8}}, color={0,127,255}));
    connect(pipe3.port_b, valve2.port_a)               annotation (Line(points=
            {{-10,-40},{8,-40}}, color={0,127,255}));
    connect(pipe2.port_b, valve1.port_a)              annotation (Line(points={
            {-14,46},{10,46}}, color={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition)              annotation (Line(
          points={{1,82},{20,82},{20,54}}, color={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition)               annotation (Line(
          points={{1,0},{18,0},{18,-32}}, color={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (Line(points={{-58,6},
            {-53.5,6},{-53.5,6},{-48,6}}, color={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (Line(points={{-34,
            46},{-38,46},{-38,16}}, color={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (Line(points={{-38,
            -4},{-38,-40},{-30,-40}}, color={0,127,255}));
  end BranchingPipes15;

  model BranchingPipes16
    replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    // replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-78,-4},{
              -58,16}}, rotation=0)));

    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=5)    annotation (Placement(transformation(extent={{10,36},{30,56}},
            rotation=0)));
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=5)    annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      nPorts=2,
      p=100000,
      T=300)   annotation (Placement(transformation(extent={{74,-20},{62,-8}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-34,36},{
              -14,56}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-30,-50},
              {-10,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          2,0; 100,0]) annotation (Placement(transformation(extent={{-20,72},{0,
              92}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 3,1;
          4,0; 100,0]) annotation (Placement(transformation(extent={{-20,-12},{
              0,8}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Junctions.JunctionVolume junctionIdeal(
                                          redeclare package Medium = Medium,
      V=1e-3,
      p_start=5.0e5,
      T_start=300,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
      annotation (Placement(transformation(
          origin={-38,6},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-78,6}},
          color={0,127,255}));
    connect(valve2.port_b, sink.ports[2])               annotation (Line(points={{28,-40},
            {46,-40},{46,-15.2},{62,-15.2}},     color={0,127,255}));
    connect(valve1.port_b, sink.ports[1])              annotation (Line(points={{30,46},
            {46,46},{46,-12.8},{62,-12.8}}, color={0,127,255}));
    connect(pipe3.port_b, valve2.port_a)               annotation (Line(points=
            {{-10,-40},{8,-40}}, color={0,127,255}));
    connect(pipe2.port_b, valve1.port_a)              annotation (Line(points={
            {-14,46},{10,46}}, color={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition)              annotation (Line(
          points={{1,82},{20,82},{20,54}}, color={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition)               annotation (Line(
          points={{1,-2},{18,-2},{18,-32}}, color={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (Line(points={{-58,6},
            {-53.5,6},{-53.5,6},{-48,6}}, color={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (Line(points={{-34,
            46},{-38,46},{-38,16}}, color={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (Line(points={{-38,
            -4},{-38,-40},{-30,-40}}, color={0,127,255}));
  end BranchingPipes16;

  model BranchingPipes17
    replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    // replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-78,-4},{
              -58,16}}, rotation=0)));

    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=5,
      dp(start=10)) 
                  annotation (Placement(transformation(extent={{10,36},{30,56}},
            rotation=0)));
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=5)    annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(graphics),
      experiment(StopTime=5, Tolerance=1e-007),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (Placement(transformation(extent={{94,-18},{82,-6}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-34,36},{
              -14,56}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-30,-50},
              {-10,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          2,0; 100,0]) annotation (Placement(transformation(extent={{-20,72},{0,
              92}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 3,1;
          4,0; 100,0]) annotation (Placement(transformation(extent={{-18,-12},{
              2,8}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Junctions.JunctionVolume junctionIdeal(
                                          redeclare package Medium = Medium,
      V=1e-3,
      p_start=5.0e5,
      T_start=300,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
      annotation (Placement(transformation(
          origin={-40,6},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Junctions.JunctionVolume junctionVolume(redeclare package Medium = Medium,
        V=1e-3) annotation (Placement(transformation(
          origin={56,-12},
          extent={{-10,10},{10,-10}},
          rotation=90)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-78,6}},
          color={0,127,255}));
    connect(pipe3.port_b, valve2.port_a)               annotation (Line(points=
            {{-10,-40},{8,-40}}, color={0,127,255}));
    connect(pipe2.port_b, valve1.port_a)              annotation (Line(points={
            {-14,46},{10,46}}, color={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition)              annotation (Line(
          points={{1,82},{20,82},{20,55}}, color={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition)               annotation (Line(
          points={{3,-2},{18,-2},{18,-31}}, color={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (Line(points={{-58,6},
            {-54.5,6},{-54.5,6},{-50,6}}, color={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (Line(points={{-34,
            46},{-40,46},{-40,16}}, color={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (Line(points={{-40,
            -4},{-40,-40},{-30,-40}}, color={0,127,255}));
    connect(junctionVolume.port_3, sink.ports[1]) annotation (Line(points={{66,-12},
            {82,-12}}, color={0,127,255}));
    connect(valve2.port_b, junctionVolume.port_1) annotation (Line(points={{28,
            -40},{56,-40},{56,-22}}, color={0,127,255}));
    connect(valve1.port_b, junctionVolume.port_2) annotation (Line(points={{30,
            46},{56,46},{56,-2}}, color={0,127,255}));
  end BranchingPipes17;

  model BranchingPipes18
    // replaceable package Medium = Modelica.Media.Air.DryAirNasa;
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-78,-4},{
              -58,16}}, rotation=0)));

    ControlValves.ValveIncompressible valve1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=5)    annotation (Placement(transformation(extent={{10,36},{30,56}},
            rotation=0)));
    ControlValves.ValveIncompressible valve2(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=4.0e5,
      m_flow_nominal=1,
      d_nominal=5)    annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (Placement(transformation(extent={{94,-18},{82,-6}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe2(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-34,36},{
              -14,56}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5,
      dp_small=10)      annotation (Placement(transformation(extent={{-30,-50},
              {-10,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1; 2,
          0; 100,0])   annotation (Placement(transformation(extent={{-20,72},{0,
              92}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 3,1; 4,
          0; 100,0])   annotation (Placement(transformation(extent={{-18,-12},{
              2,8}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Junctions.JunctionVolume junctionIdeal(
                                          redeclare package Medium = Medium,
      V=1e-3,
      p_start=5.0e5,
      T_start=300,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
      annotation (Placement(transformation(
          origin={-38,6},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Junctions.JunctionVolume junctionVolume(redeclare package Medium = Medium, V=
          1e-3) annotation (Placement(transformation(
          origin={56,-12},
          extent={{-10,10},{10,-10}},
          rotation=90)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-78,6}},
          color={0,127,255}));
    connect(pipe3.port_b, valve2.port_a)               annotation (Line(points=
            {{-10,-40},{8,-40}}, color={0,127,255}));
    connect(pipe2.port_b, valve1.port_a)              annotation (Line(points={
            {-14,46},{10,46}}, color={0,127,255}));
    connect(valveOpening1.y, valve1.stemPosition)              annotation (Line(
          points={{1,82},{20,82},{20,55}}, color={0,0,127}));
    connect(valveOpening2.y, valve2.stemPosition)               annotation (Line(
          points={{3,-2},{18,-2},{18,-31}}, color={0,0,127}));
    connect(pipe1.port_b, junctionIdeal.port_3) annotation (Line(points={{-58,6},
            {-53.5,6},{-53.5,6},{-48,6}}, color={0,127,255}));
    connect(pipe2.port_a, junctionIdeal.port_2) annotation (Line(points={{-34,
            46},{-38,46},{-38,16}}, color={0,127,255}));
    connect(junctionIdeal.port_1, pipe3.port_a) annotation (Line(points={{-38,
            -4},{-38,-40},{-30,-40}}, color={0,127,255}));
    connect(junctionVolume.port_3, sink.ports[1]) annotation (Line(points={{66,-12},
            {82,-12}}, color={0,127,255}));
    connect(valve2.port_b, junctionVolume.port_1) annotation (Line(points={{28,
            -40},{56,-40},{56,-22}}, color={0,127,255}));
    connect(valve1.port_b, junctionVolume.port_2) annotation (Line(points={{30,
            46},{56,46},{56,-2}}, color={0,127,255}));
  end BranchingPipes18;

  model BranchingPipes131
    // replaceable package Medium = Modelica.Media.Air.SimpleAir;
    replaceable package Medium = Modelica.Media.Water.StandardWater;

    Sources.FixedBoundary_pTX source(
      redeclare package Medium = Medium,
      p=5.0e5,
      T=300) annotation (Placement(transformation(extent={{-100,0},{-88,12}},
            rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe1(
      redeclare package Medium = Medium,
      p_a_start=5.0e5,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_b_start=4.95e5) annotation (Placement(transformation(extent={{-78,-4},{
              -58,16}}, rotation=0)));

    ControlValves.ValveIncompressible valveIncompressible1(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      m_flow_nominal=1,
      d_nominal=5,
      minStemPosition=0.001,
      dp(start=400000),
      dp_nominal=400000) 
                  annotation (Placement(transformation(extent={{8,-50},{28,-30}},
            rotation=0)));
    annotation (
      Diagram(graphics),
      experiment(StopTime=5),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
    Sources.FixedBoundary_pTX sink(
      redeclare package Medium = Medium,
      T=300,
      p=1.0e5) annotation (Placement(transformation(extent={{74,-20},{62,-8}},
            rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe pipe3(
      redeclare package Medium = Medium,
      use_T_start=true,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      length=10,
      diameter=2.54e-2,
      p_a_start=4.95e5,
      p_b_start=4.90e5) annotation (Placement(transformation(extent={{-30,-50},
              {-10,-30}}, rotation=0)));

    Modelica.Blocks.Sources.TimeTable valveOpening1(offset=0, table=[0,1; 1,1;
          1,0; 100,0]) annotation (Placement(transformation(extent={{-20,70},{0,
              90}}, rotation=0)));
    Modelica.Blocks.Sources.TimeTable valveOpening2(offset=0, table=[0,1; 2,1;
          2.01,1e-6; 100,0]) 
                       annotation (Placement(transformation(extent={{-20,-8},{0,
              12}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
  equation
    connect(source.ports[1], pipe1.port_a) annotation (Line(points={{-88,6},{-78,6}},
          color={0,127,255}));
    connect(valveIncompressible1.port_b, sink.ports[1]) annotation (Line(points={{
            28,-40},{46,-40},{46,-14},{62,-14}}, color={0,127,255}));
    connect(pipe3.port_b, valveIncompressible1.port_a) annotation (Line(points=
            {{-10,-40},{8,-40}}, color={0,127,255}));
    connect(valveOpening2.y, valveIncompressible1.stemPosition) annotation (Line(
          points={{1,2},{18,2},{18,-31}}, color={0,0,127}));
    connect(pipe1.port_b, pipe3.port_a) annotation (Line(points={{-58,6},{-44,6},
            {-44,-40},{-30,-40}}, color={0,127,255}));
  end BranchingPipes131;

  model DistributedPipeClosingValve "This test demonstrates the importance of smooth regularization of fluid properties for reversing flow.
 A DistributedPipe model with switching port densities and viscosities generates tons of events as the valve closes at time 2."

    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics),
      experiment(StopTime=3));
    Modelica_Fluid.Sources.FixedBoundary source(redeclare package Medium = 
          Modelica.Media.Water.StandardWater, p=200000) 
      annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
    Modelica_Fluid.Pipes.DistributedPipe pipe(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      length=1,
      diameter=0.32,
      initType=Modelica_Fluid.Types.Init.SteadyState,
      use_T_start=false,
      p_a_start=200000,
      p_b_start=200000) 
      annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
    Modelica_Fluid.ControlValves.ValveIncompressible valve(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      m_flow_nominal=10,
      Av=1e-3,
      dp_nominal=100000) 
      annotation (Placement(transformation(extent={{0,-10},{20,10}})));
    Modelica_Fluid.Sources.FixedBoundary sink(redeclare package Medium = 
          Modelica.Media.Water.StandardWater, p=100000) 
                annotation (Placement(transformation(extent={{60,-10},{40,10}})));
    Modelica.Blocks.Sources.Ramp ramp(
      height=-1,
      offset=1,
      duration=1,
      startTime=1) 
                annotation (Placement(transformation(extent={{46,30},{26,50}})));
    inner System system 
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  equation
    connect(source.ports[1], pipe.port_a)         annotation (Line(
        points={{-60,0},{-40,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipe.port_b, valve.port_a)               annotation (Line(
        points={{-20,0},{0,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(valve.port_b, sink.ports[1])                          annotation (Line(
        points={{20,0},{40,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(ramp.y, valve.stemPosition)               annotation (Line(
        points={{25,40},{10,40},{10,8}},
        color={0,0,127},
        smooth=Smooth.None));
  end DistributedPipeClosingValve;

  model DistributedPipeInitialization
    "Steady-state initialization of a distributed pipe"

    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics={
          Text(
            extent={{-43,-52},{43,-60}},
            lineColor={0,0,255},
            textString="or b) Change to ValveIncompressible and it works!?"),
          Text(
            extent={{-40,-34},{62,-54}},
            lineColor={0,0,255},
            textString=
                "a) Select WallFriction = \"No pipe wall friction\" and it works!?"),
          Text(
            extent={{-52,-24},{50,-38}},
            lineColor={0,0,255},
            textString="How to make steady-state initialization work:")}),
      experiment(StopTime=4),
      experimentSetupOutput,
      uses(Modelica_Fluid(version="1.0 Streams Beta 3"), Modelica(version="3.0")));

    Modelica_Fluid.Sources.FixedBoundary source(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      use_T=false,
      p=10000000,
      h=2e6) 
      annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
    Pipes.DistributedPipe pipe(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      nNodes=5,
      h_start=2e6,
      diameter=0.05,
      length=200,
      use_T_start=false,
      modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb,
      p_a_start=10000000,
      p_b_start=9900000) 
      annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
    ControlValves.ValveCompressible valve(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      Av=1e-3,
      dp_nominal=10000000,
      m_flow_nominal=10) 
      annotation (Placement(transformation(extent={{0,-10},{20,10}})));
    Modelica_Fluid.Sources.FixedBoundary sink(redeclare package Medium = 
          Modelica.Media.Water.StandardWaterOnePhase, p=9500000) 
                annotation (Placement(transformation(extent={{60,-10},{40,10}})));
    Modelica.Blocks.Sources.Ramp ramp(
      offset=1,
      duration=0.1,
      height=-0.5,
      startTime=2) 
                annotation (Placement(transformation(extent={{46,30},{26,50}})));
    inner Modelica_Fluid.System system(initType=Modelica_Fluid.Types.Init.SteadyState) 
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
    discrete Modelica.SIunits.MassFlowRate m_flow_initial;
  equation
    when time > 0.1 then
      m_flow_initial = valve.port_a.m_flow;
    end when;
    if pipe.initType == Modelica_Fluid.Types.Init.SteadyState or 
       pipe.dynamicsType == Modelica_Fluid.Types.Dynamics.SteadyState then
      when time > 1 then
        assert(abs(valve.port_a.m_flow - m_flow_initial) < 1e-3, "!!!THE SIMULATION DID NOT START IN STEADY-STATE!!!");
      end when;
    end if;
    connect(source.ports[1], pipe.port_a)         annotation (Line(
        points={{-60,0},{-40,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipe.port_b, valve.port_a)               annotation (Line(
        points={{-20,0},{0,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(valve.port_b, sink.ports[1])                          annotation (Line(
        points={{20,0},{40,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(ramp.y, valve.stemPosition)               annotation (Line(
        points={{25,40},{10,40},{10,8}},
        color={0,0,127},
        smooth=Smooth.None));
  end DistributedPipeInitialization;

  model LumpedPipeInitialization "Steady-state initialization of a lumped pipe"

    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics),
      experiment(StopTime=4),
      experimentSetupOutput,
      uses(Modelica_Fluid(version="1.0 Streams Beta 3"), Modelica(version="3.0")));
    Modelica_Fluid.Sources.FixedBoundary source(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      use_T=false,
      p=10000000,
      h=2e6) 
      annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
    Pipes.LumpedPipe pipe(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      h_start=2e6,
      diameter=0.05,
      length=200,
      use_T_start=false,
      p_a_start=10000000,
      p_b_start=9900000) 
      annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
    Modelica_Fluid.ControlValves.ValveCompressible valve(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      Av=1e-3,
      dp_nominal=10000000,
      m_flow_nominal=10) 
      annotation (Placement(transformation(extent={{0,-10},{20,10}})));
    Modelica_Fluid.Sources.FixedBoundary sink(redeclare package Medium = 
          Modelica.Media.Water.StandardWaterOnePhase, p=9500000) 
                annotation (Placement(transformation(extent={{60,-10},{40,10}})));
    Modelica.Blocks.Sources.Ramp ramp(
      offset=1,
      duration=0.1,
      height=-0.5,
      startTime=2) 
                annotation (Placement(transformation(extent={{46,30},{26,50}})));
    inner Modelica_Fluid.System system(initType=Modelica_Fluid.Types.Init.SteadyState) 
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
    discrete Modelica.SIunits.MassFlowRate m_flow_initial;
  equation
    when time > 0.1 then
      m_flow_initial = valve.port_a.m_flow;
    end when;
    if pipe.initType == Modelica_Fluid.Types.Init.SteadyState or 
       pipe.dynamicsType == Modelica_Fluid.Types.Dynamics.SteadyState then
      when time > 1 then
        assert(abs(valve.port_a.m_flow - m_flow_initial) < 1e-3, "!!!THE SIMULATION DID NOT START IN STEADY-STATE!!!");
      end when;
    end if;
    connect(source.ports[1], pipe.port_a)         annotation (Line(
        points={{-60,0},{-40,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipe.port_b, valve.port_a)               annotation (Line(
        points={{-20,0},{0,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(valve.port_b, sink.ports[1])                          annotation (Line(
        points={{20,0},{40,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(ramp.y, valve.stemPosition)               annotation (Line(
        points={{25,40},{10,40},{10,8}},
        color={0,0,127},
        smooth=Smooth.None));
  end LumpedPipeInitialization;
end TestCriticalCases;
