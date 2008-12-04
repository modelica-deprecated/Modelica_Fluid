within Modelica_Fluid.Examples;
model BranchingDistributedPipes
  import Modelica_Fluid;
extends Modelica.Icons.Example;
//replaceable package Medium=Modelica.Media.Water.StandardWater;
replaceable package Medium=Modelica.Media.Air.DryAirNasa;  //

 Modelica_Fluid.Pipes.DistributedPipe pipe2(
    redeclare package Medium = Medium,
    use_T_start=true,
    from_dp=true,
    diameter=0.01,
    nNodes=5,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    length=2,
    use_approxPortProperties=true,
    initType=Modelica_Fluid.Types.Init.NoInit,
    heatTransfer(alpha0=500),
    p_a_start=100000,
    p_b_start=100000,
    T_start=280) 
            annotation (Placement(transformation(extent={{-34,38},{-14,58}},
          rotation=0)));

  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics),
                       experiment(StopTime=20, Tolerance=1e-005),
    experimentSetupOutput,
    Documentation(info="<html>
<p>This model demonstrates the difference in using an explicit junction model for branching flows compared to direct connections. The user is encouraged to always use explicit junctions since otherwise errors are easily introduced. In this model two identical systems are composed independent from one another, except for the additional junction elements in the top case. Heat is added to the fluid in <tt><b>pipe2</b></tt> and <tt><b>pipe5</b></tt>, respectively. After five seconds of simulation time the flow is reversed by raising the pressure in the reservoirs in one end of each system. When looking at the port enthalpies in <tt><b>pipe2</b></tt> and <tt><b>pipe5</b></tt>, it becomes clear that in the top case they resemble the specific enthalpies of the flow just inside the component boundary, while in the second case the port properties represent the mixing properties of the infinitesimal volume at the connection point, therefore being located just outside the component boundary.</p>
<p>In order to safely use approach two without explicit junctions the component model must not use any downstream properties, since they may not belong to the component itself. This is the case for the used pipe model but cannot be guaranteed for all fluid models. Additionally it may be very confusing to find mixing properties at ports that are graphically located some distance away from the other ports in the connection set.</p>
<p>Although the introduced dynamic junction volume may break undesired non-linear systems, it must be noted that if chosen too small it may also lead to a stiff system and therefore slow down the simulation.</p>
</html>"));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX boundary2(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    usePressureInput=true,
    useTemperatureInput=false)                                      annotation (Placement(
        transformation(extent={{70,24},{50,44}}, rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe3(
    redeclare package Medium=Medium,
    length=1,
    use_T_start=true,
    from_dp=true,
    diameter=0.01,
    nNodes=5,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    use_approxPortProperties=true,
    initType=Modelica_Fluid.Types.Init.NoInit,
    p_a_start=100000,
    p_b_start=100000,
    T_start=340,
    modelStructure=Modelica_Fluid.Types.ModelStructure.av_b) 
            annotation (Placement(transformation(extent={{20,24},{40,44}},
          rotation=0)));

  Modelica_Fluid.Pipes.DistributedPipe pipe1(
    redeclare package Medium=Medium,
    use_T_start=true,
    from_dp=true,
    diameter=0.01,
    nNodes=5,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    length=0.4,
    use_approxPortProperties=true,
    initType=Modelica_Fluid.Types.Init.NoInit,
    p_a_start=100000,
    p_b_start=100000,
    T_start=300,
    modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb) 
            annotation (Placement(transformation(extent={{-88,24},{-68,44}},
          rotation=0)));

  Modelica_Fluid.Sources.FixedBoundary_pTX boundary1(
    T=280,
    redeclare package Medium = Medium,
    p=1.2e5)                                                        annotation (Placement(
        transformation(extent={{-114,24},{-94,44}}, rotation=0)));

    annotation (extent=[-90,-86; -70,-66]);
                                     annotation (points=[-22,42; -2,42; -2,32;
        7.8,32],
      style(color=69, rgbcolor={0,127,255}));
  Modelica.Blocks.Sources.Ramp ramp(
    offset=1e5,
    startTime=5,
    duration=0,
    height=0.4e5) 
                annotation (Placement(transformation(extent={{98,66},{78,86}},
          rotation=0)));

  Modelica_Fluid.Pipes.DistributedPipe pipe4(
    redeclare package Medium=Medium,
    length=1,
    use_T_start=true,
    from_dp=true,
    diameter=0.01,
    T_start=360,
    nNodes=5,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    p_b_start=1.0e5,
    heatTransfer(alpha0=1000),
    use_approxPortProperties=true,
    initType=Modelica_Fluid.Types.Init.NoInit,
    p_a_start=100000) 
            annotation (Placement(transformation(extent={{-34,8},{-14,28}},
          rotation=0)));

  inner Modelica_Fluid.System system 
    annotation (Placement(transformation(extent={{70,-100},{90,-80}}, rotation=
            0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow[
                                              pipe2.n] heat(each Q_flow=200,
      each alpha=10000) 
    annotation (Placement(transformation(extent={{-60,72},{-40,92}}, rotation=0)));
 Modelica_Fluid.Pipes.DistributedPipe pipe5(
    redeclare package Medium = Medium,
    use_T_start=true,
    from_dp=true,
    diameter=0.01,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    length=2,
    use_approxPortProperties=true,
    p_a_start=100000,
    p_b_start=100000,
    T_start=280,
    heatTransfer(alpha0=500)) 
            annotation (Placement(transformation(extent={{-28,-50},{-8,-30}},
          rotation=0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX boundary4(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    usePressureInput=true,
    useTemperatureInput=false)                                      annotation (Placement(
        transformation(extent={{70,-70},{50,-50}}, rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe6(
    redeclare package Medium=Medium,
    T_start=340,
    length=1,
    use_T_start=true,
    from_dp=true,
    p_b_start=1e5,
    diameter=0.01,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    use_approxPortProperties=true,
    p_a_start=100000) 
            annotation (Placement(transformation(extent={{16,-70},{36,-50}},
          rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe7(
    redeclare package Medium=Medium,
    use_T_start=true,
    from_dp=true,
    T_start=300,
    diameter=0.01,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    p_b_start=1.0e5,
    length=0.4,
    use_approxPortProperties=true,
    p_a_start=100000) 
            annotation (Placement(transformation(extent={{-66,-70},{-46,-50}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX boundary3(
    T=280,
    redeclare package Medium = Medium,
    p=1.2e5)                                                        annotation (Placement(
        transformation(extent={{-110,-70},{-90,-50}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp ramp1(
    offset=1e5,
    startTime=5,
    duration=0,
    height=0.4e5) 
                annotation (Placement(transformation(extent={{100,-20},{80,0}},
          rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe8(
    redeclare package Medium=Medium,
    length=1,
    use_T_start=true,
    from_dp=true,
    diameter=0.01,
    T_start=360,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    p_b_start=1.0e5,
    heatTransfer(alpha0=1000),
    use_approxPortProperties=true,
    p_a_start=100000) 
            annotation (Placement(transformation(extent={{-28,-90},{-8,-70}},
          rotation=0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow[
                                              pipe2.n] heat1(each Q_flow=200,
      each alpha=10000) 
    annotation (Placement(transformation(extent={{-62,-24},{-42,-4}}, rotation=
            0)));
  Modelica_Fluid.Junctions.MultiPort multiPort1(redeclare package Medium = 
        Medium, nPorts_b=2) 
    annotation (Placement(transformation(extent={{-62,24},{-54,44}})));
  Modelica_Fluid.Junctions.MultiPort multiPort3(redeclare package Medium = 
        Medium, nPorts_b=2) 
    annotation (Placement(transformation(extent={{14,24},{6,44}})));
equation
  connect(boundary1.ports[1], pipe1.port_a) annotation (Line(
      points={{-94,34},{-88,34}},
      color={0,127,255},
      thickness=0.5));
  connect(ramp.y, boundary2.p_in) annotation (Line(
      points={{77,76},{76,76},{76,40},{72,40}},
      color={0,0,127},
      thickness=0.5));
  connect(pipe3.port_b, boundary2.ports[1]) annotation (Line(
      points={{40,34},{50,34}},
      color={0,127,255},
      thickness=0.5));
  connect(heat.port, pipe2.heatPorts) 
                                     annotation (Line(
      points={{-40,82},{-23.9,82},{-23.9,53.2}},
      color={191,0,0},
      thickness=0.5));
  connect(ramp1.y, boundary4.p_in) annotation (Line(
      points={{79,-10},{78,-10},{78,-54},{72,-54}},
      color={0,0,127},
      thickness=0.5));
  connect(boundary3.ports[1], pipe7.port_a) annotation (Line(
      points={{-90,-60},{-80,-60},{-66,-60}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe7.port_b, pipe5.port_a) annotation (Line(
      points={{-46,-60},{-40,-60},{-40,-40},{-28,-40}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe7.port_b, pipe8.port_a) annotation (Line(
      points={{-46,-60},{-40,-60},{-40,-80},{-28,-80}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe5.port_b, pipe6.port_a) annotation (Line(
      points={{-8,-40},{4,-40},{4,-60},{16,-60}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe8.port_b, pipe6.port_a) annotation (Line(
      points={{-8,-80},{4,-80},{4,-60},{16,-60}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe6.port_b, boundary4.ports[1]) annotation (Line(
      points={{36,-60},{36,-60},{50,-60}},
      color={0,127,255},
      thickness=0.5));
  connect(heat1.port, pipe5.heatPorts) 
                                      annotation (Line(
      points={{-42,-14},{-17.9,-14},{-17.9,-34.8}},
      color={191,0,0},
      thickness=0.5));
  connect(pipe1.port_b, multiPort1.port_a) annotation (Line(
      points={{-68,34},{-62,34}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(multiPort1.ports_b[1], pipe2.port_a) annotation (Line(
      points={{-54,36},{-46,36},{-46,48},{-34,48}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(multiPort1.ports_b[2], pipe4.port_a) annotation (Line(
      points={{-54,32},{-46,32},{-46,18},{-34,18}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(multiPort3.port_a, pipe3.port_a) annotation (Line(
      points={{14,34},{20,34}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(pipe2.port_b, multiPort3.ports_b[1]) annotation (Line(
      points={{-14,48},{0,48},{0,36},{6,36}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(pipe4.port_b, multiPort3.ports_b[2]) annotation (Line(
      points={{-14,18},{0,18},{0,32},{6,32}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
end BranchingDistributedPipes;
