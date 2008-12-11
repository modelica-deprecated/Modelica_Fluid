within Modelica_Fluid.Examples;
model BranchingDistributedPipes
  "Multi-way connections of pipes with and without Fittings.MultiPort"
  import Modelica_Fluid;
extends Modelica.Icons.Example;
//replaceable package Medium=Modelica.Media.Water.StandardWater;
replaceable package Medium=Modelica.Media.Air.DryAirNasa;  //

 Modelica_Fluid.Pipes.DistributedPipe pipe2(
    redeclare package Medium = Medium,
    use_T_start=true,
    diameter=0.01,
    nNodes=5,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedFlow(from_dp=true),
    m_flow_start=0.1,
    length=2,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare model HeatTransfer = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_localLamTurb,
    p_a_start=100000,
    p_b_start=100000,
    T_start=280,
    modelStructure=Modelica_Fluid.Types.ModelStructure.a_v_b) 
            annotation (Placement(transformation(extent={{-14,50},{6,70}},
          rotation=0)));

  Modelica_Fluid.Sources.PrescribedBoundary_pTX boundary2(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    usePressureInput=true,
    useTemperatureInput=false)                                      annotation (Placement(
        transformation(extent={{90,30},{70,50}}, rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe3(
    redeclare package Medium=Medium,
    length=1,
    use_T_start=true,
    diameter=0.01,
    nNodes=5,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedFlow(from_dp=true),
    m_flow_start=0.1,
    initType=Modelica_Fluid.Types.Init.NoInit,
    p_a_start=100000,
    p_b_start=100000,
    T_start=340,
    modelStructure=Modelica_Fluid.Types.ModelStructure.av_b) 
            annotation (Placement(transformation(extent={{40,30},{60,50}},
          rotation=0)));

  Modelica_Fluid.Pipes.DistributedPipe pipe1(
    redeclare package Medium=Medium,
    use_T_start=true,
    diameter=0.01,
    nNodes=5,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedFlow(from_dp=true),
    m_flow_start=0.1,
    length=0.4,
    initType=Modelica_Fluid.Types.Init.NoInit,
    p_a_start=100000,
    p_b_start=100000,
    T_start=300,
    modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb) 
            annotation (Placement(transformation(extent={{-68,30},{-48,50}},
          rotation=0)));

  Modelica_Fluid.Sources.FixedBoundary_pTX boundary1(
    T=280,
    redeclare package Medium = Medium,
    p=1.2e5)                                                        annotation (Placement(
        transformation(extent={{-100,30},{-80,50}}, rotation=0)));

    annotation (extent=[-90,-86; -70,-66], Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
                                     annotation (points=[-22,42; -2,42; -2,32;
        7.8,32],
      style(color=69, rgbcolor={0,127,255}));
  Modelica.Blocks.Sources.Ramp ramp(
    offset=1e5,
    startTime=5,
    duration=0,
    height=0.4e5) 
                annotation (Placement(transformation(extent={{70,60},{90,80}},
          rotation=0)));

  Modelica_Fluid.Pipes.DistributedPipe pipe4(
    redeclare package Medium=Medium,
    length=1,
    use_T_start=true,
    diameter=0.01,
    nNodes=5,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedFlow(from_dp=true),
    m_flow_start=0.1,
    heatTransfer(alpha0=1000),
    initType=Modelica_Fluid.Types.Init.NoInit,
    p_a_start=100000,
    p_b_start=100000,
    T_start=360,
    modelStructure=Modelica_Fluid.Types.ModelStructure.a_v_b) 
            annotation (Placement(transformation(extent={{-14,10},{6,30}},
          rotation=0)));

  inner Modelica_Fluid.System system 
    annotation (Placement(transformation(extent={{-90,72},{-70,92}},  rotation=
            0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow[pipe2.nNodes] heat2(each
      Q_flow=200, each alpha=10000) 
    annotation (Placement(transformation(extent={{-34,70},{-14,90}}, rotation=0)));
 Modelica_Fluid.Pipes.DistributedPipe pipe5(
    redeclare package Medium = Medium,
    use_T_start=true,
    diameter=0.01,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedFlow(from_dp=true),
    m_flow_start=0.1,
    length=2,
    redeclare model HeatTransfer = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_localLamTurb,
    p_a_start=100000,
    p_b_start=100000,
    T_start=280) 
            annotation (Placement(transformation(extent={{-14,-40},{6,-20}},
          rotation=0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX boundary4(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    usePressureInput=true,
    useTemperatureInput=false)                                      annotation (Placement(
        transformation(extent={{88,-60},{68,-40}}, rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe6(
    redeclare package Medium=Medium,
    length=1,
    use_T_start=true,
    diameter=0.01,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedFlow(from_dp=true),
    m_flow_start=0.1,
    p_a_start=100000,
    p_b_start=100000,
    T_start=340,
    modelStructure=Modelica_Fluid.Types.ModelStructure.av_b) 
            annotation (Placement(transformation(extent={{40,-60},{60,-40}},
          rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe7(
    redeclare package Medium=Medium,
    use_T_start=true,
    T_start=300,
    diameter=0.01,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedFlow(from_dp=true),
    m_flow_start=0.1,
    p_b_start=1.0e5,
    length=0.4,
    p_a_start=100000) 
            annotation (Placement(transformation(extent={{-68,-60},{-48,-40}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX boundary3(
    T=280,
    redeclare package Medium = Medium,
    p=1.2e5)                                                        annotation (Placement(
        transformation(extent={{-100,-60},{-80,-40}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp ramp1(
    offset=1e5,
    startTime=5,
    duration=0,
    height=0.4e5) 
                annotation (Placement(transformation(extent={{68,-30},{88,-10}},
          rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe8(
    redeclare package Medium=Medium,
    length=1,
    use_T_start=true,
    diameter=0.01,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedFlow(from_dp=true),
    m_flow_start=0.1,
    heatTransfer(alpha0=1000),
    p_a_start=100000,
    p_b_start=100000,
    T_start=360) 
            annotation (Placement(transformation(extent={{-14,-80},{6,-60}},
          rotation=0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow[
                                              pipe5.nNodes] heat5(each Q_flow=200,
      each alpha=10000) 
    annotation (Placement(transformation(extent={{-34,-20},{-14,0}},  rotation=
            0)));
  Modelica_Fluid.Fittings.MultiPort multiPort1(redeclare package Medium = 
        Medium, nPorts_b=2) 
    annotation (Placement(transformation(extent={{-42,30},{-34,50}})));
  Modelica_Fluid.Fittings.MultiPort multiPort3(redeclare package Medium = 
        Medium, nPorts_b=2) 
    annotation (Placement(transformation(extent={{34,30},{26,50}})));
equation
  connect(boundary1.ports[1], pipe1.port_a) annotation (Line(
      points={{-80,40},{-68,40}},
      color={0,127,255},
      thickness=0.5));
  connect(ramp.y, boundary2.p_in) annotation (Line(
      points={{91,70},{96,70},{96,48},{92,48}},
      color={0,0,127},
      thickness=0.5));
  connect(pipe3.port_b, boundary2.ports[1]) annotation (Line(
      points={{60,40},{70,40}},
      color={0,127,255},
      thickness=0.5));
  connect(heat2.port, pipe2.heatPorts) 
                                     annotation (Line(
      points={{-14,80},{-3.9,80},{-3.9,65.2}},
      color={191,0,0},
      thickness=0.5));
  connect(ramp1.y, boundary4.p_in) annotation (Line(
      points={{89,-20},{96,-20},{96,-42},{90,-42}},
      color={0,0,127},
      thickness=0.5));
  connect(boundary3.ports[1], pipe7.port_a) annotation (Line(
      points={{-80,-50},{-68,-50}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe7.port_b, pipe5.port_a) annotation (Line(
      points={{-48,-50},{-30,-50},{-30,-30},{-14,-30}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe7.port_b, pipe8.port_a) annotation (Line(
      points={{-48,-50},{-30,-50},{-30,-70},{-14,-70}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe5.port_b, pipe6.port_a) annotation (Line(
      points={{6,-30},{22,-30},{22,-50},{40,-50}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe8.port_b, pipe6.port_a) annotation (Line(
      points={{6,-70},{22,-70},{22,-50},{40,-50}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe6.port_b, boundary4.ports[1]) annotation (Line(
      points={{60,-50},{60,-50},{68,-50}},
      color={0,127,255},
      thickness=0.5));
  connect(heat5.port, pipe5.heatPorts) 
                                      annotation (Line(
      points={{-14,-10},{-3.9,-10},{-3.9,-24.8}},
      color={191,0,0},
      thickness=0.5));
  connect(pipe1.port_b, multiPort1.port_a) annotation (Line(
      points={{-48,40},{-42,40}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(multiPort1.ports_b[1], pipe2.port_a) annotation (Line(
      points={{-34,42},{-30,42},{-30,60},{-14,60}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(multiPort1.ports_b[2], pipe4.port_a) annotation (Line(
      points={{-34,38},{-30,38},{-30,20},{-14,20}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(multiPort3.port_a, pipe3.port_a) annotation (Line(
      points={{34,40},{40,40}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(pipe2.port_b, multiPort3.ports_b[1]) annotation (Line(
      points={{6,60},{22,60},{22,42},{26,42}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(pipe4.port_b, multiPort3.ports_b[2]) annotation (Line(
      points={{6,20},{22,20},{22,38},{26,38}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));

  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics),
                       experiment(StopTime=20, Tolerance=1e-005),
    experimentSetupOutput,
    Documentation(info="<html>
<p>
This model demonstrates the difference in using an explicit junction model for branching flows compared to direct connections. 
</p>
<p>
Specifying appropriate model structures, the states of the outermost pipe segments get exposed to the fluid ports.
The mass balances are lumped together in connection sets. 
Energy and substance blances remain independent in the connected pipe segments with the stream concept in fluid ports.
One needs to be careful when fixing pressure states by connecting them to a prescribed boundary with 
non-differentiable boundary conditions.
</p>
<p>
Moreover, algebraic energy and substance balances result from the interconnection of multiple flow models, 
like pipe2 and pipe4 with modelStructure=a_v_b. This is normally not intended. 
Connecting to models with multiple ports, e.g. a Source, or introducing a 
<a href=\"Modelica:Modelica_Fluid.Fittings.MultiPort\">MultiPort</a>, the algebraic 
energy and substance balances are avoided; the mixing takes place in the volume model. 
</p>
</html>"));

end BranchingDistributedPipes;
