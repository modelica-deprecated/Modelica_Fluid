within Modelica_Fluid.Test.TestComponents.Pipes;
model TestDistributedPipe01
  import Modelica_Fluid;
extends Modelica.Icons.Example;
replaceable package Medium=Modelica.Media.Water.StandardWater;
//replaceable package Medium=Modelica.Media.Air.DryAirNasa;  //

 Modelica_Fluid.Pipes.DistributedPipe pipe2(
    redeclare package Medium = Medium,
    use_T_start=true,
    from_dp=true,
    T_start=280,
    diameter=0.01,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    p_b_start=1.0e5,
    length=2,
    use_approxPortProperties=true,
    heatTransfer(each alpha0=500),
    p_a_start=100000) 
            annotation (Placement(transformation(extent={{-28,68},{-8,88}},
          rotation=0)));

  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics),
                       experiment(StopTime=20, Tolerance=1e-005),
    experimentSetupOutput,
    Documentation(info="<html>
Test of different distributed pipe models. The first system uses explicit junctions, in the third system some of the pipe models are replaced by non-symmetric components.
</html>"));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX boundary2(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    usePressureInput=true,
    useTemperatureInput=false)                                      annotation (Placement(
        transformation(extent={{84,54},{64,74}}, rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe3(
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
            annotation (Placement(transformation(extent={{34,54},{54,74}},
          rotation=0)));

  Modelica_Fluid.Pipes.DistributedPipe pipe1(
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
            annotation (Placement(transformation(extent={{-82,54},{-62,74}},
          rotation=0)));

  Modelica_Fluid.Sources.FixedBoundary_pTX boundary1(
    T=280,
    redeclare package Medium = Medium,
    p=1.5e5)                                                        annotation (Placement(
        transformation(extent={{-108,54},{-88,74}}, rotation=0)));

    annotation (extent=[-90,-86; -70,-66]);
                                     annotation (points=[-22,42; -2,42; -2,32;
        7.8,32],
      style(color=69, rgbcolor={0,127,255}));
  Modelica.Blocks.Sources.Ramp ramp(
    offset=1e5,
    startTime=5,
    duration=0,
    height=1.0e5) 
                annotation (Placement(transformation(extent={{104,64},{92,76}},
          rotation=0)));

  Modelica_Fluid.Pipes.DistributedPipe pipe4(
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
    heatTransfer(each alpha0=1000),
    use_approxPortProperties=true,
    p_a_start=100000) 
            annotation (Placement(transformation(extent={{-28,38},{-8,58}},
          rotation=0)));

  inner Modelica_Fluid.System system 
    annotation (Placement(transformation(extent={{72,-94},{92,-74}}, rotation=0)));
  Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.GenericJunction junction1(
    nPorts_b=2,
    redeclare package Medium = Medium,
    V=0.0001,
    dp_nominal=100000,
    m_flow_nominal=0.01,
    p_start=100000,
    T_start=300)                       annotation (Placement(transformation(
          extent={{-58,54},{-38,74}}, rotation=0)));
  Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.GenericJunction junction2(
    nPorts_a=2,
    redeclare package Medium = Medium,
    V=0.00001,
    dp_nominal=100000,
    m_flow_nominal=0.01,
    p_start=100000,
    T_start=300)                       annotation (Placement(transformation(
          extent={{6,54},{26,74}}, rotation=0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow[
                                              pipe2.n] heat(each Q_flow=200,
      each alpha=10000,
      each T_ref=350) 
    annotation (Placement(transformation(extent={{-54,80},{-34,100}}, rotation=
            0)));
 Modelica_Fluid.Pipes.DistributedPipe pipe5(
    redeclare package Medium = Medium,
    use_T_start=true,
    from_dp=true,
    T_start=280,
    diameter=0.01,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    p_b_start=1.0e5,
    length=2,
    use_approxPortProperties=true,
    heatTransfer(each alpha0=500),
    p_a_start=100000) 
            annotation (Placement(transformation(extent={{-30,10},{-10,30}},
          rotation=0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX boundary4(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    usePressureInput=true,
    useTemperatureInput=false)                                      annotation (Placement(
        transformation(extent={{68,-2},{48,18}}, rotation=0)));
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
            annotation (Placement(transformation(extent={{10,-2},{30,18}},
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
            annotation (Placement(transformation(extent={{-68,-2},{-48,18}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX boundary3(
    T=280,
    redeclare package Medium = Medium,
    p=1.5e5)                                                        annotation (Placement(
        transformation(extent={{-108,-2},{-88,18}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp ramp1(
    offset=1e5,
    startTime=5,
    duration=0,
    height=1.0e5) 
                annotation (Placement(transformation(extent={{96,8},{84,20}},
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
    heatTransfer(each alpha0=1000),
    use_approxPortProperties=true,
    p_a_start=100000) 
            annotation (Placement(transformation(extent={{-30,-12},{-10,8}},
          rotation=0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow[
                                              pipe2.n] heat1(each Q_flow=200,
      each alpha=10000,
      each T_ref=350) 
    annotation (Placement(transformation(extent={{-72,22},{-52,42}}, rotation=0)));
 Modelica_Fluid.Pipes.DistributedPipe pipe9(
    redeclare package Medium = Medium,
    use_T_start=true,
    from_dp=true,
    T_start=280,
    diameter=0.01,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    p_b_start=1.0e5,
    length=2,
    use_approxPortProperties=true,
    heatTransfer(each alpha0=500),
    p_a_start=100000) 
            annotation (Placement(transformation(extent={{-29,-50},{-9,-30}},
          rotation=0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX boundary5(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    usePressureInput=true,
    useTemperatureInput=false)                                      annotation (Placement(
        transformation(extent={{70,-62},{50,-42}}, rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe10(
    redeclare package Medium=Medium,
    length=1,
    use_T_start=true,
    from_dp=true,
    diameter=0.01,
    nNodes=5,
    modelStructure=Modelica_Fluid.Types.ModelStructure.av_b,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    use_approxPortProperties=true,
    p_a_start=100000,
    p_b_start=100000,
    T_start=340) 
            annotation (Placement(transformation(extent={{12,-62},{32,-42}},
          rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe11(
    redeclare package Medium=Medium,
    use_T_start=true,
    from_dp=true,
    diameter=0.01,
    nNodes=5,
    modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    length=0.4,
    use_approxPortProperties=true,
    p_a_start=100000,
    p_b_start=100000,
    T_start=300) 
            annotation (Placement(transformation(extent={{-70,-62},{-50,-42}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX boundary6(
    T=280,
    redeclare package Medium = Medium,
    p=1.5e5)                                                        annotation (Placement(
        transformation(extent={{-110,-62},{-90,-42}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp ramp2(
    offset=1e5,
    startTime=5,
    duration=0,
    height=1.0e5) 
                annotation (Placement(transformation(extent={{96,-52},{84,-40}},
          rotation=0)));
  Modelica_Fluid.Pipes.DistributedPipe pipe12(
    redeclare package Medium=Medium,
    length=1,
    use_T_start=true,
    from_dp=true,
    diameter=0.01,
    nNodes=5,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    m_flow_start=0.1,
    heatTransfer(each alpha0=1000),
    use_approxPortProperties=true,
    p_a_start=100000,
    p_b_start=100000,
    T_start=360) 
            annotation (Placement(transformation(extent={{-29,-72},{-9,-52}},
          rotation=0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow[
                                              pipe2.n] heat2(each Q_flow=200,
      each alpha=10000,
      each T_ref=350) 
    annotation (Placement(transformation(extent={{-72,-38},{-52,-18}}, rotation=
           0)));
  Modelica_Fluid.Junctions.MultiPort muliPort11(nPorts_b=2, redeclare package
      Medium = Medium) 
    annotation (Placement(transformation(extent={{-46,-62},{-38,-42}})));
  Modelica_Fluid.Junctions.MultiPort multiPort10(nPorts_b=2, redeclare package
      Medium = Medium) 
    annotation (Placement(transformation(extent={{8,-62},{0,-42}})));
equation
  connect(boundary1.ports[1], pipe1.port_a) 
                                        annotation (Line(
      points={{-88,64},{-82,64}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe1.port_b, junction1.ports_a[1]) annotation (Line(
      points={{-62,64},{-59.5,64},{-58,64}},
      color={0,127,255},
      thickness=0.5));
  connect(junction1.ports_b[1], pipe2.port_a) annotation (Line(
      points={{-38,66},{-34,66},{-34,78},{-28,78}},
      color={0,127,255},
      thickness=0.5));
  connect(junction1.ports_b[2], pipe4.port_a) annotation (Line(
      points={{-38,62},{-34,62},{-34,48},{-28,48}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe2.port_b, junction2.ports_a[1]) annotation (Line(
      points={{-8,78},{-2,78},{-2,66},{6,66}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe4.port_b, junction2.ports_a[2]) annotation (Line(
      points={{-8,48},{-2,48},{-2,62},{6,62}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe3.port_b, boundary2.ports[1]) 
                                        annotation (Line(
      points={{54,64},{64,64}},
      color={0,127,255},
      thickness=0.5));
  connect(junction2.ports_b[1], pipe3.port_a) annotation (Line(
      points={{26,64},{29.5,64},{34,64}},
      color={0,127,255},
      thickness=0.5));
  connect(heat.port, pipe2.heatPorts) 
                                     annotation (Line(
      points={{-34,90},{-17.9,90},{-17.9,83.2}},
      color={191,0,0},
      thickness=0.5));
  connect(ramp1.y, boundary4.p_in) annotation (Line(
      points={{83.4,14},{70,14}},
      color={0,0,127},
      thickness=0.5));
  connect(boundary3.ports[1], pipe7.port_a) 
                                        annotation (Line(
      points={{-88,8},{-68,8}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe7.port_b, pipe5.port_a) annotation (Line(
      points={{-48,8},{-42,8},{-42,20},{-30,20}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe7.port_b, pipe8.port_a) annotation (Line(
      points={{-48,8},{-42,8},{-42,-2},{-30,-2}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe5.port_b, pipe6.port_a) annotation (Line(
      points={{-10,20},{2,20},{2,8},{10,8}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe8.port_b, pipe6.port_a) annotation (Line(
      points={{-10,-2},{2,-2},{2,8},{10,8}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe6.port_b, boundary4.ports[1]) 
                                        annotation (Line(
      points={{30,8},{48,8}},
      color={0,127,255},
      thickness=0.5));
  connect(heat1.port, pipe5.heatPorts) 
                                      annotation (Line(
      points={{-52,32},{-19.9,32},{-19.9,25.2}},
      color={191,0,0},
      thickness=0.5));
  connect(boundary2.p_in, ramp.y) annotation (Line(
      points={{86,70},{91.4,70}},
      color={0,0,127},
      thickness=0.5));
  connect(ramp2.y, boundary5.p_in) annotation (Line(
      points={{83.4,-46},{72,-46}},
      color={0,0,127},
      thickness=0.5));
  connect(boundary6.ports[1], pipe11.port_a) 
                                         annotation (Line(
      points={{-90,-52},{-70,-52}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe10.port_b, boundary5.ports[1]) 
                                         annotation (Line(
      points={{32,-52},{50,-52}},
      color={0,127,255},
      thickness=0.5));
  connect(heat2.port, pipe9.heatPorts) 
                                      annotation (Line(
      points={{-52,-28},{-18.9,-28},{-18.9,-34.8}},
      color={191,0,0},
      thickness=0.5));
  connect(pipe11.port_b, muliPort11.port_a) annotation (Line(
      points={{-50,-52},{-46,-52}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(muliPort11.ports_b[2], pipe12.port_a) annotation (Line(
      points={{-38,-54},{-32,-54},{-32,-62},{-29,-62}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(muliPort11.ports_b[1], pipe9.port_a) annotation (Line(
      points={{-38,-50},{-32,-50},{-32,-40},{-29,-40}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(multiPort10.port_a, pipe10.port_a) annotation (Line(
      points={{8,-52},{12,-52}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(pipe9.port_b, multiPort10.ports_b[1]) annotation (Line(
      points={{-9,-40},{-6,-40},{-6,-50},{0,-50}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
  connect(pipe12.port_b, multiPort10.ports_b[2]) annotation (Line(
      points={{-9,-62},{-6,-62},{-6,-54},{0,-54}},
      color={0,127,255},
      thickness=0.5,
      smooth=Smooth.None));
end TestDistributedPipe01;
