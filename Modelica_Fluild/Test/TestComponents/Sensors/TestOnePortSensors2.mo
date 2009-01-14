within Modelica_Fluid.Test.TestComponents.Sensors;
model TestOnePortSensors2
  import Modelica_Fluid;
  package Medium = Modelica.Media.Water.StandardWater;
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics),
                       experiment(StopTime=15, Algorithm="Euler"),
    experimentSetupOutput);
  Modelica_Fluid.Vessels.Volume MixingVolume1(
    V=1e-3,
    use_T_start=false,
    redeclare package Medium = Medium,
    h_start=1e5,
    nPorts=2,
    use_portDiameters=false) 
                 annotation (Placement(transformation(extent={{-34,30},{-14,50}},
          rotation=0)));

  Modelica_Fluid.Sources.MassFlowSource_h FlowSource2(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium,
    useFlowRateInput=true) 
                   annotation (Placement(transformation(extent={{-68,30},{-48,
            50}}, rotation=0)));
  Modelica_Fluid.Vessels.Volume MixingVolume2(
    V=1e-3,
    use_T_start=false,
    redeclare package Medium = Medium,
    h_start=1.5e5,
    nPorts=2,
    use_portDiameters=false) 
                 annotation (Placement(transformation(extent={{32,30},{52,50}},
          rotation=0)));
  Modelica_Fluid.Sources.Boundary_ph Sink2(             redeclare package
      Medium = Medium,
    p=101325,
    h=5e4) 
    annotation (Placement(transformation(extent={{100,30},{80,50}}, rotation=0)));
  inner Modelica_Fluid.System system 
    annotation (Placement(transformation(extent={{-100,-100},{-80,-80}},
          rotation=0)));
  Modelica.Blocks.Sources.Ramp ramp(
    height=2,
    offset=-1,
    duration=10) 
               annotation (Placement(transformation(extent={{-100,30},{-80,50}},
          rotation=0)));
  Modelica_Fluid.Vessels.Volume MixingVolume3(
    V=1e-3,
    use_T_start=false,
    redeclare package Medium = Medium,
    h_start=1e5,
    nPorts=2,
    use_portDiameters=false) 
                 annotation (Placement(transformation(extent={{-34,-30},{-14,
            -10}}, rotation=0)));
  Modelica_Fluid.Sources.MassFlowSource_h FlowSource1(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium,
    useFlowRateInput=true) 
                   annotation (Placement(transformation(extent={{-68,-30},{-48,
            -10}}, rotation=0)));
  Modelica_Fluid.Vessels.Volume MixingVolume4(
    V=1e-3,
    use_T_start=false,
    redeclare package Medium = Medium,
    h_start=1.5e5,
    nPorts=2,
    use_portDiameters=false) 
                 annotation (Placement(transformation(extent={{32,-30},{52,-10}},
          rotation=0)));
  Modelica_Fluid.Sources.Boundary_ph Sink1(             redeclare package
      Medium = Medium,
    p=101325,
    h=5e4) 
    annotation (Placement(transformation(extent={{100,-30},{80,-10}}, rotation=
            0)));
  Modelica_Fluid.Sensors.TemperatureTwoPort Tmix2(redeclare package Medium = 
        Medium) annotation (Placement(transformation(extent={{0,-30},{20,-10}},
          rotation=0)));
equation
  connect(ramp.y, FlowSource2.m_flow_in) annotation (Line(
      points={{-79,40},{-74,40},{-74,48},{-68,48}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp.y, FlowSource1.m_flow_in) annotation (Line(
      points={{-79,40},{-76,40},{-76,-12},{-68,-12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(FlowSource2.ports[1], MixingVolume1.ports[1]) annotation (Line(
      points={{-48,40},{-36,40},{-36,32},{-24,32}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(MixingVolume1.ports[2], MixingVolume2.ports[2]) annotation (Line(
      points={{-24,28},{42,28}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(MixingVolume2.ports[1], Sink2.ports[1]) annotation (Line(
      points={{42,32},{61,32},{61,40},{80,40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(FlowSource1.ports[1], MixingVolume3.ports[1]) annotation (Line(
      points={{-48,-20},{-36,-20},{-36,-28},{-24,-28}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(MixingVolume3.ports[2], Tmix2.port_a) annotation (Line(
      points={{-24,-32},{-12,-32},{-12,-20},{0,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(Tmix2.port_b, MixingVolume4.ports[2]) annotation (Line(
      points={{20,-20},{32,-20},{32,-32},{42,-32}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(MixingVolume4.ports[1], Sink1.ports[1]) annotation (Line(
      points={{42,-28},{61,-28},{61,-20},{80,-20}},
      color={0,127,255},
      smooth=Smooth.None));
end TestOnePortSensors2;
