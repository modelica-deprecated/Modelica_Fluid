within Modelica_Fluid.Test.TestComponents.Sensors;
model TestOnePortSensors2
  import Modelica_Fluid;
  package Medium = Modelica.Media.Water.StandardWater;
  annotation (Diagram(graphics),
                       experiment(StopTime=15, Algorithm="Euler"),
    experimentSetupOutput);
  Modelica_Fluid.Volumes.MixingVolume MixingVolume1(
    V=1e-3,
    use_T_start=false,
    redeclare package Medium = Medium,
    h_start=1e5) annotation (Placement(transformation(extent={{-34,30},{-14,50}}, 
          rotation=0)));

  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource2(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium,
    useFlowRateInput=true) 
                   annotation (Placement(transformation(extent={{-68,30},{-48,
            50}}, rotation=0)));
  Modelica_Fluid.Volumes.MixingVolume MixingVolume2(
    V=1e-3,
    use_T_start=false,
    redeclare package Medium = Medium,
    h_start=1.5e5) 
                 annotation (Placement(transformation(extent={{32,30},{52,50}}, 
          rotation=0)));
  Modelica_Fluid.Sensors.TemperatureOnePort Tmix1(redeclare package Medium = 
        Medium) 
    annotation (Placement(transformation(extent={{0,60},{20,80}}, rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_phX Sink2(             redeclare package
      Medium = Medium,
    p=101325,
    h=5e4) 
    annotation (Placement(transformation(extent={{100,30},{80,50}}, rotation=0)));
  inner Modelica_Fluid.Ambient ambient 
    annotation (Placement(transformation(extent={{-100,-100},{-80,-80}},
          rotation=0)));
  Modelica.Blocks.Sources.Ramp ramp(
    height=2,
    offset=-1,
    duration=10) 
               annotation (Placement(transformation(extent={{-100,30},{-80,50}}, 
          rotation=0)));
  Modelica_Fluid.Volumes.MixingVolume MixingVolume3(
    V=1e-3,
    use_T_start=false,
    redeclare package Medium = Medium,
    h_start=1e5) annotation (Placement(transformation(extent={{-34,-30},{-14,
            -10}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource1(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium,
    useFlowRateInput=true) 
                   annotation (Placement(transformation(extent={{-68,-30},{-48,
            -10}}, rotation=0)));
  Modelica_Fluid.Volumes.MixingVolume MixingVolume4(
    V=1e-3,
    use_T_start=false,
    redeclare package Medium = Medium,
    h_start=1.5e5) 
                 annotation (Placement(transformation(extent={{32,-30},{52,-10}}, 
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_phX Sink1(             redeclare package
      Medium = Medium,
    p=101325,
    h=5e4) 
    annotation (Placement(transformation(extent={{100,-30},{80,-10}}, rotation=
            0)));
  Modelica_Fluid.Sensors.TemperatureTwoPort Tmix2(redeclare package Medium = 
        Medium) annotation (Placement(transformation(extent={{0,-30},{20,-10}}, 
          rotation=0)));
equation
  connect(FlowSource2.port, MixingVolume1.port_a) annotation (Line(points={{-48,
          40},{-34.2,40}}, color={0,127,255}));
  connect(MixingVolume2.port_b, Sink2.port) annotation (Line(points={{52,40},{
          80,40}}, color={0,127,255}));
  connect(MixingVolume1.port_b, MixingVolume2.port_a) annotation (Line(
      points={{-14,40},{31.8,40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(Tmix1.port, MixingVolume1.port_b) annotation (Line(
      points={{10,60},{10,40},{-14,40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(ramp.y, FlowSource2.m_flow_in) annotation (Line(
      points={{-79,40},{-74,40},{-74,46},{-67.3,46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(FlowSource1.port,MixingVolume3. port_a) annotation (Line(points={{-48,
          -20},{-34.2,-20}}, color={0,127,255}));
  connect(MixingVolume4.port_b,Sink1. port) annotation (Line(points={{52,-20},{
          80,-20}}, color={0,127,255}));
  connect(ramp.y, FlowSource1.m_flow_in) annotation (Line(
      points={{-79,40},{-76,40},{-76,-14},{-67.3,-14}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Tmix2.port_a, MixingVolume3.port_b) annotation (Line(
      points={{0,-20},{-14,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(Tmix2.port_b, MixingVolume4.port_a) annotation (Line(
      points={{20,-20},{31.8,-20}},
      color={0,127,255},
      smooth=Smooth.None));
end TestOnePortSensors2;
