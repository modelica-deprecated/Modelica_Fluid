within Modelica_Fluid.Test.TestComponents.Volumes;
model TestMixingVolumes
  extends Modelica.Icons.Example;
  // package Medium = Modelica.Media.Water.StandardWater;
  package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater;
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                      graphics),
                       experiment(StopTime=10));
  Modelica_Fluid.Volumes.MixingVolume MixingVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (Placement(transformation(extent={{-30,30},{-10,50}},
          rotation=0)));

  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource2(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (Placement(transformation(extent={{-100,30},{-80,
            50}}, rotation=0)));
  Modelica_Fluid.Volumes.MixingVolume MixingVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (Placement(transformation(extent={{10,30},{30,50}},
          rotation=0)));
  Modelica_Fluid.Sensors.TemperatureOnePort Tmix_in(
                                         redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{-60,50},{-40,70}}, rotation=0)));
  Modelica_Fluid.Sensors.TemperatureOnePort Tmix_out(
                                          redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{40,50},{60,70}}, rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_phX Sink2(   p=101325, redeclare package
      Medium = Medium,
    h=Medium.h_default) 
    annotation (Placement(transformation(extent={{100,30},{80,50}}, rotation=0)));
  inner Modelica_Fluid.System system 
    annotation (Placement(transformation(extent={{-100,-100},{-80,-80}},
          rotation=0)));
equation
  connect(MixingVolume1.port_b, MixingVolume2.port_a) annotation (Line(points={{-10,40},
          {9.8,40}},          color={0,127,255}));
  connect(FlowSource2.port, MixingVolume1.port_a) annotation (Line(points={{-80,40},
          {-30.2,40}},     color={0,127,255}));
  connect(FlowSource2.port, Tmix_in.port) annotation (Line(points={{-80,40},{
          -50,40},{-50,50}}, color={0,127,255}));
  connect(MixingVolume2.port_b, Sink2.port) annotation (Line(points={{30,40},{
          80,40}}, color={0,127,255}));
  connect(MixingVolume2.port_b, Tmix_out.port) annotation (Line(points={{30,40},
          {50,40},{50,50}}, color={0,127,255}));
end TestMixingVolumes;
