within Modelica_Fluid.Test.TestComponents.Sensors;
model TestFlowRate
  import Modelica_Fluid;
  Modelica_Fluid.PressureLosses.SimpleGenericOrifice simpleGenericOrifice(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    zeta=2,
    diameter=0.1) annotation (Placement(transformation(extent={{-20,-10},{0,10}},
          rotation=0)));
  Modelica.Blocks.Sources.Sine sine annotation (Placement(transformation(extent=
           {{-100,0},{-80,20}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate1(
    useFlowRateInput=true,
    T=SI.Conversions.from_degC(50),
    redeclare package Medium = Modelica.Media.Water.StandardWater) 
                                    annotation (Placement(transformation(extent=
           {{-60,0},{-40,20}}, rotation=0)));
  annotation (
    Diagram(graphics),
    experiment(Tolerance=1e-006),
    experimentSetupOutput);
  Modelica_Fluid.Sources.FixedBoundary_phX boundary_fixed(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=system.p_ambient,
    h=3000e3) annotation (Placement(transformation(extent={{100,-10},{80,10}},
          rotation=0)));
  inner Modelica_Fluid.System system  annotation (Placement(transformation(
          extent={{-100,-100},{-80,-80}}, rotation=0)));
  Modelica_Fluid.Sensors.MassFlowRate massFlowRate(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (Placement(
        transformation(extent={{10,-10},{30,10}}, rotation=0)));
  Modelica_Fluid.Sensors.VolumeFlowRate volumeFlowRate(redeclare package Medium
      = Modelica.Media.Water.StandardWater) annotation (Placement(
        transformation(extent={{50,-10},{70,10}}, rotation=0)));
equation
  connect(sine.y, massFlowRate1.m_flow_in) annotation (Line(
      points={{-79,10},{-70,10},{-70,16},{-59.3,16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowRate1.port, simpleGenericOrifice.port_a) annotation (Line(
      points={{-40,10},{-32,10},{-32,0},{-20,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(volumeFlowRate.port_b, boundary_fixed.port) annotation (Line(
      points={{70,0},{80,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(massFlowRate.port_b, volumeFlowRate.port_a) annotation (Line(
      points={{30,0},{50,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(simpleGenericOrifice.port_b, massFlowRate.port_a) annotation (Line(
      points={{0,0},{10,0}},
      color={0,127,255},
      smooth=Smooth.None));
end TestFlowRate;
