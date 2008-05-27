within Modelica_Fluid.Test.TestComponents.Sensors;
model TestDensity
  import Modelica_Fluid;
  annotation (
    Diagram(graphics),
    experiment(Tolerance=1e-006),
    experimentSetupOutput);
  inner Modelica_Fluid.Ambient ambient annotation (Placement(transformation(
          extent={{-100,-100},{-80,-80}}, rotation=0)));
  Modelica_Fluid.Sensors.DensityTwoPort density2_1(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (Placement(
        transformation(extent={{-20,-30},{0,-10}}, rotation=0)));
  Modelica_Fluid.PressureLosses.SimpleGenericOrifice simpleGenericOrifice1(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    zeta=2,
    diameter=0.1) annotation (Placement(transformation(extent={{20,-30},{40,-10}}, 
          rotation=0)));
  Modelica.Blocks.Sources.Sine sine1 
                                    annotation (Placement(transformation(extent
          ={{-100,-20},{-80,0}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX massFlowRate2(
    useFlowRateInput=true,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    h=3200e3)                       annotation (Placement(transformation(extent
          ={{-60,-20},{-40,0}}, rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_phX boundary_fixed1(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=ambient.default_p_ambient,
    h=3000e3) annotation (Placement(transformation(extent={{100,-30},{80,-10}}, 
          rotation=0)));
  Modelica_Fluid.Sensors.DensityTwoPort density2_2(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (Placement(
        transformation(extent={{50,-30},{70,-10}}, rotation=0)));
  Modelica_Fluid.Sensors.DensityOnePort density1_1(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (Placement(
        transformation(extent={{-20,50},{0,70}}, rotation=0)));
  Modelica_Fluid.PressureLosses.SimpleGenericOrifice simpleGenericOrifice(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    zeta=2,
    diameter=0.1) annotation (Placement(transformation(extent={{20,40},{40,60}}, 
          rotation=0)));
  Modelica.Blocks.Sources.Sine sine annotation (Placement(transformation(extent
          ={{-100,50},{-80,70}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX massFlowRate1(
    useFlowRateInput=true,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    h=3200e3)                       annotation (Placement(transformation(extent
          ={{-60,50},{-40,70}}, rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_phX boundary_fixed(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=ambient.default_p_ambient,
    h=3000e3) annotation (Placement(transformation(extent={{100,40},{80,60}},
          rotation=0)));
  Modelica_Fluid.Sensors.DensityOnePort density1_2(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (Placement(
        transformation(extent={{50,50},{70,70}}, rotation=0)));
equation
  connect(sine1.y, massFlowRate2.m_flow_in) annotation (Line(
      points={{-79,-10},{-70,-10},{-70,-4},{-59.3,-4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowRate2.port, density2_1.port_a) annotation (Line(
      points={{-40,-10},{-30,-10},{-30,-20},{-20,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(density2_1.port_b, simpleGenericOrifice1.port_a) annotation (Line(
      points={{0,-20},{20,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(simpleGenericOrifice1.port_b, density2_2.port_a) annotation (Line(
      points={{40,-20},{50,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(density2_2.port_b, boundary_fixed1.port) annotation (Line(
      points={{70,-20},{80,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(sine.y, massFlowRate1.m_flow_in) annotation (Line(
      points={{-79,60},{-70,60},{-70,66},{-59.3,66}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowRate1.port, density1_1.port) annotation (Line(
      points={{-40,60},{-26,60},{-26,50},{-10,50}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(density1_1.port, simpleGenericOrifice.port_a) annotation (Line(
      points={{-10,50},{20,50}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(density1_2.port, boundary_fixed.port) annotation (Line(
      points={{60,50},{80,50}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(simpleGenericOrifice.port_b, density1_2.port) annotation (Line(
      points={{40,50},{60,50}},
      color={0,127,255},
      smooth=Smooth.None));
end TestDensity;
