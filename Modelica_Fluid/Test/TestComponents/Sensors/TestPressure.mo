within Modelica_Fluid.Test.TestComponents.Sensors;
model TestPressure
  import Modelica_Fluid;
  Modelica_Fluid.Sensors.Pressure pressure1(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (Placement(
        transformation(extent={{-20,0},{0,20}}, rotation=0)));
  Modelica_Fluid.PressureLosses.SimpleGenericOrifice simpleGenericOrifice(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    zeta=2,
    diameter=0.1) annotation (Placement(transformation(extent={{20,-10},{40,10}},
          rotation=0)));
  Modelica_Fluid.Sensors.RelativePressure relativePressure(redeclare package
      Medium = Modelica.Media.Water.StandardWater) 
    annotation (Placement(transformation(extent={{20,34},{40,54}}, rotation=0)));
  Modelica.Blocks.Sources.Sine sine annotation (Placement(transformation(extent=
           {{-100,0},{-80,20}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate1(
    useFlowRateInput=true,
    T=SI.Conversions.from_degC(50),
    redeclare package Medium = Modelica.Media.Water.StandardWater) 
                                    annotation (Placement(transformation(extent=
           {{-60,0},{-40,20}}, rotation=0)));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}),
            graphics),
    experiment(Tolerance=1e-006),
    experimentSetupOutput);
  Modelica_Fluid.Sources.FixedBoundary_phX boundary_fixed(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=system.p_ambient,
    h=3000e3) annotation (Placement(transformation(extent={{100,-10},{80,10}},
          rotation=0)));
  Modelica_Fluid.Sensors.Pressure pressure2(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (Placement(
        transformation(extent={{50,0},{70,20}}, rotation=0)));
  inner Modelica_Fluid.System system  annotation (Placement(transformation(
          extent={{-100,-100},{-80,-80}}, rotation=0)));
equation
  connect(sine.y, massFlowRate1.m_flow_in) annotation (Line(
      points={{-79,10},{-70,10},{-70,16},{-59.3,16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowRate1.ports[1], pressure1.port) 
                                              annotation (Line(
      points={{-40,10},{-26,10},{-26,0},{-10,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pressure1.port, simpleGenericOrifice.port_a) annotation (Line(
      points={{-10,0},{20,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(simpleGenericOrifice.port_a, relativePressure.port_a) annotation (Line(
      points={{20,0},{20,44},{20,44}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(relativePressure.port_b, simpleGenericOrifice.port_b) annotation (Line(
      points={{40,44},{40,23},{40,23},{40,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pressure2.port, boundary_fixed.ports[1]) 
                                               annotation (Line(
      points={{60,0},{80,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(simpleGenericOrifice.port_b, pressure2.port) annotation (Line(
      points={{40,0},{60,0}},
      color={0,127,255},
      smooth=Smooth.None));
end TestPressure;
