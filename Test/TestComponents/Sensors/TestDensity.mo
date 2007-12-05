within Modelica_Fluid.Test.TestComponents.Sensors;
model TestDensity 
  import Modelica_Fluid;
  annotation (
    Diagram,
    experiment(Tolerance=1e-006),
    experimentSetupOutput);
  inner Modelica_Fluid.Ambient ambient annotation (extent=[-100,-100; -80,-80]);
  Modelica_Fluid.Sensors.DensityTwoPort density2_1(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (extent=[-20,-30; 0,-10]);
  Modelica_Fluid.PressureLosses.SimpleGenericOrifice simpleGenericOrifice1(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    zeta=2,
    p_a_start=ambient.default_p_ambient,
    p_b_start=ambient.default_p_ambient,
    T_start=ambient.default_T_ambient,
    diameter=0.1,
    use_T_start=false,
    h_start=3200e3) annotation (extent=[20,-30; 40,-10]);
  Modelica.Blocks.Sources.Sine sine1 
                                    annotation (extent=[-100,-20; -80,0]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX massFlowRate2(
    useFlowRateInput=true,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    h=3200e3)                       annotation (extent=[-60,-20; -40,0]);
  Modelica_Fluid.Sources.FixedBoundary_phX boundary_fixed1(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=ambient.default_p_ambient,
    h=3000e3) annotation (extent=[100,-30; 80,-10]);
  Modelica_Fluid.Sensors.DensityTwoPort density2_2(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (extent=[50,-30; 70,-10]);
  Modelica_Fluid.Sensors.DensityOnePort density1_1(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (extent=[-20,50; 0,70]);
  Modelica_Fluid.PressureLosses.SimpleGenericOrifice simpleGenericOrifice(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    zeta=2,
    p_a_start=ambient.default_p_ambient,
    p_b_start=ambient.default_p_ambient,
    T_start=ambient.default_T_ambient,
    diameter=0.1,
    use_T_start=false,
    h_start=3200e3) annotation (extent=[20,40; 40,60]);
  Modelica.Blocks.Sources.Sine sine annotation (extent=[-100,50; -80,70]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX massFlowRate1(
    useFlowRateInput=true,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    h=3200e3)                       annotation (extent=[-60,50; -40,70]);
  Modelica_Fluid.Sources.FixedBoundary_phX boundary_fixed(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=ambient.default_p_ambient,
    h=3000e3) annotation (extent=[100,40; 80,60]);
  Modelica_Fluid.Sensors.DensityOnePort density1_2(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (extent=[50,50; 70,70]);
equation 
  connect(sine1.y, massFlowRate2.m_flow_in) annotation (points=[-79,-10; -70,
        -10; -70,-4; -59.3,-4], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(massFlowRate2.port, density2_1.port_a) annotation (points=[-40,-10;
        -30,-10; -30,-20; -20,-20], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(density2_1.port_b, simpleGenericOrifice1.port_a) annotation (points=[
        0,-20; 20,-20], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(simpleGenericOrifice1.port_b, density2_2.port_a) annotation (points=[
        40,-20; 50,-20], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(density2_2.port_b, boundary_fixed1.port) annotation (points=[70,-20;
        80,-20], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(sine.y, massFlowRate1.m_flow_in) annotation (points=[-79,60; -70,60;
        -70,66; -59.3,66], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(massFlowRate1.port, density1_1.port) annotation (points=[-40,60; -26,
        60; -26,50; -10,50], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(density1_1.port, simpleGenericOrifice.port_a) annotation (points=[-10,
        50; 20,50], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(density1_2.port, boundary_fixed.port) annotation (points=[60,50; 80,
        50], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(simpleGenericOrifice.port_b, density1_2.port) annotation (points=[40,
        50; 60,50], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
end TestDensity;
