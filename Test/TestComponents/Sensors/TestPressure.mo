within Modelica_Fluid.Test.TestComponents.Sensors;
model TestPressure 
  import Modelica_Fluid;
  Modelica_Fluid.Sensors.Pressure pressure1(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (extent=[-20,0; 0,20]);
  Modelica_Fluid.PressureLosses.SimpleGenericOrifice simpleGenericOrifice(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    zeta=2,
    diameter=0.1) annotation (extent=[20,-10; 40,10]);
  Modelica_Fluid.Sensors.RelativePressure relativePressure(redeclare package 
      Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[20,34; 40,54]);
  Modelica.Blocks.Sources.Sine sine annotation (extent=[-100,0; -80,20]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate1(
    useFlowRateInput=true,
    T=SI.Conversions.from_degC(50),
    redeclare package Medium = Modelica.Media.Water.StandardWater) 
                                    annotation (extent=[-60,0; -40,20]);
  annotation (
    Diagram,
    experiment(Tolerance=1e-006),
    experimentSetupOutput);
  Modelica_Fluid.Sources.FixedBoundary_phX boundary_fixed(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=ambient.default_p_ambient,
    h=3000e3) annotation (extent=[100,-10; 80,10]);
  Modelica_Fluid.Sensors.Pressure pressure2(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (extent=[50,0; 70,20]);
  inner Modelica_Fluid.Ambient ambient annotation (extent=[-100,-100; -80,-80]);
equation 
  connect(sine.y, massFlowRate1.m_flow_in) annotation (points=[-79,10; -70,10;
        -70,16; -59.3,16], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(massFlowRate1.port, pressure1.port) annotation (points=[-40,10; -26,
        10; -26,0; -10,0], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(pressure1.port, simpleGenericOrifice.port_a) annotation (points=[-10,
        0; 20,0], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(simpleGenericOrifice.port_a, relativePressure.port_a) annotation (
      points=[20,0; 20,44; 20,44], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(relativePressure.port_b, simpleGenericOrifice.port_b) annotation (
      points=[40,43.8; 40,23; 40,23; 40,0],
                                          style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(pressure2.port, boundary_fixed.port) annotation (points=[60,0; 80,0],
      style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(simpleGenericOrifice.port_b, pressure2.port) annotation (points=[40,0;
        60,0], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
end TestPressure;
