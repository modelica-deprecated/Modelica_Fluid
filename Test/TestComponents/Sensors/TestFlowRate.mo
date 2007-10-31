within Modelica_Fluid.Test.TestComponents.Sensors;
model TestFlowRate 
  import Modelica_Fluid;
  Modelica_Fluid.PressureLosses.SimpleGenericOrifice simpleGenericOrifice(
    redeclare package Medium = Modelica.Media.Water.StandardWater, 
    zeta=2, 
    p_a_start=ambient.default_p_ambient, 
    p_b_start=ambient.default_p_ambient, 
    T_start=ambient.default_T_ambient, 
    diameter=0.1) annotation (extent=[-20,-10; 0,10]);
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
  inner Modelica_Fluid.Ambient ambient annotation (extent=[-100,-100; -80,-80]);
  Modelica_Fluid.Sensors.MassFlowRate massFlowRate(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (extent=[10,-10; 30,10]);
  Modelica_Fluid.Sensors.VolumeFlowRate volumeFlowRate(redeclare package Medium
      = Modelica.Media.Water.StandardWater) annotation (extent=[50,-10; 70,10]);
equation 
  connect(sine.y, massFlowRate1.m_flow_in) annotation (points=[-79,10; -70,10; 
        -70,16; -59.3,16], style(
      color=74, 
      rgbcolor={0,0,127}, 
      smooth=0));
  connect(massFlowRate1.port, simpleGenericOrifice.port_a) annotation (points=[
        -40,10; -32,10; -32,0; -20,0], style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(volumeFlowRate.port_b, boundary_fixed.port) annotation (points=[70,0; 
        80,0], style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(massFlowRate.port_b, volumeFlowRate.port_a) annotation (points=[30,0; 
        50,0], style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(simpleGenericOrifice.port_b, massFlowRate.port_a) annotation (points=
        [0,0; 10,0], style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
end TestFlowRate;
