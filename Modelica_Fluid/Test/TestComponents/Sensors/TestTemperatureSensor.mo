within Modelica_Fluid.Test.TestComponents.Sensors;
model TestTemperatureSensor "Test and compare case for the difference between using one port with
   and without explicit junction model and two port sensor for fluid temperature meassuring" 
  import Modelica_Fluid;
  Modelica_Fluid.Sensors.TemperatureOnePort temperatureOnePort(redeclare 
      package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[-20,40; 0,60]);
  Modelica_Fluid.Sensors.TemperatureTwoPort temperatureTwoPort(redeclare 
      package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[-20,-20; 0,0]);
  inner Modelica_Fluid.Ambient ambient annotation (extent=[-100,-100; -80,-80]);
  Modelica_Fluid.Volumes.OpenTank openTankCold2(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    pipe_diameters={0.05}) annotation (extent=[20,0; 40,20]);
  Modelica_Fluid.Volumes.OpenTank openTankCold1(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    pipe_diameters={0.05}) annotation (extent=[20,60; 40,80]);
  Modelica_Fluid.Volumes.OpenTank openTankHot1(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    T_start=SI.Conversions.from_degC(80),
    pipe_diameters={0.05}) annotation (extent=[60,40; 80,60]);
  Modelica_Fluid.Volumes.OpenTank openTankHot2(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    T_start=SI.Conversions.from_degC(80),
    pipe_diameters={0.05}) annotation (extent=[60,-20; 80,0]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate1(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    useFlowRateInput=true,
    T=SI.Conversions.from_degC(50)) annotation (extent=[-60,30; -40,50]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate2(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    useFlowRateInput=true,
    T=SI.Conversions.from_degC(50)) annotation (extent=[-60,-20; -40,0]);
  annotation (Diagram, Documentation(info="<html>
<p align = justify>In that test model the behaviour of one port temperature sensors with and without explicit junction models and two port temperature sensor are compared. Therefor each sensor is connected to two tanks with different temperatures and a flow source with changing flow direction.<p>
<p align = justify>The one port sensor with explicit junction model and the two port sensor are showing the same expected results. The one port sensor without explicit junction model shows a different and unexpected behaviour. That test case shows that in all case where more are the sensor is between more than two components connected the one port sensor with explicit junction model or the two port sensor should be used!</p>
</html>"));
  Modelica.Blocks.Sources.Sine sine annotation (extent=[-100,10; -80,30]);
  Modelica_Fluid.Sensors.TemperatureOnePort temperatureOnePortJunction(
      redeclare package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[-20,-80; 0,-60]);
  Modelica_Fluid.Volumes.OpenTank openTankCold3(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    pipe_diameters={0.05}) annotation (extent=[20,-60; 40,-40]);
  Modelica_Fluid.Volumes.OpenTank openTankHot3(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    T_start=SI.Conversions.from_degC(80),
    pipe_diameters={0.05}) annotation (extent=[60,-80; 80,-60]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate3(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    useFlowRateInput=true,
    T=SI.Conversions.from_degC(50)) annotation (extent=[-60,-90; -40,-70]);
  Modelica_Fluid.Junctions.JunctionIdeal junctionIdeal(redeclare package Medium
      = Modelica.Media.Water.StandardWater, p_start=ambient.default_p_ambient) 
    annotation (extent=[20,-90; 40,-70]);
equation 
  connect(massFlowRate2.port, temperatureTwoPort.port_a) annotation (points=[
        -40,-10; -20,-10], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(temperatureTwoPort.port_b, openTankCold2.port[1]) annotation (points=[0,-10; 
        16,-10; 16,0.1; 29.8,0.1],         style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(temperatureTwoPort.port_b, openTankHot2.port[1]) annotation (points=[
        0,-10; 36,-10; 36,-19.9; 69.8,-19.9], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(massFlowRate1.port, temperatureOnePort.port) annotation (points=[-40,
        40; -10,40], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(temperatureOnePort.port, openTankHot1.port[1]) annotation (points=[
        -10,40; 29.9,40; 29.9,40.1; 69.8,40.1], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(temperatureOnePort.port, openTankCold1.port[1]) annotation (points=[
        -10,40; 10,40; 10,60.1; 29.8,60.1], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(sine.y, massFlowRate1.m_flow_in) annotation (points=[-79,20; -70,20;
        -70,46; -59.3,46], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(sine.y, massFlowRate2.m_flow_in) annotation (points=[-79,20; -70,20;
        -70,-4; -59.3,-4], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(massFlowRate3.port, temperatureOnePortJunction.port) annotation (
      points=[-40,-80; -10,-80], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(sine.y, massFlowRate3.m_flow_in) annotation (points=[-79,20; -70,20;
        -70,-74; -59.3,-74], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(temperatureOnePortJunction.port, junctionIdeal.port_1) annotation (
      points=[-10,-80; 19,-80], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(junctionIdeal.port_2, openTankHot3.port[1]) annotation (points=[41,
        -80; 55.4,-80; 55.4,-79.9; 69.8,-79.9], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(junctionIdeal.port_3, openTankCold3.port[1]) annotation (points=[30,
        -69; 30,-59.9; 29.8,-59.9], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
end TestTemperatureSensor;
