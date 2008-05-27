within Modelica_Fluid.Test.TestComponents.Sensors;
model TestTemperatureSensor "Test and compare case for the difference between using one port with
   and without explicit junction model and two port sensor for fluid temperature meassuring"
  import Modelica_Fluid;
  Modelica_Fluid.Sensors.TemperatureOnePort temperatureOnePort(redeclare
      package Medium = Modelica.Media.Water.StandardWater) 
    annotation (Placement(transformation(extent={{-20,40},{0,60}}, rotation=0)));
  Modelica_Fluid.Sensors.TemperatureTwoPort temperatureTwoPort(redeclare
      package Medium = Modelica.Media.Water.StandardWater) 
    annotation (Placement(transformation(extent={{-20,-20},{0,0}}, rotation=0)));
  inner Modelica_Fluid.Ambient ambient annotation (Placement(transformation(
          extent={{-100,-100},{-80,-80}}, rotation=0)));
  Modelica_Fluid.Volumes.OpenTank openTankCold2(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    pipe_diameters={0.05}) annotation (Placement(transformation(extent={{20,0},
            {40,20}}, rotation=0)));
  Modelica_Fluid.Volumes.OpenTank openTankCold1(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    pipe_diameters={0.05}) annotation (Placement(transformation(extent={{20,60},
            {40,80}}, rotation=0)));
  Modelica_Fluid.Volumes.OpenTank openTankHot1(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    T_start=SI.Conversions.from_degC(80),
    pipe_diameters={0.05}) annotation (Placement(transformation(extent={{60,40},
            {80,60}}, rotation=0)));
  Modelica_Fluid.Volumes.OpenTank openTankHot2(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    T_start=SI.Conversions.from_degC(80),
    pipe_diameters={0.05}) annotation (Placement(transformation(extent={{60,-20},
            {80,0}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate1(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    useFlowRateInput=true,
    T=SI.Conversions.from_degC(50)) annotation (Placement(transformation(extent
          ={{-60,30},{-40,50}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate2(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    useFlowRateInput=true,
    T=SI.Conversions.from_degC(50)) annotation (Placement(transformation(extent
          ={{-60,-20},{-40,0}}, rotation=0)));
  annotation (Diagram(graphics),
                       Documentation(info="<html>
<p align = justify>In that test model the behaviour of one port temperature sensors with and without explicit junction models and two port temperature sensor are compared. Therefor each sensor is connected to two tanks with different temperatures and a flow source with changing flow direction.<p>
<p align = justify>The one port sensor with explicit junction model and the two port sensor are showing the same expected results. The one port sensor without explicit junction model shows a different and unexpected behaviour. That test case shows that in all case where more are the sensor is between more than two components connected the one port sensor with explicit junction model or the two port sensor should be used!</p>
</html>"));
  Modelica.Blocks.Sources.Sine sine annotation (Placement(transformation(extent
          ={{-100,10},{-80,30}}, rotation=0)));
  Modelica_Fluid.Sensors.TemperatureOnePort temperatureOnePortJunction(
      redeclare package Medium = Modelica.Media.Water.StandardWater) 
    annotation (Placement(transformation(extent={{-20,-80},{0,-60}}, rotation=0)));
  Modelica_Fluid.Volumes.OpenTank openTankCold3(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    pipe_diameters={0.05}) annotation (Placement(transformation(extent={{20,-60},
            {40,-40}}, rotation=0)));
  Modelica_Fluid.Volumes.OpenTank openTankHot3(
    level_start=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    height=2,
    area=2,
    T_start=SI.Conversions.from_degC(80),
    pipe_diameters={0.05}) annotation (Placement(transformation(extent={{60,-80},
            {80,-60}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate3(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    useFlowRateInput=true,
    T=SI.Conversions.from_degC(50)) annotation (Placement(transformation(extent
          ={{-60,-90},{-40,-70}}, rotation=0)));
  Modelica_Fluid.Junctions.JunctionIdeal junctionIdeal(redeclare package Medium
      = Modelica.Media.Water.StandardWater) 
    annotation (Placement(transformation(extent={{20,-90},{40,-70}}, rotation=0)));
equation
  connect(massFlowRate2.port, temperatureTwoPort.port_a) annotation (Line(
      points={{-40,-10},{-20,-10}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(massFlowRate1.port, temperatureOnePort.port) annotation (Line(
      points={{-40,40},{-10,40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(sine.y, massFlowRate1.m_flow_in) annotation (Line(
      points={{-79,20},{-70,20},{-70,46},{-59.3,46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sine.y, massFlowRate2.m_flow_in) annotation (Line(
      points={{-79,20},{-70,20},{-70,-4},{-59.3,-4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowRate3.port, temperatureOnePortJunction.port) annotation (Line(
      points={{-40,-80},{-10,-80}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(sine.y, massFlowRate3.m_flow_in) annotation (Line(
      points={{-79,20},{-70,20},{-70,-74},{-59.3,-74}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(temperatureOnePortJunction.port, junctionIdeal.port_1) annotation (Line(
      points={{-10,-80},{20,-80}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(temperatureTwoPort.port_b, openTankCold2.ports[1]) annotation (Line(
        points={{0,-10},{0,-6},{30,-6},{30,-0.5}}, color={0,127,255}));
  connect(temperatureOnePort.port, openTankCold1.ports[1]) annotation (Line(
        points={{-10,40},{6,40},{6,56},{30,56},{30,59.5}}, color={0,127,255}));
  connect(temperatureOnePort.port, openTankHot1.ports[1]) annotation (Line(
        points={{-10,40},{30,40},{30,39.5},{70,39.5}}, color={0,127,255}));
  connect(temperatureTwoPort.port_b, openTankHot2.ports[1]) annotation (Line(
        points={{0,-10},{34,-10},{34,-20.5},{70,-20.5}}, color={0,127,255}));
  connect(junctionIdeal.port_3, openTankCold3.ports[1]) annotation (Line(points
        ={{30,-70},{30,-60.5}}, color={0,127,255}));
  connect(junctionIdeal.port_2, openTankHot3.ports[1]) annotation (Line(points=
          {{40,-80},{55.5,-80},{55.5,-80.5},{70,-80.5}}, color={0,127,255}));
end TestTemperatureSensor;
