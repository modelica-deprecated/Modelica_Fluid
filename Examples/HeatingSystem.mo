within Modelica_Fluid.Examples;
model HeatingSystem "Simple model of a heating system"
   replaceable package Medium = Modelica.Media.Water.StandardWater 
     constrainedby Modelica.Media.Interfaces.PartialMedium;

  Modelica_Fluid.Vessels.OpenTank tank(
    redeclare package Medium = Medium,
    neglectPortDiameters=true,
    crossArea=0.01,
    V0=0.01,
    height=2,
    level_start=1,
    nPorts=2,
    portDiameters={0.025,0.025}) 
              annotation (Placement(transformation(extent={{-76,6},{-54,28}},
          rotation=0)));
  Modelica_Fluid.Machines.ControlledPump pump(
    redeclare package Medium = Medium,
    N_nominal=1500,
    N_const=1500,
    p_a_start=1.1e5,
    p_b_start=4.0e5,
    use_T_start=true,
    T_start=Modelica.SIunits.Conversions.from_degC(40),
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Machines.BaseClasses.PumpCharacteristics.linearFlow (
          head_nominal={60.0,0}, q_nominal={0.0,0.02e-3}),
    m_flow_start=0.01,
    checkValve=false) 
    annotation (Placement(transformation(extent={{-58,-16},{-38,4}}, rotation=
           0)));
  Modelica_Fluid.Valves.ValveIncompressible valve(
    redeclare package Medium = Medium,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    m_flow_nominal=0.01,
    compute_T=true,
    dp_nominal=300000) 
    annotation (Placement(transformation(extent={{42,-12},{58,4}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput circuitFlowRate 
    annotation (Placement(transformation(extent={{88,12},{108,32}}, rotation=
            0)));
  Sensors.MassFlowRate massFlowRate(redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{-34,6},{-14,-14}}, rotation=
           0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature ambientTemperature(
                                                                    T=system.T_ambient) 
    annotation (Placement(transformation(extent={{-12,-40},{2,-26}}, rotation=
           0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1(
                                                                   G=1.6e3/20) 
    annotation (Placement(transformation(
        origin={18,-48},
        extent={{8,-10},{-8,10}},
        rotation=90)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow burner(
                                                     Q_flow=1.6e3) 
    annotation (Placement(transformation(extent={{-2,12},{18,32}}, rotation=0)));
  inner Modelica_Fluid.System system(initType=Modelica_Fluid.Types.Init.SteadyState) 
                        annotation (Placement(transformation(extent={{-88,70},{
            -68,90}},   rotation=0)));
  Modelica_Fluid.Pipes.LumpedPipe pipe(
    redeclare package Medium = Medium,
    use_T_start=true,
    diameter=0.03,
    T_start=Modelica.SIunits.Conversions.from_degC(80),
    length=2,
    redeclare model HeatTransfer = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_ideal,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.LinearPressureLoss (
          m_flow_nominal=1, dp_nominal=100),
    p_a_start=400000,
    p_b_start=390000) 
    annotation (Placement(transformation(extent={{12,-14},{32,6}}, rotation=0)));

  Modelica_Fluid.Pipes.LumpedPipe radiator(
    use_T_start=true,
    redeclare package Medium = Medium,
    length=10,
    diameter=0.05,
    T_start=Modelica.SIunits.Conversions.from_degC(40),
    redeclare model HeatTransfer = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_ideal,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.LinearPressureLoss (
          m_flow_nominal=1, dp_nominal=100),
    p_a_start=110000,
    p_b_start=105000) 
    annotation (Placement(transformation(extent={{28,-76},{8,-56}}, rotation=
            0)));

  Modelica.Blocks.Interfaces.RealOutput hotWaterTemperature 
    annotation (Placement(transformation(extent={{88,-28},{108,-8}}, rotation=
           0)));
  Modelica.Blocks.Interfaces.RealOutput coldWaterTemperature 
    annotation (Placement(transformation(extent={{88,-78},{108,-58}},
          rotation=0)));
  Modelica_Fluid.Sensors.TemperatureOnePort sensor_T_1(
                                 redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{56,-56},{36,-36}}, rotation=
           0)));
  Modelica_Fluid.Sensors.TemperatureOnePort sensor_T_2(
                                 redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{-16,-56},{-36,-36}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput tankLevel 
                                 annotation (Placement(transformation(extent=
            {{90,60},{110,80}}, rotation=0)));
  Modelica.Blocks.Sources.RealExpression realExpression 
    annotation (Placement(transformation(extent={{-74,22},{-54,42}}, rotation=
           0)));
  Modelica.Blocks.Sources.Step valveOpening(
    height=0.1,
    startTime=2000,
    offset=0.9)   annotation (Placement(transformation(extent={{20,60},{40,80}},
                  rotation=0)));
equation
tankLevel = tank.level;
  connect(pump.port_b, massFlowRate.port_a) annotation (Line(points={{-38,-6},
          {-38,-3.4},{-34,-3.4},{-34,-4}},       color={0,127,255}));
  connect(massFlowRate.m_flow, circuitFlowRate) annotation (Line(points={{-24,
          7},{-24,38},{36,38},{36,22},{98,22}}, color={0,0,127}));
  connect(massFlowRate.port_b, pipe.port_a) annotation (Line(points={{-14,-4},
          {12,-4}}, color={0,127,255}));
  connect(pipe.port_b, valve.port_a) annotation (Line(points={{32,-4},{42,
          -4}},
        color={0,127,255}));
  connect(thermalConductor1.port_b, radiator.heatPort) annotation (Line(
        points={{18,-56},{18,-60.6}}, color={191,0,0}));
  connect(burner.port, pipe.heatPort) annotation (Line(points={{18,22},{22,22},
          {22,1.4}},     color={191,0,0}));
  connect(ambientTemperature.port, thermalConductor1.port_a) annotation (Line(
        points={{2,-33},{18,-33},{18,-40}}, color={191,0,0}));
  connect(sensor_T_1.T, hotWaterTemperature) annotation (Line(points={{39,-46},
          {39,-80},{72,-80},{72,-18},{98,-18}}, color={0,0,127}));
  connect(sensor_T_2.T, coldWaterTemperature) annotation (Line(points={{-33,
          -46},{-33,-88},{76,-88},{76,-68},{98,-68}}, color={0,0,127}));
  connect(radiator.port_a, valve.port_b) annotation (Line(points={{28,-66},
          {68,-66},{68,-4},{58,-4}},color={0,127,255}));
  connect(sensor_T_2.port, radiator.port_b) annotation (Line(points={{-26,-56},
          {-26,-66},{8,-66}}, color={0,127,255}));
  connect(radiator.port_a, sensor_T_1.port) annotation (Line(points={{28,-66},
          {46,-66},{46,-56}}, color={0,127,255}));
  connect(radiator.port_b, tank.ports[1]) annotation (Line(
      points={{8,-66},{-65,-66},{-65,8.2}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(tank.ports[2], pump.port_a) annotation (Line(
      points={{-65,3.8},{-65,-6},{-58,-6}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(valveOpening.y, valve.opening) annotation (Line(
      points={{41,70},{50,70},{50,2.4}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics),
                       Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={Rectangle(extent={{-100,100},{100,-100}},
            lineColor={0,0,255}), Text(
          extent={{-60,60},{60,-60}},
          lineColor={0,0,255},
          textString="H")}), Documentation(info="<html>
Simple heating system with a closed flow cycle. 
It is initialized in steady-state, achieved with the global system.initType. 
After 2000s of simulation time the valve fully opens.
</html>
"), experiment(StopTime=6000));

end HeatingSystem;
