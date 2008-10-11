within Modelica_Fluid.Test.TestComponents.Pumps;
model TestWaterPumpPowerCharacteristic
  "Test pump with power consumption characteristic"
  extends Modelica.Icons.Example;
  Modelica_Fluid.Sources.FixedBoundary_pTX Source(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater, p=1e5,
    T=system.T_ambient) 
  annotation (Placement(transformation(extent={{-100,20},{-80,40}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX Sink(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=5e5,
    T=system.T_ambient,
    usePressureInput=true) 
  annotation (Placement(transformation(extent={{34,26},{14,46}}, rotation=0)));
  Modelica_Fluid.Pumps.Pump pump(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow_start=1,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticFlow (
          q_nom={0,0.001,0.0015}, head_nom={100,50,0}),
    usePowerCharacteristic=true,
    redeclare function powerCharacteristic =
        Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticPower (
          q_nom={0,0.001,0.0015}, W_nom={550,650,800}),
    M=0.1,
    pin_start=100000,
    pout_start=700000)     annotation (Placement(transformation(extent={{-66,20},
            {-34,50}}, rotation=0)));
  Modelica.Blocks.Sources.Constant valveOpening(k=1) 
  annotation (Placement(transformation(extent={{-60,60},{-40,80}}, rotation=0)));
  Modelica_Fluid.ControlValves.ValveLinear Valve(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater,
    m_flow_nom=1,
    dp_nom=20000) 
  annotation (Placement(transformation(extent={{-16,26},{2,46}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp downstreamPressure(
    startTime=1,
    duration=5,
    offset=1e5,
    height=10e5) 
                annotation (Placement(transformation(extent={{4,74},{24,94}},
          rotation=0)));
  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent={{64,-4},
            {84,16}},          rotation=0)));
equation
  connect(valveOpening.y, Valve.opening)  annotation (Line(points={{-39,70},{-7,
          70},{-7,45}}, color={0,0,127}));
  connect(Valve.port_b,Sink. port)       annotation (Line(points={{2,36},{14,36}},
        color={0,127,255}));
  connect(Valve.port_a,pump. outlet)      annotation (Line(points={{-16,36},{
          -26,36},{-26,39.8},{-40.4,39.8}}, color={0,127,255}));
  connect(pump.inlet,Source. port)  annotation (Line(points={{-62.8,32},{-70,32},
          {-70,30},{-80,30}}, color={0,127,255}));
  connect(downstreamPressure.y, Sink.p_in) 
                                annotation (Line(points={{25,84},{58,84},{58,42},
          {36,42}}, color={0,0,127}));
  annotation (experiment(StopTime=10));
end TestWaterPumpPowerCharacteristic;
