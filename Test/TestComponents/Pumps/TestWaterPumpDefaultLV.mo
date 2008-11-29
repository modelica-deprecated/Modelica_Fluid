within Modelica_Fluid.Test.TestComponents.Pumps;
model TestWaterPumpDefaultLV "Test pump with default options an linear valve"
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
annotation (
  Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
            100}}),
          graphics),
  experiment(StopTime=10, Tolerance=1e-006),
  Documentation(info=""));
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
          q_nominal={0,0.001,0.0015}, head_nominal={100,50,0}),
    p_a_start=100000,
    p_b_start=700000)     annotation (Placement(transformation(extent={{-66,20},
            {-34,50}}, rotation=0)));
  Modelica.Blocks.Sources.Constant valveOpening(k=1) 
  annotation (Placement(transformation(extent={{-60,60},{-40,80}}, rotation=0)));
  Modelica_Fluid.ControlValves.ValveLinear Valve(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater,
    m_flow_nominal=1,
    dp_nominal=20000) 
  annotation (Placement(transformation(extent={{-16,26},{2,46}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp downstreamPressure(
    startTime=1,
    duration=5,
    offset=1e5,
    height=10e5) 
                annotation (Placement(transformation(extent={{4,74},{24,94}},
          rotation=0)));

  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent=
            {{64,-4},{84,16}}, rotation=0)));
equation
  connect(Valve.port_b, Sink.ports[1])       annotation (Line(points={{2,36},{14,36}},
        color={0,127,255}));
  connect(Valve.port_a, pump.port_b)      annotation (Line(points={{-16,36},{
          -26,36},{-26,35},{-34,35}},       color={0,127,255}));
  connect(pump.port_a, Source.ports[1])  annotation (Line(points={{-66,35},{-70,35},
          {-70,30},{-80,30}}, color={0,127,255}));
  connect(downstreamPressure.y, Sink.p_in) 
                                annotation (Line(points={{25,84},{58,84},{58,42},
          {36,42}}, color={0,0,127}));
  connect(valveOpening.y, Valve.opening) annotation (Line(
      points={{-39,70},{-7,70},{-7,44}},
      color={0,0,127},
      smooth=Smooth.None));
end TestWaterPumpDefaultLV;
