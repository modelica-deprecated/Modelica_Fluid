within Modelica_Fluid.Test.TestComponents.Pumps;
model TestWaterPumpNPSH "Test case for WaterPump"
  extends Modelica.Icons.Example;
  import PC = Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics;
  replaceable function pumpFlowChar = PC.quadraticFlow(q_nom={0,0.001,0.0015}, head_nom={100,50,0});
annotation (
  Diagram(graphics),
  experiment(StopTime=10, Tolerance=1e-006),
  Documentation(info=""));
  Modelica_Fluid.Sources.FixedBoundary_pTX Source(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater, p=1e5,
    T=ambient.default_T_ambient) 
  annotation (Placement(transformation(extent={{-100,20},{-80,40}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX Sink(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=5e5,
    T=ambient.default_T_ambient,
    usePressureInput=true) 
  annotation (Placement(transformation(extent={{34,26},{14,46}}, rotation=0)));
  Modelica_Fluid.Pumps.PumpNPSH Pump1(
    pin_start=1e5,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    redeclare function flowCharacteristic = pumpFlowChar,
    m_flow_start=1,
    pout_start=7e5)        annotation (Placement(transformation(extent={{-66,18},
            {-34,48}}, rotation=0)));
  Modelica.Blocks.Sources.Constant Constant1 
  annotation (Placement(transformation(extent={{-60,60},{-40,80}}, rotation=0)));
  Modelica_Fluid.ControlValves.ValveLinear Valve(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater, Kv=1e-5) 
  annotation (Placement(transformation(extent={{-16,26},{2,46}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp Ramp1(
    offset=5e5,
    startTime=1,
    duration=5,
    height=6e5) annotation (Placement(transformation(extent={{4,74},{24,94}},
          rotation=0)));

  inner Modelica_Fluid.Ambient ambient 
                                   annotation (Placement(transformation(extent=
            {{64,-4},{84,16}}, rotation=0)));
equation
  connect(Constant1.y, Valve.opening)     annotation (Line(points={{-39,70},{-7,
          70},{-7,45}}, color={0,0,127}));
  connect(Valve.port_b, Sink.port)       annotation (Line(points={{2,36},{14,36}},
        color={0,127,255}));
  connect(Valve.port_a, Pump1.outlet)     annotation (Line(points={{-16,36},{
          -26,36},{-26,37.8},{-40.4,37.8}}, color={0,127,255}));
  connect(Pump1.inlet, Source.port) annotation (Line(points={{-62.8,30},{-80,30}},
        color={0,127,255}));
  connect(Ramp1.y, Sink.p_in)   annotation (Line(points={{25,84},{58,84},{58,42},
          {36,42}}, color={0,0,127}));
end TestWaterPumpNPSH;
