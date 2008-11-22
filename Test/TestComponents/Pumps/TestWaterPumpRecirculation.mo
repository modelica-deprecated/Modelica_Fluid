within Modelica_Fluid.Test.TestComponents.Pumps;
model TestWaterPumpRecirculation
  "Test pump with variable speed and recirculating flow"
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
annotation (
  Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
            100}}),
          graphics),
  experiment(StopTime=10, Tolerance=1e-006),
  Documentation(info=""));

  Modelica.Blocks.Sources.Ramp N_pump(
    startTime=1,
    height=1500,
    offset=0,
    duration=1) annotation (Placement(transformation(extent={{-100,-8},{-80,12}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX Source(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=100000) 
  annotation (Placement(transformation(extent={{-100,-38},{-80,-18}},
                                                                    rotation=0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX Sink(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    usePressureInput=false,
    p=100000) 
  annotation (Placement(transformation(extent={{62,-30},{42,-10}},
                                                                 rotation=0)));
  Modelica_Fluid.Pumps.Pump pump(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow_start=1,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticFlow (
          q_nominal={0,0.001,0.0015}, head_nominal={100,50,0}),
    N_nominal=1500,
    use_N_input=true,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    M=0.1,
    pin_start=100000,
    pout_start=700000)     annotation (Placement(transformation(extent={{-46,-40},
            {-14,-10}},rotation=0)));
  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent={{80,60},
            {100,80}},         rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible V1(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    m_flow_nominal=1,
    dp_nominal=800000) 
    annotation (Placement(transformation(extent={{-28,6},{-50,26}})));
  Modelica_Fluid.ControlValves.ValveIncompressible V2(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    m_flow_nominal=1,
    dp_nominal=800000) 
    annotation (Placement(transformation(extent={{0,-30},{22,-10}})));
  Modelica.Blocks.Sources.Ramp V1_Opening(
    duration=1,
    height=-1,
    offset=1,
    startTime=5) 
                annotation (Placement(transformation(extent={{-100,40},{-80,60}},
          rotation=0)));
  Modelica.Blocks.Sources.Ramp V2_Opening(
    offset=0,
    height=1,
    duration=1,
    startTime=3) 
                annotation (Placement(transformation(extent={{-100,70},{-80,90}},
          rotation=0)));
equation
  connect(pump.inlet,Source. port)  annotation (Line(points={{-42.8,-28},{-42.8,
          -28},{-80,-28}},    color={0,127,255}));
  connect(N_pump.y, pump.N_in) annotation (Line(
      points={{-79,2},{-34.16,2},{-34.16,-18.4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pump.outlet, V1.port_a) annotation (Line(
      points={{-20.4,-20.2},{-20,14},{-20,16},{-28,16}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(V1.port_b, pump.inlet) annotation (Line(
      points={{-50,16},{-60,16},{-60,-28},{-42.8,-28}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pump.outlet, V2.port_a) annotation (Line(
      points={{-20.4,-20.2},{-10.2,-20.2},{-10.2,-20},{0,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(V2.port_b, Sink.port) annotation (Line(
      points={{22,-20},{42,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(V1_Opening.y, V1.stemPosition) annotation (Line(
      points={{-79,50},{-39,50},{-39,25}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(V2_Opening.y, V2.stemPosition) annotation (Line(
      points={{-79,80},{11,80},{11,-11}},
      color={0,0,127},
      smooth=Smooth.None));
end TestWaterPumpRecirculation;
