within Modelica_Fluid.Test.TestComponents.Pumps;
model TestWaterPumpCheckValve "Test pump with check valve behaviour"
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
annotation (
  Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
            100}}),
          graphics),
  experiment(StopTime=3, Tolerance=1e-006),
  Documentation(info=""));

  Modelica.Blocks.Sources.Ramp N_pump(
    startTime=1,
    duration=1,
    offset=1500,
    height=-1500) 
                annotation (Placement(transformation(extent={{-96,6},{-76,26}},
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
    p=500000) 
  annotation (Placement(transformation(extent={{62,-30},{42,-10}},
                                                                 rotation=0)));
  Modelica_Fluid.Pumps.Pump pump(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow_start=1,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticFlow (
          q_nom={0,0.001,0.0015}, head_nom={100,50,0}),
    N_nom=1500,
    use_N_input=true,
    checkValve=true,
    pin_start=100000,
    pout_start=700000)     annotation (Placement(transformation(extent={{-46,-40},
            {-14,-10}},rotation=0)));
  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent={{80,60},
            {100,80}},         rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible valve(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    m_flow_nom=1,
    dp_nom=800000) 
    annotation (Placement(transformation(extent={{0,-30},{22,-10}})));
  Modelica.Blocks.Sources.Constant valveOpening(k=1) 
                annotation (Placement(transformation(extent={{-38,32},{-18,52}},
          rotation=0)));
equation
  connect(pump.inlet,Source. port)  annotation (Line(points={{-42.8,-28},{-42.8,
          -28},{-80,-28}},    color={0,127,255}));
  connect(N_pump.y, pump.N_in) annotation (Line(
      points={{-75,16},{-34.16,16},{-34.16,-18.4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pump.outlet, valve.port_a) annotation (Line(
      points={{-20.4,-20.2},{-10.2,-20.2},{-10.2,-20},{0,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(valve.port_b, Sink.port) annotation (Line(
      points={{22,-20},{42,-20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(valveOpening.y, valve.stemPosition) annotation (Line(
      points={{-17,42},{11,42},{11,-11}},
      color={0,0,127},
      smooth=Smooth.None));
end TestWaterPumpCheckValve;
