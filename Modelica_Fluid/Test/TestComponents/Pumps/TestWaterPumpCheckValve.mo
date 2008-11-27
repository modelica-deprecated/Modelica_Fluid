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
                annotation (Placement(transformation(extent={{-80,0},{-60,20}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX Source(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=100000) 
  annotation (Placement(transformation(extent={{-80,-40},{-60,-20}},rotation=0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX Sink(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    usePressureInput=false,
    p=500000) 
  annotation (Placement(transformation(extent={{60,-40},{40,-20}},
                                                                 rotation=0)));
  Modelica_Fluid.Pumps.Pump pump(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow_start=1,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticFlow (
          q_nominal={0,0.001,0.0015}, head_nominal={100,50,0}),
    N_nominal=1500,
    checkValve=true,
    M=0.1,
    initType=Modelica_Fluid.Types.Init.SteadyState,
    p_a_start=100000,
    p_b_start=700000,
    use_N_input=true)      annotation (Placement(transformation(extent={{-40,-40},
            {-20,-20}},rotation=0)));
  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent={{80,60},
            {100,80}},         rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible valve(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    m_flow_nominal=1,
    dp_nominal=800000) 
    annotation (Placement(transformation(extent={{0,-40},{22,-20}})));
  Modelica.Blocks.Sources.Constant valveOpening(k=1) 
                annotation (Placement(transformation(extent={{-38,32},{-18,52}},
          rotation=0)));
equation
  connect(pump.port_a,Source. port)  annotation (Line(points={{-40,-30},{-60,
          -30}},              color={0,127,255}));
  connect(N_pump.y, pump.N_in) annotation (Line(
      points={{-59,10},{-30,10},{-30,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pump.port_b, valve.port_a) annotation (Line(
      points={{-20,-30},{0,-30}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(valve.port_b, Sink.port) annotation (Line(
      points={{22,-30},{40,-30}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(valveOpening.y, valve.stemPosition) annotation (Line(
      points={{-17,42},{11,42},{11,-22}},
      color={0,0,127},
      smooth=Smooth.None));
end TestWaterPumpCheckValve;
