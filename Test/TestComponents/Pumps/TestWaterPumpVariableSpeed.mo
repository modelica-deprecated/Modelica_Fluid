within Modelica_Fluid.Test.TestComponents.Pumps;
model TestWaterPumpVariableSpeed
  "Test pump with variable speed (starting from zero)"
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
    duration=5,
    height=1500,
    offset=0)   annotation (Placement(transformation(extent={{-100,62},{-80,82}},
          rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX Source(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=100000) 
  annotation (Placement(transformation(extent={{-100,20},{-80,40}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX Sink(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    usePressureInput=false,
    p=100000) 
  annotation (Placement(transformation(extent={{34,26},{14,46}}, rotation=0)));
  Modelica_Fluid.Pumps.Pump pump(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow_start=1,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticFlow (
          q_nominal={0,0.001,0.0015}, head_nominal={100,50,0}),
    N_nominal=1500,
    pin_start=100000,
    pout_start=700000,
    use_N_input=true)      annotation (Placement(transformation(extent={{-66,20},
            {-34,50}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp valveOpening(
    height=-1,
    duration=1,
    offset=1,
    startTime=8) 
  annotation (Placement(transformation(extent={{-40,64},{-20,84}}, rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible Valve(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater,
    m_flow_nominal=1,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    dp_nominal=1000000) 
  annotation (Placement(transformation(extent={{-16,26},{2,46}}, rotation=0)));
  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent={{64,-4},
            {84,16}},          rotation=0)));
equation
  connect(Valve.port_b,Sink. port)       annotation (Line(points={{2,36},{14,36}},
        color={0,127,255}));
  connect(Valve.port_a,pump. outlet)      annotation (Line(points={{-16,36},{
          -26,36},{-26,39.8},{-40.4,39.8}}, color={0,127,255}));
  connect(pump.inlet,Source. port)  annotation (Line(points={{-62.8,32},{-70,32},
          {-70,30},{-80,30}}, color={0,127,255}));
  connect(N_pump.y, pump.N_in) annotation (Line(
      points={{-79,72},{-54.16,72},{-54.16,41.6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(valveOpening.y, Valve.stemPosition) annotation (Line(
      points={{-19,74},{-7,74},{-7,45}},
      color={0,0,127},
      smooth=Smooth.None));
end TestWaterPumpVariableSpeed;
