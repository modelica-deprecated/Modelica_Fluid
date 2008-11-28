within Modelica_Fluid.Test.TestComponents.Pumps;
model TestWaterPumpDCMotor "Test pump with dc motor (startup transient)"
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
annotation (
  Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,
            100}}),
          graphics),
  experiment(StopTime=8, Tolerance=1e-006),
  Documentation(info=""));

  Modelica_Fluid.Sources.FixedBoundary_pTX Source(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    p=100000) 
  annotation (Placement(transformation(extent={{-100,2},{-80,22}},  rotation=0)));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX Sink(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    T=system.T_ambient,
    usePressureInput=false,
    p=100000) 
  annotation (Placement(transformation(extent={{34,50},{14,70}}, rotation=0)));
  Modelica_Fluid.Pumps.PumpShaft pump(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow_start=1,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticFlow (
          q_nominal={0,0.001,0.0015}, head_nominal={100,50,0}),
    N_nominal=1500,
    p_a_start=100000,
    p_b_start=700000)     annotation (Placement(transformation(extent={{-66,0},
            {-34,30}}, rotation=0)));
  Modelica.Blocks.Sources.Ramp valveOpening(
    height=-0.5,
    duration=-0.1,
    offset=1,
    startTime=5) 
  annotation (Placement(transformation(extent={{-40,74},{-20,94}}, rotation=0)));
  Modelica_Fluid.ControlValves.ValveIncompressible Valve(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater,
    m_flow_nominal=1,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    dp_nominal=1000000) 
  annotation (Placement(transformation(extent={{-16,50},{2,70}}, rotation=0)));
  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent={{64,-4},
            {84,16}},          rotation=0)));
  Modelica.Electrical.Machines.BasicMachines.DCMachines.DC_PermanentMagnet
    motor(
    La=1e-3,
    Jr=0.1,
    inertiaRotor(w(
        fixed=true,
        displayUnit="1/min",
        start=0.10471975511966)),
    VaNominal=400,
    wNominal(displayUnit="1/min") = 157.07963267949,
    IaNominal=10,
    Ra=10) annotation (Placement(transformation(extent={{2,0},{-22,24}})));
  Modelica.Electrical.Analog.Sources.StepVoltage stepVoltage(
    startTime=1,
    offset=0.1,
    V=400) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=-90,
        origin={32,8})));
  Modelica.Electrical.Analog.Basic.Ground ground 
    annotation (Placement(transformation(extent={{22,-30},{42,-8}})));
equation
  connect(Valve.port_b,Sink.ports[1])    annotation (Line(points={{2,60},{14,60}},
        color={0,127,255}));
  connect(Valve.port_a,pump.port_b)      annotation (Line(points={{-16,60},{-34,
          60},{-34,15}},                    color={0,127,255}));
  connect(pump.port_a,Source.ports[1]) 
                                     annotation (Line(points={{-66,15},{-70,15},
          {-80,12}},          color={0,127,255}));
  connect(valveOpening.y, Valve.stemPosition) annotation (Line(
      points={{-19,84},{-7,84},{-7,68}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pump.shaft, motor.flange) annotation (Line(
      points={{-34.8,12},{-22,12}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(motor.pin_ap, stepVoltage.p) annotation (Line(
      points={{-17.2,24},{-18,24},{-18,30},{32,30},{32,18}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(stepVoltage.n, motor.pin_an) annotation (Line(
      points={{32,-2},{16,-2},{16,24},{-2.8,24}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(stepVoltage.n, ground.p) annotation (Line(
      points={{32,-2},{32,-5},{32,-5},{32,-8}},
      color={0,0,255},
      smooth=Smooth.None));
end TestWaterPumpDCMotor;
