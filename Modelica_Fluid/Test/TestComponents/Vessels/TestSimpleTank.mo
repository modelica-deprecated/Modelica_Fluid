within Modelica_Fluid.Test.TestComponents.Vessels;
model TestSimpleTank
  extends Modelica.Icons.Example;
  import Modelica_Fluid;
  Modelica_Fluid.Vessels.SimpleTank upperTank(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    nPorts=2,
    height=20,
    portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=0.1),
      Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=0.1)},
    level_start=2,
    crossArea=0.2) 
    annotation (Placement(transformation(extent={{-40,20},{0,60}}, rotation=0)));
  Modelica_Fluid.Sources.MassFlowSource_T massFlowRate(nPorts=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow=0.2,
    useFlowRateInput=true) 
    annotation (Placement(transformation(extent={{-60,-40},{-40,-20}}, rotation=
           0)));
  annotation (Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-160,-120},{100,120}},
        initialScale=0.1), graphics),
    experiment(StopTime=20000, Tolerance=1e-005),
    experimentSetupOutput(equdistant=false),
    Documentation(info="<html>
<p><b>Test case for open tank</b></p>
<p align=justify>The mass flow rate to the upper tank is controlled by the static pressure at the bootom of the upper tank. The fluid flows from the upper to the lower tank forced by pressure difference between the bootoms of both tanks. Increasing the simulation time leads to an error message, due to a full lower tank.</p>
</html>"));
  inner Modelica_Fluid.System system(energyDynamics=Modelica_Fluid.Types.Dynamics.FixedInitial) 
                                      annotation (Placement(transformation(
          extent={{-150,-112},{-130,-92}},  rotation=0)));
  Modelica_Fluid.Sensors.Pressure pressure(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (Placement(
        transformation(extent={{40,16},{60,36}}, rotation=0)));
  Modelica_Fluid.Pipes.StaticPipe pipe(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    diameter=0.02,
    length=200,
    height_ab=-20) 
                annotation (Placement(transformation(
        origin={0,-20},
        extent={{10,-10},{-10,10}},
        rotation=90)));

  Modelica_Fluid.Vessels.SimpleTank lowerTank(
    nPorts=1,
    height=20,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    level_start=2,
    crossArea=1,
    portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=0.1,
        height=6)}) 
    annotation (Placement(transformation(extent={{40,-60},{80,-20}}, rotation=0)));
  Modelica.Blocks.Logical.Hysteresis hysteresis(
    uLow=1.1e5,
    uHigh=2.5e5,
    pre_y_start=true) "mass flow rate signal by pressure control" 
    annotation (Placement(transformation(extent={{-140,-30},{-120,-10}},
          rotation=0)));
  Modelica.Blocks.Logical.Switch switch1 annotation (Placement(transformation(
          extent={{-100,-30},{-80,-10}}, rotation=0)));
  Modelica.Blocks.Sources.Constant m_flow_off(k=0) 
    annotation (Placement(transformation(extent={{-140,10},{-120,30}}, rotation=
           0)));
  Modelica.Blocks.Sources.Constant m_flow_on(k=2) 
    annotation (Placement(transformation(extent={{-140,-60},{-120,-40}},
          rotation=0)));
equation
  connect(massFlowRate.ports[1], upperTank.ports[1]) 
                                                 annotation (Line(
      points={{-40,-30},{-24,-30},{-24,20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pressure.p, hysteresis.u) annotation (Line(
      points={{61,26},{80,26},{80,70},{-160,70},{-160,-20},{-142,-20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(hysteresis.y, switch1.u2) annotation (Line(
      points={{-119,-20},{-102,-20}},
      color={255,0,255},
      smooth=Smooth.None));
  connect(m_flow_off.y, switch1.u1) annotation (Line(
      points={{-119,20},{-119,5},{-102,5},{-102,-12}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(m_flow_on.y, switch1.u3) annotation (Line(
      points={{-119,-50},{-110,-50},{-110,-28},{-102,-28}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(switch1.y, massFlowRate.m_flow_in) annotation (Line(
      points={{-79,-20},{-70,-20},{-70,-22},{-60,-22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(upperTank.ports[2], pipe.port_a) annotation (Line(
      points={{-16,20},{-18,20},{-18,10},{6.12323e-016,10},{6.12323e-016,-10}},
      color={0,127,255},
      smooth=Smooth.None));

  connect(pressure.port, pipe.port_a) annotation (Line(
      points={{50,16},{50,10},{6.12323e-016,10},{6.12323e-016,-10}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe.port_b, lowerTank.ports[1]) annotation (Line(
      points={{-6.12323e-016,-30},{0,-30},{0,-48},{60,-48},{60,-60}},
      color={0,127,255},
      smooth=Smooth.None));
end TestSimpleTank;
