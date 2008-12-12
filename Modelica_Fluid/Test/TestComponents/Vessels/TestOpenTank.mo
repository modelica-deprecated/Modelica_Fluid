within Modelica_Fluid.Test.TestComponents.Vessels;
model TestOpenTank
  extends Modelica.Icons.Example;
  import Modelica_Fluid;
  Modelica_Fluid.Vessels.OpenTank upperTank(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    nPorts=2,
    height=20,
    portDiameters={0.1,0.1},
    neglectPortDiameters=true,
    level_start=2,
    crossArea=0.2,
    V0=0.1) 
    annotation (Placement(transformation(extent={{-40,20},{0,60}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate(
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
  inner Modelica_Fluid.System system(initType=Modelica_Fluid.Types.Init.InitialValues) 
                                      annotation (Placement(transformation(
          extent={{-150,-112},{-130,-92}},  rotation=0)));
  Modelica_Fluid.Sensors.Pressure pressure(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (Placement(
        transformation(extent={{40,16},{60,36}}, rotation=0)));
  Modelica_Fluid.Pipes.StaticPipe pipe(
    height_ab=20,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.LaminarAndQuadraticTurbulentFlow,
    diameter=0.02,
    length=200) annotation (Placement(transformation(
        origin={0,-10},
        extent={{-10,-10},{10,10}},
        rotation=90)));

  Modelica_Fluid.Vessels.OpenTank lowerTank(
    nPorts=1,
    height=20,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    portDiameters={0.1},
    neglectPortDiameters=true,
    level_start=2,
    crossArea=1,
    V0=0.1) 
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
      points={{-40,-30},{-20,-30},{-20,24}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe.port_a, lowerTank.ports[1]) annotation (Line(
      points={{-6.12323e-016,-20},{-6.12323e-016,-70},{60,-70},{60,-60}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe.port_b, pressure.port) annotation (Line(
      points={{6.12323e-016,0},{50,0},{50,16}},
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
  connect(pipe.port_b, upperTank.ports[2]) annotation (Line(
      points={{6.12323e-016,0},{-18,0},{-18,16},{-20,16}},
      color={0,127,255},
      smooth=Smooth.None));
end TestOpenTank;
