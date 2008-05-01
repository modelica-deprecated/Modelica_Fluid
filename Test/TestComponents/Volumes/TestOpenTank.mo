within Modelica_Fluid.Test.TestComponents.Volumes;
model TestOpenTank 
  extends Modelica.Icons.Example;
  import Modelica_Fluid;
  Modelica_Fluid.Volumes.OpenTank upperTank(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    n_ports=2,
    height=20,
    pipe_diameters={0.1,0.1},
    p_static_at_port=true,
    level_start=2,
    area=0.2,
    V0=0.1) 
    annotation (extent=[-40,20; 0,60]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow=0.2,
    useFlowRateInput=true) 
    annotation (extent=[-60,-40; -40,-20]);
  annotation (Diagram,
    Coordsys(extent=[-160,-120; 100,120], scale=0.1),
    experiment(StopTime=20000, Tolerance=1e-005),
    experimentSetupOutput(equdistant=false),
    Documentation(info="<html>
<p><b>Test case for open tank</b></p>
<p align=justify>The mass flow rate to the upper tank is controlled by the static pressure at the bootom of the upper tank. The fluid flows from the upper to the lower tank forced by pressure difference between the bootoms of both tanks. Increasing the simulation time leads to an error message, due to a full lower tank.</p>
</html>"));
  inner Modelica_Fluid.Ambient ambient annotation (extent=[-160,-120; -140,-100]);
  Modelica_Fluid.Sensors.Pressure pressure(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) annotation (extent=[40,16; 60,36]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe(
    height_ab=20,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.LaminarAndQuadraticTurbulent,
    diameter=0.02,
    length=200) annotation (extent=[-10,-20; 10,0], rotation=90);
  
  Modelica_Fluid.Volumes.OpenTank lowerTank(
    n_ports=1,
    height=20,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    pipe_diameters={0.1},
    p_static_at_port=true,
    level_start=2,
    area=1,
    V0=0.1) 
    annotation (extent=[40,-60; 80,-20]);
  Modelica.Blocks.Logical.Hysteresis hysteresis(
    uLow=1.1e5,
    uHigh=2.5e5,
    pre_y_start=true) "mass flow rate signal by pressure control" 
    annotation (extent=[-140,-30; -120,-10]);
  Modelica.Blocks.Logical.Switch switch1 annotation (extent=[-100,-30; -80,-10]);
  Modelica.Blocks.Sources.Constant m_flow_off(k=0) 
    annotation (extent=[-140,10; -120,30]);
  Modelica.Blocks.Sources.Constant m_flow_on(k=2) 
    annotation (extent=[-140,-60; -120,-40]);
equation 
  connect(massFlowRate.port, upperTank.ports[1]) annotation (points=[-40,-30;
        -20,-30; -20,17], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(upperTank.ports[2], pipe.port_b) annotation (points=[-20,21; -18,21; 
        -18,0; 6.12323e-016,0], style(
      color=3,
      rgbcolor={0,0,255},
      smooth=0));
  connect(pipe.port_a, lowerTank.ports[1]) annotation (points=[-6.12323e-016,
        -20; -6.12323e-016,-70; 60,-70; 60,-61], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(pipe.port_b, pressure.port) annotation (points=[6.12323e-016,0; 50,0; 
        50,16], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(pressure.p, hysteresis.u) annotation (points=[61,26; 80,26; 80,70;
        -160,70; -160,-20; -142,-20], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(hysteresis.y, switch1.u2) annotation (points=[-119,-20; -102,-20],
      style(
      color=5,
      rgbcolor={255,0,255},
      smooth=0));
  connect(m_flow_off.y, switch1.u1) annotation (points=[-119,20; -119,5; -102,5;
        -102,-12], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(m_flow_on.y, switch1.u3) annotation (points=[-119,-50; -110,-50; -110,
        -28; -102,-28], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(switch1.y, massFlowRate.m_flow_in) annotation (points=[-79,-20; -70,
        -20; -70,-24; -59.3,-24], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
end TestOpenTank;
