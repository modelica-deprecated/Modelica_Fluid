model TestTank 
  extends Modelica.Icons.Example;
  import Modelica_Fluid;
  Modelica_Fluid.Volumes.Tank tank(
    nTopPorts=2,
    levelMax=10,
    V0=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    level_start=1,
    portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.01,
        portLevel=9),Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=
        0.01, portLevel=6),Modelica_Fluid.Volumes.BaseClasses.TankPortData(
        diameter=0.01, portLevel=4),
        Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.01,
        portLevel=2)},
    area=0.2) 
    annotation (extent=[-40,40; 0,80]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate[2](redeclare 
      package Medium = Modelica.Media.Water.StandardWater, each m_flow=0.75) 
    annotation (extent=[-82,70; -62,90]);
  annotation (Diagram,
    experiment(StopTime=5000, Tolerance=1e-005),
    experimentSetupOutput(equdistant=false));
  inner Modelica_Fluid.Ambient ambient annotation (extent=[-100,-100; -80,-80]);
  Modelica_Fluid.Sources.FixedBoundary Boundary_fixed(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=ambient.default_p_ambient,
    T=ambient.default_T_ambient) 
    annotation (extent=[100,-90; 80,-70], rotation=0);
  Modelica_Fluid.Sensors.MassFlowRate mFlow_9m(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) 
    "Mass flow rate out of the port at a lever of 9 m" 
    annotation (extent=[20,30; 40,50]);
  Modelica_Fluid.Sensors.MassFlowRate mFlow_6m(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) 
    "Mass flow rate out of the port at a lever of 6 m" 
    annotation (extent=[20,-10; 40,10]);
  Modelica_Fluid.Sensors.MassFlowRate mFlow_4m(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) 
    "Mass flow rate out of the port at a lever of 4 m" 
    annotation (extent=[20,-50; 40,-30]);
  Modelica_Fluid.Sensors.MassFlowRate mFlow_2m(redeclare package Medium = 
        Modelica.Media.Water.StandardWater) 
    "Mass flow rate out of the port at a lever of 2 m" 
    annotation (extent=[20,-90; 40,-70]);
equation 
  connect(massFlowRate.port, tank.topPorts) annotation (points=[-62,80; -40,80;
        -40,81; -20,81], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(tank.ports[1], mFlow_9m.port_a) annotation (points=[-20,39; 0,39; 0,
        40; 20,40], style(
      color=3,
      rgbcolor={0,0,255},
      smooth=0));
  connect(mFlow_9m.port_b, Boundary_fixed.port) annotation (points=[40,40; 60,
        40; 60,-80; 80,-80], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(tank.ports[2], mFlow_6m.port_a) annotation (points=[-20,39; 0,39; 0,0;
        20,0], style(
      color=3,
      rgbcolor={0,0,255},
      smooth=0));
  connect(mFlow_6m.port_b, Boundary_fixed.port) annotation (points=[40,0; 60,0;
        60,-80; 80,-80], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(tank.ports[3], mFlow_4m.port_a) annotation (points=[-20,39; 0,39; 0,
        -40; 20,-40], style(
      color=3,
      rgbcolor={0,0,255},
      smooth=0));
  connect(mFlow_4m.port_b, Boundary_fixed.port) annotation (points=[40,-40; 60,
        -40; 60,-80; 80,-80], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(tank.ports[4], mFlow_2m.port_a) annotation (points=[-20,39; 0,39; 0,
        -80; 20,-80], style(
      color=3,
      rgbcolor={0,0,255},
      smooth=0));
  connect(mFlow_2m.port_b, Boundary_fixed.port) annotation (points=[40,-80; 80,
        -80], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
end TestTank;
