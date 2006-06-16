package Tanks "Examples with Tanks" 
  extends Modelica.Icons.Library;
  
  model TwoTanksSimpleWater 
    import Modelica.SIunits.Conversions.*;
    import Modelica_Fluid;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=50),
      Documentation(info="<html>
<p>
Two tanks with different initial levels are connected by a 
horizontal pipe. After about 30 s, the levels of the
two tanks are identical and the system is in steady state.
The two tanks have different initial temperatures. Since the
water is flowing from tank1 to tank2, the temperature of tank2
changes until water from tank1 is flowing into tank2.
</p>
</html>"),
      experimentSetupOutput);
    Modelica_Fluid.WorkInProgress.Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V0=0.1,
      pipeArea=0.01) 
      annotation (extent=[-70,20; -50,40]);
    Modelica_Fluid.WorkInProgress.Components.Tank Tank2(
      area=1,
      T_start=from_degC(100),
      level_start=1,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V0=0.1,
      pipeArea=0.01) 
      annotation (extent=[10,20; 30,40]);
    Modelica_Fluid.Components.PressureLosses.PressureDropPipe shortPipe1(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[-29,0; -9,20]);
    
    inner Modelica_Fluid.Components.Ambient ambient 
                                     annotation (extent=[46,64; 66,84]);
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-60,19; -60,10; -29,10],      style(color=69));
    connect(shortPipe1.port_b, Tank2.port) annotation (points=[-9,10; 20,10; 20,
          19], style(color=69, rgbcolor={0,127,255}));
  end TwoTanksSimpleWater;
  
  model ThreeTanksSimpleWater 
    import Modelica.SIunits.Conversions.*;
    import Modelica_Fluid;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=50),
      Documentation(info="<html>
<p>
This example demonstrates the mixing of simple (constant
property) liquid water model
between three tanks with different temperature.
The water model is incompressible.
The same tank system with a compressible water model
is provided in ThreeTanksIF97.
</p>
</html>"),
      experimentSetupOutput);
    Modelica_Fluid.WorkInProgress.Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V0=0.1,
      pipeArea=0.01) 
      annotation (extent=[-90,20; -70,40]);
    Modelica_Fluid.WorkInProgress.Components.Tank Tank2(
      area=1,
      T_start=from_degC(100),
      level_start=1,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      H0=0.5,
      V0=0.1,
      pipeArea=0.01) 
      annotation (extent=[-10, 20; 10, 40]);
    Modelica_Fluid.Components.PressureLosses.PressureDropPipe shortPipe1(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[-50, -30; -30, -10]);
    Modelica_Fluid.WorkInProgress.Components.Tank Tank3(
      area=1,
      T_start=from_degC(20),
      level_start=2,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      V0=0.1,
      pipeArea=0.01) 
      annotation (extent=[70,20; 90,40]);
    Modelica_Fluid.Components.PressureLosses.PressureDropPipe shortPipe3(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar) 
      annotation (extent=[30, -30; 50, -10]);
    Modelica_Fluid.Components.PressureLosses.PressureDropPipe shortPipe2(
      m_flow_nominal=1000,
      dp_nominal=from_bar(0.1),
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      medium_b(T(stateSelect=StateSelect.avoid)),
      medium_a(T(stateSelect=StateSelect.avoid))) 
      annotation (extent=[-10, -10; 10, 10], rotation=-90);
    
    inner Modelica_Fluid.Components.Ambient ambient 
                                     annotation (extent=[49,65; 69,85]);
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-80,19; -80,-20; -50,-20],    style(color=69));
    connect(shortPipe3.port_b, Tank3.port) 
      annotation (points=[50,-20; 80,-20; 80,19],    style(color=69));
    connect(Tank2.port, shortPipe2.port_a) 
      annotation (points=[0,19; 0,10; -6.12303e-016,10],    style(color=69));
    connect(shortPipe1.port_b, shortPipe3.port_a) 
      annotation (points=[-30,-20; 30,-20],   style(color=69));
    connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
          6.12303e-016,-10; 0,-10; 0,-20; 30,-20],     style(color=69));
  end ThreeTanksSimpleWater;
  
  model ThreeTanksIF97 
    import Modelica.SIunits.Conversions.*;
    import Modelica_Fluid;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=50),
      Documentation(info="<html>
<p>
This example is the same as the \"ThreeTanksSimpleWater\"
model. The only difference is that the very detailed,
compressible medium model WaterIF97 is used.
</p>
</html>"),
      experimentSetupOutput);
    Modelica_Fluid.WorkInProgress.Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      pipeArea=0.01) 
      annotation (extent=[-90, 20; -70, 40]);
    Modelica_Fluid.WorkInProgress.Components.Tank Tank2(
      area=1,
      level_start=1,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      T_start=from_degC(90),
      pipeArea=0.01) 
      annotation (extent=[-10,20; 10,40]);
    Modelica_Fluid.Components.PressureLosses.PressureDropPipe shortPipe1(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[-50, -30; -30, -10]);
    Modelica_Fluid.WorkInProgress.Components.Tank Tank3(
      area=1,
      T_start=from_degC(20),
      level_start=2,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      pipeArea=0.01) 
      annotation (extent=[70, 20; 90, 40]);
    Modelica_Fluid.Components.PressureLosses.PressureDropPipe shortPipe3(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[30, -30; 50, -10]);
    Modelica_Fluid.Components.PressureLosses.PressureDropPipe shortPipe2(
      m_flow_nominal=1000,
      dp_nominal=from_bar(0.1),
      frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[-10,-10; 10,10],   rotation=-90);
    
    inner Modelica_Fluid.Components.Ambient ambient 
                                     annotation (extent=[50,61; 70,81]);
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-80,19; -80,-20; -50,-20],    style(color=69));
    connect(shortPipe3.port_b, Tank3.port) 
      annotation (points=[50,-20; 80,-20; 80,19],    style(color=69));
    connect(Tank2.port, shortPipe2.port_a) 
      annotation (points=[0,19; 0,10; -6.12303e-016,10],    style(color=69));
    connect(shortPipe1.port_b, shortPipe3.port_a) 
      annotation (points=[-30,-20; 30,-20],   style(color=69));
    connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
          6.12303e-016,-10; 0,-10; 0,-20; 30,-20],     style(color=69));
  end ThreeTanksIF97;
  
  model ThreeTanksIF97SteadyState 
    import Modelica.SIunits.Conversions.*;
    import Modelica_Fluid;
    extends Modelica.Icons.Example;
    annotation (
      Diagram,
      Coordsys(grid=[1, 1], component=[20, 20]),
      experiment(StopTime=50),
      Documentation(info="<html>
<p>
This example is the same as ThreeTanks97. The only difference
is that the system starts in steady state, i.e., with constant
(mixing) temperature.
</p>
</html>"),
      experimentSetupOutput);
    Modelica_Fluid.WorkInProgress.Components.Tank Tank1(
      area=1,
      T_start=from_degC(50),
      level_start=3,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      pipeArea=0.01) 
      annotation (extent=[-90, 20; -70, 40]);
    Modelica_Fluid.WorkInProgress.Components.Tank Tank2(
      area=1,
      level_start=1,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      T_start=from_degC(90),
      pipeArea=0.01) 
      annotation (extent=[-10,20; 10,40]);
    Modelica_Fluid.Components.PressureLosses.PressureDropPipe shortPipe1(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[-50, -30; -30, -10]);
    Modelica_Fluid.WorkInProgress.Components.Tank Tank3(
      area=1,
      T_start=from_degC(20),
      level_start=2,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph,
      pipeArea=0.01) 
      annotation (extent=[70, 20; 90, 40]);
    Modelica_Fluid.Components.PressureLosses.PressureDropPipe shortPipe3(
      m_flow_nominal=2000,
      dp_nominal=from_bar(0.1),
      frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[30, -30; 50, -10]);
    Modelica_Fluid.Components.PressureLosses.PressureDropPipe shortPipe2(
      m_flow_nominal=1000,
      dp_nominal=from_bar(0.1),
      frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
      redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) 
      annotation (extent=[-10, -10; 10, 10], rotation=-90);
    inner Modelica_Fluid.Components.Ambient ambient 
                                     annotation (extent=[60,65; 80,85]);
  equation 
    connect(Tank1.port, shortPipe1.port_a) 
      annotation (points=[-80,19; -80,-20; -50,-20],    style(color=69));
    connect(shortPipe3.port_b, Tank3.port) 
      annotation (points=[50,-20; 80,-20; 80,19],    style(color=69));
    connect(Tank2.port, shortPipe2.port_a) 
      annotation (points=[0,19; 0,10; -6.12303e-016,10],    style(color=69));
    connect(shortPipe1.port_b, shortPipe3.port_a) 
      annotation (points=[-30,-20; 30,-20],   style(color=69));
    connect(shortPipe2.port_b, shortPipe3.port_a) annotation (points=[
          6.12303e-016,-10; 0,-10; 0,-20; 30,-20],     style(color=69));
  end ThreeTanksIF97SteadyState;
  
  model TwoOpenTanks "Demonstrating the usage of OpenTank" 
    import Modelica_Fluid;
    extends Modelica.Icons.Example;
     // replaceable package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater extends 
    // replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase extends 
    // replaceable package Medium = Modelica.Media.Incompressible.Examples.Glycol47 extends
     replaceable package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater                           extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
    
    Modelica_Fluid.WorkInProgress.Components.OpenTank tank1(
      height=10,
      area=1,
      n_bottomPorts=1,
      bottom_diameters={0.1},
      level_start=6,
      redeclare package Medium = Medium) 
                     annotation (extent=[-80,20; -40,60]);
    Modelica_Fluid.WorkInProgress.Components.OpenTank tank2(
      height=10,
      area=1,
      n_bottomPorts=1,
      bottom_diameters={0.1},
      level_start=3,
      redeclare package Medium = Medium) 
                     annotation (extent=[20,20; 60,60]);
    Modelica_Fluid.Components.PressureLosses.WallFrictionAndGravity pipe1(
      length=1,
      p_start=ambient.default_p_ambient,
      T_start=ambient.default_T_ambient,
      diameter=0.1,
      height_ab=0,
      redeclare package Medium = Medium,
      redeclare package WallFriction = 
          Modelica_Fluid.SubClasses.PressureLosses.WallFriction.Detailed) 
            annotation (extent=[-20,-10; 0,10],  rotation=0);
    annotation (Diagram,
      experiment(StopTime=40),
      experimentSetupOutput,
      Documentation(info="<html>
  
</html>"));
    
    inner Modelica_Fluid.Components.Ambient ambient 
                                     annotation (extent=[60,-34; 80,-14]);
  equation 
    connect(pipe1.port_a, tank1.bottomPorts[1]) 
                                               annotation (points=[-20,0; -60,0;
          -60,19.2], style(color=69, rgbcolor={0,127,255}));
    connect(pipe1.port_b,tank2. bottomPorts[1]) annotation (points=[0,0; 40,0; 40,
          19.2], style(color=69, rgbcolor={0,127,255}));
  end TwoOpenTanks;
  
  model ThreeOpenTanks "Demonstrating the usage of OpenTank" 
    
    annotation (
      Diagram,
      experiment(StopTime=150),
      Coordsys(grid=[1, 1], component=[20, 20]),
      uses(Modelica_Fluid(version="0.952")),
      experimentSetupOutput,
      Documentation(info="<html> 
  
</html>"));
    
    extends Modelica.Icons.Example;
    // replaceable package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater extends 
    // replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase extends 
    // replaceable package Medium = Modelica.Media.Incompressible.Examples.Glycol47 extends
     replaceable package Medium = 
       Modelica.Media.Water.ConstantPropertyLiquidWater                    extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
    
    Modelica_Fluid.WorkInProgress.Components.OpenTank tank1(
      area=1,
      V0=0,
      bottom_heights={0},
      redeclare package Medium = Medium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      height=10,
      bottom_diameters={0.05},
      n_sidePorts=2,
      side_heights={3,2},
      n_topPorts=0,
      n_bottomPorts=1,
      side_diameters={0.1,0.05},
      level_start=1) 
                    annotation (extent=[-40,0; 0,40]);
    
    Modelica_Fluid.WorkInProgress.Components.OpenTank tank2(
      area=1,
      V0=0,
      bottom_heights={0},
      bottom_diameters={0.1},
      redeclare package Medium = Medium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      height=10,
      level_start=9,
      n_bottomPorts=1) 
                     annotation (extent=[30,60; 70,100]);
    
    Modelica_Fluid.WorkInProgress.Components.OpenTank tank3(
      area=1,
      V0=0,
      level_start=6,
      top_heights={10},
      redeclare package Medium = Medium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      n_sidePorts=1,
      side_diameters={0.05},
      n_topPorts=1,
      side_heights={6.5},
      height=20)                                            annotation (extent=[-40,-90;
          0,-50]);
    
    Modelica_Fluid.Components.PressureLosses.WallFrictionAndGravity pipe1(
      redeclare package Medium = Medium,
      length=1,
      height_ab=2,
      diameter=0.05,
      p_start=ambient.default_p_ambient,
      T_start=ambient.default_T_ambient,
      redeclare package WallFriction = 
          Modelica_Fluid.SubClasses.PressureLosses.WallFriction.NoFriction) 
            annotation (extent=[-30,-40; -10,-20],
                                                 rotation=90);
    Modelica_Fluid.Components.PressureLosses.WallFrictionAndGravity pipe2(
      redeclare package Medium = Medium,
      length=1,
      diameter=0.1,
      height_ab=2,
      p_start=ambient.default_p_ambient,
      T_start=ambient.default_T_ambient,
      redeclare package WallFriction = 
          Modelica_Fluid.SubClasses.PressureLosses.WallFriction.NoFriction) 
            annotation (extent=[40,30; 60,50],   rotation=90);
    inner Modelica_Fluid.Components.Ambient ambient 
      annotation (extent=[-90,-90; -70,-70]);
    Modelica_Fluid.Components.PressureLosses.WallFrictionAndGravity pipe3(
      redeclare package Medium = Medium,
      length=1,
      height_ab=2,
      diameter=0.05,
      p_start=ambient.default_p_ambient,
      T_start=ambient.default_T_ambient,
      redeclare package WallFriction = 
          Modelica_Fluid.SubClasses.PressureLosses.WallFriction.NoFriction) 
            annotation (extent=[20,-40; 40,-20], rotation=90);
  equation 
    connect(tank2.bottomPorts[1], pipe2.port_b) 
                                      annotation (points=[50,59.2; 50,50],
        style(color=3, rgbcolor={0,0,255}));
    connect(tank1.bottomPorts[1], pipe1.port_b) 
                                      annotation (points=[-20,-0.8; -20,-20],
        style(color=3, rgbcolor={0,0,255}));
    connect(pipe1.port_a, tank3.topPorts[1]) 
                       annotation (points=[-20,-40; -20,-49.2],
        style(color=69, rgbcolor={0,127,255}));
    connect(pipe2.port_a, tank1.sidePorts[1]) annotation (points=[50,30; 50,23;
          0.8,23], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_b, tank1.sidePorts[2]) annotation (points=[30,-20; 30,16;
          0.8,16; 0.8,17], style(color=69, rgbcolor={0,127,255}));
    connect(pipe3.port_a, tank3.sidePorts[1]) annotation (points=[30,-40; 30,
          -70; 0.8,-70],
                    style(color=69, rgbcolor={0,127,255}));
  end ThreeOpenTanks;
end Tanks;
