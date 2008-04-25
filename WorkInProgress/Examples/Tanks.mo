within Modelica_Fluid.WorkInProgress.Examples;
package Tanks "Examples with Tanks" 
  extends Modelica.Icons.Library;
  
  
  
  
  
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
    Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe1(
      length=1,
      p_a_start=ambient.default_p_ambient,
      p_b_start=ambient.default_p_ambient,
      T_start=ambient.default_T_ambient,
      diameter=0.1,
      height_ab=0,
      redeclare package Medium = Medium,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed) 
            annotation (extent=[-20,-10; 0,10],  rotation=0);
    annotation (Diagram,
      experiment(StopTime=40),
      experimentSetupOutput,
      Documentation(info="<html>
  
</html>"));
    
    inner Modelica_Fluid.Ambient ambient 
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
    
    Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe1(
      redeclare package Medium = Medium,
      length=1,
      height_ab=2,
      diameter=0.05,
      p_a_start=ambient.default_p_ambient,
      p_b_start=ambient.default_p_ambient,
      T_start=ambient.default_T_ambient,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction) 
            annotation (extent=[-30,-40; -10,-20],
                                                 rotation=90);
    Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe2(
      redeclare package Medium = Medium,
      length=1,
      diameter=0.1,
      height_ab=2,
      p_a_start=ambient.default_p_ambient,
      p_b_start=ambient.default_p_ambient,
      T_start=ambient.default_T_ambient,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction) 
            annotation (extent=[40,30; 60,50],   rotation=90);
    inner Modelica_Fluid.Ambient ambient 
      annotation (extent=[-90,-90; -70,-70]);
    Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe3(
      redeclare package Medium = Medium,
      length=1,
      height_ab=2,
      diameter=0.05,
      p_a_start=ambient.default_p_ambient,
      p_b_start=ambient.default_p_ambient,
      T_start=ambient.default_T_ambient,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction) 
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
