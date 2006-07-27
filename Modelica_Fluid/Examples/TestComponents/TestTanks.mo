package TestTanks "Test tank components" 
  extends Modelica.Icons.Library;
  model TestOneTank 
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;
    
    Volumes.OpenTank2 tank(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      area=1,
      levelMax=1,
      bottomPortData={Modelica_Fluid.Volumes.BaseClasses.TankBottomPortData(
          diameter=0.1, portLevel=0)},
      topPortDiameter={0.1}, 
      level_start=0.9) 
      annotation (extent=[-40,20; 0,60]);
    
    Sources.PrescribedMassFlowRate_TX flowSource(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      m_flow=1,
      flowDirection=Modelica_Fluid.Types.SourceFlowDirection.OutOfPort) 
      annotation (extent=[-60,70; -40,90]);
    annotation (Diagram, 
      experiment(StopTime=30), 
      experimentSetupOutput);
    inner Ambient ambient annotation (extent=[40,60; 60,80]);
    Sources.FixedAmbient_pTX ambient_fixed(redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      flowDirection=Modelica_Fluid.Types.SourceFlowDirection.InToPort) 
      annotation (extent=[-60,-10; -40,10]);
  equation 
    connect(flowSource.port, tank.topPort[1])   annotation (points=[-40,80; -20,
          80; -20,60], style(color=69, rgbcolor={0,127,255}));
    connect(tank.bottomPort[1], ambient_fixed.port) annotation (points=[-20,20; 
          -20,0; -40,0],     style(color=69, rgbcolor={0,127,255}));
  end TestOneTank;
  
  model ThreeOpenTanks "Demonstrating the usage of OpenTank" 
    import Modelica_Fluid;
    extends Modelica.Icons.Example;
     // replaceable package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater extends 
    // replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase extends 
    // replaceable package Medium = Modelica.Media.Incompressible.Examples.Glycol47 extends
     replaceable package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater                           extends 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
    
    Modelica_Fluid.Volumes.OpenTank tank1(
      area=1,
      redeclare package Medium = Medium,
      p_static_at_port=false,
      height=12,
      level_start=8,
      zeta_in={1.05},
      n_ports=1,
      pipe_diameters={0.1}) 
                     annotation (extent=[-80,20; -40,60]);
    Modelica_Fluid.Volumes.OpenTank tank2(
      area=1,
      redeclare package Medium = Medium,
      p_static_at_port=false,
      height=12,
      level_start=3,
      zeta_in={1.05},
      n_ports=1,
      pipe_diameters={0.1}) 
                     annotation (extent=[-20,20; 20,60]);
    annotation (Diagram,
      experiment(StopTime=100),
      experimentSetupOutput,
      Documentation(info="<html>
  
</html>"));
    
    inner Modelica_Fluid.Ambient ambient 
                                     annotation (extent=[76,-96; 96,-76]);
    Modelica_Fluid.Volumes.OpenTank tank3(
      area=1,
      redeclare package Medium = Medium,
      p_static_at_port=false,
      height=12,
      level_start=3,
      zeta_in={1.05},
      n_ports=1,
      pipe_diameters={0.1}) 
                     annotation (extent=[40,20; 80,60]);
    Modelica_Fluid.PressureLosses.StaticHead pipe1(           redeclare package
        Medium =                                                                       Medium,
      flowDirection=Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional,
      height_ab=2) annotation (extent=[-70,-20; -50,0], rotation=90);
    Modelica_Fluid.PressureLosses.StaticHead pipe2(           redeclare package
        Medium =                                                                       Medium,
      flowDirection=Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional,
      height_ab=2) annotation (extent=[-10,-20; 10,0], rotation=90);
    Modelica_Fluid.PressureLosses.StaticHead pipe3(           redeclare package
        Medium =                                                                       Medium,
      flowDirection=Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional,
      height_ab=-1) annotation (extent=[50,-20; 70,0], rotation=90);
  equation 
    connect(tank1.port[1], pipe1.port_b) annotation (points=[-60.4,20.2; -60.4,
          10.1; -60,10.1; -60,0], style(color=69, rgbcolor={0,127,255}));
    connect(tank2.port[1], pipe2.port_b) annotation (points=[-0.4,20.2; -0.4,
          10.1; 6.12303e-016,10.1; 6.12303e-016,0],
                                              style(color=69, rgbcolor={0,127,255}));
    connect(tank3.port[1], pipe3.port_b) annotation (points=[59.6,20.2; 59.6,9.1;
          60,9.1; 60,0], style(color=69, rgbcolor={0,127,255}));
    connect(pipe1.port_a, pipe2.port_a) annotation (points=[-60,-20; -62,-20; 
          -62,-42; -6.12303e-016,-42; -6.12303e-016,-20],
                                                      style(color=69, rgbcolor={0,
            127,255}));
    connect(pipe2.port_a, pipe3.port_a) annotation (points=[-6.12303e-016,-20; 
          0,-20; 0,-42; 60,-42; 60,-20],
                                       style(color=69, rgbcolor={0,127,255}));
  end ThreeOpenTanks;

  model TestEmptyOpenTank "Test whether an empty tank is properly handeled" 
    extends Modelica.Icons.Example;
    Modelica_Fluid.Volumes.OpenTank tank1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      height=1,
      area=1,
      level_start=1,
      n_bottomPorts=1,
      bottom_diameters={0.1}) annotation (extent=[-20,20; 20,60]);
    
    Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.LaminarAndQuadraticTurbulent,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (extent=[-10,0; 10,-20], rotation=-90);
    
    annotation (
      Diagram,
      experiment(StopTime=50),
      experimentSetupOutput);
    Modelica_Fluid.Volumes.OpenTank tank2(
      area=1,
      height=2,
      level_start=0,
      n_topPorts=1,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater) 
      annotation (extent=[-20,-80; 20,-40]);
    inner Modelica_Fluid.Ambient ambient 
                                     annotation (extent=[56,58; 76,78]);
  equation 
    connect(pipe.port_b, tank1.bottomPorts[1]) annotation (points=[6.12303e-016,0;
          0,0; 0,19.2], style(color=69, rgbcolor={0,127,255}));
    connect(pipe.port_a, tank2.topPorts[1]) annotation (points=[-6.12303e-016,-20;
          0,-20; 0,-39.2], style(color=69, rgbcolor={0,127,255}));
  end TestEmptyOpenTank;
  
end TestTanks;
