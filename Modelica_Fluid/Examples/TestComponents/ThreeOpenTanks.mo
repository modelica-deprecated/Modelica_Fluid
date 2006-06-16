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
  
  Modelica_Fluid.Components.FluidStorage.OpenTank tank1(
    area=1,
    redeclare package Medium = Medium,
    p_static_at_port=false,
    height=12,
    level_start=8,
    zeta_in={1.05},
    n_ports=1,
    pipe_diameters={0.1}) 
                   annotation (extent=[-80,20; -40,60]);
  Modelica_Fluid.Components.FluidStorage.OpenTank tank2(
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
  
  inner Modelica_Fluid.Components.Ambient ambient 
                                   annotation (extent=[76,-96; 96,-76]);
  Modelica_Fluid.Components.FluidStorage.OpenTank tank3(
    area=1,
    redeclare package Medium = Medium,
    p_static_at_port=false,
    height=12,
    level_start=3,
    zeta_in={1.05},
    n_ports=1,
    pipe_diameters={0.1}) 
                   annotation (extent=[40,20; 80,60]);
  Modelica_Fluid.Components.PressureLosses.StaticHead pipe1(redeclare package 
      Medium =                                                                       Medium,
    flowDirection=Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional,
    height_ab=2) annotation (extent=[-70,-20; -50,0], rotation=90);
  Modelica_Fluid.Components.PressureLosses.StaticHead pipe2(redeclare package 
      Medium =                                                                       Medium,
    flowDirection=Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional,
    height_ab=2) annotation (extent=[-10,-20; 10,0], rotation=90);
  Modelica_Fluid.Components.PressureLosses.StaticHead pipe3(redeclare package 
      Medium =                                                                       Medium,
    flowDirection=Modelica_Fluid.Types.FlowDirectionWithGlobalDefault.Bidirectional,
    height_ab=-1) annotation (extent=[50,-20; 70,0], rotation=90);
equation 
  connect(tank1.port[1], pipe1.port_b) annotation (points=[-60.4,20.2; -60.4,
        10.1; -60,10.1; -60,0], style(color=69, rgbcolor={0,127,255}));
  connect(tank2.port[1], pipe2.port_b) annotation (points=[-0.4,20.2; -0.4,10.1;
        6.12303e-016,10.1; 6.12303e-016,0], style(color=69, rgbcolor={0,127,255}));
  connect(tank3.port[1], pipe3.port_b) annotation (points=[59.6,20.2; 59.6,9.1;
        60,9.1; 60,0], style(color=69, rgbcolor={0,127,255}));
  connect(pipe1.port_a, pipe2.port_a) annotation (points=[-60,-20; -62,-20; -62,
        -42; -6.12303e-016,-42; -6.12303e-016,-20], style(color=69, rgbcolor={0,
          127,255}));
  connect(pipe2.port_a, pipe3.port_a) annotation (points=[-6.12303e-016,-20; 0,
        -20; 0,-42; 60,-42; 60,-20], style(color=69, rgbcolor={0,127,255}));
end ThreeOpenTanks;
