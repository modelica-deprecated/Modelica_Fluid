model TestEmptyOpenTank "Test whether an empty tank is properly handeled" 
  extends Modelica.Icons.Example;
  Components.OpenTank tank1(
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    height=1,
    area=1,
    level_start=1,
    n_bottomPorts=1,
    bottom_diameters={0.1}) annotation (extent=[-20,20; 20,60]);
  
  inner Components.FluidOptions fluidOptions 
    annotation (extent=[-100,0; -80,20]);
  PressureLosses.WallFrictionAndGravity pipe(
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.Utilities.WallFriction.LaminarAndQuadraticTurbulent,
    length=1,
    diameter=0.1,
    height_ab=1) annotation (extent=[-10,0; 10,-20], rotation=-90);
  
  annotation (
    Diagram,
    experiment(StopTime=50),
    experimentSetupOutput);
  Components.OpenTank tank2(
    area=1,
    height=2,
    level_start=0,
    n_topPorts=1,
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater) 
    annotation (extent=[-20,-80; 20,-40]);
equation 
  connect(pipe.port_b, tank1.bottomPorts[1]) annotation (points=[
        6.12303e-016,0; 0,0; 0,19.2],
                      style(color=69, rgbcolor={0,127,255}));
  connect(pipe.port_a, tank2.topPorts[1]) annotation (points=[
        -6.12303e-016,-20; 0,-20; 0,-39.2],
                         style(color=69, rgbcolor={0,127,255}));
end TestEmptyOpenTank;
