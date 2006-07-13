model ThreeTanksIF97WithShortPipe2 
  import Modelica.SIunits.Conversions;
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  
  replaceable package Medium = 
      Modelica.Media.Water.ConstantPropertyLiquidWater 
      extends Modelica.Media.Interfaces.PartialMedium "Medium in the component"
      annotation (choicesAllMatching = true);
  // replaceable package Medium = Modelica.Media.Water.WaterIF97_ph 
  replaceable package WallFriction = 
      PressureLosses.Utilities.WallFriction.QuadraticTurbulent;
  
  parameter Real D=0.1 "pipe diameters";
  final parameter Real A = Modelica.Constants.pi*(D/2)^2 "pipe area (for tank)";
  
  annotation (
    Diagram,
    Coordsys(grid=[1, 1], component=[20, 20]),
    experiment(StopTime=50),
    Documentation(info="<html>
<p>
Demonstrate usage of ShortPipe model with three tanks and three short pipes.
Results when using ConstantPropertyLiquidWater:
</p>
 
<img src=\"../Images/Examples/ThreeTanksResult3.png\">
 
<p>
Results when using Modelica.Media.Water.WaterIF97_ph:
</p>
 
<img src=\"../Images/Examples/ThreeTanksResult3.png\">
 
</html>"),
    experimentSetupOutput);
  Modelica_Fluid.Components.Tank Tank1(
    area=1,
    T_start=Conversions.from_degC(50),
    level_start=3,
    redeclare package Medium = Medium,
    pipeArea=A) 
    annotation (extent=[-90, 20; -70, 40]);
  Modelica_Fluid.Components.Tank Tank2(
    area=1,
    level_start=1,
    initOption=Modelica_Fluid.Types.Init.InitialValues,
    T_start=Conversions.from_degC(90),
    redeclare package Medium = Medium,
    pipeArea=A) 
    annotation (extent=[-10,20; 10,40]);
  Modelica_Fluid.Components.Tank Tank3(
    area=1,
    T_start=Conversions.from_degC(20),
    level_start=2,
    redeclare package Medium = Medium,
    pipeArea=A) 
    annotation (extent=[70, 20; 90, 40]);
  Modelica_Fluid.WorkInProgress.Components.ShortPipe2 shortPipe1(
    redeclare package Medium = Medium,
    length=1,
    height_ab=-1,
    T_start=Conversions.from_degC(20),
    diameter=D,
    redeclare package WallFriction = WallFriction,
    initVolume2=Modelica_Fluid.Types.Init.NoInit) 
    annotation (extent=[-50, -30; -30, -10]);
  Modelica_Fluid.WorkInProgress.Components.ShortPipe2 shortPipe2(
    redeclare package Medium = Medium,
    T_start=Conversions.from_degC(20),
    length=1,
    diameter=D,
    height_ab=-1,
    redeclare package WallFriction = WallFriction,
    initVolume2=Modelica_Fluid.Types.Init.NoInit) 
    annotation (extent=[-10,-10; 10,10],   rotation=-90);
  Modelica_Fluid.WorkInProgress.Components.ShortPipe2 shortPipe3(
    redeclare package Medium = Medium,
    T_start=Conversions.from_degC(20),
    length=1,
    height_ab=1,
    diameter=D,
    redeclare package WallFriction = WallFriction) 
    annotation (extent=[30, -30; 50, -10]);
  inner Modelica_Fluid.Components.FluidOptions fluidOptions(default_initOption=
        Modelica_Fluid.Types.Init.InitialValues) 
    annotation (extent=[-100,-100; -80,-80]);
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
end ThreeTanksIF97WithShortPipe2;
