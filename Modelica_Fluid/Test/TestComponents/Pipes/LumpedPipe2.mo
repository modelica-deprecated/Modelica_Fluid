within Modelica_Fluid.Test.TestComponents.Pipes;
model LumpedPipe2 
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  replaceable package Medium = 
      Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater;
  //Modelica.Media.Water.StandardWater;
  
  Modelica_Fluid.Sources.FixedBoundary_pTX source(
    redeclare package Medium = Medium,
    p=5.0e5,
    T=300) annotation (extent=[-76,4; -64,16]);
  
  annotation (
    Diagram,
    experiment(StopTime=5),
    experimentSetupOutput(equdistant=false),
    Documentation(info="<html>
Simulation starts with both valves open. At t=1, valve 1 closes; at t=2 valve 2 closes, and the simulation fails.
</html>"));
  Modelica_Fluid.Sources.FixedBoundary_pTX sink(
    redeclare package Medium = Medium,
    p=100000,
    T=300)   annotation (extent=[76,4; 64,16]);
  
  inner Ambient ambient annotation (extent=[-100,60; -80,80]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity wallFriction1(
    redeclare package Medium = Medium, 
    length=5, 
    diameter=2.54e-2, 
    height_ab=0, 
    roughness=0, 
    use_nominal=false, 
    eta_nominal=1, 
    d_nominal=1, 
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed, 
    p_a_start=100000, 
    p_b_start=100000, 
    use_T_start=true, 
    T_start=293.15)    annotation (extent=[-46,0; -26,20]);
  Modelica_Fluid.Volumes.MixingVolume volume(
    redeclare package Medium = Medium, 
    V=Modelica.Constants.pi*(2.54e-2/2)^2*10, 
    p_start=100000, 
    T_start=293.15) 
    annotation (extent=[-8,0; 12,20]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity wallFriction2(
    redeclare package Medium = Medium, 
    length=5, 
    diameter=2.54e-2, 
    height_ab=0, 
    roughness=0, 
    use_nominal=false, 
    eta_nominal=1, 
    d_nominal=1, 
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed, 
    p_a_start=100000, 
    p_b_start=100000, 
    T_start=293.15)    annotation (extent=[30,0; 50,20]);
equation 
  connect(volume.port_a, wallFriction1.port_b)       annotation (points=[-8.2,10; 
        -26,10],         style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(volume.port_b, wallFriction2.port_a)       annotation (points=[12,10; 
        30,10],style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(wallFriction2.port_b, sink.port) annotation (points=[50,10; 64,10], 
      style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(wallFriction1.port_a, source.port) annotation (points=[-46,10; -64,10], 
      style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
end LumpedPipe2;
