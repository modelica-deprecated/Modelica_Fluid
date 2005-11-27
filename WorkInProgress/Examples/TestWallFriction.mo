model TestWallFriction 
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  replaceable package Medium = 
      Modelica.Media.Water.ConstantPropertyLiquidWater 
    extends Modelica.Media.Interfaces.PartialMedium "Medium in all components" 
                                                      annotation (
    choicesAllMatching =                                                                            true);
  parameter Modelica.SIunits.Length roughness = 0.025e-3;
  
  annotation (Diagram,
    experiment(StopTime=10, NumberOfIntervals=10000),
    experimentSetupOutput,
    Coordsys(extent=[-100,-100; 100,100]),
    Documentation(info="<html>
</html>"));
  Sources.PrescribedAmbient_pTX ambient(redeclare package Medium = Medium) 
    annotation (extent=[-40,40; -20,60]);
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.9999e5; 10,1.0001e5]) 
    annotation (extent=[-80,40; -60,60]);
  
  Sources.FixedAmbient_pTX ambient_p1(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[62,40; 42,60]);
  Sources.FixedAmbient_pTX ambient_p2(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,10; 40,30]);
  Components.WallFriction pipe1(
    length=1,
    diameter=0.1,
    redeclare package Medium = Medium,
    roughness=roughness) annotation (extent=[0,40; 20,60]);
  Components.WallFriction pipe2(
    length=1,
    diameter=0.1,
    from_dp=false,
    redeclare package Medium = Medium,
    roughness=roughness) annotation (extent=[0,10; 20,30]);
  Modelica_Fluid.Components.PressureDropPipe pressureDropPipe(
    redeclare package Medium = Medium,
    frictionType=Modelica_Fluid.Types.FrictionTypes.DetailedFriction,
    length=1,
    diameter=0.1,
    roughness=roughness) annotation (extent=[0,-20; 20,0]);
  Sources.FixedAmbient_pTX ambient_p3(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,-20; 40,0]);
equation 
  connect(p_table.y, ambient.p_in)  annotation (points=[-59,50; -52,50; -52,56;
        -42,56], style(color=74, rgbcolor={0,0,127}));
  connect(ambient.port, pipe1.port_a) annotation (points=[-19,50; -1,50],
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe1.port_b, ambient_p1.port) 
    annotation (points=[21,50; 41,50],style(color=69, rgbcolor={0,127,255}));
  connect(pipe2.port_b, ambient_p2.port) 
    annotation (points=[21,20; 39,20],style(color=69, rgbcolor={0,127,255}));
  connect(ambient.port, pipe2.port_a) annotation (points=[-19,50; -12,50; -12,
        20; -1,20],  style(color=69, rgbcolor={0,127,255}));
  connect(pressureDropPipe.port_b, ambient_p3.port) 
    annotation (points=[21,-10; 39,-10],
                                      style(color=69, rgbcolor={0,127,255}));
  connect(pressureDropPipe.port_a, ambient.port)  annotation (points=[-1,-10;
        -12,-10; -12,50; -19,50],style(color=69, rgbcolor={0,127,255}));
end TestWallFriction;
