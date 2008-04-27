within Modelica_Fluid.Test.TestComponents.PressureLosses;
model TestWallFriction 
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  // Modelica.Media.Water.WaterIF97_ph
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
<p>
5 different wall friction models are compared.
</p>
 
<p>
pipe1 and pipe2 should be identical to pipe3 and pipe4 for
pipe1.WallFriction = WallFriction.QuadraticTurbulent (since the same equations).
</p>
 
<p>
pipe1 and pipe2 should be identical to pressureDropPipe for
pipe1.WallFriction = WallFriction.Detailed (since the same equations).
</p>
</html>"));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX ambient_a(
                                                     redeclare package Medium 
      = Medium,
    p=ambient.default_p_ambient,
    T=ambient.default_T_ambient,
    usePressureInput=true) 
    annotation (extent=[-38,40; -18,60]);
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.99999e5; 10,1.00001e5]) 
    annotation (extent=[-80,40; -60,60]);
  
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient_p1(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[62,40; 42,60]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient_p2(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,10; 40,30]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe1(
    length=1,
    diameter=0.1,
    redeclare package Medium = Medium,
    roughness=roughness,
    port_a(m_flow(start=-0.6)),
    dp_small=0.1,
    show_Re=true,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed) 
                      annotation (extent=[0,40; 20,60]);
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe2(
    length=1,
    diameter=0.1,
    from_dp=false,
    redeclare package Medium = Medium,
    roughness=roughness,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    port_a(m_flow(start=-0.6)),
    show_Re=true)     annotation (extent=[0,10; 20,30]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient_p3(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,-20; 40,0]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient_p4(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,-50; 40,-30]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.TestWallFriction
    pipe3(
    length=1,
    diameter=0.1,
    redeclare package Medium = Medium,
    roughness=roughness,
    dp_small=0.1)        annotation (extent=[0,-20; 20,0]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.TestWallFriction
    pipe4(
    length=1,
    diameter=0.1,
    from_dp=false,
    redeclare package Medium = Medium,
    roughness=roughness) annotation (extent=[0,-50; 20,-30]);
  
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[68,74; 88,94]);
equation 
  connect(p_table.y, ambient_a.p_in) 
                                    annotation (points=[-59,50; -52,50; -52,56;
        -40,56], style(color=74, rgbcolor={0,0,127}));
  connect(ambient_a.port, pipe1.port_a) 
                                      annotation (points=[-18,50; 0,50],
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe1.port_b, ambient_p1.port) 
    annotation (points=[20,50; 42,50],style(color=69, rgbcolor={0,127,255}));
  connect(pipe2.port_b, ambient_p2.port) 
    annotation (points=[20,20; 40,20],style(color=69, rgbcolor={0,127,255}));
  connect(ambient_a.port, pipe2.port_a) 
                                      annotation (points=[-18,50; -12,50; -12,
        20; 0,20],   style(color=69, rgbcolor={0,127,255}));
  connect(pipe3.port_b,ambient_p3. port) 
    annotation (points=[20,-10; 40,-10],
                                      style(color=69, rgbcolor={0,127,255}));
  connect(pipe4.port_b,ambient_p4. port) 
    annotation (points=[20,-40; 40,-40],
                                      style(color=69, rgbcolor={0,127,255}));
  connect(ambient_a.port, pipe3.port_a) 
                                      annotation (points=[-18,50; -12,50; -12,
        -10; 0,-10],
                style(color=69, rgbcolor={0,127,255}));
  connect(ambient_a.port, pipe4.port_a) 
                                      annotation (points=[-18,50; -12,50; -12,
        -40; 0,-40], style(color=69, rgbcolor={0,127,255}));
end TestWallFriction;
