model TestThreeLongPipesPressure "Test ShortPipe componet" 
  
  import Modelica.SIunits.Conversions.*;
  
  extends Modelica.Icons.Example;
  parameter Real pressurePulsHeight=1E4;
  parameter Real pressurePulsWidth=1E-3;
  parameter Real pressurePulsStart=0.1E-3;
  parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient=300;
  parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient1=310;
  
annotation (
  Diagram,
  experiment(
        StopTime=0.001,
        NumberOfIntervals=5000,
        Tolerance=1e-006,
        fixedstepsize=1e-006,
        Algorithm="Dassl"),
  experimentSetupOutput,
  Commands(file="Simulate and plot pressure.mos", file=
     "Simulate and Plot Temperature.mos"));
  
  FiniteVolume.Components.SimpleLongPipe LongPipe1(
    m_dot_nominal=1, 
    A_a=0.01, 
    dp_nominal=50, 
    dynamicMomentumBalance=true, 
    n=100, 
    L=0.5, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
           annotation (extent=[12,70; 32,90]);
  FiniteVolume.Sources.VaryingAmbientPressure VaryingAmbientPressure1(
    h_ambient=1.e4, 
    T_ambient=T_ambient1, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
                  annotation (extent=[-22,70; -2,90]);
  Modelica.Blocks.Sources.Ramp Ramp1(
      duration=scalar({pressurePulsWidth/2}),
      height=scalar({pressurePulsHeight}),
      offset=1E5,
      startTime=scalar({pressurePulsStart})) 
                 annotation (extent=[-100,40; -80,60]);
  Modelica.Blocks.Sources.Ramp Ramp2(
      duration=scalar({pressurePulsWidth/2}),
      height=scalar({-pressurePulsHeight}),
      offset=0,
      startTime=scalar({pressurePulsStart + pressurePulsWidth/2})) 
                  annotation (extent=[-100,100; -80,80]);
  Modelica.Blocks.Math.Add Add1 
                              annotation (extent=[-60,70; -40,90]);
  Modelica_Fluid.Sources.FixedAmbient ambient1(
    p_ambient=1E5,
    port(h(start=10000)),
    T_ambient=T_ambient,
      redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
  annotation (extent=[100,70; 80,90]);
  FiniteVolume.Components.SimpleLongPipe LongPipe2(
    m_dot_nominal=1, 
    dp_nominal=50, 
    A_a=0.01/2, 
    dynamicMomentumBalance=true, 
    n=100, 
    L=0.5, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
           annotation (extent=[48,70; 68,90]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot1(
    maxY=1.1e5,
    minY=0.9e5,
      x=linspace(0, 1, LongPipe1.n),
      y=LongPipe1.shortPipe.medium.p) 
                        annotation(extent=[-100,-40; 0,20]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot2(
    maxY=1.1e5,
    minY=0.9e5,
      x=linspace(0, 1, LongPipe2.n),
      y=LongPipe2.shortPipe.medium.p) 
                        annotation(extent=[0,-40; 100,20]);
  Modelica_Fluid.Sources.FixedAmbient ambient2(
    p_ambient=1E5,
    port(h(start=10000)),
    T_ambient=T_ambient,
      redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
  annotation (extent=[100,30; 80,50]);
  FiniteVolume.Components.SimpleLongPipe LongPipe3(
    m_dot_nominal=1, 
    dp_nominal=25, 
    A_a=0.01/2, 
    dynamicMomentumBalance=true, 
    n=50, 
    L=0.25, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
           annotation (extent=[48,30; 68,50]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot3(
    maxY=1.1e5,
    minY=0.9e5,
      x=linspace(0, 1, LongPipe3.n),
      y=LongPipe3.shortPipe.medium.p) 
                        annotation(extent=[6,-100; 54,-40]);
equation 
  
  connect(VaryingAmbientPressure1.port, LongPipe1.port_a) 
  annotation (points=[-1,80; 11,80], style(color=69));
  connect(Ramp1.y,Add1.u2)           annotation (points=[-79,50; -72,50;
        -72,74; -62,74],  style(color=3));
  connect(Ramp2.y,Add1.u1) 
  annotation (points=[-79,90; -70,90; -70,86; -62,86],   style(color=3));
  connect(Add1.y,       VaryingAmbientPressure1.p_in) 
  annotation (points=[-39,80; -24,80], style(color=3));
  connect(LongPipe2.port_a, LongPipe1.port_b) 
                                            annotation(points=[47,80; 33,80],
      style(color=69, rgbcolor={0,127,255}));
  connect(LongPipe2.port_b, ambient1.port) 
                                         annotation(points=[69,80; 79,80],
      style(color=69, rgbcolor={0,127,255}));
  connect(LongPipe3.port_a, LongPipe1.port_b) 
                                            annotation(points=[47,40; 40,40;
        40,80; 33,80], style(color=69, rgbcolor={0,127,255}));
  connect(LongPipe3.port_b,ambient2. port) 
                                         annotation(points=[69,40; 79,40],
      style(color=69, rgbcolor={0,127,255}));
  
end TestThreeLongPipesPressure;
