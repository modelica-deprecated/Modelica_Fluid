model TestThreeLongPipesTemperature "Test ShortPipe componet" 
  
  import Modelica.SIunits.Conversions.*;
  
  extends Modelica.Icons.Example;
  
  parameter Integer n=20;
  parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient=300;
  parameter Modelica_Media.Interfaces.PartialMedium.Temperature T_ambient1=310;
  
  FiniteVolume.Sources.FixedAmbient ambient(
    T_ambient=T_ambient1, 
    p_ambient=1E5, 
    port(h(start=10000)), 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
  annotation (extent=[-74,70; -54,90]);
  
annotation (
  Diagram,
  experiment(
        StopTime=100,
        NumberOfIntervals=5000,
        Tolerance=1e-006,
        fixedstepsize=1e-006,
        Algorithm=""),
  experimentSetupOutput,
  Commands(file="Simulate and plot pressure.mos", file=
     "Simulate and Plot Temperature.mos"));
  
  FiniteVolume.Components.SimpleLongPipe LongPipe1(
    m_dot_nominal=1, 
    A_a=0.01, 
    dp_nominal=50, 
    L=0.05, 
    n=20, 
    includeThermalConductance=true, 
    k=0.1, 
    dynamicMomentumBalance=false, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
           annotation (extent=[-40,70; -20,90]);
  FiniteVolume.Sources.FixedAmbient ambient1(
    p_ambient=1E5, 
    port(h(start=10000)), 
    T_ambient=T_ambient1, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
  annotation (extent=[48,70; 28,90]);
  FiniteVolume.Components.SimpleLongPipe LongPipe2(
    m_dot_nominal=1, 
    dp_nominal=50, 
    A_a=0.01/2, 
    L=0.05, 
    n=20, 
    includeThermalConductance=true, 
    k=0.1, 
    dynamicMomentumBalance=false, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
           annotation (extent=[-4,70; 16,90]);
  FiniteVolume.Sources.FixedAmbient ambient2(
    p_ambient=1E5, 
    port(h(start=10000)), 
    T_ambient=T_ambient, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
  annotation (extent=[48,30; 28,50]);
  FiniteVolume.Components.SimpleLongPipe LongPipe3(
    m_dot_nominal=1, 
    dp_nominal=25, 
    A_a=0.01/2, 
    L=0.025, 
    n=20, 
    includeThermalConductance=true, 
    k=0.1, 
    dynamicMomentumBalance=false, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
           annotation (extent=[-4,30; 16,50]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot1(
    minY=295,
    maxY=315,
      maxX=LongPipe1.L,
      y=LongPipe1.shortPipe.medium.T,
      x=linspace(0, LongPipe1.L, LongPipe1.n)) 
                                   annotation(extent=[-100,-40; 0,20]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot2(
    minY=295,
    maxY=315,
      maxX=LongPipe2.L,
      x=linspace(0, LongPipe2.L, LongPipe2.n),
      y=LongPipe2.shortPipe.medium.T) 
                                   annotation(extent=[0,-40; 100,20]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot3(
    minY=295,
    maxY=315,
      maxX=LongPipe3.L,
      x=linspace(0, LongPipe3.L, LongPipe3.n),
      y=LongPipe3.shortPipe.medium.T) 
                                   annotation(extent=[6,-100; 54,-40]);
equation 
  
  connect(LongPipe2.port_a,LongPipe1. port_b) 
                                            annotation(points=[-5,80; -19,80],
             style(color=69, rgbcolor={0,127,255}));
  
  connect(LongPipe2.port_b,ambient1. port) 
                                         annotation(points=[17,80; 27,80],
      style(color=69, rgbcolor={0,127,255}));
  
  connect(LongPipe3.port_a,LongPipe1. port_b) 
                                            annotation(points=[-5,40; -12,40;
          -12,80; -19,80],
                        style(color=69, rgbcolor={0,127,255}));
  
  connect(LongPipe3.port_b,ambient2. port) 
                                         annotation(points=[17,40; 27,40],
      style(color=69, rgbcolor={0,127,255}));
  
  connect(ambient.port, LongPipe1.port_a) 
                                        annotation(points=[-53,80; -41,80],
      style(
      color=69,
      rgbcolor={0,127,255},
      gradient=2,
      fillColor=76,
      rgbfillColor={170,170,255}));
  
end TestThreeLongPipesTemperature;
