model TestLongPipeTemperature "Test ShortPipe componet" 
  
  import Modelica.SIunits.Conversions.*;
  
  extends Modelica.Icons.Example;
  FiniteVolume.Sources.FixedAmbient ambient(
    T_ambient=300, 
    p_ambient=1E5, 
    port(h(start=10000)), 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
  annotation (extent=[40,60; 20,80]);
  
annotation (
  Diagram,
  experiment(
        StopTime=10,
        NumberOfIntervals=5000,
        Tolerance=1e-006,
        fixedstepsize=1e-006,
        Algorithm=""),
  experimentSetupOutput,
  Commands(file="Simulate and plot pressure.mos", file=
     "Simulate and Plot Temperature.mos"));
  
  FiniteVolume.Components.SimpleLongPipe LongPipe1(
    dp_nominal=500, 
    m_dot_nominal=1, 
    A_a=0.01, 
    L=0.025, 
    n=20, 
    k=0.001, 
    includeThermalConductance=true, 
    dynamicMomentumBalance=false, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
           annotation (extent=[-20,60; 0,80]);
  
//  Real T[LongPipe1.n](start=cat(1,fill(310, div(LongPipe1.n,2)),fill(300, LongPipe1.n-div(LongPipe1.n,2))),
//    fixed=fill(true,LongPipe1.n));
  
  FiniteVolume.Sources.FixedAmbient ambient1(
    p_ambient=1E5, 
    port(h(start=10000)), 
    T_ambient=310, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
  annotation (extent=[-40,80; -60,60], rotation=180);
  UserInteraction.Outputs.SpatialPlot SpatialPlot1(
    minY=295,
    maxY=315,
    x=linspace(0, 1, LongPipe1.n),
    y=LongPipe1.shortPipe.medium.T) 
            annotation(extent=[-100,-100; 100,0]);
equation 
  
  connect(LongPipe1.port_b, ambient.port) 
  annotation (points=[1,70; 19,70],  style(color=69));
  
  connect(ambient1.port, LongPipe1.port_a) 
                                         annotation(points=[-39,70; -21,70],
      style(color=69, rgbcolor={0,127,255}));
  
end TestLongPipeTemperature;
