model ShockTube 
  parameter Modelica.SIunits.Pressure p0=5E5;
  parameter Modelica.SIunits.Temperature T0=1200;
  FiniteVolume.Components.SimpleLongPipe LongPipe1(
    dynamicMomentumBalance=true, 
    L=1, 
    A_a=0.01, 
    k=0, 
    includeThermalConductance=false, 
    shortPipe(p_start=cat(1, fill(p0, div(LongPipe1.n, 2)), fill(1E5, div(
          LongPipe1.n, 2))), T_start=cat(1, fill(T0, div(LongPipe1.n, 2)), fill(
          300, div(LongPipe1.n, 2)))), 
    m_dot_nominal=1E-3, 
    includeKineticTerm=true, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa, 
    n=100, 
    dp_nominal=0.001, 
    viscosityFactor2=0.9, 
    viscosityFactor1=0.5) 
           annotation (extent=[-20,60; 0,80]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot1(
    x=linspace(0, 1, LongPipe1.n),
      y=LongPipe1.shortPipe.medium.T,
      minY=300,
      maxY=1200)           annotation(extent=[-100,-100; 100,-20]);
    annotation(experiment(StopTime=0.0005, Tolerance=1e-006),
        experimentSetupOutput,
      Diagram);
  UserInteraction.Outputs.SpatialPlot SpatialPlot2(
    y=LongPipe1.shortPipe.medium.p,
    x=linspace(0, 1, LongPipe1.n),
      minY=1E5,
      maxY=5E5)            annotation(extent=[-100,-20; 100,60]);
end ShockTube;
