model TestSimpleLongPipePressure "Test ShortPipe componet" 
  
  import Modelica.SIunits.Conversions.*;
  
  extends Modelica.Icons.Example;
  
  parameter Real pressurePulseHeight=1E5;
  parameter Real pressurePulseWidth=50E-6;
  parameter Real pressurePulseStart=0.1E-3;
  parameter Real pressurePulseBase=2E5;
  
  Modelica_Fluid.Sources.FixedAmbient ambient(
    port(h(start=10000)),
      p_ambient=pressurePulseBase,
      T_ambient=320,
      redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
  annotation (extent=[80,50; 60,70]);
  
annotation (
  Diagram,
  experiment(
        StopTime=0.0003,
        NumberOfIntervals=5000,
        fixedstepsize=1e-006,
        Algorithm=""),
  experimentSetupOutput,
  Commands(file="Simulate and plot pressure.mos", file=
     "Simulate and Plot Temperature.mos"));
  
  FiniteVolume.Components.SimpleLongPipe LongPipe1(
    dp_nominal=50, 
    L=0.025, 
    dynamicMomentumBalance=true, 
    shortPipe(each p_start=2e5, each T_start=320), 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa, 
    includeKineticTerm=true, 
    n=100) annotation (extent=[20,50; 40,70]);
  FiniteVolume.Sources.VaryingAmbientPressure VaryingAmbientPressure1(
    h_ambient=1.e4, 
    T_ambient=320, 
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
                  annotation (extent=[-22,50; -2,70]);
  Modelica.Blocks.Math.Add Add1 
                              annotation (extent=[-60,50; -40,70]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot1(
    y=LongPipe1.shortPipe.medium.p,
    x=linspace(0, 1, LongPipe1.n),
      minY=pressurePulseBase - pressurePulseHeight,
      maxY=pressurePulseBase + pressurePulseHeight) 
                           annotation(extent=[-100,-100; 100,20]);
  Modelica.Blocks.Sources.Ramp Ramp3(
      duration=scalar({pressurePulseWidth/2}),
      height=scalar({pressurePulseHeight}),
      offset=scalar({pressurePulseBase}),
      startTime=scalar({pressurePulseStart})) 
                 annotation (extent=[-100,30; -80,50]);
  Modelica.Blocks.Sources.Ramp Ramp4(
      duration=scalar({pressurePulseWidth/2}),
      height=scalar({-pressurePulseHeight}),
      offset=0,
      startTime=scalar({pressurePulseStart + pressurePulseWidth/2})) 
                  annotation (extent=[-100,90; -80,70]);
equation 
  
  connect(LongPipe1.port_b, ambient.port) 
  annotation (points=[41,60; 59,60], style(color=69));
    connect(Add1.y,       VaryingAmbientPressure1.p_in) 
                                                      annotation(points=[-39,60;
        -24,60], style(color=3, rgbcolor={0,0,255}));
    connect(VaryingAmbientPressure1.port, LongPipe1.port_a) 
                                                          annotation(points=
       [-1,60; 19,60], style(color=69, rgbcolor={0,127,255}));
  connect(Ramp3.y,Add1.u2)           annotation (points=[-79,40; -72,40; -72,
          54; -62,54],    style(color=3));
  connect(Ramp4.y,Add1.u1) 
  annotation (points=[-79,80; -70,80; -70,66; -62,66],   style(color=3));
  
end TestSimpleLongPipePressure;
