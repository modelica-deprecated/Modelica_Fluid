model TestSimpleIsolatedPipePressure "Test IsolatedPipe component" 
  
  import Modelica.SIunits.Conversions.*;
  
  extends Modelica.Icons.Example;
  
  parameter Real pressurePulseHeight=1E5;
  parameter Real pressurePulseWidth=50E-6;
  parameter Real pressurePulseStart=0.1E-3;
  parameter Real pressurePulseBase=2E5;
  
  Modelica_Fluid.Sources.FixedAmbient_pTX ambient(
    port(h(start=10000)),
      p_ambient=pressurePulseBase,
      T_ambient=320,
      redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
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
  
  Modelica_Fluid.Components.IsolatedPipe isolatedPipe(
    redeclare package Medium = Modelica.Media.Air.DryAirNasa,
    dp_nominal=50,
    L=0.025,
    dynamicMomentumBalance=true,
    includeKineticTerm=true,
    nVolumes=25,
    p_start=2e5,
    T_start=320,
    A_a=1,
    includeViscosity=true,
    initType=Modelica_Fluid.Types.InitTypes.InitialStates) 
           annotation (extent=[20,50; 40,70]);
  Modelica_Fluid.Sources.PrescribedAmbient_pT prescribedAmbient(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa) 
                  annotation (extent=[-22,50; -2,70]);
  Modelica.Blocks.Math.Add Add1 
                              annotation (extent=[-58,56; -38,76]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot1(
    y=isolatedPipe.pipeSegment.medium.p,
    x=linspace(0, 1, isolatedPipe.nVolumes),
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
  Modelica.Blocks.Sources.Constant Constant1(k=320) 
    annotation (extent=[-60,20; -40,40]);
equation 
  
  connect(isolatedPipe.port_b, ambient.port) 
  annotation (points=[41,60; 59,60], style(color=69));
  connect(prescribedAmbient.port, isolatedPipe.port_a)    annotation(points=
       [-1,60; 19,60], style(color=69, rgbcolor={0,127,255}));
  connect(Ramp3.y,Add1.u2)           annotation (points=[-79,40; -72,40; -72,60;
        -60,60],          style(color=3));
  connect(Ramp4.y,Add1.u1) 
  annotation (points=[-79,80; -70,80; -70,72; -60,72],   style(color=3));
  
  connect(Add1.y, prescribedAmbient.p_ambient) 
    annotation (points=[-37,66; -24,66], style(color=3, rgbcolor={0,0,255}));
  connect(Constant1.y, prescribedAmbient.T_ambient) annotation (points=[-39,30;
        -32,30; -32,54; -24,54], style(color=3, rgbcolor={0,0,255}));
end TestSimpleIsolatedPipePressure;
