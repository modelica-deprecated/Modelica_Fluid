package ShockWaves "Examples demonstrating shock waves in pipes" 
extends Modelica.Icons.Library;
  
model TestSimpleIsolatedPipePressure "Test IsolatedPipe component" 
    
    import Modelica.SIunits.Conversions.*;
    
  extends Modelica.Icons.Example;
    
  parameter Real pressurePulseHeight=1E5;
  parameter Real pressurePulseWidth=50E-6;
  parameter Real pressurePulseStart=0.1E-3;
  parameter Real pressurePulseBase=2E5;
    
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient(
    port(h(start=10000)),
      p=pressurePulseBase,
      T=320,
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
    
  Modelica_Fluid.WorkInProgress.Components.IsolatedPipe isolatedPipe(
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
    initType=Modelica_Fluid.Types.Init.InitialValues) 
           annotation (extent=[20,50; 40,70]);
  Modelica_Fluid.Sources.PrescribedBoundary_pTX prescribedAmbient(
                                                                 redeclare 
        package Medium=Modelica.Media.Air.DryAirNasa) 
                  annotation (extent=[-22,50; -2,70]);
  Modelica.Blocks.Math.Add Add1 
                              annotation (extent=[-58,56; -38,76]);
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
    
  UserInteraction.Outputs.SpatialPlot SpatialPlot1(
    y=isolatedPipe.pipeSegment.medium.p,
    x=linspace(0, 1, isolatedPipe.nVolumes),
      minY=pressurePulseBase - pressurePulseHeight,
      maxY=pressurePulseBase + pressurePulseHeight) 
                           annotation(extent=[-100,-100; 100,20]);
    
equation 
  connect(isolatedPipe.port_b, ambient.port) 
  annotation (points=[41,60; 60,60], style(color=69));
  connect(prescribedAmbient.port, isolatedPipe.port_a)    annotation(points=[-2,60;
          19,60],      style(color=69, rgbcolor={0,127,255}));
  connect(Ramp3.y,Add1.u2)           annotation (points=[-79,40; -72,40; -72,60;
        -60,60],          style(color=3));
  connect(Ramp4.y,Add1.u1) 
  annotation (points=[-79,80; -70,80; -70,72; -60,72],   style(color=3));
    
  connect(Constant1.y, prescribedAmbient.T_in) annotation (points=[-39,30; -32,30;
        -32,60; -24,60], style(color=74, rgbcolor={0,0,127}));
    connect(Add1.y, prescribedAmbient.p_in) annotation (points=[-37,66; -24,66],
        style(color=74, rgbcolor={0,0,127}));
end TestSimpleIsolatedPipePressure;
  
model TestThreeIsolatedPipesPressure "Test ShortPipe componet" 
    
    import Modelica.SIunits.Conversions.*;
    
  extends Modelica.Icons.Example;
  parameter Real pressurePulsHeight=1E4;
  parameter Real pressurePulsWidth=1E-3;
  parameter Real pressurePulsStart=0.1E-3;
  parameter Modelica.Media.Interfaces.PartialMedium.Temperature T_ambient=300;
  parameter Modelica.Media.Interfaces.PartialMedium.Temperature T_ambient1=310;
    
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
    
  Modelica_Fluid.WorkInProgress.Components.IsolatedPipe IsolatedPipe1(
    m_flow_nominal=1,
    A_a=0.01,
    dp_nominal=50,
    dynamicMomentumBalance=true,
    L=0.5,
    redeclare package Medium = Modelica.Media.Air.DryAirNasa,
    includeKineticTerm=true,
    includeViscosity=true,
    nVolumes=25,
      h_start=1e4,
      T_start=293.15,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
           annotation (extent=[12,70; 32,90]);
  Modelica_Fluid.Sources.PrescribedBoundary_pTX prescribedAmbient(
      redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                  annotation (extent=[-22,70; -2,90]);
  Modelica.Blocks.Sources.Ramp Ramp1(
      duration=scalar({pressurePulsWidth/2}),
      height=scalar({pressurePulsHeight}),
      offset=1E5,
      startTime=scalar({pressurePulsStart})) 
                 annotation (extent=[-100,52; -80,72]);
  Modelica.Blocks.Sources.Ramp Ramp2(
      duration=scalar({pressurePulsWidth/2}),
      height=scalar({-pressurePulsHeight}),
      offset=0,
      startTime=scalar({pressurePulsStart + pressurePulsWidth/2})) 
                  annotation (extent=[-100,106; -80,86]);
  Modelica.Blocks.Math.Add Add1 
                              annotation (extent=[-60,76; -40,96]);
  Modelica_Fluid.Sources.FixedBoundary ambient1(
    p=1E5,
    port(h(start=10000)),
    T=T_ambient,
      redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
  annotation (extent=[100,70; 80,90]);
  Modelica_Fluid.WorkInProgress.Components.IsolatedPipe IsolatedPipe2(
    m_flow_nominal=1,
    dp_nominal=50,
    A_a=0.01/2,
    dynamicMomentumBalance=true,
    L=0.5,
    redeclare package Medium = Modelica.Media.Air.DryAirNasa,
    includeKineticTerm=true,
    includeViscosity=true,
    nVolumes=25,
      T_start=293.15,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
           annotation (extent=[48,70; 68,90]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot1(
    maxY=1.1e5,
    minY=0.9e5,
      x=linspace(0, 1, IsolatedPipe1.nVolumes),
    y=IsolatedPipe1.pipeSegment.medium.p) 
                        annotation(extent=[-100,-40; 0,20]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot2(
    maxY=1.1e5,
    minY=0.9e5,
      x=linspace(0, 1, IsolatedPipe2.nVolumes),
    y=IsolatedPipe2.pipeSegment.medium.p) 
                        annotation(extent=[0,-40; 100,20]);
  Modelica_Fluid.Sources.FixedBoundary ambient2(
    p=1E5,
    port(h(start=10000)),
    T=T_ambient,
      redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
  annotation (extent=[100,30; 80,50]);
  Modelica_Fluid.WorkInProgress.Components.IsolatedPipe IsolatedPipe3(
    m_flow_nominal=1,
    dp_nominal=25,
    A_a=0.01/2,
    dynamicMomentumBalance=true,
    L=0.25,
    redeclare package Medium = Modelica.Media.Air.DryAirNasa,
    includeKineticTerm=true,
    includeViscosity=true,
    nVolumes=20,
      T_start=293.15,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
           annotation (extent=[48,30; 68,50]);
  UserInteraction.Outputs.SpatialPlot SpatialPlot3(
    maxY=1.1e5,
    minY=0.9e5,
      x=linspace(0, 1, IsolatedPipe3.nVolumes),
    y=IsolatedPipe3.pipeSegment.medium.p) 
                        annotation(extent=[6,-100; 54,-40]);
  Modelica.Blocks.Sources.Constant Constant1(k=T_ambient1) 
    annotation (extent=[-60,46; -40,66]);
equation 
    
  connect(prescribedAmbient.port, IsolatedPipe1.port_a) 
  annotation (points=[-2,80; 11,80], style(color=69));
  connect(Ramp1.y,Add1.u2)           annotation (points=[-79,62; -72,62; -72,80;
        -62,80],          style(color=3));
  connect(Ramp2.y,Add1.u1) 
  annotation (points=[-79,96; -70,96; -70,92; -62,92],   style(color=3));
  connect(IsolatedPipe2.port_a, IsolatedPipe1.port_b) 
                                            annotation(points=[47,80; 33,80],
      style(color=69, rgbcolor={0,127,255}));
  connect(IsolatedPipe2.port_b, ambient1.port) 
                                         annotation(points=[69,80; 80,80],
      style(color=69, rgbcolor={0,127,255}));
  connect(IsolatedPipe3.port_a, IsolatedPipe1.port_b) 
                                            annotation(points=[47,40; 40,40;
        40,80; 33,80], style(color=69, rgbcolor={0,127,255}));
  connect(IsolatedPipe3.port_b,ambient2. port) 
                                         annotation(points=[69,40; 80,40],
      style(color=69, rgbcolor={0,127,255}));
    
  connect(Add1.y, prescribedAmbient.p_in) annotation (points=[-39,86; -32,86;
        -32,86; -24,86], style(color=74, rgbcolor={0,0,127}));
  connect(Constant1.y, prescribedAmbient.T_in) annotation (points=[-39,56; -32,
        56; -32,80; -24,80], style(color=74, rgbcolor={0,0,127}));
end TestThreeIsolatedPipesPressure;
  
end ShockWaves;