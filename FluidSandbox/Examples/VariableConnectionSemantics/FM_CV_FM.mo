within FluidSandbox.Examples.VariableConnectionSemantics;
model FM_CV_FM 
  "Boundary condition (pressure) - flow model - control volume - flow model - boundary condition (pressure)" 
  extends Icons.Example;
  FluidSandbox.Sources.PrescribedBoundary_pTX_A source2(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[70,-30; 90,-10],
                                              rotation=90);
  
  FluidSandbox.PressureLosses.WallFriction pipeFriction1(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package FluidInterface = FluidInterface,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent)
                                              annotation (extent=[-60,0; -40,20],
             rotation=0);
  
  FluidSandbox.Volumes.Volume mixingVolume(
    V=0.1,
    redeclare package Medium = Medium,
    n_a=1,
    n_b=1,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    p_start=1e5,
    use_T_start=true,
    T_start=340,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[-10,0; 10,20]);
  FluidSandbox.PressureLosses.WallFriction pipeFriction2(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package FluidInterface = FluidInterface,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent)
                                              annotation (extent=[40,20; 60,0],
      rotation=180);
  
  FluidSandbox.Sources.PrescribedBoundary_pTX_A prescribedMassFlowRate_TX_A(
    redeclare package Medium = Medium,
    T=320,
    usePressureInput=true,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[-90,30; -70,50], rotation=270);
  Modelica.Blocks.Sources.Sine sine1(
    phase=0,
    freqHz=2,
    amplitude=0.05e5,
    offset=1e5) 
              annotation (extent=[-84,70; -64,90], rotation=270);
  annotation (Diagram);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics(redeclare 
      package Medium = Medium, redeclare package FluidInterface = 
        FluidInterface) if 
    FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[-65,7; -75,9]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics1(redeclare 
      package Medium = Medium, redeclare package FluidInterface = 
        FluidInterface) if 
    FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[-28,7; -18,9]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics2(redeclare 
      package Medium = Medium, redeclare package FluidInterface = 
        FluidInterface) if 
    FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[28,7; 18,9]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics3(redeclare 
      package Medium = Medium, redeclare package FluidInterface = 
        FluidInterface) if 
    FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[66,7; 76,9]);
equation 
  connect(sine1.y, prescribedMassFlowRate_TX_A.p_in) annotation (points=[-74,69;
        -74,52],             style(color=74, rgbcolor={0,0,127}));
  connect(prescribedMassFlowRate_TX_A.port, semantics.port_b) annotation (
      points=[-80,30; -80,8; -75,8],   style(color=69, rgbcolor={0,127,255}));
  connect(semantics.port_a, pipeFriction1.port_a) annotation (points=[-65,8;
        -62,8; -62,10; -60,10],
                 style(color=69, rgbcolor={0,127,255}));
  connect(pipeFriction1.port_b, semantics1.port_a) annotation (points=[-40,10;
        -34,10; -34,8; -28,8],
                 style(color=69, rgbcolor={0,127,255}));
  connect(semantics1.port_b, mixingVolume.port_a[1]) annotation (points=[-18,8;
        -14,8; -14,10; -10,10],
                 style(color=69, rgbcolor={0,127,255}));
  connect(semantics2.port_a, pipeFriction2.port_b) annotation (points=[28,8; 34,
        8; 34,10; 40,10],     style(color=69, rgbcolor={0,127,255}));
  connect(semantics2.port_b, mixingVolume.port_b[1]) 
    annotation (points=[18,8; 14,8; 14,10; 10,10],
                                       style(color=69, rgbcolor={0,127,255}));
  connect(semantics3.port_a, pipeFriction2.port_a) annotation (points=[66,8; 64,
        8; 64,10; 60,10],     style(color=69, rgbcolor={0,127,255}));
  connect(semantics3.port_b, source2.port) annotation (points=[76,8; 80,8; 80,
        -10], style(color=69, rgbcolor={0,127,255}));
  
  // Use plain connections is no new semantics are required for this approach
  if not FluidInterface.usesNewConnectionSemantics then
    connect(pipeFriction1.port_a, prescribedMassFlowRate_TX_A.port) 
                                                                  annotation (
      points=[-60,10; -62,10; -62,12; -80,12; -80,30], style(color=69, rgbcolor=
         {0,127,255}));
    connect(pipeFriction1.port_b, mixingVolume.port_a[1]) 
                                                        annotation (points=[-40,
        10; -34,10; -34,12; -14,12; -14,10; -10,10], style(color=69, rgbcolor={
          0,127,255}));
    connect(mixingVolume.port_b[1], pipeFriction2.port_b) 
                                                        annotation (points=[10,10;
          14,10; 14,12; 34,12; 34,10; 40,10],   style(color=69, rgbcolor={0,127,
          255}));
    connect(pipeFriction2.port_a, source2.port) 
                                              annotation (points=[60,10;
          64,10; 64,12; 80,12; 80,-10],
                               style(color=69, rgbcolor={0,127,255}));
  end if;
end FM_CV_FM;
