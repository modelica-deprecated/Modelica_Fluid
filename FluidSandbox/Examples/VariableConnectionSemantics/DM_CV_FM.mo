model DM_CV_FM 
  "Boundary condition (pressure) - distributed model (1D discretized) - control volume - flow model - boundary condition (pressure)" 
  extends Icons.Example;
                          FluidSandbox.Sources.PrescribedBoundary_pTX_A source2
    (
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    redeclare package FluidInterface = FluidInterface,
    redeclare package FluidDiscretization = FluidDiscretization) 
           annotation (extent=[70,-30; 90,-10],
                                              rotation=90);
  
  FluidSandbox.Volumes.Volume mixingVolume(
    V=0.1,
    redeclare package Medium = Medium,
    n_a=1,
    n_b=1,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    p_start=1e5,
    use_T_start=true,
    T_start=340,
    redeclare package FluidInterface = FluidInterface,
    redeclare package FluidDiscretization = FluidDiscretization) 
           annotation (extent=[-10,0; 10,20]);
  FluidSandbox.PressureLosses.WallFriction pipeFriction2(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package FluidInterface = FluidInterface,
    redeclare package WallFriction = 
        PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent,
    redeclare package FluidDiscretization = FluidDiscretization) 
                                              annotation (extent=[40,20; 60,0],
      rotation=180);
  
  FluidSandbox.Sources.PrescribedBoundary_pTX_B prescribedMassFlowRate_TX_A(
    redeclare package Medium = Medium,
    T=320,
    usePressureInput=true,
    redeclare package FluidInterface = FluidInterface,
    redeclare package FluidDiscretization = FluidDiscretization) 
           annotation (extent=[-90,30; -70,50], rotation=270);
  annotation (Diagram);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics(redeclare 
      package Medium = Medium, redeclare package FluidInterface = 
        FluidInterface,
    redeclare package FluidDiscretization = FluidDiscretization) if 
    FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[-75,7; -65,9]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics1(redeclare 
      package Medium = Medium, redeclare package FluidInterface = 
        FluidInterface,
    redeclare package FluidDiscretization = FluidDiscretization) if 
    FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[-28,7; -18,9]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics2(redeclare 
      package Medium = Medium, redeclare package FluidInterface = 
        FluidInterface,
    redeclare package FluidDiscretization = FluidDiscretization) if 
    FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[28,7; 18,9]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics3(redeclare 
      package Medium = Medium, redeclare package FluidInterface = 
        FluidInterface,
    redeclare package FluidDiscretization = FluidDiscretization) if 
    FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[66,7; 76,9]);
  FluidSandbox.Pipes.AsymmetricDistributedPipe distributedPipe(
    redeclare package Medium = Medium,
    isCircular=true,
    length=2,
    diameter=0.1,
    height_ab=0,
    n=10,
    redeclare package WallFriction = 
        PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent,
    from_dp=true,
    redeclare package FluidInterface = FluidInterface,
    redeclare package FluidDiscretization = FluidDiscretization) 
                  annotation (extent=[-60,20; -40,0]);
  
  Modelica.Blocks.Sources.Sine sine1(
    freqHz=2,
    amplitude=0.05e5,
    offset=1e5,
    phase=0.1) 
              annotation (extent=[-84,70; -64,90], rotation=270);
equation 
  connect(semantics.port_a, prescribedMassFlowRate_TX_A.port) annotation (
      points=[-75,8; -80,8; -80,30], style(color=69, rgbcolor={0,127,255}));
  connect(semantics.port_b, distributedPipe.port_a) annotation (points=[-65,8;
        -62,8; -62,10; -60,10], style(color=69, rgbcolor={0,127,255}));
  connect(distributedPipe.port_b, semantics1.port_a) 
                                                   annotation (points=[-40,10;
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
    connect(distributedPipe.port_a, prescribedMassFlowRate_TX_A.port) 
                                                                  annotation (
      points=[-60,10; -62,10; -62,12; -80,12; -80,30], style(color=69, rgbcolor=
         {0,127,255}));
    connect(distributedPipe.port_b, mixingVolume.port_a[1]) 
                                                        annotation (points=[-40,10;
          -34,10; -34,12; -14,12; -14,10; -10,10],   style(color=69, rgbcolor={
          0,127,255}));
    connect(mixingVolume.port_b[1], pipeFriction2.port_b) 
                                                        annotation (points=[10,10; 
          14,10; 14,12; 34,12; 34,10; 40,10],   style(color=69, rgbcolor={0,127,
          255}));
    connect(pipeFriction2.port_a, source2.port) 
                                              annotation (points=[60,10; 64,10; 
          64,12; 80,12; 80,-10],
                               style(color=69, rgbcolor={0,127,255}));
  end if;
  
  connect(sine1.y, prescribedMassFlowRate_TX_A.p_in) annotation (points=[
        -74,69; -74,52], style(color=74, rgbcolor={0,0,127}));
end DM_CV_FM;
