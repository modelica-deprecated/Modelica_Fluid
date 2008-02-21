within FluidSandbox.Examples.VariableConnectionSemantics;
model SimpleTest4 
  
  extends Icons.Example;
  
  FluidSandbox.Sources.PrescribedBoundary_pTX_A source2(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[80,-80; 100,-60],
                                              rotation=90);
  
  FluidSandbox.PressureLosses.WallFriction pipeFriction2(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package FluidInterface = FluidInterface,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent)
                                              annotation (extent=[80,-20; 100,
        -40],
      rotation=270);
  
  FluidSandbox.Volumes.Volume mixingVolume2(
    V=0.1,
    redeclare package Medium = Medium,
    n_b=1,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    p_start=1e5,
    use_T_start=true,
    T_start=340,
    n_a=2,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[60,10; 80,-10]);
  FluidSandbox.Volumes.Volume mixingVolume3(
    V=0.1,
    redeclare package Medium = Medium,
    n_a=1,
    p_start=1e5,
    use_T_start=true,
    T_start=340,
    initType=Modelica_Fluid.Types.Init.NoInit,
    n_b=2,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[-80,10; -60,-10]);
  FluidSandbox.PressureLosses.WallFriction pipeFriction3(
    redeclare package Medium = Medium,
    length=2,
    diameter=0.1,
    redeclare package FluidInterface = FluidInterface,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent)
                                              annotation (extent=[-80,40; -100,
        20],
      rotation=90);
  FluidSandbox.Sources.PrescribedBoundary_pTX_A source3(
    redeclare package Medium = Medium,
    p=1e5,
    usePressureInput=true,
    T=340,
    redeclare package FluidInterface = FluidInterface) 
           annotation (extent=[-100,60; -80,80],
                                              rotation=270);
  Modelica.Blocks.Sources.Sine sine1(
    phase=0,
    freqHz=3,
    amplitude=0.02e5,
    offset=1.02e5) 
              annotation (extent=[-48,82; -60,94], rotation=0);
  annotation (Diagram);
  FluidSandbox.Pipes.AsymmetricDistributedPipe distributedPipe(
    redeclare package Medium = Medium,
    n=4,
    isCircular=true,
    length=2,
    diameter=0.1,
    height_ab=0,
    redeclare package WallFriction = 
        FluidSandbox.Interfaces.PartialModellingApproaches.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent,
    redeclare package FluidInterface = FluidInterface) 
    annotation (extent=[10,10; 30,30]);
  
  FluidSandbox.PressureLosses.WallFriction pipeFriction1(
    redeclare package Medium = Medium,
    diameter=0.1,
    length=1,
    redeclare package FluidInterface = FluidInterface,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent)
                                              annotation (extent=[-10,30; -30,
        10],
      rotation=180);
  
  FluidSandbox.PressureLosses.WallFriction pipeFriction4(
    redeclare package Medium = Medium,
    length=1,
    diameter=0.05,
    redeclare package FluidInterface = FluidInterface,
    redeclare package WallFriction = 
        FluidSandbox.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent)
                                              annotation (extent=[-10,-10; -30,
        -30],
      rotation=180);
  
  FluidSandbox.Pipes.AsymmetricDistributedPipe distributedPipe1(
    redeclare package Medium = Medium,
    n=4,
    isCircular=true,
    length=2,
    height_ab=0,
    redeclare package WallFriction = 
        FluidSandbox.Interfaces.PartialModellingApproaches.PressureLosses.WallFrictionCorrelations.LaminarAndQuadraticTurbulent,
    diameter=0.05,
    redeclare package FluidInterface = FluidInterface) 
                   annotation (extent=[10,-30; 30,-10]);
  
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if        FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[-97,48; -87,50], rotation=90);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics1(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if        FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[-97,8; -87,10], rotation=270);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics2(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if        FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[-36,21; -46,23]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics3(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if        FluidInterface.usesNewConnectionSemantics annotation (extent=[-5,21; 5,
        23]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics4(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if        FluidInterface.usesNewConnectionSemantics annotation (extent=[36,21;
        46,23]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics5(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if        FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[-36,-23; -46,-21]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics6(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if        FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[-5,-23; 5,-21]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics7(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if        FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[36,-23; 46,-21]);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics8(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if        FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[87,-10; 97,-8], rotation=90);
  FluidSandbox.ConnectionSemantics.ConnectionSemantics semantics9(redeclare 
      package Medium = Modelica.Media.Air.DryAirNasa, redeclare package 
      FluidInterface = FluidInterface) if        FluidInterface.usesNewConnectionSemantics 
    annotation (extent=[87,-51; 97,-49], rotation=270);
equation 
  connect(sine1.y, source3.p_in) annotation (points=[-60.6,88; -84,88; -84,82],
                                                                        style(
        color=74, rgbcolor={0,0,127}));
  connect(semantics.port_b, source3.port) annotation (points=[-92,54; -92,58;
        -90,58; -90,60],
      style(color=69, rgbcolor={0,127,255}));
  connect(semantics.port_a, pipeFriction3.port_a) annotation (points=[-92,44;
        -92,42; -90,42; -90,40],
                 style(color=69, rgbcolor={0,127,255}));
  connect(semantics1.port_a, pipeFriction3.port_b) annotation (points=[-92,14;
        -92,18; -90,18; -90,20],
                 style(color=69, rgbcolor={0,127,255}));
  connect(semantics1.port_b, mixingVolume3.port_a[1]) annotation (points=[-92,4;
        -92,0; -80,0], style(color=69, rgbcolor={0,127,255}));
  connect(mixingVolume3.port_b[1], semantics2.port_b) annotation (points=[-60,0.5;
        -50,0.5; -50,22; -46,22],      style(color=69, rgbcolor={0,127,255}));
  connect(semantics2.port_a, pipeFriction1.port_a) annotation (points=[-36,22;
        -34,22; -34,20; -30,20],
                 style(color=69, rgbcolor={0,127,255}));
  connect(semantics3.port_a, pipeFriction1.port_b) 
    annotation (points=[-5,22; -8,22; -8,20; -10,20],
                                        style(color=69, rgbcolor={0,127,255}));
  connect(semantics3.port_b, distributedPipe.port_a) 
    annotation (points=[5,22; 8,22; 8,20; 10,20],
                                      style(color=69, rgbcolor={0,127,255}));
  connect(semantics4.port_a, distributedPipe.port_b) 
    annotation (points=[36,22; 34,22; 34,20; 30,20],
                                       style(color=69, rgbcolor={0,127,255}));
  connect(semantics4.port_b, mixingVolume2.port_a[1]) annotation (points=[46,22;
        50,22; 50,0.5; 60,0.5], style(color=69, rgbcolor={0,127,255}));
  connect(semantics5.port_b, mixingVolume3.port_b[2]) annotation (points=[-46,-22;
        -50,-22; -50,-0.5; -60,-0.5],      style(color=69, rgbcolor={0,127,255}));
  connect(semantics5.port_a, pipeFriction4.port_a) annotation (points=[-36,-22;
        -34,-22; -34,-20; -30,-20],
                  style(color=69, rgbcolor={0,127,255}));
  connect(semantics6.port_a, pipeFriction4.port_b) annotation (points=[-5,-22;
        -8,-22; -8,-20; -10,-20],
                  style(color=69, rgbcolor={0,127,255}));
  connect(semantics6.port_b, distributedPipe1.port_a) 
    annotation (points=[5,-22; 8,-22; 8,-20; 10,-20],
                                        style(color=69, rgbcolor={0,127,255}));
  connect(semantics7.port_a, distributedPipe1.port_b) annotation (points=[36,-22;
        34,-22; 34,-20; 30,-20],
                      style(color=69, rgbcolor={0,127,255}));
  connect(semantics7.port_b, mixingVolume2.port_a[2]) annotation (points=[46,-22;
        50,-22; 50,-0.5; 60,-0.5],      style(color=69, rgbcolor={0,127,255}));
  connect(semantics8.port_b, mixingVolume2.port_b[1]) annotation (points=[92,-4;
        92,0; 80,0], style(color=69, rgbcolor={0,127,255}));
  connect(semantics8.port_a, pipeFriction2.port_b) annotation (points=[92,-14;
        92,-18; 90,-18; 90,-20],
                 style(color=69, rgbcolor={0,127,255}));
  connect(semantics9.port_a, pipeFriction2.port_a) annotation (points=[92,-45;
        92,-42; 90,-42; 90,-40],
                 style(color=69, rgbcolor={0,127,255}));
  connect(semantics9.port_b, source2.port) annotation (points=[92,-55; 92,-58;
        90,-58; 90,-60],
      style(color=69, rgbcolor={0,127,255}));
  // Use plain connections is no new semantics are required for this approach
  if not FluidInterface.usesNewConnectionSemantics then
    connect(source3.port, pipeFriction3.port_a) 
                                              annotation (points=[-90,60; -90,
        58; -88,58; -88,42; -90,42; -90,40], style(color=69, rgbcolor={0,127,
          255}));
    connect(pipeFriction3.port_b, mixingVolume3.port_a[1]) 
                                                         annotation (points=[
        -90,20; -90,18; -88,18; -88,0; -80,0], style(color=69, rgbcolor={0,127,
          255}));
    connect(mixingVolume3.port_b[1], pipeFriction1.port_a) 
                                                         annotation (points=[
        -60,0.5; -50,0.5; -50,18; -34,18; -34,20; -30,20], style(color=69,
        rgbcolor={0,127,255}));
    connect(mixingVolume3.port_b[2], pipeFriction4.port_a) 
                                                         annotation (points=[
        -60,-0.5; -50,-0.5; -50,-18; -34,-18; -34,-20; -30,-20], style(color=69,
        rgbcolor={0,127,255}));
    connect(pipeFriction1.port_b, distributedPipe.port_a) 
                                                        annotation (points=[-10,
        20; -8,20; -8,18; 8,18; 8,20; 10,20], style(color=69, rgbcolor={0,127,
          255}));
    connect(pipeFriction4.port_b, distributedPipe1.port_a) 
                                                         annotation (points=[
        -10,-20; -8,-20; -8,-18; 8,-18; 8,-20; 10,-20], style(color=69,
        rgbcolor={0,127,255}));
    connect(distributedPipe1.port_b, mixingVolume2.port_a[2]) 
                                                            annotation (points=
        [30,-20; 34,-20; 34,-18; 50,-18; 50,-0.5; 60,-0.5], style(color=69,
        rgbcolor={0,127,255}));
    connect(distributedPipe.port_b, mixingVolume2.port_a[1]) 
                                                           annotation (points=[
        30,20; 34,20; 34,18; 50,18; 50,0.5; 60,0.5], style(color=69, rgbcolor={
          0,127,255}));
    connect(mixingVolume2.port_b[1], pipeFriction2.port_b) 
                                                         annotation (points=[80,
        0; 88,0; 88,-18; 90,-18; 90,-20], style(color=69, rgbcolor={0,127,255}));
    connect(pipeFriction2.port_a, source2.port) 
                                              annotation (points=[90,-40; 90,-42; 88,
        -42; 88,-58; 90,-58; 90,-60], style(color=69, rgbcolor={0,127,255}));
  end if;
end SimpleTest4;
