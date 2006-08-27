model TestDistributedPipe 
  import Modelica_Fluid;
extends Modelica.Icons.Example;
  
replaceable package Medium=Modelica.Media.Water.StandardWater;
  
  Modelica_Fluid.Pipes.DistributedPipe pipe2(
    redeclare package Medium = Medium,
    allowFlowReversal=true,
    length=1,
    use_eta_nominal=true,
    use_T_start=true,
    from_dp=true,
    T_start=280,
    mflow_start=0,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Laminar,
    diameter=0.01,
    p_a_start=1.2e5,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    kineticTerm=false,
    use_d_nominal=false,
    n=5,
    static=false,
    p_b_start=1.1e5) 
            annotation (extent=[-38,42; -18,62]);
  
  annotation (Diagram, experiment(StopTime=20, Tolerance=1e-005),
    experimentSetupOutput);
  Modelica_Fluid.Sources.PrescribedBoundary_pTX ambient6(
    redeclare package Medium = Medium,
    p=1e5,
    T=300)                                                          annotation (extent=[68,20; 
        48,40]);
  Modelica_Fluid.Pipes.DistributedPipe pipe5(
    redeclare package Medium=Medium,
    allowFlowReversal=true,
    T_start=340,
    length=1,
    use_eta_nominal=true,
    use_T_start=true,
    from_dp=true,
    p_b_start=1e5,
    mflow_start=0,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Laminar,
    diameter=0.01,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    kineticTerm=false,
    use_d_nominal=false,
    n=5,
    static=false,
    p_a_start=1.1e5) 
            annotation (extent=[12,22; 32,42]);
  
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient1(
                                   redeclare package Medium=Medium,
    T=300,
    p=1.2e5)                                                        annotation (extent=[-88,16;
        -68,36]);
  Modelica_Fluid.Pipes.DistributedPipe pipe1(
    redeclare package Medium=Medium,
    allowFlowReversal=true,
    length=1,
    use_eta_nominal=true,
    use_T_start=true,
    from_dp=true,
    p_a_start=2e5,
    T_start=300,
    mflow_start=0,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Laminar,
    diameter=0.01,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    kineticTerm=false,
    use_d_nominal=false,
    n=5,
    static=false,
    p_b_start=1.1e5) 
            annotation (extent=[-38,16; -18,36]);
  
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient2(
    T=280,
    redeclare package Medium = Medium,
    p=1.2e5)                                                        annotation (extent=[-88,42;
        -68,62]);
  
    annotation (extent=[-90,-86; -70,-66]);
                                     annotation (points=[-22,42; -2,42; -2,32;
        7.8,32],
      style(color=69, rgbcolor={0,127,255}));
  Modelica.Blocks.Sources.Ramp ramp(
    offset=1e5,
    duration=10,
    startTime=5,
    height=0.4e5) 
                annotation (extent=[102,62; 82,82], rotation=0);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient3(
                                   redeclare package Medium=Medium,
    T=330,
    p=1.2e5)                                                        annotation (extent=[-88,-12;
        -68,8]);
  Modelica_Fluid.Pipes.DistributedPipe pipe3(
    redeclare package Medium=Medium,
    allowFlowReversal=true,
    length=1,
    use_eta_nominal=true,
    use_T_start=true,
    from_dp=true,
    p_a_start=2e5,
    T_start=340,
    mflow_start=0,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Laminar,
    diameter=0.01,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    kineticTerm=false,
    use_d_nominal=false,
    n=5,
    static=false,
    p_b_start=1.1e5) 
            annotation (extent=[-38,-12; -18,8]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient4(
                                   redeclare package Medium=Medium,
    T=360,
    p=1.2e5)                                                        annotation (extent=[-88,-40;
        -68,-20]);
  Modelica_Fluid.Pipes.DistributedPipe pipe4(
    redeclare package Medium=Medium,
    allowFlowReversal=true,
    length=1,
    use_eta_nominal=true,
    use_T_start=true,
    from_dp=true,
    p_a_start=2e5,
    mflow_start=0,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Laminar,
    diameter=0.01,
    T_start=360,
    initType=Modelica_Fluid.Types.Init.InitialValues,
    kineticTerm=false,
    use_d_nominal=false,
    n=5,
    static=false,
    p_b_start=1.1e5) 
            annotation (extent=[-38,-40; -18,-20]);
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[62,-78; 82,-58]);
equation 
  connect(ambient2.port,pipe2. port_a) annotation (points=[-68,52; -38,52],
                       style(color=69, rgbcolor={0,127,255}));
  connect(ambient1.port, pipe1.port_a) annotation (points=[-68,26; -38,26],
                       style(color=69, rgbcolor={0,127,255}));
  connect(pipe5.port_b, ambient6.port) 
                                      annotation (points=[32,32; 40,32; 40,30; 
        48,30], style(color=69, rgbcolor={0,127,255}));
  connect(ramp.y, ambient6.p_in) 
                                annotation (points=[81,72; 76,72; 76,36; 70,36],
      style(color=74, rgbcolor={0,0,127}));
  connect(ambient3.port, pipe3.port_a) annotation (points=[-68,-2; -38,-2],
                   style(color=69, rgbcolor={0,127,255}));
  connect(ambient4.port, pipe4.port_a) annotation (points=[-68,-30; -38,-30],
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe2.port_b, pipe5.port_a) annotation (points=[-18,52; 0,52; 0,32;
        12,32], style(color=69, rgbcolor={0,127,255}));
  connect(pipe1.port_b, pipe5.port_a) annotation (points=[-18,26; 0,26; 0,32;
        12,32], style(color=69, rgbcolor={0,127,255}));
  connect(pipe3.port_b, pipe5.port_a) annotation (points=[-18,-2; 0,-2; 0,32;
        12,32], style(color=69, rgbcolor={0,127,255}));
  connect(pipe4.port_b, pipe5.port_a) annotation (points=[-18,-30; 0,-30; 0,32;
        12,32], style(color=69, rgbcolor={0,127,255}));
end TestDistributedPipe;
