model TestDistributedPipe_01 
  import Modelica_Fluid;
extends Modelica.Icons.Example;
replaceable package Medium=Modelica.Media.Water.StandardWater;
//replaceable package Medium=Modelica.Media.Air.DryAirNasa;  //
  
 Modelica_Fluid.Pipes.DistributedPipe_a_v_b pipe2(
    redeclare package Medium = Medium,
    allowFlowReversal=true,
    use_T_start=true,
    from_dp=true,
    T_start=280,
    diameter=0.01,
    use_d_nominal=false,
    n=5,
    static=false,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    mflow_start=0.1,
    p_a_start=1.0e5,
    p_b_start=1.0e5,
    length=2,
    use_eta_nominal=false,
    use_approxPortProperties=false) 
            annotation (extent=[-32,32; -4,60]);
  
  annotation (Diagram, experiment(StopTime=20, Tolerance=1e-005),
    experimentSetupOutput);
  Modelica_Fluid.Sources.PrescribedBoundary_pTX ambient6(
    redeclare package Medium = Medium,
    p=1e5,
    T=300,
    usePressureInput=true,
    useTemperatureInput=false)                                      annotation (extent=[88,18;
        68,38]);
  Modelica_Fluid.Pipes.DistributedPipe_a_v_b pipe3(
    redeclare package Medium=Medium,
    allowFlowReversal=true,
    T_start=340,
    length=1,
    use_T_start=true,
    from_dp=true,
    p_b_start=1e5,
    diameter=0.01,
    use_d_nominal=false,
    n=5,
    static=false,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    mflow_start=0.1,
    p_a_start=1.0e5,
    use_eta_nominal=false,
    use_approxPortProperties=false) 
            annotation (extent=[22,6; 48,34]);
  
  Modelica_Fluid.Pipes.DistributedPipe_a_v_b pipe1(
    redeclare package Medium=Medium,
    allowFlowReversal=true,
    use_T_start=true,
    from_dp=true,
    T_start=300,
    diameter=0.01,
    use_d_nominal=false,
    n=5,
    static=false,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    mflow_start=0.1,
    p_a_start=1.0e5,
    p_b_start=1.0e5,
    length=0.4,
    use_eta_nominal=false,
    use_approxPortProperties=false) 
            annotation (extent=[-90,6; -62,34]);
  
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient2(
    T=280,
    redeclare package Medium = Medium,
    p=1.5e5)                                                        annotation (extent=[-114,10;
        -94,30]);
  
    annotation (extent=[-90,-86; -70,-66]);
                                     annotation (points=[-22,42; -2,42; -2,32;
        7.8,32],
      style(color=69, rgbcolor={0,127,255}));
  Modelica.Blocks.Sources.Ramp ramp(
    offset=1e5,
    startTime=5,
    duration=0,
    height=1.0e5) 
                annotation (extent=[102,62; 82,82], rotation=0);
  
  Modelica_Fluid.Pipes.DistributedPipe_a_v_b pipe4(
    redeclare package Medium=Medium,
    allowFlowReversal=true,
    length=1,
    use_T_start=true,
    from_dp=true,
    diameter=0.01,
    T_start=360,
    use_d_nominal=false,
    n=5,
    static=false,
    initType=Modelica_Fluid.Types.Init.NoInit,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    mflow_start=0.1,
    p_a_start=1.0e5,
    p_b_start=1.0e5,
    heatTransfer(alpha0=1000),
    use_eta_nominal=false,
    use_approxPortProperties=false) 
            annotation (extent=[-32,-20; -4,8]);
  
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[62,-78; 82,-58]);
  Modelica_Fluid.Junctions.GenericJunction genericJunction(
    n_b=2,
    V=0.0001,
    p_start=1e5,
    T_start=300,
    redeclare package Medium = Medium) annotation (extent=[-54,10; -34,30]);
  Modelica_Fluid.Junctions.GenericJunction genericJunction1(
    n_a=2,
    V=0.0001,
    p_start=1e5,
    T_start=300,
    redeclare package Medium = Medium) annotation (extent=[-4,8; 16,28]);
equation 
  connect(ramp.y, ambient6.p_in) annotation (points=[81,72; 80,72; 80,38; 90,38;
        90,34], style(color=74, rgbcolor={0,0,127}));
  connect(ambient2.port, pipe1.port_a) annotation (points=[-94,20; -90,20],
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe3.port_b, ambient6.port) annotation (points=[48,20; 64,20; 64,28;
        68,28], style(color=69, rgbcolor={0,127,255}));
  connect(pipe1.port_b, genericJunction.ports_a[1]) annotation (points=[-62,20;
        -53,20], style(color=69, rgbcolor={0,127,255}));
  connect(genericJunction.ports_b[1], pipe2.port_a) annotation (points=[-35,18; 
        -35,47; -32,47; -32,46], style(color=3, rgbcolor={0,0,255}));
  connect(genericJunction.ports_b[2], pipe4.port_a) annotation (points=[-35,22; 
        -36,22; -36,-6; -32,-6], style(color=3, rgbcolor={0,0,255}));
  connect(pipe2.port_b, genericJunction1.ports_a[1]) annotation (points=[-4,46; 
        -3,46; -3,20], style(color=69, rgbcolor={0,127,255}));
  connect(genericJunction1.ports_a[2], pipe4.port_b) 
    annotation (points=[-3,16; -4,-6], style(color=3, rgbcolor={0,0,255}));
  connect(genericJunction1.ports_b[1], pipe3.port_a) annotation (points=[15,18;
        18,18; 18,20; 22,20], style(color=3, rgbcolor={0,0,255}));
end TestDistributedPipe_01;
