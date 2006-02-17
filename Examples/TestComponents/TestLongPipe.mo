model TestLongPipe 
  import Modelica_Fluid;
extends Modelica.Icons.Example;
 // replaceable package Medium=Modelica.Media.IdealGases.SingleGases.Air;
replaceable package Medium=Modelica.Media.Air.DryAirNasa;
  
  Modelica_Fluid.Components.LongPipe pipe2(
    redeclare package Medium = Medium,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
    allowFlowReversal=true,
    redeclare Modelica_Fluid.HeatTransfer.PipeHT_constAlpha heat(
           alpha0=2000),
    length=1,
    use_eta_nominal=true,
    use_T_start=true,
    dynamicTerm=false,
    from_dp=true,
    lumped_dp=false,
    p_b_start=1.5e5,
    T_start=280,
    n=10,
    mflow_start=0,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.Utilities.WallFriction.Laminar,
    d_inner=0.01,
    redeclare model Wall = 
        Modelica_Fluid.Components.Wall_constProps (
        d_wall=1000,
        c_wall=800,
        T_start=350,
        initOption=Modelica_Fluid.Types.Init.InitialValues),
    use_wall=false,
    p_a_start=1.2e5,
    kineticTerm=false,
    static=false,
    initOption=Modelica_Fluid.Types.Init.SteadyState) 
            annotation (extent=[-40,42; -20,62]);
  
  annotation (Diagram, experiment(StopTime=20, Tolerance=1e-005),
    experimentSetupOutput);
  Modelica_Fluid.Sources.PrescribedAmbient_pTX ambient(
                                   redeclare package Medium=Medium,
    p=1e5,
    T=300)                                                          annotation (extent=[70,20;
        50,40]);
  Modelica_Fluid.Components.LongPipe pipe5(
    redeclare package Medium=Medium,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
    allowFlowReversal=true,
    redeclare Modelica_Fluid.HeatTransfer.PipeHT_constAlpha heat(
           alpha0=2000),
    T_start=340,
    length=1,
    use_eta_nominal=true,
    use_T_start=true,
    dynamicTerm=false,
    from_dp=true,
    lumped_dp=false,
    p_a_start=1.5e5,
    p_b_start=1e5,
    n=10,
    mflow_start=0,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.Utilities.WallFriction.Laminar,
    d_inner=0.01,
    redeclare model Wall = 
        Modelica_Fluid.Components.Wall_constProps (
        d_wall=1000,
        c_wall=800,
        T_start=350,
        initOption=Modelica_Fluid.Types.Init.InitialValues),
    use_wall=false,
    kineticTerm=false,
    static=false,
    initOption=Modelica_Fluid.Types.Init.SteadyState) 
            annotation (extent=[10,22; 30,42]);
  
  Sources.FixedAmbient_pTX ambient1(
                                   redeclare package Medium=Medium,
    p=2e5,
    T=300)                                                          annotation (extent=[-88,16;
        -68,36]);
  Modelica_Fluid.Components.LongPipe pipe1(
    redeclare package Medium=Medium,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
    allowFlowReversal=true,
    redeclare Modelica_Fluid.HeatTransfer.PipeHT_constAlpha heat(
           alpha0=2000),
    length=1,
    use_eta_nominal=true,
    use_T_start=true,
    dynamicTerm=false,
    from_dp=true,
    lumped_dp=false,
    p_a_start=2e5,
    p_b_start=1.5e5,
    T_start=300,
    n=10,
    mflow_start=0,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.Utilities.WallFriction.Laminar,
    d_inner=0.01,
    redeclare model Wall = 
        Modelica_Fluid.Components.Wall_constProps (
        d_wall=1000,
        c_wall=800,
        T_start=350,
        initOption=Modelica_Fluid.Types.Init.InitialValues),
    use_wall=false,
    kineticTerm=false,
    static=false,
    initOption=Modelica_Fluid.Types.Init.SteadyState) 
            annotation (extent=[-40,16; -20,36]);
  
  Sources.FixedAmbient_pTX ambient2(
                                   redeclare package Medium=Medium,
    T=280,
    p=1.2e5)                                                        annotation (extent=[-88,42;
        -68,62]);
  inner Modelica_Fluid.Components.FluidOptions fluidOptions(default_initOption=
        Modelica_Fluid.Types.Init.SteadyState) 
    annotation (extent=[-90,-86; -70,-66]);
                                     annotation (points=[-22,42; -2,42; -2,32;
        7.8,32],
      style(color=69, rgbcolor={0,127,255}));
  Modelica.Blocks.Sources.Ramp ramp(
    offset=1e5,
    duration=1,
    height=2e5,
      startTime=10) 
                annotation (extent=[102,62; 82,82], rotation=0);
  Sources.FixedAmbient_pTX ambient3(
                                   redeclare package Medium=Medium,
    p=2e5,
    T=330)                                                          annotation (extent=[-88,-12;
        -68,8]);
  Modelica_Fluid.Components.LongPipe pipe3(
    redeclare package Medium=Medium,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
    allowFlowReversal=true,
    redeclare Modelica_Fluid.HeatTransfer.PipeHT_constAlpha heat(
           alpha0=2000),
    length=1,
    use_eta_nominal=true,
    use_T_start=true,
    dynamicTerm=false,
    from_dp=true,
    lumped_dp=false,
    p_a_start=2e5,
    p_b_start=1.5e5,
    T_start=340,
    n=10,
    mflow_start=0,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.Utilities.WallFriction.Laminar,
    d_inner=0.01,
    redeclare model Wall = 
        Modelica_Fluid.Components.Wall_constProps (
        d_wall=1000,
        c_wall=800,
        T_start=350,
        initOption=Modelica_Fluid.Types.Init.InitialValues),
    use_wall=false,
    kineticTerm=false,
    static=false,
    initOption=Modelica_Fluid.Types.Init.SteadyState) 
            annotation (extent=[-40,-12; -20,8]);
  Sources.FixedAmbient_pTX ambient4(
                                   redeclare package Medium=Medium,
    p=2e5,
    T=360)                                                          annotation (extent=[-88,-40;
        -68,-20]);
  Modelica_Fluid.Components.LongPipe pipe4(
    redeclare package Medium=Medium,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
    allowFlowReversal=true,
    redeclare Modelica_Fluid.HeatTransfer.PipeHT_constAlpha heat(
           alpha0=2000),
    length=1,
    use_eta_nominal=true,
    use_T_start=true,
    dynamicTerm=false,
    from_dp=true,
    lumped_dp=false,
    p_a_start=2e5,
    p_b_start=1.5e5,
    n=10,
    mflow_start=0,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.Utilities.WallFriction.Laminar,
    d_inner=0.01,
    redeclare model Wall = 
        Modelica_Fluid.Components.Wall_constProps (
        d_wall=1000,
        c_wall=800,
        T_start=350,
        initOption=Modelica_Fluid.Types.Init.InitialValues),
    use_wall=false,
    kineticTerm=false,
    T_start=360,
    static=false,
    initOption=Modelica_Fluid.Types.Init.SteadyState) 
            annotation (extent=[-40,-40; -20,-20]);
equation 
  connect(ambient2.port,pipe2. port_a) annotation (points=[-68,52; -40.2,52],
                       style(color=69, rgbcolor={0,127,255}));
  connect(ambient1.port, pipe1.port_a) annotation (points=[-68,26; -40.2,26],
                       style(color=69, rgbcolor={0,127,255}));
  connect(pipe2.port_b,pipe5. port_a) annotation (points=[-20,52; -10,52; -10,
          32; 9.8,32],
                     style(color=69, rgbcolor={0,127,255}));
  connect(pipe1.port_b,pipe5. port_a) annotation (points=[-20,26; -4,26; -4,32;
          9.8,32],   style(color=69, rgbcolor={0,127,255}));
  connect(pipe5.port_b, ambient.port) annotation (points=[30,32; 40,32; 40,30;
        50,30], style(color=69, rgbcolor={0,127,255}));
  connect(ramp.y, ambient.p_in) annotation (points=[81,72; 76,72; 76,36; 72,36],
      style(color=74, rgbcolor={0,0,127}));
  connect(ambient3.port, pipe3.port_a) annotation (points=[-68,-2; -40.2,-2],
                   style(color=69, rgbcolor={0,127,255}));
  connect(pipe5.port_a, pipe3.port_b) annotation (points=[9.8,32; -4,32; -4,-2;
          -20,-2],
                 style(color=69, rgbcolor={0,127,255}));
  connect(ambient4.port, pipe4.port_a) annotation (points=[-68,-30; -40.2,-30],
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe5.port_a, pipe4.port_b) annotation (points=[9.8,32; 6,32; 6,-30;
          -20,-30],
                  style(color=69, rgbcolor={0,127,255}));
end TestLongPipe;
