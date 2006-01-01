model TestPipe 
  
extends Modelica.Icons.Example;
  replaceable package Medium=Modelica.Media.IdealGases.SingleGases.Air;
  Modelica_Fluid.WorkInProgress.Components.Pipe pipe3(
    redeclare package Medium = Medium,
    use_T_start=true,
    length=5,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
    allowFlowReversal=true,
    d_inner=0.2,
    H_ab=0,
    mflow_start=0.5,
    gravityTerm=false,
    static=false,
    lumped_dp=false,
    redeclare 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PipeHT_constAlpha
      heat(alpha0=2000),
    n=5,
    kineticTerm=false,
    initOption=Modelica_Fluid.Types.InitTypes.SteadyState,
    T_start=400,
    p_start=1.5e5,
    redeclare 
      Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PipeFriction_SimpleLinear
      friction(
      m_flow_nominal=0.5,
      dp_nominal=1000,
      d_nominal=1.2)) 
            annotation (extent=[-42,32; -22,52]);
  annotation (Diagram, experiment(Tolerance=1e-005));
  Sources.FixedAmbient_pTX ambient(redeclare package Medium=Medium,
    p=1e5,
    T=300)                                                          annotation (extent=[70,20;
        50,40]);
  Modelica_Fluid.WorkInProgress.Components.Pipe pipe2(
    redeclare package Medium=Medium,
    use_T_start=true,
    length=5,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
    allowFlowReversal=true,
    d_inner=0.2,
    H_ab=0,
    mflow_start=0.5,
    gravityTerm=false,
    static=false,
    lumped_dp=false,
    redeclare 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PipeHT_constAlpha
      heat(alpha0=2000),
    n=5,
    kineticTerm=false,
    initOption=Modelica_Fluid.Types.InitTypes.SteadyState,
    p_start=1e5,
    T_start=350,
    redeclare 
      Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PipeFriction_SimpleLinear
      friction(
      m_flow_nominal=0.5,
      dp_nominal=1000,
      d_nominal=1.2)) 
            annotation (extent=[8,22; 28,42]);
  Sources.FixedAmbient_pTX ambient1(
                                   redeclare package Medium=Medium,
    T=300,
    p=1.5e5)                                                        annotation (extent=[-88,-14;
        -68,6]);
  Modelica_Fluid.WorkInProgress.Components.Pipe pipe1(
    redeclare package Medium=Medium,
    use_T_start=true,
    length=5,
    crossSectionType=Modelica_Fluid.Types.CrossSectionTypes.Circular,
    allowFlowReversal=true,
    d_inner=0.2,
    H_ab=0,
    mflow_start=0.5,
    gravityTerm=false,
    static=false,
    lumped_dp=false,
    redeclare 
      Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PipeHT_constAlpha
      heat(alpha0=2000),
    n=5,
    kineticTerm=false,
    initOption=Modelica_Fluid.Types.InitTypes.SteadyState,
    p_start=1.5e5,
    T_start=300,
    redeclare 
      Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PipeFriction_SimpleLinear
      friction(
      m_flow_nominal=0.5,
      dp_nominal=1000,
      d_nominal=1.2)) 
            annotation (extent=[-42,8; -22,28]);
  Sources.FixedAmbient_pTX ambient2(
                                   redeclare package Medium=Medium,
    T=400,
    p=1.5e5)                                                        annotation (extent=[-88,54;
        -68,74]);
equation 
  connect(pipe3.port_b, pipe2.port_a) 
                                     annotation (points=[-22,42; -2,42; -2,32;
        7.8,32],
      style(color=69, rgbcolor={0,127,255}));
  connect(pipe1.port_b, pipe2.port_a) annotation (points=[-22,18; -2,18; -2,32;
        7.8,32],  style(color=69, rgbcolor={0,127,255}));
  connect(pipe2.port_b, ambient.port) annotation (points=[28,32; 40,32; 40,30;
        49,30], style(color=69, rgbcolor={0,127,255}));
  connect(ambient2.port, pipe3.port_a) annotation (points=[-67,64; -62,64; -62,
        42; -42.2,42], style(color=69, rgbcolor={0,127,255}));
  connect(ambient1.port, pipe1.port_a) annotation (points=[-67,-4; -62,-4; -62,
        20; -42.2,20; -42.2,18], style(color=69, rgbcolor={0,127,255}));
end TestPipe;
