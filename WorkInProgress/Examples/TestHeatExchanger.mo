model TestHeatExchanger 
  
extends Modelica.Icons.Example;
replaceable package Medium = Modelica.Media.Water.StandardWater;
  Modelica_Fluid.WorkInProgress.Components.HeatExchanger HEX(
    d_wall=200,
    c_wall=500,
    T_start_w=300,
    use_T_start_1=true,
    p_start_1=1.5e5,
    T_start_1=300,
    use_T_start_2=true,
    p_start_2=1.5e5,
    di_1=0.024,
    da_1=0.03,
    di_2=0.054,
    T_start_2=300,
    n=20,
    length=2,
    lumped_dp=false,
    redeclare package Medium_1 = 
        Medium,
    redeclare package Medium_2 = 
        Medium,
    kineticTerm=false,
    mflow_start_1=0.2,
    static=false,
    redeclare model HeatTransfer = 
        Modelica_Fluid.WorkInProgress.Utilities.PipeHeatTransfer.PipeHT_constAlpha
        (alpha0=1000),
    initOption_1=Modelica_Fluid.Types.Init.InitialValues,
    initOption_2=Modelica_Fluid.Types.Init.InitialValues,
    mflow_start_2=0.2,
      K1=3,
    redeclare model PipeFriction = 
        Modelica_Fluid.WorkInProgress.Utilities.PipeFriction.PipeFriction_SimpleLinear
        (
        m_flow_nominal=0.2,
        dp_nominal=600,
        d_nominal=1000))       annotation (extent=[-26,-12; 34,48]);
  
  Sources.FixedAmbient_pTX ambient2(
    redeclare package Medium = Medium,
    p=1e5,
    T=280)                                                          annotation (extent=[82,-28;
        62,-8]);
  Sources.FixedAmbient_pTX ambient1(
                                   redeclare package Medium=Medium,
    p=1e5,
    T=320)                                                          annotation (extent=[82,24;
        62,44]);
  Sources.PrescribedMassFlowRate_TX massFlowRate2(redeclare package Medium = Medium,
    m_flow=0.2,
    T=360)      annotation (extent=[-66,24; -46,44]);
  Sources.PrescribedMassFlowRate_TX massFlowRate1(redeclare package Medium = Medium,
    m_flow=0.2,
    T=320)       annotation (extent=[-66,-10; -46,10]);
  annotation (Diagram);
  Modelica.Blocks.Sources.Ramp Ramp1(
    height=-0.4,
    duration=10,
    offset=0.2,
    startTime=50) annotation (extent=[-100,24; -80,44]);
equation 
  connect(massFlowRate2.port, HEX.port_a2)            annotation (points=[-45,34;
        -40,34; -40,31.8; -29,31.8],     style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=70,
      rgbfillColor={0,63,125},
      fillPattern=1));
  connect(massFlowRate1.port, HEX.port_a1)            annotation (points=[-45,0;
        -40,0; -40,17.4; -29,17.4], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=70,
      rgbfillColor={0,63,125},
      fillPattern=1));
  connect(HEX.port_b1, ambient1.port)            annotation (points=[37,17.4;
        48.5,17.4; 48.5,34; 61,34], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=70,
      rgbfillColor={0,63,125},
      fillPattern=1));
  connect(HEX.port_b2, ambient2.port)            annotation (points=[37,4.2;
        49.5,4.2; 49.5,-18; 61,-18], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=70,
      rgbfillColor={0,63,125},
      fillPattern=1));
  connect(Ramp1.y, massFlowRate2.m_flow_in) annotation (points=[-79,34; -74,34;
        -74,40; -65.3,40], style(color=74, rgbcolor={0,0,127}));
end TestHeatExchanger;
