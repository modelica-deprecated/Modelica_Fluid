model TestValvesReverse "Test case for valves with reverse and zero flow" 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP1(
                                               p=10e5,
  redeclare package Medium = Medium, 
    T=ambient.default_T_ambient) 
  annotation (extent=[-100,26; -80,46]);
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP2(
                                               p=8e5,
  redeclare package Medium = Medium, 
    T=ambient.default_T_ambient) 
  annotation (extent=[-100, -50; -80, -30]);
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP1(
                                             p=1e5,
  redeclare package Medium = Medium, 
    T=ambient.default_T_ambient) 
  annotation (extent=[82,-4; 62,16]);
  Modelica_Fluid.ControlValves.ValveIncompressible V1(
    dp_nom=9e5,
    m_flow_nom=1.5,
  redeclare package Medium = Medium,
    p_nom=10e5,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    Av=0.1) annotation (extent=[-50, 58; -30, 78]);
  Modelica_Fluid.ControlValves.ValveIncompressible V2(
    dp_nom=5e5,
    m_flow_nom=1.2,
  redeclare package Medium = Medium,
    p_nom=10e5,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    Av=0.1) annotation (extent=[-38, 26; -18, 46]);
  Modelica_Fluid.ControlValves.ValveIncompressible V3(
    dp_nom=3e5,
    m_flow_nom=1.1,
  redeclare package Medium = Medium,
    p_nom=8e5,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint) 
            annotation (extent=[-38, -38; -18, -18]);
  Modelica_Fluid.ControlValves.ValveIncompressible V4(
    dp_nom=8e5,
    m_flow_nom=1.3,
  redeclare package Medium = Medium,
    p_nom=8e5,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint) 
            annotation (extent=[-40,-78; -20,-58]);
  Modelica_Fluid.ControlValves.ValveIncompressible V5(
    dp_nom=4e5,
    m_flow_nom=2,
  redeclare package Medium = Medium,
    p_nom=6e5,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint) 
            annotation (extent=[30,-4; 50,16]);
  
annotation (
  Diagram,
  experiment(StopTime=4, Tolerance=1e-006),
  Documentation(info=""));
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP2(
                                             p=1e5,
  redeclare package Medium = Medium, 
    T=ambient.default_T_ambient) 
  annotation (extent=[4,58; -16,78]);
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP3(
                                             p=1e5, redeclare package Medium = Medium, 
    T=ambient.default_T_ambient) 
  annotation (extent=[26,-78; 6,-58]);
  Modelica.Blocks.Sources.Ramp CloseLoad(
    duration=1,
    offset=1,
    startTime=1,
    height=-0.99) annotation (extent=[8,26; 28,46]);
  Modelica.Blocks.Sources.Ramp OpenRelief(
    duration=2,
    height=1,
    offset=0,
    startTime=1) 
              annotation (extent=[-92,70; -72,90]);
  Modelica.Blocks.Sources.Ramp CloseValves(
    duration=2,
    offset=1,
    startTime=1,
    height=-1) 
              annotation (extent=[-96, -12; -76, 8]);
  
  inner Modelica_Fluid.Ambient ambient 
                                   annotation (extent=[60,68; 80,88]);
equation 
  connect(V1.port_b, SinkP2.port) annotation (points=[-30,68; -16,68]);
  connect(V4.port_b, SinkP3.port) annotation (points=[-20,-68; 6,-68]);
  connect(SourceP1.port, V1.port_a) 
                                   annotation (points=[-80,36; -68,36; -68,68; 
        -50,68],
               style(color=69, rgbcolor={0,127,255}));
  connect(SourceP1.port, V2.port_a) 
                                   annotation (points=[-80,36; -38,36],
               style(color=69, rgbcolor={0,127,255}));
  connect(V2.port_b, V5.port_a) 
                             annotation (points=[-18,36; 5,36; 5,6; 30,6],
    style(color=69, rgbcolor={0,127,255}));
  connect(V3.port_b, V5.port_a) 
                             annotation (points=[-18,-28; 6,-28; 6,6; 30,
        6],
    style(color=69, rgbcolor={0,127,255}));
  connect(SourceP2.port, V4.port_a) 
                                   annotation (points=[-80,-40; -60,-40;
        -60,-68; -40,-68],
                     style(color=69, rgbcolor={0,127,255}));
  connect(SourceP2.port, V3.port_a) 
                                   annotation (points=[-80,-40; -60,-40;
        -60,-28; -38,-28],
                     style(color=69, rgbcolor={0,127,255}));
  connect(OpenRelief.y, V1.stemPosition) annotation (points=[-71,80; -40,
        80; -40,77],
                 style(color=74, rgbcolor={0,0,127}));
  connect(OpenRelief.y, V4.stemPosition) annotation (points=[-71,80; -64,
        80; -64,-52; -30,-52; -30,-59],
                                    style(color=74, rgbcolor={0,0,127}));
  connect(CloseValves.y, V2.stemPosition) annotation (points=[-75,-2; -46,
        -2; -46,54; -28,54; -28,45],
                                 style(color=74, rgbcolor={0,0,127}));
  connect(CloseValves.y, V3.stemPosition) annotation (points=[-75,-2; -28,
        -2; -28,-19],
                  style(color=74, rgbcolor={0,0,127}));
  connect(CloseLoad.y, V5.stemPosition) annotation (points=[29,36; 40,36;
        40,15],
      style(color=74, rgbcolor={0,0,127}));
  connect(V5.port_b, SinkP1.port) 
    annotation (points=[50,6; 62,6], style(color=69, rgbcolor={0,127,255}));
end TestValvesReverse;
