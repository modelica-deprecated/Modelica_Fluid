model TestValvesReverse "Test case for valves with reverse and zero flow" 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  Sources.FixedAmbient_pTX SourceP1(p=10e5,
  redeclare package Medium = Medium) 
  annotation (extent=[-100,30; -80,50]);
  Sources.FixedAmbient_pTX SourceP2(p=8e5,
  redeclare package Medium = Medium) 
  annotation (extent=[-100, -50; -80, -30]);
  Sources.FixedAmbient_pTX SinkP1(p=1e5,
  redeclare package Medium = Medium) 
  annotation (extent=[82,-4; 62,16]);
  Components.ValveIncompressible V1(
    dpnom=9e5,
    m_flow_nom=1.5,
  redeclare package Medium = Medium,
    pnom=10e5,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    Av=0.1) annotation (extent=[-50, 58; -30, 78]);
  Components.ValveIncompressible V2(
    dpnom=5e5,
    m_flow_nom=1.2,
  redeclare package Medium = Medium,
    pnom=10e5,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    Av=0.1) annotation (extent=[-38, 26; -18, 46]);
  Components.ValveIncompressible V3(
    dpnom=3e5,
    m_flow_nom=1.1,
  redeclare package Medium = Medium,
    pnom=8e5,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint) 
            annotation (extent=[-38, -38; -18, -18]);
  Components.ValveIncompressible V4(
    dpnom=8e5,
    m_flow_nom=1.3,
  redeclare package Medium = Medium,
    pnom=8e5,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint) 
            annotation (extent=[-40,-78; -20,-58]);
  Components.ValveIncompressible V5(
    dpnom=4e5,
    m_flow_nom=2,
  redeclare package Medium = Medium,
    pnom=6e5,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint) 
            annotation (extent=[30,-4; 50,16]);
  
annotation (
  Diagram,
  experiment(StopTime=4, Tolerance=1e-006),
  Documentation(info="<HTML>
<p>This model tests the <tt>ValveLiq</tt> model zero or reverse flow conditions.
<p>Simulate the model for 4 s. At t = 1 s the V5 valve closes in 1 s, the V2 and V3 valves close in 2 s and the V1 and V4 valves open in 2 s. The flow in valve V3 reverses between t = 1.83 and t = 1.93.
<p><b>Revision history:</b></p>
<ul>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco
Casella</a>:<br>
       First release.</li>
</ul>
</HTML>"));
  Sources.FixedAmbient_pTX SinkP2(p=1e5,
  redeclare package Medium = Medium) 
  annotation (extent=[4,58; -16,78]);
  Sources.FixedAmbient_pTX SinkP3(p=1e5, redeclare package Medium = Medium) 
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
equation 
  connect(V1.port_b, SinkP2.port) annotation (points=[-29,68; -17,68]);
  connect(V4.port_b, SinkP3.port) annotation (points=[-19,-68; 5,-68]);
  connect(SourceP1.port, V1.port_a) 
                                   annotation (points=[-79,40; -68,40; -68,68;
        -51,68],
               style(color=69, rgbcolor={0,127,255}));
  connect(SourceP1.port, V2.port_a) 
                                   annotation (points=[-79,40; -60,40; -60,36;
        -39,36],
               style(color=69, rgbcolor={0,127,255}));
  connect(V2.port_b, V5.port_a) 
                             annotation (points=[-17,36; 5,36; 5,6; 29,6],
    style(color=69, rgbcolor={0,127,255}));
  connect(V3.port_b, V5.port_a) 
                             annotation (points=[-17,-28; 6,-28; 6,6; 29,6],
    style(color=69, rgbcolor={0,127,255}));
  connect(SourceP2.port, V4.port_a) 
                                   annotation (points=[-79,-40; -60,-40; -60,
        -68; -41,-68],
                     style(color=69, rgbcolor={0,127,255}));
  connect(SourceP2.port, V3.port_a) 
                                   annotation (points=[-79,-40; -60,-40; -60,
        -28; -39,-28],
                     style(color=69, rgbcolor={0,127,255}));
  connect(OpenRelief.y, V1.stemPosition) annotation (points=[-71,80; -40,80;
        -40,76], style(color=74, rgbcolor={0,0,127}));
  connect(OpenRelief.y, V4.stemPosition) annotation (points=[-71,80; -64,80;
        -64,-52; -30,-52; -30,-60], style(color=74, rgbcolor={0,0,127}));
  connect(CloseValves.y, V2.stemPosition) annotation (points=[-75,-2; -46,-2;
        -46,54; -28,54; -28,44], style(color=74, rgbcolor={0,0,127}));
  connect(CloseValves.y, V3.stemPosition) annotation (points=[-75,-2; -28,-2;
        -28,-20], style(color=74, rgbcolor={0,0,127}));
  connect(CloseLoad.y, V5.stemPosition) annotation (points=[29,36; 40,36; 40,14],
      style(color=74, rgbcolor={0,0,127}));
  connect(V5.port_b, SinkP1.port) 
    annotation (points=[51,6; 61,6], style(color=69, rgbcolor={0,127,255}));
end TestValvesReverse;
