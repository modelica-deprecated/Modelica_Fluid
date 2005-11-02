model TestValveIncompressible "Test case for valves" 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  Sources.FixedAmbient_pTX SourceP1(p=10e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater) 
  annotation (extent=[-100,30; -80,50]);
  Components.ValveIncompressible V1(
    d_nom=1000,
    dp_nom=9e5,
    m_flow_nom=1.5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.Cv,
    Cv=10,
    p_nom=10e5) 
            annotation (extent=[-50,30; -30,50]);
  
annotation (
  Diagram,
  experiment(StopTime=4, Tolerance=1e-006),
  Documentation(info=""));
  Sources.FixedAmbient_pTX SinkP2(p=1e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater) 
  annotation (extent=[22,30; 2,50]);
  Modelica.Blocks.Sources.Ramp Opening(
    duration=2,
    height=1,
    offset=0,
    startTime=1) 
              annotation (extent=[-92, 74; -72, 94]);
  Sources.FixedAmbient_pTX SourceP2(p=10e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater) 
  annotation (extent=[-100,-10; -80,10]);
  Components.ValveIncompressible V2(
    d_nom=1000,
    dp_nom=9e5,
    m_flow_nom=1.5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.Cv,
    Cv=10,
    p_nom=10e5,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Components.ValveCharacteristics.equalPercentage) 
            annotation (extent=[-50,-10; -30,10]);
  Sources.FixedAmbient_pTX SinkP1(p=1e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater) 
  annotation (extent=[22,-10; 2,10]);
  Sources.FixedAmbient_pTX SourceP3(p=10e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater) 
  annotation (extent=[-100,-50; -80,-30]);
  Components.ValveIncompressible V3(
    d_nom=1000,
    dp_nom=9e5,
    m_flow_nom=1.5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.Cv,
    Cv=10,
    p_nom=10e5,
    redeclare function flowCharacteristic = 
        Components.ValveCharacteristics.equalPercentage (rangeability=10)) 
            annotation (extent=[-50,-50; -30,-30]);
  Sources.FixedAmbient_pTX SinkP3(p=1e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater) 
  annotation (extent=[22,-50; 2,-30]);
equation 
  connect(V1.port_b, SinkP2.port) annotation (points=[-29,40; 1,40]);
  connect(Opening.y, V1.stemPosition) 
  annotation (points=[-71,84; -40,84; -40,48],    style(color=3));
  connect(SourceP1.port, V1.port_a) 
                                   annotation (points=[-79,40; -51,40],
               style(color=69, rgbcolor={0,127,255}));
  connect(Opening.y, V2.stemPosition) annotation (points=[-71,84; -64,84; -64,
        20; -40,20; -40,8], style(color=74, rgbcolor={0,0,127}));
  connect(Opening.y, V3.stemPosition) annotation (points=[-71,84; -64,84; -64,
        -22; -40,-22; -40,-32], style(color=74, rgbcolor={0,0,127}));
  connect(SourceP2.port, V2.port_a) 
    annotation (points=[-79,0; -51,0], style(color=69, rgbcolor={0,127,255}));
  connect(V2.port_b, SinkP1.port) annotation (points=[-29,0; -14,0; -14,0; 1,0],
      style(color=69, rgbcolor={0,127,255}));
  connect(SourceP3.port, V3.port_a) annotation (points=[-79,-40; -51,-40],
      style(color=69, rgbcolor={0,127,255}));
  connect(V3.port_b, SinkP3.port) annotation (points=[-29,-40; 1,-40], style(
        color=69, rgbcolor={0,127,255}));
end TestValveIncompressible;