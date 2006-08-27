model TestValveIncompressible "Test case for valves" 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP1(
                                               p=10e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater, 
    T=ambient.default_T_ambient) 
  annotation (extent=[-100,30; -80,50]);
  Modelica_Fluid.ControlValves.ValveIncompressible V1(
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
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP2(
                                             p=1e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater, 
    T=ambient.default_T_ambient) 
  annotation (extent=[22,30; 2,50]);
  Modelica.Blocks.Sources.Ramp Opening(
    duration=2,
    height=1,
    offset=0,
    startTime=1) 
              annotation (extent=[-92, 74; -72, 94]);
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP2(
                                               p=10e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater, 
    T=ambient.default_T_ambient) 
  annotation (extent=[-100,-10; -80,10]);
  Modelica_Fluid.ControlValves.ValveIncompressible V2(
    d_nom=1000,
    dp_nom=9e5,
    m_flow_nom=1.5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.Cv,
    Cv=10,
    p_nom=10e5,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.equalPercentage)
            annotation (extent=[-50,-10; -30,10]);
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP1(
                                             p=1e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater, 
    T=ambient.default_T_ambient) 
  annotation (extent=[22,-10; 2,10]);
  Modelica_Fluid.Sources.FixedBoundary_pTX SourceP3(
                                               p=10e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater, 
    T=ambient.default_T_ambient) 
  annotation (extent=[-100,-50; -80,-30]);
  Modelica_Fluid.ControlValves.ValveIncompressible V3(
    d_nom=1000,
    dp_nom=9e5,
    m_flow_nom=1.5,
  redeclare package Medium = Modelica.Media.Water.StandardWater,
    CvData=Modelica_Fluid.Types.CvTypes.Cv,
    Cv=10,
    p_nom=10e5,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.equalPercentage
        (                                                                              rangeability=10)) 
            annotation (extent=[-50,-50; -30,-30]);
  Modelica_Fluid.Sources.FixedBoundary_pTX SinkP3(
                                             p=1e5,
  redeclare package Medium = Modelica.Media.Water.StandardWater, 
    T=ambient.default_T_ambient) 
  annotation (extent=[22,-50; 2,-30]);
  
  inner Modelica_Fluid.Ambient ambient 
                                   annotation (extent=[58,72; 78,92]);
equation 
  connect(V1.port_b, SinkP2.port) annotation (points=[-30,40; 2,40]);
  connect(Opening.y, V1.stemPosition) 
  annotation (points=[-71,84; -40,84; -40,49],    style(color=3));
  connect(SourceP1.port, V1.port_a) 
                                   annotation (points=[-80,40; -50,40],
               style(color=69, rgbcolor={0,127,255}));
  connect(Opening.y, V2.stemPosition) annotation (points=[-71,84; -64,84;
        -64,20; -40,20; -40,9],
                            style(color=74, rgbcolor={0,0,127}));
  connect(Opening.y, V3.stemPosition) annotation (points=[-71,84; -64,84;
        -64,-22; -40,-22; -40,-31],
                                style(color=74, rgbcolor={0,0,127}));
  connect(SourceP2.port, V2.port_a) 
    annotation (points=[-80,0; -50,0], style(color=69, rgbcolor={0,127,255}));
  connect(V2.port_b, SinkP1.port) annotation (points=[-30,0; 2,0],
      style(color=69, rgbcolor={0,127,255}));
  connect(SourceP3.port, V3.port_a) annotation (points=[-80,-40; -50,-40],
      style(color=69, rgbcolor={0,127,255}));
  connect(V3.port_b, SinkP3.port) annotation (points=[-30,-40; 2,-40], style(
        color=69, rgbcolor={0,127,255}));
end TestValveIncompressible;
