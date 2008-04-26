model TestPortVolumes 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[-30,-10; -10,10]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource1(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-100,-10; -80,10]);
  Modelica_Fluid.Sources.FixedBoundary_phX Sink1(
                                         p=101325, redeclare package Medium = 
               Medium,
    h=Medium.h_default) 
    annotation (extent=[100,-10; 80,10]);
  Modelica_Fluid.Pipes.BaseClasses.PortVolume PortVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[10,-10; 30,10]);
  Modelica_Fluid.Sensors.TemperatureOnePort Tin(
                                     redeclare package Medium = Medium) 
    annotation (extent=[-60,10; -40,30]);
  Modelica_Fluid.Sensors.TemperatureOnePort Tout(
                                      redeclare package Medium = Medium) 
    annotation (extent=[40,10; 60,30]);
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[-100,-100; -80,-80]);
equation 
  connect(PortVolume1.port, PortVolume2.port) 
    annotation (points=[-20,0; 20,0],style(color=69, rgbcolor={0,127,255}));
  annotation (Diagram);
  connect(PortVolume2.port, Sink1.port) 
    annotation (points=[20,0; 80,0], style(color=69, rgbcolor={0,127,255}));
  connect(PortVolume2.port, Tout.port) annotation (points=[20,0; 50,0; 50,10],
      style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource1.port, PortVolume1.port) 
    annotation (points=[-80,0; -20,0], style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource1.port, Tin.port) annotation (points=[-80,0; -50,0; -50,10],
      style(color=69, rgbcolor={0,127,255}));
end TestPortVolumes;
