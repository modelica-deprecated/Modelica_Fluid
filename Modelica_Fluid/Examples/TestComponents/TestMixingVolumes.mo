model TestMixingVolumes 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  annotation (Diagram, experiment(StopTime=10));
  Modelica_Fluid.Volumes.MixingVolume MixingVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initType=Modelica_Fluid.Types.Init.InitialValues) 
                 annotation (extent=[-30,-40; -10,-20]);
  
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource2(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-100,-40; -80,-20]);
  Modelica_Fluid.Volumes.MixingVolume MixingVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initType=Modelica_Fluid.Types.Init.InitialValues) 
                 annotation (extent=[10,-40; 30,-20]);
  Modelica_Fluid.Sensors.Temperature Tmix_in(
                                         redeclare package Medium = Medium) 
    annotation (extent=[-60,-40; -40,-20]);
  Modelica_Fluid.Sensors.Temperature Tmix_out(
                                          redeclare package Medium = Medium) 
    annotation (extent=[40,-40; 60,-20]);
  BaseClasses.Pipes.PortVolume PortVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initType=Modelica_Fluid.Types.Init.InitialValues) 
                 annotation (extent=[-30,10; -10,30]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource1(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-100,10; -80,30]);
  BaseClasses.Pipes.PortVolume PortVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initType=Modelica_Fluid.Types.Init.InitialValues) 
                 annotation (extent=[10,10; 30,30]);
  Modelica_Fluid.Sensors.Temperature Tport_in(
                                          redeclare package Medium = Medium) 
    annotation (extent=[-60,10; -40,30]);
  Modelica_Fluid.Sensors.Temperature Tport_out(
                                           redeclare package Medium = Medium) 
    annotation (extent=[40,10; 60,30]);
  Modelica_Fluid.Sources.FixedAmbient_phX Sink1(    p=101325, redeclare package
      Medium = Medium) 
    annotation (extent=[100,10; 80,30]);
  Modelica_Fluid.Sources.FixedAmbient_phX Sink2(    p=101325, redeclare package
      Medium = Medium) 
    annotation (extent=[100,-40; 80,-20]);
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[-100,-100; -80,-80]);
equation 
  connect(MixingVolume1.port_b, MixingVolume2.port_a) annotation (points=[-10,-30;
        9.8,-30],                                     style(color=69, rgbcolor=
          {0,127,255}));
  connect(PortVolume1.port,PortVolume2. port) 
    annotation (points=[-20,20; 20,20],
                                     style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource1.port, Tport_in.port_a) annotation (points=[-80,20;
        -60,20],
      style(color=69, rgbcolor={0,127,255}));
  connect(Tport_in.port_b, PortVolume1.port) annotation (points=[-40,20;
        -20,20],
      style(color=69, rgbcolor={0,127,255}));
  connect(PortVolume2.port, Tport_out.port_a) 
    annotation (points=[20,20; 40,20], style(color=69, rgbcolor={0,127,255}));
  connect(Tport_out.port_b, Sink1.port) 
    annotation (points=[60,20; 80,20], style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource2.port, Tmix_in.port_a) annotation (points=[-80,-30;
        -60,-30],
              style(color=69, rgbcolor={0,127,255}));
  connect(Tmix_in.port_b, MixingVolume1.port_a) annotation (points=[-40,-30;
        -30.2,-30], style(color=69, rgbcolor={0,127,255}));
  connect(MixingVolume2.port_b, Tmix_out.port_a) annotation (points=[30,-30;
        40,-30],
              style(color=69, rgbcolor={0,127,255}));
  connect(Tmix_out.port_b, Sink2.port) annotation (points=[60,-30; 80,-30],
      style(color=69, rgbcolor={0,127,255}));
end TestMixingVolumes;
