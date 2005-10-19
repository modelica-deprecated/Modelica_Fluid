model TestMixingVolumes 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  annotation (Diagram, experiment(StopTime=10));
  Modelica_Fluid.Components.MixingVolume MixingVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initOption=Modelica_Fluid.Types.InitTypes.InitialStates) 
                 annotation (extent=[-40,-70; -20,-50]);
  
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource2(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-86,-70; -66,-50]);
  Modelica_Fluid.Components.MixingVolume MixingVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[0,-70; 20,-50]);
  Sensors.Temperature Tmix_in(redeclare package Medium = Medium) 
    annotation (extent=[-64,-36; -44,-16]);
  Sensors.Temperature Tmix_out(redeclare package Medium = Medium) 
    annotation (extent=[34,-38; 54,-18]);
  Components.PortVolume PortVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium,
    initOption=Modelica_Fluid.Types.InitTypes.InitialStates) 
                 annotation (extent=[-40,10; -20,30]);
  Sources.PrescribedMassFlowRate_hX FlowSource1(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-88,10; -68,30]);
  Components.PortVolume PortVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[0,10; 20,30]);
  Sensors.Temperature Tport_in(redeclare package Medium = Medium) 
    annotation (extent=[-64,50; -44,70]);
  Sensors.Temperature Tport_out(redeclare package Medium = Medium) 
    annotation (extent=[34,50; 54,70]);
equation 
  connect(FlowSource2.port, MixingVolume1.port_a) 
                                                annotation (points=[-65,-60;
        -40.2,-60],                           style(color=69, rgbcolor={0,127,
          255}));
  connect(MixingVolume1.port_b, MixingVolume2.port_a) annotation (points=[-20,-60;
        -0.2,-60],                                    style(color=69, rgbcolor=
          {0,127,255}));
  connect(Tmix_in.port, MixingVolume1.port_a) annotation (points=[-54,-37; -54,
        -60; -40.2,-60],     style(color=69, rgbcolor={0,127,255}));
  connect(Tmix_out.port, MixingVolume2.port_b) annotation (points=[44,-39;
        44,-60; 20,-60],    style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource1.port,PortVolume1. port) 
    annotation (points=[-67,20; -30,20],
                                       style(color=69, rgbcolor={0,127,255}));
  connect(PortVolume1.port,PortVolume2. port) 
    annotation (points=[-30,20; 10,20],
                                     style(color=69, rgbcolor={0,127,255}));
  connect(Tport_out.port, PortVolume2.port) annotation (points=[44,49; 44,20; 10,
        20], style(color=69, rgbcolor={0,127,255}));
  connect(Tport_in.port, FlowSource1.port) annotation (points=[-54,49; -54,20;
        -67,20], style(color=69, rgbcolor={0,127,255}));
end TestMixingVolumes;
