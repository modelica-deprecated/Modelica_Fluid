model TestNewMixingVolume 
  extends Modelica.Icons.Example;
  Components.MixingVolume Volume(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    initOption=Modelica_Fluid.Types.InitTypes.SteadyState,
    p_start=2e5,
    V=1,
    use_T_start=false,
    h_start=3e6) 
         annotation (extent=[-42,0; -22,20]);
  Sources.PrescribedMassFlowRate_hX FlowSource(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow=1,
    h=3e6) annotation (extent=[-82,0; -62,20]);
  Sources.FixedAmbient_pTX Sink(redeclare package Medium = 
        Modelica.Media.Water.StandardWater, p=101325) 
    annotation (extent=[60,0; 40,20]);
  Components.ValveLinear Valve(redeclare package Medium = 
        Modelica.Media.Water.StandardWater, Kv=1) 
                                            annotation (extent=[2,0; 22,20]);
  annotation (Diagram, experiment(StopTime=5));
  Modelica.Blocks.Sources.Step Step1(
    startTime=1,
    height=-0.5,
    offset=1) annotation (extent=[-36,48; -16,68]);
equation 
  connect(FlowSource.port, Volume.port_a) annotation (points=[-61,10; -42.2,10],
      style(color=69, rgbcolor={0,127,255}));
  connect(Volume.port_b, Valve.port_a) 
    annotation (points=[-22,10; 1,10], style(color=69, rgbcolor={0,127,255}));
  connect(Valve.port_b, Sink.port) 
    annotation (points=[23,10; 39,10], style(color=69, rgbcolor={0,127,255}));
  connect(Step1.y, Valve.opening) annotation (points=[-15,58; 12,58; 12,18],
      style(color=74, rgbcolor={0,0,127}));
end TestNewMixingVolume;
