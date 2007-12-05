within Modelica_Fluid.Test.TestComponents.Volumes;
model TestNewMixingVolume 
  extends Modelica.Icons.Example;
  Modelica_Fluid.Volumes.MixingVolume Volume(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    initType=Modelica_Fluid.Types.Init.SteadyState,
    V=1,
    use_T_start=false,
    h_start=3e6) 
         annotation (extent=[-42,0; -22,20]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    m_flow=1,
    h=3e6) annotation (extent=[-82,0; -62,20]);
  Modelica_Fluid.Sources.FixedBoundary_pTX Sink(
                                           redeclare package Medium = 
        Modelica.Media.Water.StandardWater, p=101325,
    T=ambient.default_T_ambient) 
    annotation (extent=[60,0; 40,20]);
  Modelica_Fluid.ControlValves.ValveLinear Valve(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater, Kv=1,
    p_a_start=2e5,
    p_b_start=101325,
    T_start=ambient.default_T_ambient)      annotation (extent=[2,0; 22,20]);
  annotation (Diagram, experiment(StopTime=5));
  Modelica.Blocks.Sources.Step Step1(
    startTime=1,
    height=-0.5,
    offset=1) annotation (extent=[-40,48; -20,68]);
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[-100,-100; -80,-80]);
equation 
  connect(FlowSource.port, Volume.port_a) annotation (points=[-62,10;
        -42.2,10],
      style(color=69, rgbcolor={0,127,255}));
  connect(Volume.port_b, Valve.port_a) 
    annotation (points=[-22,10; 2,10], style(color=69, rgbcolor={0,127,255}));
  connect(Valve.port_b, Sink.port) 
    annotation (points=[22,10; 40,10], style(color=69, rgbcolor={0,127,255}));
  connect(Step1.y, Valve.opening) annotation (points=[-19,58; 12,58; 12,19],
      style(color=74, rgbcolor={0,0,127}));
end TestNewMixingVolume;
