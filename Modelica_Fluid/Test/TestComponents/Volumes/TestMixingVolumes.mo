within Modelica_Fluid.Test.TestComponents.Volumes;
model TestMixingVolumes 
  extends Modelica.Icons.Example;
  package Medium = Modelica.Media.Water.StandardWater;
  annotation (Diagram, experiment(StopTime=10));
  Modelica_Fluid.Volumes.MixingVolume MixingVolume1(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[-30,30; -10,50]);
  
  Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource2(
    m_flow=1,
    h=2e5,
    redeclare package Medium = Medium) 
                   annotation (extent=[-100,30; -80,50]);
  Modelica_Fluid.Volumes.MixingVolume MixingVolume2(
    V=1e-3,
    use_T_start=false,
    h_start=1e5,
    redeclare package Medium = Medium) 
                 annotation (extent=[10,30; 30,50]);
  Modelica_Fluid.Sensors.Temperature Tmix_in(
                                         redeclare package Medium = Medium) 
    annotation (extent=[-60,50; -40,70]);
  Modelica_Fluid.Sensors.Temperature Tmix_out(
                                          redeclare package Medium = Medium) 
    annotation (extent=[40,50; 60,70]);
  Modelica_Fluid.Sources.FixedBoundary_phX Sink2(   p=101325, redeclare package
      Medium = Medium,
    h=Medium.h_default) 
    annotation (extent=[100,30; 80,50]);
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[-100,-100; -80,-80]);
equation 
  connect(MixingVolume1.port_b, MixingVolume2.port_a) annotation (points=[-10,40;
        9.8,40],                                      style(color=69, rgbcolor=
          {0,127,255}));
  connect(FlowSource2.port, MixingVolume1.port_a) annotation (points=[-80,40;
        -30.2,40],  style(color=69, rgbcolor={0,127,255}));
  connect(FlowSource2.port, Tmix_in.port) annotation (points=[-80,40; -50,40;
        -50,50],  style(color=69, rgbcolor={0,127,255}));
  connect(MixingVolume2.port_b, Sink2.port) annotation (points=[30,40; 80,40],
      style(color=69, rgbcolor={0,127,255}));
  connect(MixingVolume2.port_b, Tmix_out.port) annotation (points=[30,40; 50,40;
        50,50],       style(color=69, rgbcolor={0,127,255}));
end TestMixingVolumes;
