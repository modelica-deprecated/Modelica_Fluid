model TestIsolatedPipe "Test ShortPipe component" 
  import Modelica.SIunits.Conversions.*;
  
  extends Modelica.Icons.Example;
  annotation (
    Diagram(Text(
        extent=[-89,43; -43,35], 
        string="water", 
        style(color=0, rgbcolor={0,0,0})), Text(
        extent=[-82,-49; -36,-57], 
        style(color=0, rgbcolor={0,0,0}), 
        string="DetailedAir")),
    experiment(StopTime=3),
    Coordsys(grid=[1, 1], component=[20, 20]));
  Modelica_Fluid.Components.IsolatedPipe isolatedPipe(
    L=1,
    dp_nominal=1.e4,
    A_a=0.02*0.02,
    A_b=0.02*0.02,
    T_start=Modelica.SIunits.Conversions.from_degC(20),
    initType=Modelica_Fluid.Types.InitTypes.InitialStates,
    m_flow_nominal=0.1,
    nVolumes=10,
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
    annotation (extent=[-9,0; 11,20]);
  
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX MassFlowSource(T_ambient=
        from_degC(30), redeclare package Medium = Modelica_Media.Air.DryAirNasa)
    annotation (extent=[-50,0; -30,20]);
  Modelica_Fluid.Sources.FixedAmbient_pTX ambient(T_ambient=from_degC(15),
      redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
    annotation (extent=[50,0; 30,20]);
  Modelica.Blocks.Sources.Ramp ramp(
    duration=3,
    height=3,
    offset=0.01) annotation (extent=[-90,0; -70,20]);
  Modelica_Fluid.Components.IsolatedPipe isolatedPipe1(
    L=1,
    dp_nominal=1.e4,
    A_a=0.02*0.02,
    A_b=0.02*0.02,
    T_start=Modelica.SIunits.Conversions.from_degC(20),
    initType=Modelica_Fluid.Types.InitTypes.InitialStates,
    m_flow_nominal=0.1,
    nVolumes=10,
    redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
    annotation (extent=[-9,-40; 11,-20]);
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX MassFlowSource1(
                                                                  T_ambient=
        from_degC(30), redeclare package Medium = Modelica_Media.Air.DryAirNasa)
    annotation (extent=[-50,-40; -30,-20]);
  Modelica_Fluid.Sources.FixedAmbient_pTX ambient1(
                                                  T_ambient=from_degC(15),
      redeclare package Medium = Modelica_Media.Air.DryAirNasa) 
    annotation (extent=[50,-40; 30,-20]);
equation 
  connect(MassFlowSource.port, isolatedPipe.port_a) 
    annotation (points=[-29,10; -10,10],  style(color=69));
  connect(isolatedPipe.port_b, ambient.port) 
    annotation (points=[12,10; 29,10],   style(color=69));
  connect(ramp.y, MassFlowSource.m_flow_ambient) 
    annotation (points=[-69,10; -52,10], style(color=3, rgbcolor={0,0,255}));
  connect(MassFlowSource1.port, isolatedPipe1.port_a) 
    annotation (points=[-29,-30; -10,-30],style(color=69));
  connect(isolatedPipe1.port_b, ambient1.port) 
    annotation (points=[12,-30; 29,-30], style(color=69));
  connect(ramp.y, MassFlowSource1.m_flow_ambient) annotation (points=[-69,10; 
        -62,10; -62,-30; -52,-30], style(color=3, rgbcolor={0,0,255}));
end TestIsolatedPipe;
