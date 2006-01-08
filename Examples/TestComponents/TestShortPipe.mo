model TestShortPipe "Test ShortPipe component" 
  import Modelica.SIunits.Conversions.*;
  
  extends Modelica.Icons.Example;
  annotation (
    Diagram,
    experiment(StopTime=3),
    Coordsys(grid=[1, 1], component=[20, 20]));
  Modelica_Fluid.Components.PressureDropPipe shortPipe(
    dp_nominal=from_bar(0.1),
    roughness=2e-5,
    from_dp=true,
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    frictionType=Modelica_Fluid.Types.FrictionTypes.DetailedFriction,
    diameter=0.02) 
    annotation (extent=[-10,0; 10,20]);
  
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX m_flow_source(T=from_degC(30),
      redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater) 
    annotation (extent=[-50,0; -30,20]);
  Modelica_Fluid.Sources.FixedAmbient_pTX ambient(T=from_degC(15),
      redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater) 
    annotation (extent=[50,0; 30,20]);
  Modelica.Blocks.Sources.Ramp ramp(
    duration=3,
    height=6,
    offset=-3)   annotation (extent=[-90,0; -70,20]);
  inner Components.FluidOptions fluidOptions 
    annotation (extent=[-100,-100; -80,-80]);
equation 
  connect(m_flow_source.port, shortPipe.port_a) 
    annotation (points=[-30,10; -10,10],  style(color=69));
  connect(shortPipe.port_b, ambient.port) 
    annotation (points=[10,10; 30,10],   style(color=69));
  connect(ramp.y, m_flow_source.m_flow_in) 
    annotation (points=[-69,10; -59.9,10; -59.9,16; -49.3,16],
                                         style(color=3, rgbcolor={0,0,255}));
end TestShortPipe;
