model TestShortPipeWithVolume "Test ShortPipe with PortVolume" 
  import Modelica_Fluid.Types.Init;
  import SI = Modelica.SIunits;
  extends Modelica.Icons.Example;
  replaceable package Medium = Modelica.Media.Air.SimpleAir 
                       extends Modelica.Media.Interfaces.PartialMedium 
    "Medium model" annotation (choicesAllMatching=true);
  parameter SI.AbsolutePressure p_start = 1.0e5 "Initial value of pressure";
  parameter SI.Temperature T_start = 300 "Initial value of temperature";
  parameter Real X_start[Medium.nX] = Medium.X_default 
    "Initial value of mass fractions";
  Sources.PrescribedMassFlowRate_TX pump(
    m_flow=1,
    T=1.2*T_start,
    X=X_start,
    redeclare package Medium = Medium) annotation (extent=[-80,0; -60,20]);
  Utilities.PortVolume volume(
    redeclare package Medium = Medium,
    V=0.1,
    initOption=Init.InitialValues,
    p_start=p_start,
    T_start=T_start,
    X_start=X_start) annotation (extent=[-40,0; -20,20]);
  Components.PressureDropPipe pipe(
    redeclare package Medium = Medium,
    frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantLaminar,
    dp_nominal=0.1e5,
    m_flow_nominal=1) annotation (extent=[0,0; 20,20]);
  Sources.FixedAmbient_pTX ambient(
    redeclare package Medium = Medium,
    p=p_start,
    T=T_start,
    X=X_start) annotation (extent=[60,0; 40,20]);
  inner Components.FluidOptions fluidOptions 
    annotation (extent=[-100,-100; -80,-80]);
equation 
  connect(pump.port, volume.port) annotation (points=[-60,10; -30,10], style(
        color=69, rgbcolor={0,127,255}));
  connect(volume.port, pipe.port_a) annotation (points=[-30,10; 0,10],  style(
        color=69, rgbcolor={0,127,255}));
  connect(pipe.port_b, ambient.port) annotation (points=[20,10; 40,10], style(
        color=69, rgbcolor={0,127,255}));
  annotation (Diagram);
end TestShortPipeWithVolume;
