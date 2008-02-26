within FluidSandbox.Examples.StandardConnectionSemantics;
model TurbineTest "Series connection of two transport components" 
  extends Icons.Example;
  
  annotation (Diagram, experiment(StopTime=7200));
  import Modelica.SIunits.Conversions.*;
  Turbomachinery.TurbineStage turbine(
    redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase,
    G=0,
    H=0,
    redeclare package FluidInterface = FluidInterface, 
    provide_p_a=false, 
    provide_p_b=false, 
    provide_T_a=false, 
    provide_m_flow_ab=false, 
    provide_T_b=false)                    annotation (extent=[10,-10; 30,
      10]);
  Sources.PrescribedBoundary_pTX_A source(
    p=from_bar(100),
    T=from_degC(500),
    redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase,
    redeclare package FluidInterface = FluidInterface) 
  annotation (extent=[-90,-10; -70,10]);
  Sources.PrescribedBoundary_pTX_A sink(
    p=from_bar(30),
    T=from_degC(250),
    redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase,
    redeclare package FluidInterface = FluidInterface) 
  annotation (extent=[70,-40; 50,-20]);
  annotation (Diagram, experiment(StopTime=7200));
  Modelica.Blocks.Sources.TimeTable valveTable(table=[0,0; 1800,0; 3600,1; 7210,
        1], offset=0) 
            annotation (extent=[-80,50; -60,70]);
  Valves.ValveLinearAB valve(
    Kv=1e-4,
    redeclare package Medium = Modelica.Media.Water.StandardWaterOnePhase,
    redeclare package FluidInterface = FluidInterface, 
    provide_p_a=false, 
    provide_p_b=false, 
    provide_T_a=false, 
    provide_T_b=false, 
    provide_m_flow_ab=false) 
  annotation (extent=[-30,-10; -50,10]);
  Modelica.Mechanics.Rotational.ConstantSpeed load(w_fixed=-50) 
  annotation (extent=[77,-4.5; 62,10.5]);
equation 
  connect(valveTable.y, valve.opening)     annotation (points=[-59,60;
      -40,60; -40,9],  style(
    color=74,
    rgbcolor={0,0,127},
    fillColor=46,
    rgbfillColor={216,62,1},
    fillPattern=1));
  connect(turbine.flange, load.flange) 
  annotation (points=[30,3; 62,3], style(color=0, rgbcolor={0,0,0}));
  connect(valve.port_b, source.port) annotation (points=[-50,0; -70,0], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=71,
      rgbfillColor={85,170,255},
      fillPattern=1));
  connect(valve.port_a, turbine.port_a) annotation (points=[-30,0; 10,0], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=71,
      rgbfillColor={85,170,255},
      fillPattern=1));
  connect(turbine.port_b, sink.port) annotation (points=[30,0; 40,0; 40,-30; 50,
        -30], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=71,
      rgbfillColor={85,170,255},
      fillPattern=1));
end TurbineTest;
