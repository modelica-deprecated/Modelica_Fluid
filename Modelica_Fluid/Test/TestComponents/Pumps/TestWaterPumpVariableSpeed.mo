within Modelica_Fluid.Test.TestComponents.Pumps;
model TestWaterPumpVariableSpeed "Test case for WaterPump" 
  extends Modelica.Icons.Example;
  import PC = Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics;
  replaceable function pumpFlowChar = PC.quadraticFlow(q_nom={0,0.001,0.0015}, head_nom={100,50,0});
annotation (
  Diagram,
  experiment(StopTime=10, Tolerance=1e-006),
  Documentation(info=""));
  Modelica_Fluid.Sources.FixedBoundary_pTX Source(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater, p=1e5,
    T=ambient.default_T_ambient) 
  annotation (extent=[-100,20; -80,40]);
  Modelica_Fluid.Sources.PrescribedBoundary_pTX Sink(
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=5e5,
    T=ambient.default_T_ambient) 
  annotation (extent=[34,26; 14,46]);
  Modelica_Fluid.Pumps.Pump Pump1(
    pin_start=1e5,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    redeclare function flowCharacteristic = pumpFlowChar,
    m_flow_start=1,
    pout_start=7e5,
    use_N_input=true)      annotation (extent=[-66,20; -34,50]);
  Modelica.Blocks.Sources.Constant Constant1 
  annotation (extent=[-40,80; -20,100]);
  Modelica_Fluid.ControlValves.ValveLinear Valve(
                                             redeclare package Medium = 
        Modelica.Media.Water.StandardWater, Kv=1e-5) 
  annotation (extent=[-16,26; 2,46]);
  
  inner Modelica_Fluid.Ambient ambient 
                                   annotation (extent=[64,-4; 84,16]);
  Modelica.Blocks.Sources.Ramp N_pump(
    startTime=1,
    duration=5,
    height=-500,
    offset=1500) 
                annotation (extent=[-100,62; -80,82]);
equation 
  connect(Constant1.y, Valve.opening)     annotation (points=[-19,90; -7,90; -7,
        45],        style(color=74, rgbcolor={0,0,127}));
  connect(Valve.port_b, Sink.port)       annotation (points=[2,36; 14,36],
                       style(color=69, rgbcolor={0,127,255}));
  connect(Valve.port_a, Pump1.outlet)     annotation (points=[-16,36; -26,36;
        -26,39.8; -40.4,39.8],   style(color=69, rgbcolor={0,127,255}));
  connect(Pump1.inlet, Source.port) annotation (points=[-62.8,32; -70,32; -70,
        30; -80,30],           style(color=69, rgbcolor={0,127,255}));
  connect(N_pump.y, Pump1.N_in) annotation (points=[-79,72; -54.16,72; -54.16,
        41.6], style(color=74, rgbcolor={0,0,127}));
end TestWaterPumpVariableSpeed;
