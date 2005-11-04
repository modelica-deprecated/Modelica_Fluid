model PumpingSystem "Model of a pumping system for drinking water" 
  extends Modelica.Icons.Example;
  Sources.FixedAmbient source(
    redeclare package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater,
    use_T=true,
    T=Modelica.SIunits.Conversions.from_degC(20)) 
    annotation (extent=[-100,-90; -80,-70]);
  
  Components.PressureDropPipe pipeFriction(
    redeclare package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater,
    frictionType=Modelica_Fluid.Types.FrictionTypes.ConstantTurbulent,
    m_flow_nominal=1000,
    dp_nominal=Modelica.SIunits.Conversions.from_bar(0.3),
    allowFlowReversal=false) 
    annotation (extent=[-32,-40; -16,-20]);
  
  Components.StaticHead pipeHead(H_b_a=50, redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    allowFlowReversal=false) 
    annotation (extent=[-48,-64; -30,-38], rotation=90);
  Components.Pump pumps(
    checkValve=true,
    redeclare package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater,
    N_nom=1200,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Types.PumpCharacteristics.quadraticFlow (q_nom={0,
            0.25,0.5}, head_nom={100,60,0}),
    Np_nom=4,
    M=50,
    T_start=Modelica.SIunits.Conversions.from_degC(20)) 
    annotation (extent=[-74,-88; -48,-62]);
  
  Components.Tank reservoir(
    H0=18,
    initOption=Modelica_Fluid.Types.InitTypes.InitialValues,
    redeclare package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater,
    area=500,
    level_start=1.8,
    T_start=Modelica.SIunits.Conversions.from_degC(20)) 
    annotation (extent=[-14,-16; 6,4]);
  
  Components.ValveLinear userValve(redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater, Kv=200/2e5,
    allowFlowReversal=false) 
    annotation (extent=[58,-38; 74,-22]);
  Sources.FixedAmbient ambient1(redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater) 
    annotation (extent=[100,-40; 80,-20]);
  Modelica.Blocks.Sources.Step valveOpening(height=0, offset=1) 
    annotation (extent=[64,0; 84,20]);
  Modelica.Blocks.Sources.Constant PressureSetPoint(k=2e5) 
    annotation (extent=[-100,60; -80,80]);
  Modelica.Blocks.Logical.OnOffController controller(bandwidth=1000,
      pre_y_start=true) annotation (extent=[-40,60; -20,80]);
  Modelica.Blocks.Logical.TriggeredTrapezoid PumpRPMGenerator(
    rising=3,
    falling=3,
    amplitude=1200,
    offset=0.001) annotation (extent=[0,60; 20,80]);
  Sensors.RelativePressure reservoirPressure(redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater) 
    annotation (extent=[8,-12; 28,-32]);
  Sources.FixedAmbient ambient2(redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater, use_T=false) 
    annotation (extent=[46,-28; 34,-16]);
equation 
  
  annotation (
    Diagram,
    Documentation(info="<html>
Water is pumped from a source by 4 pumps in parallel (fitted with check valves), through a pipe whose outlet is 50 m higher than the source, into a reservoir placed on a 18-m high tower. The users are represented by an equivalent valve, connected to the reservoir.
<p>
The water controller is a simple on-off controller, acting on the gauge pressure measured at the base of the tower; the output of the controller is the rotational speed of the pumps. A small but nonzero rotational speed is used to represent the standby state of the pumps, in order to avoid singularities in the flow characteristic.
<p>
Simulate for 2000 s. The pump turns on and off to keep the reservoir level around 2.5 meters, which means 20.5 meters higher than the base of the tower, corresponding to a gauge pressure of 2 bar
<p>
If using Dymola, turn off \"Equidistant time grid\" to avoid numerical errors.
</html>", revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created.</li>
</ul>
</html>"),
    experiment(StopTime=2000, NumberOfIntervals=5000),
    experimentSetupOutput(equdistant=false));
  connect(pumps.outlet, pipeHead.port_a) annotation (points=[-53.2,-70.84; -46,
        -72; -38,-72; -38,-65.3; -39,-65.3], style(color=69, rgbcolor={0,127,
          255}));
  connect(pipeFriction.port_b, reservoir.port) annotation (points=[-15.2,-30;
        -4,-30; -4,-17], style(color=69, rgbcolor={0,127,255}));
  connect(reservoir.port, userValve.port_a) annotation (points=[-4,-17; -4,-30;
        57.2,-30], style(color=69, rgbcolor={0,127,255}));
  connect(userValve.port_b, ambient1.port) annotation (points=[74.8,-30; 79,-30],
      style(color=69, rgbcolor={0,127,255}));
  connect(source.port, pumps.inlet) annotation (points=[-79,-80; -73.2,-80;
        -73.2,-77.6; -71.4,-77.6], style(color=69, rgbcolor={0,127,255}));
  connect(valveOpening.y, userValve.opening) annotation (points=[85,10; 98,10;
        98,-12; 66,-12; 66,-23.6], style(color=74, rgbcolor={0,0,127}));
  connect(pipeHead.port_b, pipeFriction.port_a) annotation (points=[-39,-36.7;
        -39,-30; -32.8,-30], style(color=69, rgbcolor={0,127,255}));
  connect(PressureSetPoint.y, controller.reference) annotation (points=[-79,70;
        -60,70; -60,76; -42,76], style(color=74, rgbcolor={0,0,127}));
  connect(controller.y, PumpRPMGenerator.u) 
    annotation (points=[-19,70; -2,70], style(color=5, rgbcolor={255,0,255}));
  connect(PumpRPMGenerator.y, pumps.N_in) annotation (points=[21,70; 46,70; 46,
        36; -64.38,36; -64.38,-69.28], style(color=74, rgbcolor={0,0,127}));
  connect(reservoir.port, reservoirPressure.port_a) annotation (points=[-4,-17;
        3,-17; 3,-22; 7,-22], style(color=69, rgbcolor={0,127,255}));
  connect(reservoirPressure.port_b, ambient2.port) annotation (points=[29,-22;
        33.4,-22], style(color=69, rgbcolor={0,127,255}));
  connect(reservoirPressure.p_rel, controller.u) annotation (points=[18,-13; 18,
        50; -52,50; -52,64; -42,64], style(color=74, rgbcolor={0,0,127}));
end PumpingSystem;
