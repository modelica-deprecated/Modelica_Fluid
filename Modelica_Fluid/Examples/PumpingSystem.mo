model PumpingSystem "Model of a pumping system for drinking water" 
  extends Modelica.Icons.Example;
  Modelica_Fluid.Sources.FixedAmbient source(
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    use_T=true,
    T=Modelica.SIunits.Conversions.from_degC(20)) 
    annotation (extent=[-100,-90; -80,-70]);
  
  Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe(
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    flowDirection=Modelica_Fluid.Types.FlowDirection.Bidirectional,
    height_ab=50,
    diameter=1,
    p_start=ambient.default_p_ambient,
    length=100,
    redeclare package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
    use_nominal=true) 
    annotation (extent=[-40,-40; -22,-20]);
  
  Modelica_Fluid.Pumps.Pump pumps(
    checkValve=true,
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    N_nom=1200,
    redeclare function flowCharacteristic = 
        Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticFlow (          q_nom={0,
            0.25,0.5}, head_nom={100,60,0}),
    Np_nom=4,
    M=50,
    T_start=Modelica.SIunits.Conversions.from_degC(20)) 
    annotation (extent=[-74,-88; -48,-62]);
  
  Modelica_Fluid.Volumes.OpenTank reservoir(
    initType=Modelica_Fluid.Types.Init.InitialValues,
    redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater,
    area=500,
    level_start=1.8,
    T_start=Modelica.SIunits.Conversions.from_degC(20),
    height=10,
    p_static_at_port=false,
    pipe_diameters={1}) 
    annotation (extent=[-14,-16; 6,4]);
  
  Modelica_Fluid.ControlValves.ValveLinear userValve(
                                                 redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater, Kv=200/2e5,
    flowDirection= Modelica_Fluid.Types.FlowDirection.Unidirectional) 
    annotation (extent=[58,-38; 74,-22]);
  Modelica_Fluid.Sources.FixedAmbient sink(
                                       redeclare package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater) 
    annotation (extent=[100,-40; 80,-20]);
  Modelica.Blocks.Sources.Step valveOpening(          offset=1) 
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
  Modelica_Fluid.Sensors.RelativePressure reservoirPressure(
                                                        redeclare package 
      Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater) 
    annotation (extent=[10,-12; 30,-32]);
  Modelica.Blocks.Continuous.FirstOrder PT1(T=50) 
    annotation (extent=[40,60; 60,80]);
  
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
  inner Modelica_Fluid.Ambient ambient 
                                   annotation (extent=[60,-96; 80,-76]);
equation 
  connect(userValve.port_b, sink.port)     annotation (points=[74,-30; 80,
        -30],
      style(color=69, rgbcolor={0,127,255}));
  connect(source.port, pumps.inlet) annotation (points=[-80,-80; -73.2,-80;
        -73.2,-77.6; -71.4,-77.6], style(color=69, rgbcolor={0,127,255}));
  connect(valveOpening.y, userValve.opening) annotation (points=[85,10; 98,
        10; 98,-12; 66,-12; 66,-22.8],
                                   style(color=74, rgbcolor={0,0,127}));
  connect(PressureSetPoint.y, controller.reference) annotation (points=[-79,70;
        -60,70; -60,76; -42,76], style(color=74, rgbcolor={0,0,127}));
  connect(controller.y, PumpRPMGenerator.u) 
    annotation (points=[-19,70; -2,70], style(color=5, rgbcolor={255,0,255}));
  connect(reservoirPressure.p_rel, controller.u) annotation (points=[20,-13; 20,
        50; -52,50; -52,64; -42,64], style(color=74, rgbcolor={0,0,127}));
  connect(reservoirPressure.port_b, sink.port)    annotation (points=[31,-22;
        44,-22; 44,-48; 80,-48; 80,-30],
                                      style(
      color=69,
      rgbcolor={0,127,255},
      pattern=3));
  connect(PumpRPMGenerator.y, PT1.u) 
    annotation (points=[21,70; 38,70], style(color=74, rgbcolor={0,0,127}));
  connect(PT1.y, pumps.N_in) annotation (points=[61,70; 74,70; 74,30; -64.38,30;
        -64.38,-69.28], style(color=74, rgbcolor={0,0,127}));
  connect(pipe.port_a, pumps.outlet)         annotation (points=[-40,-30; -53.2,
        -30; -53.2,-70.84], style(color=69, rgbcolor={0,127,255}));
  connect(pipe.port_b, reservoir.port[1]) annotation (points=[-22,-30; -4.2,-30;
        -4.2,-15.9], style(color=69, rgbcolor={0,127,255}));
  connect(userValve.port_a, reservoir.port[1]) annotation (points=[58,-30; -4.2,
        -30; -4.2,-15.9], style(color=69, rgbcolor={0,127,255}));
  connect(reservoirPressure.port_a, reservoir.port[1]) annotation (points=[9,-22;
        -2,-22; -2,-15.9; -4.2,-15.9],      style(
      color=69,
      rgbcolor={0,127,255},
      pattern=3));
end PumpingSystem;
