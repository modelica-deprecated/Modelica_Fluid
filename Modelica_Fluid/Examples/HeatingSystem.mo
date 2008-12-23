within Modelica_Fluid.Examples;
model HeatingSystem "Simple model of a heating system"
   replaceable package Medium = 
      Modelica.Media.Water.StandardWater 
     constrainedby Modelica.Media.Interfaces.PartialMedium;

  Modelica_Fluid.Vessels.OpenTank tank(
    redeclare package Medium = Medium,
    crossArea=0.01,
    V0=0.01,
    height=2,
    level_start=1,
    nPorts=2,
    massDynamics=Modelica_Fluid.Types.Dynamics.FixedInitial,
    use_HeatTransfer=true,
    redeclare model HeatTransfer = 
        Modelica_Fluid.Vessels.BaseClasses.HeatTransfer.ConstantHeatTransfer (
          alpha0=10)) 
              annotation (Placement(transformation(extent={{-80,30},{-60,50}},
          rotation=0)));
  Machines.ControlledPump pump(
    redeclare package Medium = Medium,
    N_nominal=1500,
    use_T_start=true,
    T_start=Modelica.SIunits.Conversions.from_degC(40),
    m_flow_start=0.01,
    checkValve=false,
    m_flow_nominal=0.01,
    control_m_flow=false,
    p_a_start=110000,
    p_b_start=400000,
    p_a_nominal=100000,
    p_b_nominal=130000) 
    annotation (Placement(transformation(extent={{-50,10},{-30,30}}, rotation=
           0)));
  Modelica_Fluid.Valves.ValveIncompressible valve(
    redeclare package Medium = Medium,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    m_flow_nominal=0.01,
    show_T=true,
    dp_nominal=10000) 
    annotation (Placement(transformation(extent={{70,-80},{50,-60}},
                                                                   rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput flowRate 
    annotation (Placement(transformation(extent={{-6,34},{6,46}},   rotation=
            0)));
  Sensors.MassFlowRate massFlowRate(redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{-20,10},{0,30}},   rotation=
           0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature ambientTemperature(
                                                                    T=system.T_ambient) 
    annotation (Placement(transformation(extent={{-15,-27},{-1,-13}},rotation=
           0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor wall(G=1.6e3/20) 
    annotation (Placement(transformation(
        origin={10,-48},
        extent={{8,-10},{-8,10}},
        rotation=90)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow burner(
                                                     Q_flow=1.6e3,
    T_ref=343.15,
    alpha=-0.5) 
    annotation (Placement(transformation(extent={{16,30},{36,50}}, rotation=0)));
  inner Modelica_Fluid.System system(energyDynamics=Modelica_Fluid.Types.Dynamics.SteadyStateInitial) 
                        annotation (Placement(transformation(extent={{-90,70},{
            -70,90}},   rotation=0)));
  Pipes.DistributedPipe heater(
    redeclare package Medium = Medium,
    use_T_start=true,
    T_start=Modelica.SIunits.Conversions.from_degC(80),
    length=2,
    redeclare model HeatTransfer = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.IdealFlowHeatTransfer,
    diameter=0.01,
    nNodes=1,
    modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedWallFriction,
    use_HeatTransfer=true) 
    annotation (Placement(transformation(extent={{30,10},{50,30}}, rotation=0)));

  Pipes.DistributedPipe radiator(
    use_T_start=true,
    redeclare package Medium = Medium,
    length=10,
    T_start=Modelica.SIunits.Conversions.from_degC(40),
    redeclare model HeatTransfer = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.IdealFlowHeatTransfer,
    diameter=0.01,
    nNodes=1,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedWallFriction,
    modelStructure=Modelica_Fluid.Types.ModelStructure.av_b,
    use_HeatTransfer=true) 
    annotation (Placement(transformation(extent={{20,-80},{0,-60}}, rotation=
            0)));

  Modelica.Blocks.Interfaces.RealOutput T_forward 
    annotation (Placement(transformation(extent={{74,34},{86,46}},   rotation=
           0)));
  Modelica.Blocks.Interfaces.RealOutput T_return 
    annotation (Placement(transformation(extent={{-46,-56},{-58,-44}},
          rotation=0)));
  Modelica_Fluid.Sensors.Temperature sensor_T_forward(redeclare package Medium
      = Medium) 
    annotation (Placement(transformation(extent={{50,30},{70,50}},   rotation=
           0)));
  Modelica_Fluid.Sensors.Temperature sensor_T_return(redeclare package Medium
      = Medium) 
    annotation (Placement(transformation(extent={{-20,-60},{-40,-40}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput tankLevel 
                                 annotation (Placement(transformation(extent={{-56,34},
            {-44,46}},          rotation=0)));
  Modelica.Blocks.Sources.Step valveOpening(
    startTime=2000,
    height=0.9,
    offset=0.1)   annotation (Placement(transformation(extent={{36,-27},{50,-13}},
                  rotation=0)));
  Pipes.DistributedPipe pipe(
    redeclare package Medium = Medium,
    use_T_start=true,
    T_start=Modelica.SIunits.Conversions.from_degC(80),
    redeclare model HeatTransfer = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.IdealFlowHeatTransfer,
    diameter=0.01,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedWallFriction,
    length=10) 
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                                                                   rotation=-90,
        origin={80,-20})));

  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature ambientTemperature1(
                                                                    T=system.T_ambient) 
    annotation (Placement(transformation(extent={{-95,35},{-85,45}}, rotation=
           0)));
equation
tankLevel = tank.level;
  connect(massFlowRate.m_flow, flowRate)        annotation (Line(points={{-10,31},
          {-10,40},{0,40}},                     color={0,0,127}));
  connect(massFlowRate.port_b, heater.port_a) 
                                            annotation (Line(points={{0,20},{0,
          20},{30,20}},
                    color={0,127,255}));
  connect(ambientTemperature.port, wall.port_a)              annotation (Line(
        points={{-1,-20},{10,-20},{10,-40}},color={191,0,0}));
  connect(sensor_T_forward.T, T_forward)     annotation (Line(points={{67,40},{
          80,40}},                              color={0,0,127}));
  connect(radiator.port_a, valve.port_b) annotation (Line(points={{20,-70},{20,
          -70},{50,-70}},           color={0,127,255}));
  connect(sensor_T_return.port, radiator.port_b) 
                                            annotation (Line(points={{-30,-60},
          {-30,-70},{0,-70}}, color={0,127,255}));
  connect(tank.ports[2], pump.port_a) annotation (Line(
      points={{-70,28},{-70,20},{-50,20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(valveOpening.y, valve.opening) annotation (Line(
      points={{50.7,-20},{60,-20},{60,-62}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pump.port_b, massFlowRate.port_a) annotation (Line(
      points={{-30,20},{-20,20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(sensor_T_return.T, T_return)        annotation (Line(
      points={{-37,-50},{-52,-50}},
      color={0,0,127},
      smooth=Smooth.None));

  connect(burner.port, heater.heatPorts[1]) 
                                          annotation (Line(
      points={{36,40},{40.1,40},{40.1,25.2}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(wall.port_b, radiator.heatPorts[1])              annotation (Line(
      points={{10,-56},{10,-64.8},{9.9,-64.8}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(sensor_T_forward.port, heater.port_b) 
                                              annotation (Line(
      points={{60,30},{60,20},{50,20}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(heater.port_b, pipe.port_a) annotation (Line(
      points={{50,20},{80,20},{80,-10}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe.port_b, valve.port_a) annotation (Line(
      points={{80,-30},{80,-70},{70,-70}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(radiator.port_b, tank.ports[1]) annotation (Line(
      points={{0,-70},{-70,-70},{-70,32}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(ambientTemperature1.port, tank.heatPort) annotation (Line(
      points={{-85,40},{-80,40}},
      color={191,0,0},
      smooth=Smooth.None));

  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics),
                       Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={Rectangle(extent={{-100,100},{100,-100}}, 
            lineColor={0,0,255}), Text(
          extent={{-60,60},{60,-60}},
          lineColor={0,0,255},
          textString="H")}), Documentation(info="<html>
<p>
Simple heating system with a closed flow cycle. It is set up for steady-state initial values.
After 2000s of simulation time the valve fully opens. A simple idealized control is embedded 
into the respective components, so that the heating system can be regulated with the valve:
the pump controls the pressure, the burner controls the temperature. 
</p>
<p>
One can investigate the temperatures and flows for different settings of <tt>system.energyDynamics</tt> 
(see Assumptions tab of the system object).
With <tt>system.energyDynamics==Types.Dynamics.SteadyState</tt> all but one dynamic states are eliminated.
The left state <tt>tank.m</tt> is to account for the closed flow cycle. It is constant as outflow and inflow are equal 
in a steady-state simulation.
</p>
<p>
Note that a closed flow cycle generally causes circular equalities for the mass flow rates and leaves the pressure undefined.
This is why the tank.massDynamics, i.e. the tank level determining the port pressure, is modified locally to Types.Dynamics.FixedInitial.
</p>
<p>
Also note that the tank is thermally isolated againts its ambient. This way the temperature of the tank is also
well defined for zero flow rate in the heating system. The pipe however is assumed to be perfectly isolated. 
If a steady-state simultion shall be started with the valve fully closed, then a thermal 
coupling between the pipe and its ambient should be established.
</p>
<p>
Moreover it is worth noting that the idialized direct connection between the heater and the pipe, resulting in equal port pressures,
is treated as high-index DAE, as opposed to a nonlinear equation system for connected pressure loss correlations. A pressure loss correlation 
could be additionally introduced to model the fitting between the heater and the pipe, e.g. to adapt different diameters.
</p>
</html>
"), experiment(StopTime=6000),
    Commands(file(ensureSimulated=true)=
        "Scripts/Examples/HeatingSystem/plotResults.mos" "plotResults"));
end HeatingSystem;
