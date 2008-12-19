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
    use_d_nominal=true) 
              annotation (Placement(transformation(extent={{-80,20},{-60,40}},
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
    p_b_nominal=200000) 
    annotation (Placement(transformation(extent={{-60,-10},{-40,10}},rotation=
           0)));
  Modelica_Fluid.Valves.ValveIncompressible valve(
    redeclare package Medium = Medium,
    CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
    m_flow_nominal=0.01,
    compute_T=true,
    dp_nominal=50000) 
    annotation (Placement(transformation(extent={{42,-8},{58,8}},  rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput circuitFlowRate 
    annotation (Placement(transformation(extent={{88,12},{108,32}}, rotation=
            0)));
  Sensors.MassFlowRate massFlowRate(redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{-30,10},{-10,-10}},rotation=
           0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature ambientTemperature(
                                                                    T=system.T_ambient) 
    annotation (Placement(transformation(extent={{-24,-38},{-10,-24}},
                                                                     rotation=
           0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1(
                                                                   G=1.6e3/20) 
    annotation (Placement(transformation(
        origin={10,-48},
        extent={{8,-10},{-8,10}},
        rotation=90)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow burner(
                                                     Q_flow=1.6e3,
    T_ref=343.15,
    alpha=-0.5) 
    annotation (Placement(transformation(extent={{-4,12},{16,32}}, rotation=0)));
  inner Modelica_Fluid.System system(energyDynamics=Modelica_Fluid.Types.Dynamics.SteadyStateInitial) 
                        annotation (Placement(transformation(extent={{-90,70},{
            -70,90}},   rotation=0)));
  Pipes.DistributedPipe pipe(
    redeclare package Medium = Medium,
    use_T_start=true,
    T_start=Modelica.SIunits.Conversions.from_degC(80),
    length=2,
    redeclare model HeatTransfer = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.IdealHeatTransfer,
    diameter=0.01,
    nNodes=1,
    modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb,
    p_a_start=400000,
    p_b_start=390000,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.WallFrictionPressureLoss)
    annotation (Placement(transformation(extent={{10,-10},{30,10}},rotation=0)));

  Pipes.DistributedPipe radiator(
    use_T_start=true,
    redeclare package Medium = Medium,
    length=10,
    T_start=Modelica.SIunits.Conversions.from_degC(40),
    redeclare model HeatTransfer = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.IdealHeatTransfer,
    diameter=0.01,
    modelStructure=Modelica_Fluid.Types.ModelStructure.av_b,
    p_a_start=110000,
    p_b_start=105000,
    nNodes=1,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.WallFrictionPressureLoss)
    annotation (Placement(transformation(extent={{20,-80},{0,-60}}, rotation=
            0)));

  Modelica.Blocks.Interfaces.RealOutput hotWaterTemperature 
    annotation (Placement(transformation(extent={{88,-28},{108,-8}}, rotation=
           0)));
  Modelica.Blocks.Interfaces.RealOutput coldWaterTemperature 
    annotation (Placement(transformation(extent={{88,-78},{108,-58}},
          rotation=0)));
  Modelica_Fluid.Sensors.Temperature sensor_T_1(
                                 redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{60,-56},{40,-36}}, rotation=
           0)));
  Modelica_Fluid.Sensors.Temperature sensor_T_2(
                                 redeclare package Medium = Medium) 
    annotation (Placement(transformation(extent={{-30,-56},{-50,-36}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput tankLevel 
                                 annotation (Placement(transformation(extent=
            {{90,60},{110,80}}, rotation=0)));
  Modelica.Blocks.Sources.Step valveOpening(
    startTime=2000,
    height=0.9,
    offset=0.1)   annotation (Placement(transformation(extent={{20,60},{40,80}},
                  rotation=0)));
equation
tankLevel = tank.level;
  connect(massFlowRate.m_flow, circuitFlowRate) annotation (Line(points={{-20,11},
          {-20,40},{80,40},{80,22},{98,22}},    color={0,0,127}));
  connect(massFlowRate.port_b, pipe.port_a) annotation (Line(points={{-10,0},{
          -10,0},{10,0}},
                    color={0,127,255}));
  connect(pipe.port_b, valve.port_a) annotation (Line(points={{30,0},{42,0}},
        color={0,127,255}));
  connect(ambientTemperature.port, thermalConductor1.port_a) annotation (Line(
        points={{-10,-31},{10,-31},{10,-40}},
                                            color={191,0,0}));
  connect(sensor_T_1.T, hotWaterTemperature) annotation (Line(points={{43,-46},
          {43,-80},{76,-80},{76,-18},{98,-18}}, color={0,0,127}));
  connect(radiator.port_a, valve.port_b) annotation (Line(points={{20,-70},{70,
          -70},{70,0},{58,0}},      color={0,127,255}));
  connect(sensor_T_2.port, radiator.port_b) annotation (Line(points={{-40,-56},
          {-40,-70},{0,-70}}, color={0,127,255}));
  connect(radiator.port_a, sensor_T_1.port) annotation (Line(points={{20,-70},{
          50,-70},{50,-56}},  color={0,127,255}));
  connect(radiator.port_b, tank.ports[1]) annotation (Line(
      points={{0,-70},{-70,-70},{-70,22}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(tank.ports[2], pump.port_a) annotation (Line(
      points={{-70,18},{-70,0},{-60,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(valveOpening.y, valve.opening) annotation (Line(
      points={{41,70},{50,70},{50,6.4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pump.port_b, massFlowRate.port_a) annotation (Line(
      points={{-40,0},{-30,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(sensor_T_2.T, coldWaterTemperature) annotation (Line(
      points={{-47,-46},{-50,-46},{-50,-90},{80,-90},{80,-68},{98,-68}},
      color={0,0,127},
      smooth=Smooth.None));

  connect(burner.port, pipe.heatPorts[1]) annotation (Line(
      points={{16,22},{20.1,22},{20.1,5.2}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(thermalConductor1.port_b, radiator.heatPorts[1]) annotation (Line(
      points={{10,-56},{10,-64.8},{9.9,-64.8}},
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
After 2000s of simulation time the valve fully opens. A simple idealized control is been embedded 
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
Moreover tank.use_d_nominal is configured to true, in order to decouple the mass balance 
from the energy balance when using the <tt>Modelica.Media.Water.StandardWater</tt> medium and performing a steady-state simulation. 
Alternatively a simpler medium model, such as <tt>Modelica.Media.Water.ConstantPropertyLiquidWater</tt>, 
could be used.
</p>
Furthermore note that a steady-state temperature of the tank is only well defined if fluid flows through the cycle. This is because the
tank is assumed to be perfectly isolated. If a steady-state simultion shall be started with the valve fully closed, then a thermal 
coupling between the tank and its environment should be established.
<p>
</p>
</html>
"), experiment(StopTime=6000),
    Commands(file(ensureSimulated=true)=
        "Scripts/Examples/HeatingSystem/plotResults.mos" "plotResults"));
end HeatingSystem;
