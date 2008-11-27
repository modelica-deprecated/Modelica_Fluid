within Modelica_Fluid.Test;
package TestOverdeterminedSteadyStateInit
  "Contains test cases to test overdetermined systems of initial equations"
  model HeatingSystem "Simple model of a heating system"
     replaceable package Medium = Modelica.Media.Water.StandardWater 
       constrainedby Modelica.Media.Interfaces.PartialMedium;
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                        graphics),
                         Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={Rectangle(extent={{-100,100},{100,-100}}, 
              lineColor={0,0,255}), Text(
            extent={{-60,60},{60,-60}},
            lineColor={0,0,255},
            textString="P")}));
    Volumes.OpenTank tank(
      redeclare package Medium = Medium,
      p_static_at_port=true,
      area=0.01,
      V0=0.01,
      pipe_diameters={0.025},
      initType=Modelica_Fluid.Types.Init.InitialValues,
      height=2,
      level_start=1) 
                annotation (Placement(transformation(extent={{-76,6},{-54,28}},
            rotation=0)));
    Pumps.Pump pump(
      redeclare package Medium = Medium,
      N_nominal=1500,
      N_const=1500,
      pin_start=1.1e5,
      pout_start=4.0e5,
      use_T_start=true,
      T_start=Modelica.SIunits.Conversions.from_degC(40),
      redeclare function flowCharacteristic = 
          Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.linearFlow (
            head_nominal={60.0,0}, q_nominal={0.0,0.02e-3}),
      m_flow_start=0.01,
      checkValve=false) 
      annotation (Placement(transformation(extent={{-58,-16},{-38,4}}, rotation=
             0)));
    ControlValves.ValveIncompressible valve(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      dp_nominal=3e5,
      m_flow_nominal=0.01) 
      annotation (Placement(transformation(extent={{42,-12},{58,4}}, rotation=0)));
    Modelica.Blocks.Interfaces.RealInput valvePosition 
      annotation (Placement(transformation(extent={{-128,-20},{-88,20}},
            rotation=0)));
    Modelica.Blocks.Interfaces.RealOutput circuitFlowRate 
      annotation (Placement(transformation(extent={{88,12},{108,32}}, rotation=
              0)));
    Sensors.MassFlowRate massFlowRate(redeclare package Medium = Medium) 
      annotation (Placement(transformation(extent={{-34,6},{-14,-14}}, rotation=
             0)));
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature ambientTemperature(
                                                                      T=system.T_ambient) 
      annotation (Placement(transformation(extent={{-12,-40},{2,-26}}, rotation=
             0)));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor1
      (                                                              G=1.6e3/20) 
      annotation (Placement(transformation(
          origin={18,-48},
          extent={{8,-10},{-8,10}},
          rotation=90)));
    Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow burner(
                                                       Q_flow=1.6e3) 
      annotation (Placement(transformation(extent={{-2,12},{18,32}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,80},
              {-80,100}}, rotation=0)));
    Modelica_Fluid.Pipes.LumpedPipe pipe(
      redeclare package Medium = Medium,
      use_T_start=true,
      diameter=0.03,
      T_start=Modelica.SIunits.Conversions.from_degC(80),
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      use_eta_nominal=true,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      length=2,
      p_a_start=400000,
      p_b_start=390000,
      redeclare model HeatTransfer = 
          Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_ideal) 
      annotation (Placement(transformation(extent={{12,-14},{32,6}}, rotation=0)));

    Modelica_Fluid.Pipes.LumpedPipe radiator(
      use_T_start=true,
      redeclare package Medium = Medium,
      length=10,
      diameter=0.05,
      T_start=Modelica.SIunits.Conversions.from_degC(40),
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      use_eta_nominal=true,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      p_a_start=110000,
      p_b_start=105000,
      redeclare model HeatTransfer = 
          Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_ideal) 
      annotation (Placement(transformation(extent={{28,-76},{8,-56}}, rotation=
              0)));

    Modelica.Blocks.Interfaces.RealOutput hotWaterTemperature 
      annotation (Placement(transformation(extent={{88,-28},{108,-8}}, rotation=
             0)));
    Modelica.Blocks.Interfaces.RealOutput coldWaterTemperature 
      annotation (Placement(transformation(extent={{88,-78},{108,-58}},
            rotation=0)));
    Modelica_Fluid.Sensors.TemperatureOnePort sensor_T_1(
                                   redeclare package Medium = Medium) 
      annotation (Placement(transformation(extent={{56,-56},{36,-36}}, rotation=
             0)));
    Modelica_Fluid.Sensors.TemperatureOnePort sensor_T_2(
                                   redeclare package Medium = Medium) 
      annotation (Placement(transformation(extent={{-16,-56},{-36,-36}},
            rotation=0)));
    Modelica.Blocks.Interfaces.RealOutput tankLevel 
                                   annotation (Placement(transformation(extent=
              {{90,60},{110,80}}, rotation=0)));
    Modelica.Blocks.Sources.RealExpression realExpression 
      annotation (Placement(transformation(extent={{-74,22},{-54,42}}, rotation=
             0)));
  equation
  tankLevel = tank.level;
    connect(valvePosition, valve.stemPosition) annotation (Line(points={{-108,0},
            {-86,0},{-86,66},{50,66},{50,2.4}}, color={0,0,127}));
    connect(pump.outlet, massFlowRate.port_a) annotation (Line(points={{-40,-6},{
            -40,-3.4},{-34,-3.4},{-34,-4}},        color={0,127,255}));
    connect(massFlowRate.m_flow, circuitFlowRate) annotation (Line(points={{-24,
            7},{-24,38},{36,38},{36,22},{98,22}}, color={0,0,127}));
    connect(massFlowRate.port_b, pipe.port_a) annotation (Line(points={{-14,-4},
            {12,-4}}, color={0,127,255}));
    connect(pipe.port_b, valve.port_a) annotation (Line(points={{32,-4},{42,-4}},
          color={0,127,255}));
    connect(thermalConductor1.port_b, radiator.heatPort) annotation (Line(
          points={{18,-56},{18,-60.6}}, color={191,0,0}));
    connect(burner.port, pipe.heatPort) annotation (Line(points={{18,22},{22,22},
            {22,1.4}},     color={191,0,0}));
    connect(ambientTemperature.port, thermalConductor1.port_a) annotation (Line(
          points={{2,-33},{18,-33},{18,-40}}, color={191,0,0}));
    connect(sensor_T_1.T, hotWaterTemperature) annotation (Line(points={{39,-46},
            {39,-80},{72,-80},{72,-18},{98,-18}}, color={0,0,127}));
    connect(sensor_T_2.T, coldWaterTemperature) annotation (Line(points={{-33,
            -46},{-33,-88},{76,-88},{76,-68},{98,-68}}, color={0,0,127}));
    connect(radiator.port_a, valve.port_b) annotation (Line(points={{28,-66},{
            68,-66},{68,-4},{58,-4}}, color={0,127,255}));
    connect(pump.inlet, radiator.port_b) annotation (Line(points={{-56,-6},{-56,
            -66},{8,-66}}, color={0,127,255}));
    connect(sensor_T_2.port, radiator.port_b) annotation (Line(points={{-26,-56},
            {-26,-66},{8,-66}}, color={0,127,255}));
    connect(radiator.port_a, sensor_T_1.port) annotation (Line(points={{28,-66},
            {46,-66},{46,-56}}, color={0,127,255}));
    connect(pump.inlet, tank.ports[1]) annotation (Line(points={{-56,-6},{-65.55,
            -6},{-65.55,6}}, color={0,127,255}));
  end HeatingSystem;

  model Test1 "Prescribed inputs, initial values"

    Modelica_Fluid.Test.TestOverdeterminedSteadyStateInit.HeatingSystem plant 
                annotation (Placement(transformation(extent={{0,0},{20,20}},
            rotation=0)));
    Modelica.Blocks.Sources.Step valveOpening(
      height=0.1,
      startTime=2000,
      offset=0.9)   annotation (Placement(transformation(extent={{-40,0},{-20,
              20}}, rotation=0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                        graphics),
      experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Initial equations with initial values for the states close to the steady state are selected for all components.
<p>
The simulation initializes and runs for 6000 seconds without problems.
</html>"));
  equation
    connect(valveOpening.y, plant.valvePosition) annotation (Line(points={{-19,10},
            {-0.8,10}},     color={0,0,127}));
  end Test1;

  model Test2 "Prescribed inputs, all derivatives equal to zero"
    extends Test1(plant(
        tank(initType=Modelica_Fluid.Types.Init.SteadyState),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(initType=Modelica_Fluid.Types.Init.SteadyState)));
    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for all components.
<p>
The simulation initializes and runs. However, due to the mathematical structure of the model, the initial level and temperature of the tank are not defined: any initial level and temperature in fact would satisfy the initial equations. Therefore, it is not possible for the user to actually set those values precisely. Accordingly, it turns out that when a simulation is run, the initial values of those two states in the simulation depend in an unpredictable fashion from the start values; in other words, if the start value for the level and temperature are changed, the initial values also change: this should not happen if the steady-state initialization problem had a well-defined, unique solution, at least if the start values are chosen close to the solution. 
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  end Test2;

  model Test3 "Prescribed inputs, all derivatives equal to zero"
    extends Test1(plant(
        tank(initType=Modelica_Fluid.Types.Init.SteadyState),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(initType=Modelica_Fluid.Types.Init.SteadyState)));
    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for all components. Moreover, additional initial equations are given to specify initial values for the tank temperature and tank pressure are specified at the system level (this is necessary because the tank components currently does not support that initialization option).
<p>The resulting initialization problem has more equations than unknowns states, but still has a unique solution. Unfortunately, Dymola 6.0d cannot find it. Also note that the equation Dymola suggests to remove is in fact necessary to determine the plant.valve.Av parameter.
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  initial equation
   // plant.tank.level = 1;
   // plant.tank.medium.T = Modelica.SIunits.Conversions.from_degC(20);
  end Test3;

  model Test4
    "Prescribed inputs, all derivatives equal to zero for the pipes, initial states given for the tank"
    extends Test1(plant(
        tank(initType=Modelica_Fluid.Types.Init.InitialValues),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(initType=Modelica_Fluid.Types.Init.SteadyState)));
    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for the pipe components, initial values are specified for the tank. Due to the mathematical structure of the tank model, steady state conditions can be obtained with arbitrary initial pressure and temperature, so this initial condition equations actually correspond to a steady-state condition. However, a user which does not know the exact mathematical structure of the model could not be aware of that.
<p>
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  end Test4;

  model Test5
    "Prescribed inputs, all derivatives equal to zero, zero pressure loss in the radiator"
    extends Test1(plant(
        tank(initType=Modelica_Fluid.Types.Init.SteadyState),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(initType=Modelica_Fluid.Types.Init.SteadyState,
                 redeclare package WallFriction = 
              Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction)));

    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for all components. The model of the radiator pipe has zero pressure losses.
<p>
The radiator pipe has no pressure losses in the momentum balances, so the pressure corresponding to its PortVolume component is equal to the pressure of the tank, leading to an index 2 DAE. One pressure state is removed by the dummy derivative algorithm, then Dymola complains there is one initial equation too many. In fact, the initialization problem is overdetermined, but still has infinitely many solutions. Also note that the equation Dymola suggests to remove is in fact necessary to determine the plant.valve.Av parameter.
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  end Test5;

  model Test6
    "Prescribed inputs, all derivatives equal to zero, zero pressure loss in the radiator"
    extends Test1(plant(
        tank(initType=Modelica_Fluid.Types.Init.SteadyState),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(initType=Modelica_Fluid.Types.Init.SteadyState,
                 redeclare package WallFriction = 
              Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction)));
    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for all components, plus additional initial equations to set the initial level and temperature of the tank. The model of the radiator pipe has zero pressure losses.
<p>
The radiator pipe has no pressure losses in the momentum balances, so the pressure corresponding to its PortVolume component is equal to the pressure of the tank, leading to an index 2 DAE. One pressure state is removed by the dummy derivative algorithm, then Dymola complains there are three initial equations too many. In fact, the initialization problem is overdetermined, but has one unique solution. Also note that the equation Dymola suggests to remove is in fact necessary to determine the plant.valve.Av parameter.
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  initial equation
    plant.tank.level = 1;
    plant.tank.medium.T = Modelica.SIunits.Conversions.from_degC(20);
  end Test6;

  model Test7
    "Prescribed inputs, carefully selected initial conditions, zero pressure loss in the radiator"
    extends Test1(plant(
        tank(initType=Modelica_Fluid.Types.Init.InitialValues),
        pipe(initType=Modelica_Fluid.Types.Init.SteadyState),
        radiator(redeclare package WallFriction = 
              Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction,
            initType=Modelica_Fluid.Types.Init.NoInit)));
    annotation (
      Documentation(info="<html>
Initial equations for steady-state are selected for all components, plus additional initial equations to set the initial level and temperature of the tank. The model of the radiator pipe has zero pressure losses.
<p>
The initial conditions here have been carefully crafted: 
<ul>
<li>Steady-state conditions for the tank have been replaced by initial states, because otherwise the problem has multiple solutions
<li>Steady-state conditions for the radiator have been replaced by a steady-state condition only on the thermal state (enthalpy), not on the pressure that has been removed by the dummy derivative algorithm. 
</ul>
The simulation initializes and runs correctly.
<p>
However, it is apparent that such a special treatment is out of the question for more complex models, whose mathematical structure might be much harder to understand. It is also extremely inconvenient that the user has to change initial conditions depending on index reduction: a user shouldn't need to know the inner mathematical details of a model (in this case, the fact that there is actually a state at the port of two directly connected components, due to a special selection of the WallFriction parameter) in order to be able to run it.
<p> 
On the other hand, it is quite straightforward to set up overdetermined conditions such as those of Test6, corresponding to steadyState for every model, plus additional equations to actually set the tank level and temperature. The latter can be easily added when one discovers that the steadyState option is not enough to get the desired initial values of those two variables.
</html>
"),   experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false));
  initial equation
    der(plant.radiator.volume.medium.h) = 0;
  end Test7;
end TestOverdeterminedSteadyStateInit;
