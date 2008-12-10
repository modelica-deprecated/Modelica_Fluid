within Modelica_Fluid.Test;
package TestOverdeterminedSteadyStateInit
  "Contains test cases to test overdetermined systems of initial equations"
  model DistributedPipeLumpedPressureInitialization
    "Steady-state initialization of a distributed pipe"

    Modelica_Fluid.Sources.FixedBoundary source(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      use_T=false,
      p=10000000,
      h=2e6) 
      annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
    Pipes.DistributedPipe pipe(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      h_start=2e6,
      diameter=0.05,
      length=200,
      use_T_start=false,
      modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb,
      lumpedPressure=true,
      p_a_start=10000000,
      p_b_start=9900000,
      nNodes=5) 
      annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
    Modelica_Fluid.Valves.ValveCompressible valve(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      Av=1e-3,
      dp_nominal=10000000,
      m_flow_nominal=10) 
      annotation (Placement(transformation(extent={{0,-10},{20,10}})));
    Modelica_Fluid.Sources.FixedBoundary sink(redeclare package Medium = 
          Modelica.Media.Water.StandardWaterOnePhase, p=9500000) 
                annotation (Placement(transformation(extent={{60,-10},{40,10}})));
    Modelica.Blocks.Sources.Ramp ramp(
      offset=1,
      duration=0.1,
      height=-0.5,
      startTime=2) 
                annotation (Placement(transformation(extent={{46,30},{26,50}})));
    inner Modelica_Fluid.System system(initType=Modelica_Fluid.Types.Init.SteadyState) 
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
    discrete Modelica.SIunits.MassFlowRate m_flow_initial;
  equation
    when time > 0.1 then
      m_flow_initial = valve.port_a.m_flow;
    end when;
    if pipe.initType == Modelica_Fluid.Types.Init.SteadyState or 
       pipe.dynamicsType == Modelica_Fluid.Types.Dynamics.SteadyState then
      when time > 1 then
        assert(abs(valve.port_a.m_flow - m_flow_initial) < 1e-3, "!!!THE SIMULATION DID NOT START IN STEADY-STATE!!!");
      end when;
    end if;
    connect(source.ports[1], pipe.port_a)         annotation (Line(
        points={{-60,0},{-40,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipe.port_b, valve.port_a)               annotation (Line(
        points={{-20,0},{0,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(valve.port_b, sink.ports[1])                          annotation (Line(
        points={{20,0},{40,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(ramp.y, valve.opening)               annotation (Line(
        points={{25,40},{10,40},{10,8}},
        color={0,0,127},
        smooth=Smooth.None));

    annotation (Documentation(info="<html>
All pressure states of the pipe are lumped into one. 
The steady-state initial conditions become overdetermined as they are now specified nNodes times for the same pressure state.
The initial equations are consistent however and a tool shall reduce them appropriately.
</html>"),
    Diagram(coordinateSystem(preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics),
      experiment(StopTime=4),
      experimentSetupOutput);
  end DistributedPipeLumpedPressureInitialization;


  model Test1 "Prescribed inputs, initial values"

    Modelica_Fluid.Examples.HeatingSystem plant(
        tank(initType=Modelica_Fluid.Types.Init.InitialValues),
        pipe(initType=Modelica_Fluid.Types.Init.InitialValues),
        radiator(initType=Modelica_Fluid.Types.Init.InitialValues)) 
                annotation (Placement(transformation(extent={{0,0},{20,20}},
            rotation=0)));
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
                 redeclare model PressureDrop = 
              Modelica_Fluid.Pipes.BaseClasses.PressureDrop.NominalPressureDrop(dp_nominal=0,smoothFlowReversal=true))));

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
                 redeclare model PressureDrop = 
              Modelica_Fluid.Pipes.BaseClasses.PressureDrop.NominalPressureDrop(dp_nominal=0,smoothFlowReversal=true))));
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
        radiator(redeclare model PressureDrop = 
              Modelica_Fluid.Pipes.BaseClasses.PressureDrop.NominalPressureDrop(dp_nominal=0,smoothFlowReversal=true),
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
