within Modelica_Fluid.Test;
package TestOverdeterminedSteadyStateInit 
  "Contains test cases to test overdetermined systems of initial equations" 
  model HeatingSystem "Simple model of a heating system" 
     replaceable package Medium = Modelica.Media.Water.StandardWater 
       extends Modelica.Media.Interfaces.PartialMedium;
    annotation (Diagram, Icon(Rectangle(extent=[-100,100; 100,-100], style(
              color=3, rgbcolor={0,0,255})), Text(
          extent=[-60,60; 60,-60],
          style(color=3, rgbcolor={0,0,255}),
          string="P")));
    Volumes.OpenTank tank(
      redeclare package Medium = Medium,
      p_static_at_port=true,
      area=0.01,
      V0=0.01,
      pipe_diameters={0.025},
      initType=Modelica_Fluid.Types.Init.InitialValues,
      height=2,
      level_start=1) 
                annotation (extent=[-76,6; -54,28]);
    Pumps.Pump pump(
      redeclare package Medium = Medium,
      N_nom=1500,
      N_const=1500,
      pin_start=1.1e5,
      pout_start=4.0e5,
      use_T_start=true,
      T_start=Modelica.SIunits.Conversions.from_degC(40),
      redeclare function flowCharacteristic = 
          Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.linearFlow (
            head_nom={60.0,0}, q_nom={0.0,0.02e-3}),
      m_flow_start=0.01,
      checkValve=false) 
      annotation (extent=[-58,-16; -38,4]);
    ControlValves.ValveIncompressible valve(
      redeclare package Medium = Medium,
      CvData=Modelica_Fluid.Types.CvTypes.OpPoint,
      p_nom=4e5,
      dp_nom=3e5,
      T_start=Modelica.SIunits.Conversions.from_degC(80),
      m_flow_nom=0.01) 
      annotation (extent=[42,-12; 58,4]);
    Modelica.Blocks.Interfaces.RealInput valvePosition 
      annotation (extent=[-128,-20; -88,20]);
    Modelica.Blocks.Interfaces.RealOutput circuitFlowRate(redeclare type 
        SignalType = Modelica.SIunits.MassFlowRate) 
      annotation (extent=[88,12; 108,32]);
    Sensors.MassFlowRate massFlowRate(redeclare package Medium = Medium) 
      annotation (extent=[-34,6; -14,-14]);
    Modelica.Thermal.HeatTransfer.FixedTemperature ambientTemperature(T=ambient.default_T_ambient) 
      annotation (extent=[-12,-40; 2,-26]);
    Modelica.Thermal.HeatTransfer.ThermalConductor thermalConductor1(G=1.6e3/20) 
      annotation (extent=[8,-40; 28,-56], rotation=90);
    Modelica.Thermal.HeatTransfer.FixedHeatFlow burner(Q_flow=1.6e3) 
      annotation (extent=[-2,12; 18,32]);
    inner Ambient ambient annotation (extent=[-100,80; -80,100]);
    Modelica_Fluid.Pipes.LumpedPipe pipe(
      redeclare package Medium = Medium,
      p_a_start=4.0e5,
      p_b_start=3.9e5,
      use_T_start=true,
      diameter=0.03,
      T_start=Modelica.SIunits.Conversions.from_degC(80),
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      use_nominal=true,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      length=2) 
      annotation (extent=[12,-14; 32,6]);
    
    Modelica_Fluid.Pipes.LumpedPipe radiator(
      p_a_start=1.1e5,
      p_b_start=1.05e5,
      use_T_start=true,
      redeclare package Medium = Medium,
      length=10,
      diameter=0.05,
      T_start=Modelica.SIunits.Conversions.from_degC(40),
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent,
      use_nominal=true,
      initType=Modelica_Fluid.Types.Init.InitialValues) 
      annotation (extent=[28,-76; 8,-56]);
    
    Modelica.Blocks.Interfaces.RealOutput hotWaterTemperature(redeclare type 
        SignalType = Modelica.SIunits.Temperature) 
      annotation (extent=[88,-28; 108,-8]);
    Modelica.Blocks.Interfaces.RealOutput coldWaterTemperature(redeclare type 
        SignalType = Modelica.SIunits.Temperature) 
      annotation (extent=[88,-78; 108,-58]);
    Modelica_Fluid.Sensors.TemperatureOnePort sensor_T_1(
                                   redeclare package Medium = Medium) 
      annotation (extent=[56,-56; 36,-36]);
    Modelica_Fluid.Sensors.TemperatureOnePort sensor_T_2(
                                   redeclare package Medium = Medium) 
      annotation (extent=[-16,-56; -36,-36]);
    Modelica.Blocks.Interfaces.RealOutput tankLevel(redeclare type SignalType 
        = Modelica.SIunits.Height) annotation (extent=[90,60; 110,80]);
    Modelica.Blocks.Sources.RealExpression realExpression 
      annotation (extent=[-74,22; -54,42]);
  equation 
  tankLevel = tank.level;
    connect(valvePosition, valve.stemPosition) annotation (points=[-108,0; -86,
          0; -86,66; 50,66; 50,3.2],   style(color=74, rgbcolor={0,0,127}));
    connect(pump.outlet, massFlowRate.port_a) annotation (points=[-42,-2.8; -42,
          -3.4; -34,-3.4; -34,-4],        style(color=69, rgbcolor={0,127,255}));
    connect(massFlowRate.m_flow, circuitFlowRate) annotation (points=[-24,7;
          -24,38; 36,38; 36,22; 98,22],                             style(color=
           74, rgbcolor={0,0,127}));
    connect(massFlowRate.port_b, pipe.port_a) annotation (points=[-14,-4; 12,-4],
                style(color=69, rgbcolor={0,127,255}));
    connect(pipe.port_b, valve.port_a) annotation (points=[32,-4; 42,-4],
        style(color=69, rgbcolor={0,127,255}));
    connect(thermalConductor1.port_b, radiator.thermalPort) annotation (points=[18,-56;
          18,-60.6],                  style(color=42, rgbcolor={191,0,0}));
    connect(burner.port, pipe.thermalPort) annotation (points=[18,22; 22,22; 22,
          1.4],   style(color=42, rgbcolor={191,0,0}));
    connect(ambientTemperature.port, thermalConductor1.port_a) annotation (
        points=[2,-33; 18,-33; 18,-40], style(color=42, rgbcolor={191,0,0}));
    connect(sensor_T_1.T, hotWaterTemperature) annotation (points=[39,-46; 39,
          -80; 72,-80; 72,-18; 98,-18], style(color=74, rgbcolor={0,0,127}));
    connect(sensor_T_2.T, coldWaterTemperature) annotation (points=[-33,-46;
          -33,-88; 76,-88; 76,-68; 98,-68], style(color=74, rgbcolor={0,0,127}));
    connect(radiator.port_a, valve.port_b) annotation (points=[28,-66; 68,-66; 68,-4;
          58,-4], style(color=69, rgbcolor={0,127,255}));
    connect(pump.inlet, radiator.port_b) annotation (points=[-56,-8; -56,-66; 8,
          -66], style(color=69, rgbcolor={0,127,255}));
    connect(sensor_T_2.port, radiator.port_b) annotation (points=[-26,-56; -26,
          -66; 8,-66], style(color=69, rgbcolor={0,127,255}));
    connect(radiator.port_a, sensor_T_1.port) annotation (points=[28,-66; 46,
          -66; 46,-56], style(color=69, rgbcolor={0,127,255}));
    connect(pump.inlet, tank.ports[1]) annotation (points=[-56,-8; -65,-8; -65,
          5.45], style(color=69, rgbcolor={0,127,255}));
  end HeatingSystem;
  
  model Test1 "Prescribed inputs, initial values" 
    
    Modelica_Fluid.Test.TestOverdeterminedSteadyStateInit.HeatingSystem plant 
                annotation (extent=[0,0; 20,20]);
    Modelica.Blocks.Sources.Step valveOpening(
      height=0.1,
      startTime=2000,
      offset=0.9)   annotation (extent=[-40,0; -20,20]);
    annotation (Diagram,
      experiment(StopTime=6000, Tolerance=1e-006),
      experimentSetupOutput(equdistant=false),
      Documentation(info="<html>
Initial equations with initial values for the states close to the steady state are selected for all components.
<p>
The simulation initializes and runs for 6000 seconds without problems.
</html>"));
  equation 
    connect(valveOpening.y, plant.valvePosition) annotation (points=[-19,10; 
          -0.8,10], style(color=74, rgbcolor={0,0,127}));
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
