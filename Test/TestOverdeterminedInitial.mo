within Modelica_Fluid.Test;
package TestOverdeterminedInitial
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
      useLumpedPressure=true,
      nNodes=5,
      p_a_start=10000000,
      p_b_start=9900000,
      modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb) 
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
    inner Modelica_Fluid.System system(energyDynamics=Modelica_Fluid.Types.Dynamics.SteadyStateInitial) 
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
    discrete Modelica.SIunits.MassFlowRate m_flow_initial;
  equation
    when time > 0.1 then
      m_flow_initial = valve.port_a.m_flow;
    end when;
    if pipe.energyDynamics >= Modelica_Fluid.Types.Dynamics.SteadyStateInitial and 
       pipe.massDynamics >= Modelica_Fluid.Types.Dynamics.SteadyStateInitial then
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
            extent={{-100,-100},{100,100}}), graphics={Text(
            extent={{-100,-20},{100,-40}},
            lineColor={0,0,255},
            textString=
                "Problem: pipe.medium.p[1:5] are equal and have initial equations der(medium.p)=zeros(5);"),
            Text(
            extent={{-76,-40},{80,-58}},
            lineColor={0,0,255},
            textString=
                "A translator should remove consistently overdetermined initial equations.")}),
      experiment(StopTime=4),
      experimentSetupOutput);
  end DistributedPipeLumpedPressureInitialization;

  model DistributedPipeInitialValues
    "Initialization of a distributed pipe with values"

    Modelica_Fluid.Sources.FixedBoundary source(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      use_T=false,
      h=2e6,
      p=10000000) 
      annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
    Pipes.DistributedPipe pipe(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      h_start=2e6,
      diameter=0.05,
      length=200,
      use_T_start=false,
      nNodes=5,
      modelStructure=Modelica_Fluid.Types.ModelStructure.av_vb,
      p_a_start=10000000,
      p_b_start=9900000) 
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
    inner Modelica_Fluid.System system(energyDynamics=Modelica_Fluid.Types.Dynamics.FixedInitial) 
      annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
    discrete Modelica.SIunits.MassFlowRate m_flow_initial;
  equation
    when time > 0.1 then
      m_flow_initial = valve.port_a.m_flow;
    end when;
    if pipe.energyDynamics >= Modelica_Fluid.Types.Dynamics.SteadyStateInitial and 
       pipe.massDynamics >= Modelica_Fluid.Types.Dynamics.SteadyStateInitial then
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
            extent={{-100,-100},{100,100}}), graphics={
          Text(
            extent={{-100,-20},{100,-40}},
            lineColor={0,0,255},
            textString=
                "Problem: pipe.medium[1].p is equal to source.p and  has a consistent initial value  of 100 bar;"),
          Text(
            extent={{-76,-36},{76,-54}},
            lineColor={0,0,255},
            textString=
                "A translator should remove consistently overdetermined initial equations."),
          Text(
            extent={{-100,-64},{90,-84}},
            lineColor={0,0,255},
            textString=
                "Work-around 2: change system.energyDynamics from FixedInitial to DynamicFreeInitial"),
          Text(
            extent={{-100,-54},{42,-74}},
            lineColor={0,0,255},
            textString=
                "Work-around 1: change pipe.modelStructure from av_vb to a_vb")}),
      experiment(StopTime=4),
      experimentSetupOutput);

  end DistributedPipeInitialValues;
end TestOverdeterminedInitial;
