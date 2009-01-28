within Modelica_Fluid.Examples;
model BranchingDynamicPipes
  "Multi-way connections of pipes with dynamic momentum balance, pressure wave and flow reversal"
  import Modelica_Fluid;
extends Modelica.Icons.Example;
replaceable package Medium=Modelica.Media.Air.MoistAir;
//replaceable package Medium=Modelica.Media.Water.StandardWater;

    annotation (extent=[-90,-86; -70,-66], Diagram(coordinateSystem(
          preserveAspectRatio=true,  extent={{-100,-100},{100,100}}), graphics),
    Documentation(info="<html>
<p>
This model demonstrates the use of distributed pipe models with dynamic energy, mass and momentum balances. 
At time=2s the pressure of boundary4 jumps, which causes a pressure wave and flow reversal.
</p>
<p>
Change system.momentumDynamics on the Assumptions tab of the system object from DynamicFreeInitial to SteadyState,
in order to assume a steady-state momentum balance. This is the default for all models of the library.
</p>
<p>
Change the Medium from MoistAir to StandardWater, in order to investigate a medium with significantly different density. 
Note the static head caused by the elevation of the pipes.
</p>
 
<p>
Note, pipe4.modelStructure = av_b, i.e., the pipe has no volume at port_b.
It is not possible to have a volume at port_b, since otherwise the pressure of the volume is
defined by the connected boundary source. This in turn means that the
derivative of the pressure of the boundary source is needed, since the volume
requires this derivative. It is, however, not possible to compute this
derivative because the input pressure is changing disontinuously and its
derivative would be a dirac impulse.
</p>
</html>"));
                                     annotation (points=[-22,42; -2,42; -2,32;
        7.8,32],
      style(color=69, rgbcolor={0,127,255}));

  inner Modelica_Fluid.System system(momentumDynamics=Modelica_Fluid.Types.Dynamics.DynamicFreeInitial) 
    annotation (Placement(transformation(extent={{-90,70},{-70,90}},  rotation=
            0)));
  Modelica_Fluid.Sources.Boundary_pT boundary1(nPorts=1,
    redeclare package Medium = Medium, p=150000)                    annotation (Placement(
        transformation(extent={{-10,-10},{10,10}},    rotation=90,
        origin={0,-80})));
  Modelica_Fluid.Pipes.DynamicPipe pipe1(
    redeclare package Medium=Medium,
    use_T_start=true,
    nNodes=5,
    diameter=2.54e-2,
    m_flow_start=0.02,
    p_a_start=150000,
    p_b_start=130000,
    height_ab=50,
    length=50) 
            annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={0,-50})));
 Modelica_Fluid.Pipes.DynamicPipe pipe2(
    redeclare package Medium = Medium,
    use_T_start=true,
    nNodes=5,
    redeclare model HeatTransfer = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.LocalPipeFlowHeatTransfer,
    use_HeatTransfer=true,
    diameter=2.54e-2,
    m_flow_start=0.01,
    length=50,
    height_ab=25,
    p_a_start=130000,
    p_b_start=120000) 
            annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={-20,-10})));

  Modelica_Fluid.Pipes.DynamicPipe pipe3(
    redeclare package Medium=Medium,
    use_T_start=true,
    nNodes=5,
    diameter=2.54e-2,
    m_flow_start=0.01,
    length=25,
    p_a_start=130000,
    p_b_start=120000,
    height_ab=25) 
            annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={20,-10})));
  Modelica_Fluid.Pipes.DynamicPipe pipe4(
    redeclare package Medium=Medium,
    use_T_start=true,
    nNodes=5,
    modelStructure=Modelica_Fluid.Types.ModelStructure.av_b,
    diameter=2.54e-2,
    m_flow_start=0.02,
    p_a_start=120000,
    p_b_start=100000,
    height_ab=50,
    length=50) 
            annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=90,
        origin={0,30})));
  Modelica_Fluid.Sources.Boundary_pT boundary4(nPorts=1,
    redeclare package Medium = Medium,
    use_p_in=true,
    use_T_in=false,
    p=100000)                                            annotation (Placement(
        transformation(extent={{10,-10},{-10,10}}, rotation=90,
        origin={0,60})));
  Modelica.Blocks.Sources.Ramp ramp1(
    offset=1e5,
    duration=0,
    startTime=2,
    height=1e5) annotation (Placement(transformation(extent={{-40,70},{-20,90}},
          rotation=0)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow[
                                              pipe2.nNodes] heat2(Q_flow=200*
        pipe2.dxs) 
    annotation (Placement(transformation(extent={{-60,-20},{-40,0}},  rotation=
            0)));
equation
  connect(ramp1.y, boundary4.p_in) annotation (Line(
      points={{-19,80},{-8,80},{-8,72}},
      color={0,0,127},
      thickness=0.5));
  connect(boundary1.ports[1],pipe1. port_a) annotation (Line(
      points={{6.66134e-016,-70},{0,-70},{0,-60},{-6.12323e-016,-60}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe1.port_b,pipe2. port_a) annotation (Line(
      points={{6.12323e-016,-40},{0,-40},{0,-30},{-20,-30},{-20,-20}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe1.port_b,pipe3. port_a) annotation (Line(
      points={{6.12323e-016,-40},{0,-40},{0,-30},{20,-30},{20,-20}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe2.port_b,pipe4. port_a) annotation (Line(
      points={{-20,0},{-20,0},{-20,10},{0,10},{0,16},{0,20},{-6.12323e-016,20}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe3.port_b,pipe4. port_a) annotation (Line(
      points={{20,0},{20,0},{20,10},{0,10},{0,16},{0,20},{-6.12323e-016,20}},
      color={0,127,255},
      thickness=0.5));
  connect(pipe4.port_b, boundary4.ports[1]) annotation (Line(
      points={{6.12323e-016,40},{6.12323e-016,50},{-8.88178e-016,50}},
      color={0,127,255},
      thickness=0.5));
  connect(heat2.port,pipe2. heatPorts) 
                                      annotation (Line(
      points={{-40,-10},{-25.2,-10},{-25.2,-9.9}},
      color={191,0,0},
      thickness=0.5));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics),
                       experiment(StopTime=10),
    Commands(file(ensureSimulated=true)=
        "Scripts/Examples/BranchingDynamicPipes/plotResults.mos" "plotResults"),
    experimentSetupOutput,
    Documentation(info="<html>
<p>
This model demonstrates the use of distributed pipe models with dynamic energy, mass and momentum balances. 
At time=2s the pressure of boundary4 jumps, which causes a pressure wave and flow reversal.
</p>
<p>
Change system.momentumDynamics on the Assumptions tab of the system object from DynamicFreeInitial to SteadyState,
in order to assume a steady-state momentum balance. This is the default for all models of the library.
</p>
<p>
Change the Medium from MoistAir to StandardWater, in order to investigate a medium with sigificantly different density. 
Note the static head caused by the elevation of the pipes.
</p>
</html>"));

end BranchingDynamicPipes;
