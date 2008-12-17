within Modelica_Fluid.Examples;
model InverseParameterization
  "Demonstrates the parameterization of a pump and a pipe for given nominal values"
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  replaceable package Medium = Modelica.Media.Water.StandardWater;
      //Modelica.Media.Water.ConstantPropertyLiquidWater;

  Modelica_Fluid.Sources.Boundary_pT source(
    redeclare package Medium = Medium,
    nPorts=1,
    use_p_in=false,
    p=100000) 
    annotation (Placement(transformation(extent={{-76,-6},{-64,6}},  rotation=0)));
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,1.9e5; 10,2.1e5]) 
    annotation (Placement(transformation(extent={{-60,30},{-40,50}}, rotation=0)));
  Modelica_Fluid.Machines.ControlledPump pump(
    m_flow_nominal=1,
    control_m_flow=false,
    use_p_set=true,
    redeclare package Medium = Medium,
    p_a_nominal=100000,
    p_b_nominal=200000) 
    annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
  Modelica_Fluid.Pipes.StaticPipe pipe1(
    redeclare package Medium = Medium,
    diameter=2.54e-2,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.NominalPressureLoss (
          m_flow_nominal=1, dp_nominal=100000),
    length=50)        annotation (Placement(transformation(extent={{20,-10},{40,
            10}}, rotation=0)));

  Modelica_Fluid.Sources.Boundary_pT sink1(redeclare package Medium = Medium, p=
       100000) 
             annotation (Placement(transformation(extent={{76,-6},{64,6}},
          rotation=0)));

  inner Modelica_Fluid.System system 
                        annotation (Placement(transformation(extent={{-90,70},{
            -70,90}},  rotation=0)));
  Modelica_Fluid.Pipes.StaticPipe pipe2(
    redeclare package Medium = Medium,
    diameter=2.54e-2,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.NominalPressureLoss (
        m_flow_nominal=1,
        show_Re=true,
        dp_nominal=100000),
    length=50)        annotation (Placement(transformation(extent={{20,-50},{40,
            -30}},rotation=0)));
  Modelica_Fluid.Sources.Boundary_pT sink2(
    redeclare package Medium = Medium, p=200000) 
             annotation (Placement(transformation(extent={{76,-46},{64,-34}},
          rotation=0)));
equation
  connect(pipe1.port_b, sink1.ports[1]) 
                                       annotation (Line(
      points={{40,0},{64,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(source.ports[1], pump.port_a) annotation (Line(
      points={{-64,0},{-40,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pump.port_b, pipe1.port_a) 
                                    annotation (Line(
      points={{-20,0},{20,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(p_table.y, pump.p_set) annotation (Line(
      points={{-39,40},{-25,40},{-25,8.2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pipe2.port_b, sink2.ports[1]) annotation (Line(
      points={{40,-40},{64,-40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe2.port_a, pump.port_b) annotation (Line(
      points={{20,-40},{0,-40},{0,0},{-20,0}},
      color={0,127,255},
      smooth=Smooth.None));

  annotation (
    Commands(file(ensureSimulated=true)="Scripts/Examples/InverseParameterization/plotResults.mos"
        "plotResults"),
    Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
            100,100}}),
            graphics),
    experiment(StopTime=10, NumberOfIntervals=10000),
    Documentation(info="<html>
<p>
A pump and two pipes are parameterized with simple nominal values. The pump controls a pressure ramp from 1.9 bar to 2.1 bar. 
This causes an appropriate ramp on the mass flow rate for pipe1, which has a boundary pressure of 1 bar. 
Flow reversal occurs in pipe2, which has a boundary pressure of 2 bar.
</p>
<p>
The Command plotResults can be used to see the pump speed N, which is controlled ideally to obtain the pressure ramp.
Moreover the Reynolds number as well as m_flow_turbulent and dp_turbulent are plotted for pipe2.
</p>
<p>
Next the pressureLoss models of the pipes could be investigated for <tt>length_nominal</tt>, 
which are obtained internally to fulfill the nominal pressure loss for given pipe diameter and roughness. 
Once the geometry has been designed, the NominalPressureLoss correlations can easily be replaced with 
WallFrictionPressureLoss correlations. Similarily the ControlledPump can be replaced with a PrescribedPump 
to investigate a real controller or with a Pump with rotational shaft to investigate inertia effects. 
</p>
</html>"));
end InverseParameterization;
